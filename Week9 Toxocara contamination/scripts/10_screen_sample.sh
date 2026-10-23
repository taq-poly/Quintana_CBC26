#!/bin/bash
# Per-sample screen: trim -> remove worm reads (HISAT2) -> Kraken2/Bracken, sendsketch, sourmash
# on the non-worm reads. Fractions are reported against total trimmed read pairs.
# Usage: 10_screen_sample.sh SRRxxxx
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen
S=$1; T=32
RAW=$B/reads/raw; O=$B/results/$S; mkdir -p $O
cd $O
log () { echo "[$(date +%H:%M:%S)] $S $*"; }

# 1. Trim. SRA strips Illumina headers so fastp can't auto-detect NextSeq: force poly-G trimming
#    (2-colour chemistry calls 'no signal' as G).
if [ ! -s trim_1.fq.gz ]; then
  log fastp
  fastp -i $RAW/${S}_1.fastq.gz -I $RAW/${S}_2.fastq.gz -o trim_1.fq.gz -O trim_2.fq.gz \
    --detect_adapter_for_pe --trim_poly_g --length_required 40 -w 16 \
    -j fastp.json -h fastp.html 2> fastp.log
fi

# 2. Remove worm reads. Spliced aligner so intron-spanning worm reads don't leak into 'non-target'.
#    Not using --un-conc-gz: it routes every read through hisat2's single-threaded wrapper (stalled
#    32 threads), and it would also count one-mate-mapped pairs as non-worm.
#    non-worm = both mates unmapped (-f 12); one-mate-mapped pairs kept for the chimera check.
if [ ! -s nonworm_1.fq.gz ]; then
  log hisat2
  hisat2 -p 24 -x $B/refs/idx/tcanis_hisat2 -1 trim_1.fq.gz -2 trim_2.fq.gz 2> hisat2.log \
    | samtools view -@ 4 -u -F 0x900 - | tee >(samtools view -@ 2 -b -f 1 -e '(flag.unmap && !flag.munmap) || (!flag.unmap && flag.munmap)' -o half_mapped.bam -) \
    | samtools fastq -@ 4 -f 12 -n -1 nonworm_1.fq.gz -2 nonworm_2.fq.gz -0 /dev/null -s /dev/null - 2> samtools_fastq.log
  samtools view -c half_mapped.bam > half_mapped_reads.txt
fi

# 3. Kraken2 + Bracken (species and genus)
if [ ! -s kraken.report ]; then
  log kraken2
  kraken2 --db $B/db/k2_pluspf16 --threads $T --paired --gzip-compressed --confidence 0.1 \
    --report kraken.report --output kraken.out nonworm_1.fq.gz nonworm_2.fq.gz 2> kraken.log
  gzip -f kraken.out
  for lvl in S G; do
    bracken -d $B/db/k2_pluspf16 -i kraken.report -o bracken_$lvl.tsv -r 75 -l $lvl -t 10 > /dev/null 2>&1 || true
  done
fi

# 4. sendsketch (JGI server, nt + RefSeq). Sends MinHash sketches only.
if [ ! -s sendsketch_nt.txt ]; then
  log sendsketch
  for db in nt refseq; do
    sendsketch.sh in=nonworm_1.fq.gz in2=nonworm_2.fq.gz $db records=40 minhits=3 \
      format=3 out=sendsketch_$db.txt > /dev/null 2>&1 || log "sendsketch $db failed"
  done
fi

# 5. sourmash via branchwater plugin (Rust, multithreaded; plain `sourmash sketch` took >17 min/sample).
#    Mates sketched separately then merged; k31 vs GTDB reps + NCBI viruses.
#    Threshold 3 kb (3 hashes at scaled=1000): hits <10 kb are low-confidence until confirmed by mapping.
D=$B/db/sourmash
if [ ! -s gather_euks.csv ] && [ -e $D/.complete ]; then
  log sourmash
  rm -f reads.sig.zip
  printf "name,genome_filename,protein_filename\n$S,nonworm_1.fq.gz,\n$S,nonworm_2.fq.gz,\n" > sketch.csv
  sourmash scripts manysketch sketch.csv -p dna,k=31,scaled=1000,abund -o reads.sig.zip -c 16 > /dev/null 2>&1
  for k in 31; do sourmash sig merge reads.sig.zip -k $k --name $S -o m$k.sig 2> /dev/null; done
  # NCBI eukaryote DBs are scaled=10000: a 2M-pair non-worm RNA sketch has only ~1.5k hashes at that
  # resolution, so trace eukaryotes (dog) are invisible to it -> eukaryotes are measured by mapping instead.
  for db in gtdb-reps-rs226-k31.dna.zip ncbi-viruses-2025.01.dna.k=31.sig.zip; do
    k=$(echo $db | grep -oE 'k=?[0-9]+' | grep -oE '[0-9]+'); tag=$(echo $db | sed -E 's/(-2025.01|-rs226)?[-.]+(dna|k).*//')
    sourmash scripts fastgather m$k.sig $D/$db -k $k -s 1000 -t 3000 -c 24 -o fg_$tag.csv > /dev/null 2>&1 || log "fastgather $tag failed"
  done
  python3 - <<'PYIN'
import csv, glob
rows = []
for f in sorted(glob.glob("fg_*.csv")):
    for r in csv.DictReader(open(f)):
        rows.append([f[3:-4], r["match_name"], r["intersect_bp"], r["unique_intersect_bp"], r["f_unique_weighted"],
                     r["average_abund"], r["match_containment_ani"], r["f_match"]])
with open("gather_euks.csv", "w") as o:   # combined table (name kept so the skip-check works)
    o.write("db,match_name,intersect_bp,unique_intersect_bp,f_unique_weighted,average_abund,match_containment_ani,f_match\n")
    for r in rows: o.write(",".join('"%s"' % x if "," in x else x for x in r) + "\n")
PYIN
fi

[ -s gather_euks.csv ] || log "sourmash skipped (DBs not ready) - rerun this script later"

# Summary counts
python3 - <<PY
import json
j=json.load(open("fastp.json"))
raw=j["summary"]["before_filtering"]["total_reads"]//2
trim=j["summary"]["after_filtering"]["total_reads"]//2
import gzip
nw=sum(1 for _ in gzip.open("nonworm_1.fq.gz","rt"))//4
open("counts.tsv","w").write(f"sample\traw_pairs\ttrim_pairs\tnonworm_pairs\tnonworm_frac\n$S\t{raw}\t{trim}\t{nw}\t{nw/trim:.4f}\n")
PY
log done
