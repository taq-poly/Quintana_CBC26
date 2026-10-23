#!/bin/bash
# Genome blobplot of the published T. canis assembly (GCA_000803305.1):
#   coverage  = SRR1706691 (170 bp insert, HiSeq 2x100) mapped with minimap2 -ax sr
#   taxonomy  = DIAMOND blastx of up to 10 x 100 kb chunks per contig vs UniRef50 (BTK pipeline approach)
#   extra     = Kraken2 contig class + dog/human minimap2 fraction from reference QC
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen; O=$B/results/genome_blob; mkdir -p $O; cd $O
BTK=$B/env_btk/bin/blobtools
DMD=/opt/conda/envs/env_diamond/bin/diamond   # >=2.1 needed: 2.0.4 rejects NCBI "domain" rank
log () { echo "[$(date +%H:%M:%S)] blob $*"; }
ASM=$B/refs/refqc/tcanis.fa
[ -s $ASM ] || seqkit seq -m 1 $B/refs/tcanis.fna.gz > $ASM

# 1. UniRef50 DIAMOND db with taxonomy (taxid = cluster common taxon from header TaxID=)
U=$B/db/uniref50
if [ ! -s $U/uniref50.dmnd ]; then
  log "diamond makedb"
  zcat $U/uniref50.fasta.gz | awk 'BEGIN{print "accession\taccession.version\ttaxid\tgi"} /^>/{a=substr($1,2); t=""; for(i=2;i<=NF;i++) if($i~/^TaxID=/){t=substr($i,7)}; if(t=="") t=1; sub(/^UniRef50_/,"",a); print a"\t"a"\t"t"\t0"}' > $U/taxonmap.tsv
  $DMD makedb --in $U/uniref50.fasta.gz -d $U/uniref50 -p ${MAKEDB_THREADS:-32} \
    --taxonmap $U/taxonmap.tsv --taxonnodes $B/db/taxdump/nodes.dmp --taxonnames $B/db/taxdump/names.dmp > $U/makedb.log 2>&1
fi

# 2. Chunk contigs: 100 kb windows, at most 10 evenly spaced per contig, contigs >= 300 bp
if [ ! -s chunks.fa ]; then
  log chunk
  python3 - "$ASM" <<'PY'
import sys
C=100000; M=10
def fa(p):
    n=None; s=[]
    for l in open(p):
        if l[0]==">":
            if n: yield n,"".join(s)
            n=l[1:].split()[0]; s=[]
        else: s.append(l.strip())
    if n: yield n,"".join(s)
out=open("chunks.fa","w")
for n,s in fa(sys.argv[1]):
    if len(s)<300: continue   # BTK default is 1 kb; lowered so short contaminant contigs get hits
    starts=list(range(0,len(s),C))
    if len(starts)>M:
        step=len(starts)/M; starts=[starts[int(i*step)] for i in range(M)]
    for st in starts:
        out.write(f">{n}_-_{st}\n{s[st:st+C]}\n")
PY
fi

# 3. DIAMOND blastx (BTK settings), then shift chunk coordinates back onto contigs
if [ ! -s diamond.contig.tsv ]; then
  log "diamond blastx"
  $DMD blastx -q chunks.fa -d $U/uniref50 -p 24 --fast -b 1 -c 4 \
    --evalue 1e-25 --max-target-seqs 10 --tmpdir $O \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    -o diamond.chunks.tsv > diamond.log 2>&1
  awk 'BEGIN{OFS="\t"} {split($1,a,"_-_"); $1=a[1]; $4=a[1]; $10+=a[2]; $11+=a[2]; print}' diamond.chunks.tsv > diamond.contig.tsv
fi

# 4. Map DNA reads
if [ ! -s SRR1706691.bam.bai ]; then
  log minimap2
  minimap2 -ax sr -t 24 $ASM $B/reads/dna/SRR1706691_1.fastq.gz $B/reads/dna/SRR1706691_2.fastq.gz 2> minimap2.log \
    | samtools sort -@ 4 -m 1G -T $O/sorttmp -o SRR1706691.bam - && samtools index SRR1706691.bam && samtools index -c SRR1706691.bam  # BTK 4.4.5 expects a .csi
  samtools flagstat -@ 8 SRR1706691.bam > SRR1706691.flagstat
fi

# 5. BlobToolKit dataset + plots
if [ ! -d dataset ]; then
  log "blobtools create"
  $BTK create --fasta $ASM --taxdump $B/db/taxdump --taxid 6265 dataset > btk.log 2>&1
  $BTK add --cov SRR1706691.bam --threads 16 dataset >> btk.log 2>&1
  $BTK add --hits diamond.contig.tsv --taxrule bestsumorder --taxdump $B/db/taxdump dataset >> btk.log 2>&1
fi
mkdir -p plots
for v in blob cumulative snail; do
  $BTK view --view $v --format png --out plots dataset >> btk.log 2>&1 || log "view $v failed"
done
$BTK view --view blob --format svg --out plots dataset >> btk.log 2>&1 || true
log done
