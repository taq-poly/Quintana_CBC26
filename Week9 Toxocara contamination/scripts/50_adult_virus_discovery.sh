#!/bin/bash
# Co-assemble adult non-worm reads (MEGAHIT), find viral contigs (geNomad), completeness (CheckV),
# then map each adult back to the contigs for per-sample counts.
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen; O=$B/results/virus_discovery; mkdir -p $O; cd $O
A="SRR26866479 SRR26866491 SRR26866490 SRR26866487 SRR26866489 SRR26866488"
log () { echo "[$(date +%H:%M:%S)] vdisc $*"; }
if [ ! -s megahit/final.contigs.fa ]; then
  log megahit
  R1=$(for r in $A; do printf "$B/results/$r/nonworm_1.fq.gz,"; done | sed 's/,$//')
  R2=$(for r in $A; do printf "$B/results/$r/nonworm_2.fq.gz,"; done | sed 's/,$//')
  rm -rf megahit
  megahit -1 $R1 -2 $R2 -o megahit -t 24 -m 0.4 --min-contig-len 500 > megahit.log 2>&1
fi
if [ ! -d genomad ]; then
  log genomad
  genomad end-to-end --cleanup --threads 24 megahit/final.contigs.fa genomad $B/db/genomad_db > genomad.log 2>&1
fi
V=genomad/final.contigs_summary/final.contigs_virus.fna
if [ -s $V ] && [ ! -d checkv ]; then
  log checkv
  checkv end_to_end $V checkv -t 24 -d $B/db/checkv-db-v1.5 > checkv.log 2>&1 || checkv end_to_end $V checkv -t 24 -d $(ls -d $B/db/checkv-db* | head -1) > checkv.log 2>&1
fi
log "map back"
bowtie2-build -q --threads 16 megahit/final.contigs.fa contigs
for r in $A; do
  [ -s $r.contigs.bam ] && continue
  bowtie2 -p 24 --no-unal -x contigs -1 $B/results/$r/nonworm_1.fq.gz -2 $B/results/$r/nonworm_2.fq.gz 2> $r.bt2.log \
    | samtools sort -@4 -o $r.contigs.bam - && samtools index $r.contigs.bam
  samtools coverage -q 20 $r.contigs.bam > $r.contigs.cov.tsv
done
log done
