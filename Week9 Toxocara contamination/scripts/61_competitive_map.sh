#!/bin/bash
# Split competitive mapping of one sample's non-worm reads against worm / microbes / dog / human.
# Each target gets its own minimap2 index (RAM); pairs are assigned to the best target afterwards
# by 62_assign_competitive.py using a preset-independent score (aligned bases - 2*NM).
# Usage: 61_competitive_map.sh SRRxxxx [R1 R2]
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen; M=$B/refs/mm2; S=$1
R1=$(realpath ${2:-$B/results/$S/nonworm_1.fq.gz}); R2=$(realpath ${3:-$B/results/$S/nonworm_2.fq.gz})
O=$B/results/$S/compete; mkdir -p $O; cd $O
log () { echo "[$(date +%H:%M:%S)] $S compete $*"; }
for t in worm microbes dog human; do
  [ -s $t.bam ] && continue
  case $t in dog|human) x=${EUK_PRESET:-splice:sr} ;; *) x=sr ;; esac   # EUK_PRESET=sr for DNA reads
  log "map $t ($x)"
  minimap2 -a -x $x -t 24 --secondary=no $M/$t.mmi $R1 $R2 2> $t.mm2.log \
    | samtools view -@ 4 -b -F 0x904 -o $t.unsorted.bam - && mv $t.unsorted.bam $t.bam
done
log assign
python3 $B/scripts/62_assign_competitive.py $S
for t in worm microbes dog human; do
  samtools view -@ 4 -b -N assigned_$t.names $t.bam | samtools sort -@ 4 -o assigned_$t.bam - && samtools index assigned_$t.bam
  samtools coverage -q 1 assigned_$t.bam > assigned_$t.cov.tsv
done
log done
