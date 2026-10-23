#!/bin/bash
# Download + convert the 6 adult RNA-seq runs (PRJNA1041894)
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
cd /workspace/toxocara_screen/reads/raw
for r in SRR26866479 SRR26866491 SRR26866490 SRR26866487 SRR26866489 SRR26866488; do
  [ -s ${r}_1.fastq.gz ] && continue
  start=$(date +%s)
  prefetch --max-size 20G -O . $r
  fasterq-dump --split-files -e 16 -t /workspace/toxocara_screen/reads/tmp -O . $r/$r.sra
  pigz -p 16 ${r}_1.fastq ${r}_2.fastq
  rm -rf $r
  echo "$r done in $(( $(date +%s)-start ))s"
done
