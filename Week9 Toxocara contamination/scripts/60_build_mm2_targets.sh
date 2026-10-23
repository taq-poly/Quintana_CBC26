#!/bin/bash
# minimap2 indexes for split competitive mapping (one index per target keeps peak RAM < ~12 GB)
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
R=/workspace/toxocara_screen/refs; M=$R/mm2; cd $M
# worm = cleaned assembly + full mito + co-assembled worm-false-positive kept out; contaminant contigs NOT included
[ -s worm.fa ] || cat $R/tcanis.clean.fa $R/tcanis_mito.fa > worm.fa
# microbes: candidate genomes + C. perfringens + PhiX + E. coli K-12, names prefixed with genome tag
if [ ! -s microbes.fa ]; then
  for f in $R/candidates/*.fna.gz $R/cperf.fna.gz $R/phix.fna.gz $R/ecoli_k12.fna.gz; do
    t=$(basename $f .fna.gz); zcat $f | sed "s/^>/>${t}|/" | cut -d' ' -f1
  done > microbes.fa
fi
for t in worm microbes; do [ -s $t.mmi ] || minimap2 -x sr -t 16 -d $t.mmi $t.fa 2> $t.mmi.log; echo "$t.mmi ok"; done
for t in dog human; do [ -s $t.mmi ] || minimap2 -x splice:sr -t 16 -d $t.mmi $R/$t.fna.gz 2> $t.mmi.log; echo "$t.mmi ok"; done
