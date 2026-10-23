#!/bin/bash
# Early start of 05_clean_reference.sh's minimap2 step (writes .tmp, renames when complete so 05 can skip it)
export PATH=/workspace/toxocara_screen/env/bin:$PATH
R=/workspace/toxocara_screen/refs; cd $R/refqc
for g in phix ecoli_k12 dog human; do
  [ -s $g.paf ] && continue
  minimap2 -x asm10 -t 24 --secondary=no $R/$g.fna.gz $R/tcanis.fna.gz > $g.paf.tmp 2> $g.minimap.log && mv $g.paf.tmp $g.paf && echo "refqc minimap $g ok"
done
