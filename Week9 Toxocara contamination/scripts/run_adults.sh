#!/bin/bash
# Orchestrate: wait for DBs -> clean reference + HISAT2 index -> screen each adult as its download finishes
set -uo pipefail
B=/workspace/toxocara_screen; cd $B
waitfor () { until grep -q "$1" $2; do sleep 60; done; }
waitfor "kraken ok" logs/01_refs_dbs.log
until [ -s refs/human.fna.gz ] && grep -q "ref human ok" logs/01_refs_dbs.log; do sleep 60; done
echo "[$(date)] reference QC start"; scripts/05_clean_reference.sh > logs/05_refqc.log 2>&1 || { echo "refqc failed"; exit 1; }
for r in SRR26866479 SRR26866491 SRR26866490 SRR26866487 SRR26866489 SRR26866488; do
  waitfor "$r done" logs/00_download.log
  echo "[$(date)] screen $r"; scripts/10_screen_sample.sh $r >> logs/10_screen.log 2>&1 || echo "[$(date)] $r FAILED"
done
echo "[$(date)] ALL ADULTS SCREENED (pass 1)"
# pass 2: sourmash once its DBs are complete
waitfor "plants-2025.01.dna.k=51.sig.zip ok" logs/01_refs_dbs.log && touch db/sourmash/.complete
for r in SRR26866479 SRR26866491 SRR26866490 SRR26866487 SRR26866489 SRR26866488; do
  scripts/10_screen_sample.sh $r >> logs/10_screen.log 2>&1 || echo "[$(date)] $r sourmash FAILED"
done
echo "[$(date)] ALL ADULTS SCREENED (sourmash)"
