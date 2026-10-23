#!/bin/bash
# Pass 2: sourmash on the adults (rest of 10_screen_sample.sh is skipped as complete)
cd /workspace/toxocara_screen
for r in SRR26866479 SRR26866491 SRR26866490 SRR26866487 SRR26866489 SRR26866488; do
  scripts/10_screen_sample.sh $r >> logs/10_screen.log 2>&1 || echo "[$(date)] $r sourmash FAILED" >> logs/10_screen.log
done
echo "[$(date)] ALL ADULTS SOURMASH DONE" >> logs/run_adults.log
