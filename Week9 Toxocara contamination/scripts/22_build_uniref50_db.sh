#!/bin/bash
# Build the UniRef50 DIAMOND db now (8 threads, alongside screening): runs only step 1 of 21_genome_blob.sh
set -euo pipefail
sed -n '1,/^fi$/p' /workspace/toxocara_screen/scripts/21_genome_blob.sh | sed '/^ASM=/,/seqkit seq/d' > ${TMPDIR:-/tmp}/makedb_only.sh
MAKEDB_THREADS=8 bash ${TMPDIR:-/tmp}/makedb_only.sh && echo "uniref50 dmnd ok"
