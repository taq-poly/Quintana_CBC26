#!/usr/bin/env bash
# 03_quantify.sh -- build a Salmon index for the D. immitis reference
# transcriptome, then pseudo-align + quantify each subsampled sample
# (single-end reads, as in the original GAIIx sequencing).
set -euo pipefail

THREADS=${THREADS:-4}
REF_DIR=data/reference
SUB_DIR=data/subsampled
IDX_DIR=results/salmon_index
QUANT_DIR=results/quant
mkdir -p "$QUANT_DIR"

echo ">> Building Salmon index..."
salmon index \
  -t "$REF_DIR/dimmitis_cds.fa" \
  -i "$IDX_DIR" \
  -k 25 \
  -p "$THREADS"

tail -n +2 config/samples.tsv | while IFS=$'\t' read -r sample_id tissue srr replicate; do
  echo ">> Quantifying ${sample_id} (${tissue}, rep ${replicate})..."
  salmon quant \
    -i "$IDX_DIR" \
    -l A \
    -r "$SUB_DIR/${sample_id}.fastq.gz" \
    -p "$THREADS" \
    --validateMappings \
    --seqBias \
    -o "$QUANT_DIR/${sample_id}"
done

echo ">> Salmon quant.sf files written under $QUANT_DIR/<sample_id>/"
