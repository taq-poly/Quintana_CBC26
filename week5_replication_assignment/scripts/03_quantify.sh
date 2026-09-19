#!/usr/bin/env bash
# 03_quantify.sh -- build a Salmon index for the T. canis reference
# transcriptome, then quantify each sample in PAIRED-END mode.
set -euo pipefail

THREADS=${THREADS:-4}
REF_DIR=data/reference
READS_DIR=data/raw
IDX_DIR=results/salmon_index
QUANT_DIR=results/quant
mkdir -p "$QUANT_DIR"

echo ">> Building Salmon index..."
salmon index \
  -t "$REF_DIR/tcanis_cds.fa" \
  -i "$IDX_DIR" \
  -k 31 \
  -p "$THREADS"

mapfile -t QUANT_LINES < <(tail -n +2 config/samples.tsv)
for line in "${QUANT_LINES[@]}"; do
  IFS=$'\t' read -r sample_id condition replicate srr <<< "$line"
  echo ">> Quantifying ${sample_id} (${condition}, rep ${replicate})..."
  salmon quant \
    -i "$IDX_DIR" \
    -l A \
    -1 "$READS_DIR/${sample_id}_R1.fastq.gz" \
    -2 "$READS_DIR/${sample_id}_R2.fastq.gz" \
    -p "$THREADS" \
    --validateMappings \
    --seqBias \
    -o "$QUANT_DIR/${sample_id}"
done

echo ">> Salmon quant.sf files written under $QUANT_DIR/<sample_id>/"
