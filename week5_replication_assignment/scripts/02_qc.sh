#!/usr/bin/env bash
# 02_qc.sh -- FastQC on each full-depth FASTQ (both mates), then MultiQC.
set -euo pipefail

THREADS=${THREADS:-4}
READS_DIR=data/raw
QC_DIR=results/qc
mkdir -p "$QC_DIR"

fastqc -t "$THREADS" -o "$QC_DIR" "$READS_DIR"/*.fastq.gz

multiqc "$QC_DIR" -o "$QC_DIR" -n multiqc_report

echo ">> QC reports written to $QC_DIR/ (see multiqc_report.html)"
