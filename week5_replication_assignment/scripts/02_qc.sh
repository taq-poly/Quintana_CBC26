#!/usr/bin/env bash
# 02_qc.sh -- FastQC on each subsampled FASTQ, then MultiQC summary report.
set -euo pipefail

THREADS=${THREADS:-4}
SUB_DIR=data/subsampled
QC_DIR=results/qc
mkdir -p "$QC_DIR"

fastqc -t "$THREADS" -o "$QC_DIR" "$SUB_DIR"/*.fastq.gz

multiqc "$QC_DIR" -o "$QC_DIR" -n multiqc_report

echo ">> QC reports written to $QC_DIR/ (see multiqc_report.html)"
