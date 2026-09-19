#!/usr/bin/env bash
# run_all.sh -- runs the full pipeline end to end.
# NOTE: full read depth, no subsampling (~11.5 GB total across 6
# samples). Step 1 downloads samples concurrently (see
# scripts/01_fetch_data.sh) to make better use of available bandwidth.
set -euo pipefail

echo "== [1/4] Fetching reference + reads ================================"
bash scripts/01_fetch_data.sh

echo "== [2/4] QC ==========================================================="
bash scripts/02_qc.sh

echo "== [3/4] Salmon quantification ========================================"
bash scripts/03_quantify.sh

echo "== [4/4] DESeq2 differential expression ==============================="
Rscript scripts/04_deseq2.R

echo "== Generating CHECKSUMS.txt ==========================================="
{
  echo "# SHA-256 checksums -- generated $(date -u +%FT%TZ)"
  find data/reference data/raw -type f | sort | xargs sha256sum
  find results -type f \( -name "*.csv" -o -name "*.png" -o -name "*.txt" \) | sort | xargs sha256sum
} > CHECKSUMS.txt

echo "== Done. Key outputs: =================================================="
echo "  results/deseq2/DE_ivermectin_vs_control.csv"
echo "  results/deseq2/DE_ivermectin_vs_control_significant.csv"
echo "  results/deseq2/volcano.png, pca.png, sample_distance_heatmap.png"
echo "  results/qc/multiqc_report.html"
echo "  CHECKSUMS.txt"
