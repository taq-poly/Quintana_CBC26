#!/usr/bin/env bash
# run_all.sh -- runs the full pipeline end to end.
# Expected wall-clock time: ~15-25 min on a laptop (4 cores), dominated by
# SRA download speed. See README.md for setup instructions.
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
  find data/reference data/subsampled -type f | sort | xargs sha256sum
  find results -type f \( -name "*.csv" -o -name "*.png" -o -name "*.txt" \) | sort | xargs sha256sum
} > CHECKSUMS.txt

echo "== Done. Key outputs: =================================================="
echo "  results/deseq2/DE_testis_vs_uterus.csv"
echo "  results/deseq2/DE_testis_vs_uterus_significant.csv"
echo "  results/deseq2/volcano.png, pca.png, sample_distance_heatmap.png"
echo "  results/qc/multiqc_report.html"
echo "  CHECKSUMS.txt"
