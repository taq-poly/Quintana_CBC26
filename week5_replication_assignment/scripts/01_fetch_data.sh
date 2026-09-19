#!/usr/bin/env bash
# ==============================================================================
# 01_fetch_data.sh
#
# Fetches:
#   1. The Toxocara canis reference CDS transcriptome from WormBase ParaSite
#      (genome PRJNA248777, release WBPS19 -- the same assembly used by
#      Quintana et al. 2025, matching their Tcan_XXXXX gene IDs).
#   2. FULL-DEPTH raw paired-end RNA-seq reads for 4 samples (2 control +
#      2 ivermectin biological replicates -- the two smallest/fastest
#      replicates of each condition, chosen to keep full-depth download
#      practical) from NCBI SRA, BioProject PRJNA1041894 / study
#      SRP472648 (Quintana, Brewer & Jesudoss Chelladurai 2025,
#      Int J Parasitol Drugs Drug Resist 29:100614).
#
# Comparison: CONTROL vs IVERMECTIN (n=2 per group).
#
# NO SUBSAMPLING: every read of every run is downloaded and used, exactly
# matching the source study's per-sample depth (~7.7 GB total across the
# 4 selected samples).
#
# All 4 samples are fetched CONCURRENTLY in a single batch (default
# MAX_PARALLEL=4), to make better use of available bandwidth. Each
# sample's `prefetch` writes to its own subdirectory and
# `fasterq-dump`/`gzip`/`mv` operate on distinctly-named files
# throughout, so concurrent runs do not write to shared files.
# ==============================================================================
set -euo pipefail

THREADS=${THREADS:-4}
MAX_PARALLEL=${MAX_PARALLEL:-4}
DATA_DIR=data
RAW_DIR=$DATA_DIR/raw
REF_DIR=$DATA_DIR/reference
mkdir -p "$RAW_DIR" "$REF_DIR" config

# ------------------------------------------------------------------
# 1. Reference transcriptome (CDS transcripts, T. canis PRJNA248777)
# ------------------------------------------------------------------
REF_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/toxocara_canis/PRJNA248777/toxocara_canis.PRJNA248777.WBPS19.CDS_transcripts.fa.gz"
if [[ ! -f "$REF_DIR/tcanis_cds.fa" ]]; then
  echo ">> Downloading T. canis reference CDS transcriptome..."
  wget -q --show-progress -O "$REF_DIR/tcanis_cds.fa.gz" "$REF_URL" || {
    echo "!! Download failed. Verify the URL resolves first with:"
    echo "   wget --spider \"$REF_URL\""
    exit 1
  }
  gunzip -f "$REF_DIR/tcanis_cds.fa.gz"
else
  echo ">> Reference already present, skipping download."
fi

# ------------------------------------------------------------------
# 2. Sample sheet: confirmed SRR accessions from NCBI SRA
# ------------------------------------------------------------------
cat > config/samples.tsv << 'EOF'
sample_id	condition	replicate	srr
CONTROL_T1	control	1	SRR26866493
CONTROL_T16	control	2	SRR26866492
IVM_T3	ivermectin	1	SRR26866484
IVM_T4	ivermectin	2	SRR26866483
EOF

echo ">> Sample sheet (config/samples.tsv):"
cat config/samples.tsv
echo ">> Downloading full read depth for all 6 samples (up to $MAX_PARALLEL concurrently)..."

# ------------------------------------------------------------------
# 3. Fetch one sample's full read depth. Each sample gets its own
#    prefetch output directory (avoids any shared-file contention),
#    and each writes uniquely-named final FASTQs.
# ------------------------------------------------------------------
fetch_one() {
  local sample_id="$1" condition="$2" srr="$3"
  local log="$RAW_DIR/${sample_id}.log"
  if [[ -f "$RAW_DIR/${sample_id}_R1.fastq.gz" && -f "$RAW_DIR/${sample_id}_R2.fastq.gz" ]]; then
    echo "[$(date +%T)] SKIP  ${srr} (${sample_id}) -- already downloaded" > "$log"
    return
  fi
  {
    echo "[$(date +%T)] START ${srr} (${sample_id}, ${condition})"
    prefetch "$srr" -O "$RAW_DIR" < /dev/null
    fasterq-dump "$RAW_DIR/$srr" -O "$RAW_DIR" -e "$THREADS" --split-files < /dev/null
    gzip -f "$RAW_DIR/${srr}_1.fastq" "$RAW_DIR/${srr}_2.fastq"
    mv "$RAW_DIR/${srr}_1.fastq.gz" "$RAW_DIR/${sample_id}_R1.fastq.gz"
    mv "$RAW_DIR/${srr}_2.fastq.gz" "$RAW_DIR/${sample_id}_R2.fastq.gz"
    echo "[$(date +%T)] DONE  ${srr} (${sample_id}, ${condition})"
  } > "$log" 2>&1
}
export -f fetch_one
export RAW_DIR THREADS

mapfile -t FETCH_LINES < <(tail -n +2 config/samples.tsv)

# Launch samples in fixed-size batches: start up to MAX_PARALLEL jobs,
# wait for that whole batch to finish (plain `wait`, no job-counting),
# then start the next batch. Simpler and more reliable than a sliding
# window across bash job control, which behaves inconsistently in
# non-interactive scripts.
batch=()
launch_batch() {
  for entry in "${batch[@]}"; do
    IFS=$'\t' read -r sample_id condition replicate srr <<< "$entry"
    echo ">> Launching ${sample_id} in background (log: $RAW_DIR/${sample_id}.log)"
    fetch_one "$sample_id" "$condition" "$srr" &
  done
  wait
  batch=()
}

for line in "${FETCH_LINES[@]}"; do
  batch+=("$line")
  if [[ "${#batch[@]}" -ge "$MAX_PARALLEL" ]]; then
    launch_batch
  fi
done
# Launch any remaining samples that didn't fill a full batch.
if [[ "${#batch[@]}" -gt 0 ]]; then
  launch_batch
fi

echo ">> All fetches complete. Per-sample logs are in $RAW_DIR/*.log -- check them if"
echo "   anything looks off (e.g. grep -c DONE $RAW_DIR/*.log should each show 1)."
for line in "${FETCH_LINES[@]}"; do
  IFS=$'\t' read -r sample_id condition replicate srr <<< "$line"
  cat "$RAW_DIR/${sample_id}.log"
done

echo ">> Done. Full-depth paired FASTQs are in $RAW_DIR/, reference in $REF_DIR/."
echo ">> NOTE: these are the 'datafile(s) used' -- checksum them for"
echo "   CHECKSUMS.txt, but do NOT commit them to git (see .gitignore)."
