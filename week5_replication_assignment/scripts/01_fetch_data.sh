#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-4}
N_READS=${N_READS:-1000000}
DATA_DIR=data
RAW_DIR=$DATA_DIR/raw
SUB_DIR=$DATA_DIR/subsampled
REF_DIR=$DATA_DIR/reference
mkdir -p "$RAW_DIR" "$SUB_DIR" "$REF_DIR" config

REF_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/dirofilaria_immitis/PRJEB1797/dirofilaria_immitis.PRJEB1797.WBPS19.CDS_transcripts.fa.gz"
if [[ ! -f "$REF_DIR/dimmitis_cds.fa" ]]; then
  echo ">> Downloading D. immitis reference CDS transcriptome..."
  wget -q --show-progress -O "$REF_DIR/dimmitis_cds.fa.gz" "$REF_URL"
  gunzip -f "$REF_DIR/dimmitis_cds.fa.gz"
else
  echo ">> Reference already present, skipping download."
fi

# Verified directly via NCBI SRA (GEO series GSE67894, study SRP057178):
#   GSM1657649 (FU1)  -> SRX995516 -> SRR1974180
#   GSM1657650 (FU2)  -> SRX995517 -> SRR1974181
#   GSM1657643 (FBW1) -> SRX995510 -> SRR1974174
#   GSM1657644 (FBW2) -> SRX995511 -> SRR1974175
cat > config/samples.tsv << 'EOF'
sample_id	tissue	replicate	srr
FU1	uterus	1	SRR1974180
FU2	uterus	2	SRR1974181
FBW1	body_wall	1	SRR1974174
FBW2	body_wall	2	SRR1974175
EOF

echo ">> Sample sheet (config/samples.tsv):"
cat config/samples.tsv

mapfile -t FETCH_LINES < <(tail -n +2 config/samples.tsv)
for line in "${FETCH_LINES[@]}"; do
  IFS=$'\t' read -r sample_id tissue replicate srr <<< "$line"
  echo ">> Fetching ${srr} (${sample_id}, ${tissue})..."
  prefetch "$srr" -O "$RAW_DIR" < /dev/null
  fasterq-dump "$RAW_DIR/$srr" -O "$RAW_DIR" -e "$THREADS" < /dev/null
  seqtk sample -s100 "$RAW_DIR/${srr}.fastq" "$N_READS" > "$SUB_DIR/${sample_id}.fastq"
  gzip -f "$SUB_DIR/${sample_id}.fastq"
  rm -f "$RAW_DIR/${srr}.fastq"
done

echo ">> Done. Subsampled FASTQs are in $SUB_DIR/, reference in $REF_DIR/."
