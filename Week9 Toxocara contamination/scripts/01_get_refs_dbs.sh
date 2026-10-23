#!/bin/bash
# Reference genomes + screening databases
set -uo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all
get_asm () { # accession_name dir_path outname
  wget -q -c -O $B/refs/$3.fna.gz $NCBI/$2/$1/$1_genomic.fna.gz && echo "ref $3 ok"; }
get_asm GCA_000803305.1_Toxocara_canis_adult_r1.0 GCA/000/803/305 tcanis
get_asm GCF_000819615.1_ViralProj14015            GCF/000/819/615 phix
get_asm GCF_000005845.2_ASM584v2                  GCF/000/005/845 ecoli_k12
get_asm GCF_011100685.1_UU_Cfam_GSD_1.0           GCF/011/100/685 dog
get_asm GCF_009914755.1_T2T-CHM13v2.0             GCF/009/914/755 human
# Kraken2 PlusPF-16
mkdir -p $B/db/k2_pluspf16 && cd $B/db/k2_pluspf16
wget -q -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20260626.tar.gz && tar xzf k2_pluspf_16_GB_20260626.tar.gz && rm k2_pluspf_16_GB_20260626.tar.gz && echo "kraken ok"
# sourmash: GTDB reps k31, NCBI viruses k31, NCBI eukaryotes k51
mkdir -p $B/db/sourmash && cd $B/db/sourmash
S=https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new
for f in gtdb-rs226/gtdb-reps-rs226-k31.dna.zip gtdb-rs226/gtdb-rs226.lineages.csv \
         ncbi-viruses-2025.01/ncbi-viruses-2025.01.dna.k=31.sig.zip ncbi-viruses-2025.01/ncbi-viruses.2025.01.lineages.csv \
         ncbi-euks-2025.01/ncbi-eukaryotes.2025.01.lineages.csv \
         ncbi-euks-2025.01/ncbi-euks-vertebrates-2025.01.dna.k=51.sig.zip \
         ncbi-euks-2025.01/ncbi-euks-bilateria-minus-vertebrates-2025.01.dna.k=51.sig.zip \
         ncbi-euks-2025.01/ncbi-euks-fungi-2025.01.dna.k=51.sig.zip \
         ncbi-euks-2025.01/ncbi-euks-metazoa-minus-bilateria-2025.01.dna.k=51.sig.zip \
         ncbi-euks-2025.01/ncbi-euks-other-2025.01.dna.k=51.sig.zip \
         ncbi-euks-2025.01/ncbi-euks-plants-2025.01.dna.k=51.sig.zip; do
  wget -q -c "$S/$f" && echo "sourmash $(basename $f) ok"; done
# geNomad + CheckV
cd $B/db && genomad download-database . > $B/logs/genomad_db.log 2>&1 && echo "genomad ok"
checkv download_database $B/db/ > $B/logs/checkv_db.log 2>&1 && echo "checkv ok"
echo ALL DONE
