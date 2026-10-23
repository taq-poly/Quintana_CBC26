#!/bin/bash
# Inputs for the genome blobplot: DNA reads (SRR1706691), UniRef50 + NCBI taxdump, BlobToolKit env
set -uo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen
# DNA short-insert run (HiSeq 2x100, 170 bp insert)
( cd $B/reads/dna && prefetch --max-size 30G -O . SRR1706691 && \
  fasterq-dump --split-files -e 16 -t $B/reads/tmp -O . SRR1706691/SRR1706691.sra && \
  pigz -p 16 SRR1706691_1.fastq SRR1706691_2.fastq && rm -rf SRR1706691 && echo "dna SRR1706691 ok" ) &
# NCBI taxdump
( cd $B/db/taxdump && wget -q -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz && \
  tar xzf new_taxdump.tar.gz && echo "taxdump ok" ) &
# UniRef50
( cd $B/db/uniref50 && wget -q -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz && echo "uniref50 download ok" ) &
# BlobToolKit in its own env (pins differ from the main env); firefox+geckodriver for headless PNG/SVG export
mamba create -y -p $B/env_btk -c conda-forge -c bioconda python=3.11 firefox geckodriver pip > $B/logs/btk_env.log 2>&1 && \
  $B/env_btk/bin/pip install "blobtoolkit[full]" >> $B/logs/btk_env.log 2>&1 && echo "btk env ok"
wait; echo "BLOB INPUTS DONE"
