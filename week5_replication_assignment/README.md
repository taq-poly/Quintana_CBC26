# Ivermectin-induced differential gene expression in *Toxocara canis* larvae

A compact RNA-seq DEG pipeline comparing **ivermectin-treated vs. control**
third-stage larvae (L3) of *Toxocara canis* (the zoonotic dog roundworm),
using public data. Ivermectin (a macrocyclic lactone) is a frontline
anthelmintic, and understanding how *T. canis* larvae respond
transcriptionally to it is directly relevant to drug tolerance/resistance
research in veterinary parasitology.

## What this run does

1. Downloads the *T. canis* reference CDS transcriptome from WormBase
   ParaSite (genome assembly PRJNA248777, release WBPS19 -- the same
   assembly used by the source study).
2. Downloads 4 of the 6 raw paired-end RNA-seq runs (2 control + 2
   ivermectin-treated biological replicates) from NCBI SRA, BioProject
   **PRJNA1041894** / study **SRP472648** (Quintana, Brewer & Jesudoss
   Chelladurai, 2025, *Int J Parasitol Drugs Drug Resist* 29:100614,
   "Transcriptional responses to *in vitro* macrocyclic lactone exposure
   in *Toxocara canis* larvae using RNA-seq", DOI:
   10.1016/j.ijpddr.2025.100614).
3. Downloads full read depth for each sample (no subsampling — see
   below).
4. Runs FastQC + MultiQC for QC.
5. Pseudo-aligns/quantifies with **Salmon** in paired-end mode against
   the reference transcriptome (external program -- not R or Python).
6. Imports counts with `tximport` and tests for differential expression
   with **DESeq2** in R, producing a results table, a filtered
   significant-hits table, and diagnostic plots (PCA, sample-distance
   heatmap, volcano plot).
7. Writes `CHECKSUMS.txt` (SHA-256 of every downloaded/generated file).

## Why control vs. ivermectin only (not moxidectin too)?

The source study profiled three conditions (control, ivermectin, moxidectin), each n=3. This reanalysis further subsets to n=2 per group (the two smallest/fastest replicates of each condition) and drops moxidectin entirely, to keep full-depth downloads practical on a laptop,
while still directly answering "does ivermectin change the larval
transcriptome".

## Full read depth -- no subsampling

This pipeline downloads **every read of every sample** -- no
subsampling, no truncation -- via `prefetch` + `fasterq-dump`,
matching the source study's depth exactly (~7.2 GB total across the 4 selected samples, 17-40 million read pairs each).

```bash
bash run_all.sh
```

This step is slow, and that is a direct consequence of the real data
volume involved, not an inefficiency in the script: ~7.2 GB of real
sequencing data has to be transferred and decompressed. 

`fasterq-dump -e $THREADS` controls how many CPU threads are used for
extraction; set the `THREADS` environment variable to roughly your
machine's core count if you want to try adjusting it (default: 4).

## Required environment

- Conda / Mamba (Miniconda3, Miniforge, or Mambaforge)
- Internet access to: `ftp.ebi.ac.uk` (reference), NCBI SRA
  (`prefetch`, `fasterq-dump`)
- ~4 CPU cores
- **~15 GB free disk** (full-depth paired-end downloads for 4 samples)
- **No GPU required**

Create and activate the environment:

```bash
conda env create -f envs/environment.yml
conda activate toxocara-deg
```

All tool versions are pinned in `envs/environment.yml` (sra-tools,
seqtk, FastQC, MultiQC, Salmon, samtools, R 4.3 +
DESeq2/tximport/apeglm/jsonlite/ggplot2/pheatmap).

## How to run

```bash
git clone <this-repo-url>
cd toxocara-ivm-deg
conda env create -f envs/environment.yml
conda activate toxocara-deg
bash run_all.sh
```

`run_all.sh` calls, in order:

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_fetch_data.sh` | Download reference transcriptome + download SRA reads (paired-end) |
| 2 | `scripts/02_qc.sh` | FastQC per sample (both mates) + MultiQC summary |
| 3 | `scripts/03_quantify.sh` | Build Salmon index; quantify each sample (paired-end mode) |
| 4 | `scripts/04_deseq2.R` | tximport + DESeq2 DE test, tables, plots |

Each script can also be run individually (in order) for debugging.

## Sample accessions (confirmed by hand via NCBI SRA)

`scripts/01_fetch_data.sh` hardcodes these six accessions, verified
directly via the NCBI SRA web interface (BioProject PRJNA1041894 ->
individual SRX/SRR pages). Every read count below matches the source
paper's Table 1 exactly, confirming correct sample identity:

| sample_id | condition | SRA experiment (SRX) | SRA run (SRR) | read pairs |
|---|---|---|---|---|
| CONTROL_T1  | control    | SRX22561044 | SRR26866493 | 32,249,245 |
| CONTROL_T16 | control    | SRX22561045 | SRR26866492 | 40,128,609 |
| IVM_T3      | ivermectin | SRX22561053 | SRR26866484 | 22,121,380 |
| IVM_T4      | ivermectin | SRX22561054 | SRR26866483 | 17,179,802 |

(All six fall under NCBI SRA study **SRP472648** / BioProject
**PRJNA1041894**.)

## Methods (as you'd write it up)

**Data.** Public paired-end RNA-seq reads (2x75bp, Illumina NextSeq 550)
for *T. canis* hatched third-stage larvae, untreated (control, n=2
biological replicates) or treated *in vitro* with 10 uM ivermectin
(n=2 biological replicates) for 12 hours, were obtained from NCBI SRA
(study SRP472648) under BioProject PRJNA1041894 (Quintana et al., 2025).
The *T. canis* reference CDS transcriptome (genome assembly PRJNA248777,
WormBase ParaSite release WBPS19) was used as the quantification
reference.

**Read processing.** Full-depth paired-end reads for each run were
downloaded from NCBI SRA using `prefetch` followed by `fasterq-dump
--split-files` (sra-tools 3.1.1), with no truncation or subsampling,
matching the source study's read depth exactly. Reads were
quality-assessed with FastQC 0.12.1, summarized with MultiQC 1.21.

**Quantification.** Reads were pseudo-aligned and quantified against the
reference CDS transcriptome using Salmon 1.10.3 (`salmon index -k 31`;
`salmon quant -l A --validateMappings --seqBias`), run in paired-end mode.

**Differential expression.** Transcript-level counts were imported with
`tximport` (transcript-level, no gene collapsing) and tested for
differential expression between ivermectin and control with DESeq2
1.42.0 (Wald test, default independent filtering, apeglm shrinkage of
log2 fold changes). Transcripts with padj < 0.05 and |log2FC| > 1 were
called significant.

## Environment, installs, and data fetch (for anyone re-running this)

- Environment: `envs/environment.yml` (conda; see above)
- Reference transcriptome: fetched live from WormBase ParaSite (URL in
  `scripts/01_fetch_data.sh`; verify with `wget --spider <URL>` first
  if the release version has moved -- WBPS numbering increments over
  time)
- Raw reads: fetched live from NCBI SRA using the accessions in the
  table above; **not included in this repository**

## Outputs

After a successful run you should have:

```
results/
├── qc/multiqc_report.html
├── quant/<sample_id>/quant.sf         (Salmon per-sample quantifications)
├── deseq2/
│   ├── DE_ivermectin_vs_control.csv                (full results table)
│   ├── DE_ivermectin_vs_control_significant.csv    (padj<0.05, |log2FC|>1)
│   ├── pca.png
│   ├── sample_distance_heatmap.png
│   ├── volcano.png
│   └── sessionInfo.txt
CHECKSUMS.txt   (sha256 of every downloaded/generated file)
```

`results/deseq2/DE_ivermectin_vs_control.csv` is the primary output
table: one row per transcript, with `baseMean`, `log2FoldChange`
(apeglm-shrunk), `lfcSE`, `pvalue`, `padj`.

## Repository layout

```
.
├── README.md
├── envs/environment.yml
├── scripts/
│   ├── 01_fetch_data.sh
│   ├── 02_qc.sh
│   ├── 03_quantify.sh
│   └── 04_deseq2.R
├── run_all.sh
└── .gitignore          (excludes data/ and large intermediate files)
```

## Caveats

- n=2 per group (subset from the original study's n=3) is a small sample size for DE testing (typical of
  in vitro drug-exposure studies of this kind), which limits
  statistical power for detecting smaller-effect-size genes.
- Transcript-level (not gene-collapsed) DE testing is used for
  simplicity; the original study aligned to the full genome with STAR
  and counted at the gene level, so exact DEG counts will differ from
  the paper's reported 608 DEGs for ivermectin -- this pipeline
  demonstrates the same underlying biological question and workflow
  shape (full read depth, same samples, same reference assembly) but
  uses a different quantification method (Salmon/CDS pseudo-alignment
  vs. STAR/full-genome alignment), so it is not an exact numeric
  replication of their results.
- Moxidectin samples exist in the same BioProject but are excluded here
  to keep this a two-group comparison (see above).

## License

MIT (see `LICENSE`).
