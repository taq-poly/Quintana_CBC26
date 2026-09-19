# Tissue-specific differential gene expression in the canine heartworm, *Dirofilaria immitis*

A compact RNA-seq DEG pipeline comparing **female uterus vs. female body
wall** transcriptomes of the veterinary parasitic nematode *Dirofilaria
immitis* (dog heartworm), using public data. Comparing reproductive
tissue (uterus) to somatic tissue (body wall) within the same sex
isolates the tissue effect cleanly and is expected to show large,
easily detectable differences (oogenesis/reproduction-associated genes
vs. structural/cuticle genes) -- a good demonstration of a standard
RNA-seq quantification + DEG workflow at teaching scale.

## What this run does

1. Downloads the *D. immitis* reference CDS transcriptome from WormBase
   ParaSite (genome assembly PRJEB1797, release WBPS19).
2. Downloads 4 raw RNA-seq runs (2 uterus + 2 body-wall biological
   replicates) from NCBI SRA, linked to **GEO series
   [GSE67894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67894)**
   (Luck et al. 2015, *BMC Genomics* 16:920 -- "Tissue-specific
   transcriptomics and proteomics of a filarial nematode and its
   *Wolbachia* endosymbiont", DOI: 10.1186/s12864-015-2083-2).
3. Subsamples each FASTQ to 1,000,000 reads (see "Why subsample" below).
4. Runs FastQC + MultiQC for QC.
5. Pseudo-aligns/quantifies with **Salmon** against the reference
   transcriptome (external program -- not R or Python).
6. Imports counts with `tximport` and tests for differential expression
   with **DESeq2** in R, producing a results table, a filtered
   significant-hits table, and diagnostic plots (PCA, sample-distance
   heatmap, volcano plot).
7. Writes `CHECKSUMS.txt` (SHA-256 of every downloaded/generated file).

## Why uterus vs. body wall (not testis vs. uterus)?

The source study (GSE67894) profiled 9 adult *D. immitis* tissues, most
with 2 biological replicates -- **except male testis (MT) and male
intestine (MI), which have only 1 replicate each**. A DE test needs
within-group replication to estimate variance, so testis vs. uterus is
not statistically supportable from this dataset. Female uterus (FU,
n=2) vs. female body wall (FBW, n=2) is the cleanest same-sex,
fully-replicated tissue contrast available, and also controls for sex
as a confound.

## Why subsample reads?

The original GAIIx libraries are single-end 50 bp and range from
~1.2M to ~2.3M total reads per sample (small by modern standards, but
downloading + quantifying all of them for 4 samples can still take
longer than 30 minutes depending on connection speed). Each FASTQ is
randomly subsampled (`seqtk sample -s100`, fixed seed for
reproducibility) to 1,000,000 reads before quantification -- chosen to
be below the smallest sample's total read count (FBW2 has ~1.2M reads).
Increase/remove the subsampling (`N_READS` env var in
`scripts/01_fetch_data.sh`) if you want the full read depth.

## Required environment

- Conda / Mamba (Miniconda3, Miniforge, or Mambaforge)
- Internet access to: `ftp.ebi.ac.uk` (reference), NCBI SRA
  (`prefetch`, `fasterq-dump`)
- ~4 CPU cores and ~10 GB free disk recommended
- **No GPU required**
- Tested on WSL2 (Ubuntu) on Windows; should work identically on
  native Linux or macOS

Create and activate the environment:

```bash
conda env create -f envs/environment.yml
conda activate dimmitis-deg
```

If `conda activate` doesn't work in your shell, either run
`conda init bash` once (then restart your terminal), or use
`mamba activate` after running `mamba shell init --shell bash` (note:
this may point to a *different* root prefix than your existing conda
install -- if so, just use `conda activate dimmitis-deg` instead, since
the environment lives in your normal conda installation regardless of
whether `conda` or `mamba` was used to create it).

All tool versions are pinned in `envs/environment.yml` (sra-tools,
seqtk, FastQC, MultiQC, Salmon, samtools, Entrez Direct, R 4.3 +
DESeq2/tximport/apeglm/jsonlite/ggplot2/pheatmap).

## How to run

```bash
git clone <this-repo-url>
cd dimmitis-tissue-deg
conda env create -f envs/environment.yml
conda activate dimmitis-deg
bash run_all.sh
```

`run_all.sh` calls, in order:

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_fetch_data.sh` | Download reference transcriptome + download SRA reads + subsample |
| 2 | `scripts/02_qc.sh` | FastQC per sample + MultiQC summary |
| 3 | `scripts/03_quantify.sh` | Build Salmon index; quantify each sample |
| 4 | `scripts/04_deseq2.R` | tximport + DESeq2 DE test, tables, plots |

Each script can also be run individually (in order) for debugging.

## Sample accessions (confirmed by hand via NCBI SRA)

`scripts/01_fetch_data.sh` uses these four accessions directly, rather
than resolving them live via NCBI's Entrez Direct API. In testing,
Entrez Direct's `esearch`/`elink`/`efetch` pipeline was unreliable over
some network setups (notably WSL2) -- it either returned stale/cached
results, or its internal `curl` calls threw intermittent
`SSL_ERROR_SYSCALL` errors. Hardcoding accessions that were manually
verified via the GEO/SRA web interface (each GSM page -> "Relations: SRA"
link -> run accession) is both simpler and more reliable here, given
there are only 4 samples to track down:

| sample_id | tissue | GEO (GSM) | SRA experiment (SRX) | SRA run (SRR) |
|---|---|---|---|---|
| FU1  | uterus    | GSM1657649 | SRX995516 | SRR1974180 |
| FU2  | uterus    | GSM1657650 | SRX995517 | SRR1974181 |
| FBW1 | body_wall | GSM1657643 | SRX995510 | SRR1974174 |
| FBW2 | body_wall | GSM1657644 | SRX995511 | SRR1974175 |

(All four fall under NCBI SRA study **SRP057178** / BioProject
**PRJNA281132**.)

## Methods (as you'd write it up)

**Data.** Public RNA-seq reads (single-end, 50 bp, Illumina Genome
Analyzer II) for adult *D. immitis* female uterus (n=2 biological
replicates) and female body wall (n=2 biological replicates) were
obtained from NCBI SRA (study SRP057178) via GEO series GSE67894
(Luck et al., 2015). The *D. immitis* reference CDS transcriptome
(genome assembly PRJEB1797, WormBase ParaSite release WBPS19) was used
as the quantification reference.

**Read processing.** Each FASTQ file was randomly subsampled to 1x10^6
reads (seqtk 1.4, seed 100) and quality-assessed with FastQC 0.12.1,
summarized with MultiQC 1.21.

**Quantification.** Reads were pseudo-aligned and quantified against the
reference CDS transcriptome using Salmon 1.10.3 (`salmon index -k 25`;
`salmon quant -l A --validateMappings --seqBias`), run in single-end mode.

**Differential expression.** Transcript-level counts were imported with
`tximport` (transcript-level, no gene collapsing, since *D. immitis*
annotation has minimal documented alternative splicing) and tested for
differential expression between uterus and body wall with DESeq2
1.42.0 (Wald test, default independent filtering, apeglm shrinkage of
log2 fold changes). Transcripts with padj < 0.05 and |log2FC| > 1 were
called significant.

## Environment, installs, and data fetch (for anyone re-running this)

- Environment: `envs/environment.yml` (conda; see above)
- Reference transcriptome: fetched live from WormBase ParaSite (URL in
  `scripts/01_fetch_data.sh`; if the release path has moved, the script
  prints the manual-download fallback URL)
- Raw reads: fetched live from NCBI SRA using the accessions in the
  table above; **not included in this repository**

## Outputs

After a successful run you should have:

```
results/
├── qc/multiqc_report.html
├── quant/<sample_id>/quant.sf         (Salmon per-sample quantifications)
├── deseq2/
│   ├── DE_uterus_vs_bodywall.csv                (full results table)
│   ├── DE_uterus_vs_bodywall_significant.csv    (padj<0.05, |log2FC|>1)
│   ├── pca.png
│   ├── sample_distance_heatmap.png
│   ├── volcano.png
│   └── sessionInfo.txt
CHECKSUMS.txt   (sha256 of every downloaded/generated file)
```

`results/deseq2/DE_uterus_vs_bodywall.csv` is the primary output table:
one row per transcript, with `baseMean`, `log2FoldChange` (apeglm-shrunk),
`lfcSE`, `pvalue`, `padj`.

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

- n=2 per group and read subsampling mean this is a **teaching/demo-scale**
  analysis, not a publication-grade re-analysis of Luck et al. 2015.
- Transcript-level (not gene-collapsed) DE testing is used for simplicity;
  a production analysis would map transcripts to gene models via the
  WormBase ParaSite GFF3 and collapse with `tximport`'s `tx2gene`.
- No multi-mapping-to-*Wolbachia* filtering was performed (the CDS
  reference is nematode-only, so *w*Di reads simply fail to map and
  are excluded). Mapping rates against the CDS-only reference were
  modest (~20-25%), consistent with older, single-end, whole-tissue
  libraries mapped against a coding-sequence-only (not full-genome)
  reference.

## License

MIT (see `LICENSE`).
