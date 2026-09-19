# Tissue-specific differential gene expression in the canine heartworm, *Dirofilaria immitis*

A compact RNA-seq DEG pipeline comparing **male testis vs. female uterus**
transcriptomes of the veterinary parasitic nematode *Dirofilaria immitis*
(dog heartworm), using public data. Reproductive-tissue-restricted genes
in filarial nematodes are of direct interest as anthelmintic and vaccine
targets, so this is a biologically meaningful (and statistically clean,
large-effect-size) demonstration of a standard RNA-seq quantification +
DEG workflow.

## What this run does

1. Downloads the *D. immitis* reference CDS transcriptome from WormBase
   ParaSite (genome assembly PRJEB1797).
2. Resolves and downloads 4 raw RNA-seq runs (2 testis + 2 uterus
   biological replicates) from NCBI SRA, linked to **GEO series
   [GSE67894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67894)**
   (Luck et al. 2015, *BMC Genomics* 16:920 — "Tissue-specific
   transcriptomics and proteomics of a filarial nematode and its
   *Wolbachia* endosymbiont", DOI: 10.1186/s12864-015-2083-2).
3. Subsamples each FASTQ to 2,000,000 reads (see "Why subsample" below).
4. Runs FastQC + MultiQC for QC.
5. Pseudo-aligns/quantifies with **Salmon** against the reference
   transcriptome (external program — not R or Python).
6. Imports counts with `tximport` and tests for differential expression
   with **DESeq2** in R, producing a results table, a filtered
   significant-hits table, and diagnostic plots (PCA, sample-distance
   heatmap, volcano plot).
7. Writes `CHECKSUMS.txt` (SHA-256 of every downloaded/generated file).

## Why subsample reads?

The original GAIIx libraries are single-end 50 bp and, while modest by
modern standards, downloading + quantifying all reads for 4 samples can
still take well over 30 minutes depending on connection speed. To keep
this a genuine "runs in ~30 minutes on a laptop" exercise, each FASTQ is
randomly subsampled (`seqtk sample -s100`, fixed seed for reproducibility)
to 2,000,000 reads before quantification. This is declared explicitly so
results are not mistaken for a full re-analysis of the original study —
it is a demonstration pipeline, not a publication-grade re-analysis.
Increase/remove the subsampling (`N_READS` env var in
`scripts/01_fetch_data.sh`) if you have more time or a faster connection.

## Required environment

- Conda / Mamba (Miniconda3 or Miniforge)
- Internet access to: `ftp.ebi.ac.uk` / `parasite.wormbase.org` (reference),
  NCBI `eutils`/SRA (`prefetch`, `fasterq-dump`)
- ~4 CPU cores and ~10 GB free disk recommended
- **No GPU required**

Create and activate the environment:

```bash
conda env create -f envs/environment.yml
conda activate dimmitis-deg
```

All tool versions are pinned in `envs/environment.yml` (sra-tools,
seqtk, FastQC, MultiQC, Salmon, samtools, Entrez Direct, R 4.3 +
DESeq2/tximport/apeglm/ggplot2/pheatmap).

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
| 1 | `scripts/01_fetch_data.sh` | Download reference transcriptome + resolve & download SRA reads + subsample |
| 2 | `scripts/02_qc.sh` | FastQC per sample + MultiQC summary |
| 3 | `scripts/03_quantify.sh` | Build Salmon index; quantify each sample |
| 4 | `scripts/04_deseq2.R` | tximport + DESeq2 DE test, tables, plots |

Each script can also be run individually (in order) for debugging.

### A note on step 1 (sample resolution)

`01_fetch_data.sh` uses NCBI Entrez Direct (`esearch`/`elink`/`efetch`/
`xtract`) to programmatically resolve GEO sample (GSM) accessions under
series GSE67894 to their linked SRA run (SRR) accessions, then
keyword-matches sample titles containing "testis" / "uterus" to build
`config/samples.tsv`. This avoids hard-coding SRR accession numbers that
could silently go stale; it also means **you should sanity-check
`data/gse67894_gsm_to_srr.tsv` and `config/samples.tsv` after step 1**
before trusting the tissue labels, and fill them in by hand from the
[GEO series page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67894)
/ [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) if the
automatic keyword match doesn't find exactly 2+2 runs.

## Methods (as you'd write it up)

**Data.** Public RNA-seq reads (single-end, 50 bp, Illumina GAIIx) for
adult *D. immitis* male testis (n=2 biological replicates) and female
uterus (n=2 biological replicates) were obtained from NCBI SRA via GEO
series GSE67894 (Luck et al., 2015). The *D. immitis* reference CDS
transcriptome (genome assembly PRJEB1797, WormBase ParaSite) was used as
the quantification reference.

**Read processing.** Each FASTQ file was randomly subsampled to 2×10⁶
reads (seqtk 1.4, seed 100) and quality-assessed with FastQC 0.12.1,
summarized with MultiQC 1.21.

**Quantification.** Reads were pseudo-aligned and quantified against the
reference CDS transcriptome using Salmon 1.10.3 (`salmon index -k 25`;
`salmon quant -l A --validateMappings --seqBias`), run in single-end mode.

**Differential expression.** Transcript-level counts were imported with
`tximport` (transcript-level, no gene collapsing, since *D. immitis*
annotation has minimal documented alternative splicing) and tested for
differential expression between testis and uterus with DESeq2 1.42.0
(Wald test, default independent filtering, apeglm shrinkage of log2 fold
changes). Transcripts with padj < 0.05 and |log2FC| > 1 were called
significant.

## Environment, installs, and data fetch (for anyone re-running this)

- Environment: `envs/environment.yml` (conda; see above)
- Reference transcriptome: fetched live from WormBase ParaSite (URL in
  `scripts/01_fetch_data.sh`; if the release path has moved, the script
  prints the manual-download fallback URL)
- Raw reads: fetched live from NCBI SRA (accessions resolved
  programmatically, see above); **not included in this repository**

## Outputs

After a successful run you should have:

```
results/
├── qc/multiqc_report.html
├── quant/<sample_id>/quant.sf         (Salmon per-sample quantifications)
├── deseq2/
│   ├── DE_testis_vs_uterus.csv                (full results table)
│   ├── DE_testis_vs_uterus_significant.csv    (padj<0.05, |log2FC|>1)
│   ├── pca.png
│   ├── sample_distance_heatmap.png
│   ├── volcano.png
│   └── sessionInfo.txt
CHECKSUMS.txt   (sha256 of every downloaded/generated file)
```

`results/deseq2/DE_testis_vs_uterus.csv` is the primary output table:
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
- No multi-mapping-to-Wolbachia filtering was performed (the CDS reference
  is nematode-only, so *w*Di reads simply fail to map and are excluded).

## License

MIT (see `LICENSE`).
