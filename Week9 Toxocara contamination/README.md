# Week 9: Non-target sequence screen of *Toxocara canis* DNA and RNA-seq reads

## Overview

This project screens one DNA readset and one RNA-seq readset from the dog roundworm *Toxocara canis* for non-target sequence (contaminants, associates, viruses, database artefacts). Each detection is reported with read fraction, depth and breadth of coverage. It is placed in one of four bins (noise, trace, consequential, biological) and given a likely source with supporting evidence.

**Full results:** open [`report/toxocara_screen_report.html`](report/toxocara_screen_report.html) in a browser. It is a single self-contained page with all tables and the blobplot. To view it from GitHub, download the file and open it locally.

## Data

| Readset | BioProject | Runs | Notes |
|---|---|---|---|
| DNA (genome) | PRJNA248777 | SRR1706691 | HiSeq 2×100, 170 bp insert, adult worms (2014); assembly GCA_000803305.1 |
| RNA-seq | PRJNA1041894 | SRR26866479, SRR26866491, SRR26866490, SRR26866487, SRR26866489, SRR26866488 | 6 adult worms (3 female, 3 male), NextSeq 2×73, poly-A; one flowcell |

SRA run metadata is in `tables/*.runinfo.csv`.

## Main findings

**Genome readset**
- ***Clostridium perfringens***: 0.32% of read pairs, 10.6× depth, 80% breadth. Consequential: 283 contigs (256 kb, including a multicopy plasmid) are part of the published reference assembly. Source: real gut associate.
- ***E. coli***: 0.06% of pairs, 1.3× depth, 58% breadth. Trace. A wild strain (98.7% identity to K-12, against a 99.8% library baseline), not a lab strain.
- **Human and dog**: about 30 ppm each. Trace (handling and host carry-over).

**RNA-seq (6 adults)**
- **TcLA virus** (Jingmen-like, segmented; 2 new segments): 7–9% of read pairs in two males. Biological and consequential. Several *T. canis* "antigen" entries in UniProt are actually TcLA proteins (database error).
- **RV1–RV4-like totivirus-like viruses** and one new toti-like contig. Biological. Each worm carries its own strain. Low counts in other worms are bleed-through between samples on the shared flowcell, shown by shared variants.
- **Reagent bacteria** (*Herbaspirillum*, *Delftia*, *Leifsonia*, *Ralstonia*), **human skin** signal, **gut bacteria** and **dog**: all trace (≤ 200 ppm).
- **Classifier false positives** (diatom, cloning vector pRNA8, *Salmonella*, *Goezia*): noise from database errors.

## Methods (scripts in run order)

| Step | Script | Tools |
|---|---|---|
| Download RNA reads | `00_download_adults.sh`, `run_adults.sh` | sra-tools |
| References and databases | `01_get_refs_dbs.sh` | NCBI genomes, Kraken2 PlusPF-16, sourmash GTDB/NCBI, geNomad, CheckV |
| Reference QC (contaminant contigs in the assembly) | `04_refqc_minimap.sh`, `05_clean_reference.sh` | minimap2, Kraken2 |
| Per-sample screen | `10_screen_sample.sh`, `11_/12_sourmash_*` | fastp, HISAT2, Kraken2 + Bracken, sendsketch, sourmash fastgather |
| Genome blobplot | `20_get_blob_inputs.sh`, `21_genome_blob.sh`, `22_build_uniref50_db.sh` | minimap2, DIAMOND vs UniRef50, BlobToolKit |
| Summaries | `30_summarize_pass1.py` | |
| RNA virus discovery | `50_adult_virus_discovery.sh`, `40_virus_haplotypes.py` | MEGAHIT, geNomad, CheckV, bowtie2 |
| Depth/breadth by competitive mapping | `60_build_mm2_targets.sh`, `61_competitive_map.sh`, `62_assign_competitive.py`, `63_summarize_competitive.py` | minimap2, samtools coverage |

Scripts use absolute paths under `/workspace/toxocara_screen`; replace that prefix to rerun elsewhere (e.g. `sed -i 's|/workspace/toxocara_screen|/your/path|g' scripts/*`). Software versions are in `environment/`.

## Repository contents

```
report/        toxocara_screen_report.html (full write-up)
figures/       BlobToolKit blob, snail and cumulative plots of the T. canis assembly
tables/
  genome/        competitive-mapping coverage (samtools coverage) and bacterial contigs, SRR1706691
  rna/           per-adult Kraken2/Bracken, sendsketch, competitive-mapping summaries
  reference_qc/  contaminant-flagged contigs in GCA_000803305.1
  viruses/       geNomad viral contigs, CheckV quality, per-sample coverage, matches to known T. canis viruses
scripts/       pipeline
environment/   conda environment exports
```

Raw reads, BAMs and databases (~155 GB) are not included. They can be regenerated with the scripts.

## Scope and limitations

- Only the 6 adult RNA-seq samples of PRJNA1041894 were screened; the 9 larval samples were not.
- The optional spike-in detection-limit test was not run.
