#!/usr/bin/env Rscript
# ==============================================================================
# 04_deseq2.R
#
# Imports Salmon quant.sf files (transcript-level; txOut = TRUE, no gene
# collapsing) and runs a DESeq2 differential expression test: ivermectin
# vs control, for Toxocara canis third-stage larvae (n=3 per group).
# ==============================================================================
suppressPackageStartupMessages({
  library(tximport)
  library(jsonlite)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(readr)
})

dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)

samples <- read_tsv("config/samples.tsv", show_col_types = FALSE)
samples$condition <- factor(samples$condition, levels = c("control", "ivermectin"))

files <- file.path("results/quant", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", txOut = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # basic low-count filter
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "ivermectin", "control"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "condition_ivermectin_vs_control", type = "apeglm")

res_df <- as.data.frame(res_shrunk)
res_df$transcript_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
res_df <- res_df[, c("transcript_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write_csv(res_df, "results/deseq2/DE_ivermectin_vs_control.csv")

sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write_csv(sig, "results/deseq2/DE_ivermectin_vs_control_significant.csv")

cat(sprintf(
  "Total transcripts tested: %d\nSignificant (padj<0.05, |log2FC|>1): %d\n",
  nrow(res_df), nrow(sig)
))

# ---- Plots ----
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA: T. canis L3 ivermectin vs control (VST counts)")
ggsave("results/deseq2/pca.png", p_pca, width = 6, height = 5, dpi = 150)

sample_dist <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)
png("results/deseq2/sample_distance_heatmap.png", width = 900, height = 800, res = 150)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

res_df$sig <- with(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  theme_bw() +
  labs(title = "Ivermectin vs control: volcano plot",
       color = "padj<0.05 & |log2FC|>1")
ggsave("results/deseq2/volcano.png", p_volcano, width = 6, height = 5, dpi = 150)

writeLines(capture.output(sessionInfo()), "results/deseq2/sessionInfo.txt")

cat(">> DESeq2 outputs written to results/deseq2/\n")
