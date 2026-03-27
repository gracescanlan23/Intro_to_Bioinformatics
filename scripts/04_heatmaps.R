if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

library(DESeq2)
library(pheatmap)
library(matrixStats)
library(here)

# ---------------------------------
# Load analysis objects
# ---------------------------------
source(here("scripts", "03_visualizations.R"))

# Make sure output folder exists
output_dir <- here("output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---------------------------------
# 1. Sample distance heatmap
# ---------------------------------
mat <- assay(vsd)

sample_dists <- dist(t(mat))
sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- colnames(mat)
colnames(sample_dist_matrix) <- colnames(mat)

pdf(here("output", "sample_distance_heatmap.pdf"), width = 8, height = 8)

pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample Distance Heatmap"
)

dev.off()

# ---------------------------------
# 2. Top significant genes heatmap
# ---------------------------------

# keep only significant genes
sig_only <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

# stop if there are no significant genes
if (nrow(sig_only) == 0) {
  stop("No significant genes found with padj < 0.05.")
}

# sort by adjusted p-value
sig_only <- sig_only[order(sig_only$padj), ]

# choose top genes
top_n <- 20
top_genes <- sig_only$gene_id[1:min(top_n, nrow(sig_only))]

# extract vst matrix values for those genes
mat_top <- assay(vsd)[top_genes, , drop = FALSE]

# scale by gene (row)
mat_scaled <- t(scale(t(mat_top)))

# remove rows with non-finite values after scaling
keep_rows <- rowSums(is.na(mat_scaled) | is.infinite(mat_scaled)) == 0
mat_scaled <- mat_scaled[keep_rows, , drop = FALSE]

# build row labels using gene symbols when available
gene_lookup <- res_df
rownames(gene_lookup) <- gene_lookup$gene_id

row_labels <- gene_lookup[rownames(mat_scaled), "gene_symbol"]
row_labels[is.na(row_labels) | row_labels == ""] <- rownames(mat_scaled)

rownames(mat_scaled) <- row_labels

# sample annotation from DESeq object
annotation_col <- as.data.frame(colData(vsd))[, "tumor_type", drop = FALSE]

pdf(here("output", "top_significant_genes_heatmap.pdf"), width = 8, height = 10)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  main = "Top Significant Genes Heatmap",
  scale = "none"
)

dev.off()