library(DESeq2)
library(pheatmap)
library(matrixStats)
library(here)

# Load normalized data
vsd <- readRDS(here("normalized_counts", "vst_object.rds"))
metadata <- read.table(here("config", "metadata.tsv"), header = TRUE, sep = "\t")

mat <- assay(vsd)

# Find top 50 most variable genes
rv <- rowVars(mat)
topVarGenes <- order(rv, decreasing = TRUE)[1:50]

# Sample Distance Heatmap
sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- colnames(mat)
colnames(sampleDistMatrix) <- colnames(mat)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# Hierarchical Clustering Dendrogram
hc <- hclust(sampleDists)
plot(hc, main = "Hierarchical Clustering of Samples")

# Heatmap of Top Variable Genes
rownames(metadata) <- metadata$sample_id
anno <- metadata[colnames(mat), c("condition", "tumor_type"), drop = FALSE]
stopifnot(identical(rownames(anno), colnames(mat)))

pheatmap(mat[topVarGenes, ],
         scale = "row",
         annotation_col = anno,
         main = "Top 50 Variable Genes Heatmap")

# Save heatmaps to file
pdf("normalized_counts/sample_distance_heatmap.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")
dev.off()


library(pheatmap)

# 1) keep only significant genes
sig_only <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

# 2) sort by adjusted p-value (smallest first)
sig_only <- sig_only[order(sig_only$padj), ]

# 3) choose top genes for the heatmap
top_n <- 20
top_genes <- rownames(sig_only)[1:min(top_n, nrow(sig_only))]

# 4) extract vst values for those genes
mat <- assay(vsd)[top_genes, ]

# 5) scale each gene across samples
mat_scaled <- t(scale(t(mat)))

# 6) optional: use gene symbols as row names if available
row_labels <- sig_only$gene_symbol[match(top_genes, rownames(sig_only))]
row_labels[is.na(row_labels) | row_labels == ""] <- top_genes[is.na(row_labels) | row_labels == ""]
rownames(mat_scaled) <- row_labels

# 7) make sure metadata rownames match sample names
rownames(metadata) <- metadata$sample_id
annotation_col <- metadata[colnames(mat_scaled), c("condition", "tumor_type"), drop = FALSE]

# 8) draw heatmap
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  main = "Top Significant Genes Heatmap",
  scale = "none"
)

pdf("normalized_counts/top_significant_genes_heatmap.pdf", width = 8, height = 10)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  main = "Top Significant Genes Heatmap",
  scale = "none"
)

dev.off()