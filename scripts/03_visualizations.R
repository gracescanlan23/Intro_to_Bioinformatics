if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(here)
library(AnnotationDbi)
library(org.Hs.eg.db)

source(here("scripts", "02_normalize_counts.R"))

output_dir <- here("output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------
# Differential expression
# -------------------------
dds <- DESeq(dds)

# Make sure the contrast matches reference level setup
res <- results(dds, contrast = c("tumor_type", "Advanced_Tumor", "Normal"))

# -------------------------
# Clean result table
# -------------------------
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# Add gene symbols and names
res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene_id,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# -------------------------
# DE summary counts
# -------------------------
genes_tested <- nrow(res_df)
sig_count <- sum(res_df$padj < 0.05, na.rm = TRUE)
up_count <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
down_count <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)

cat("Genes tested:", genes_tested, "\n")
cat("Significant genes (FDR < 0.05):", sig_count, "\n")
cat("Upregulated genes:", up_count, "\n")
cat("Downregulated genes:", down_count, "\n")

# -------------------------
# Regulation groups
# -------------------------
res_df$Regulation <- "NS"
res_df$Regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$Regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# -------------------------
# Top 10 genes for labeling
# -------------------------
top_genes <- res_df[order(res_df$padj), ][1:10, ]
top_genes$label <- ifelse(
  is.na(top_genes$gene_symbol) | top_genes$gene_symbol == "",
  top_genes$gene_id,
  top_genes$gene_symbol
)

# -------------------------
# Export DEG results
# -------------------------
deg_results <- res_df[, c("gene_id", "log2FoldChange", "pvalue", "padj", "Regulation", "gene_symbol", "gene_name")]
colnames(deg_results) <- c(
  "Gene_ID", "log2FC", "p_value", "adjusted_p_value",
  "Regulation", "Gene_Symbol", "Gene_Name"
)

write.csv(deg_results, here("DEG_results.csv"), row.names = FALSE)

# --------------------------
# number of significant DEGs (FDR < 0.05) and directionality 
# --------------------------

num_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
up_genes <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
down_genes <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)


# -------------------------
# Volcano plot
# -------------------------
volcano_plot <- ggplot(
  res_df,
  aes(x = log2FoldChange, y = -log10(padj), color = Regulation)
) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = label),
    size = 3,
    max.overlaps = 100
  ) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "gray", "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate(
  "text",
  x = max(res_df$log2FoldChange, na.rm = TRUE) * 0.35,
  y = max(-log10(res_df$padj), na.rm = TRUE) * 0.95,
  label = paste0(
    "Significant DEGs: ", num_sig, "\n",
    "Upregulated: ", up_genes, "\n",
    "Downregulated: ", down_genes
  ),
  hjust = 0,
  size = 4
)+
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
ggsave(
  filename = "output/volcano_plot.pdf",
  plot = volcano_plot,
  width = 8,
  height = 6
)
# -------------------------
# MA plot
# -------------------------
ma_plot <- ggplot(
  res_df,
  aes(x = baseMean, y = log2FoldChange, color = Regulation)
) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = label),
    size = 3,
    max.overlaps = 100
  ) +
  scale_x_log10() +
  scale_color_manual(values = c("Down" = "blue", "NS" = "gray", "Up" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean Expression (log scale)",
    y = "Log2 Fold Change"
  )

ggsave(
  filename = here("output", "ma_plot.pdf"),
  plot = ma_plot,
  width = 8,
  height = 6
)

# -------------------------
# PCA plot
# -------------------------
vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vsd, intgroup = "tumor_type") +
  ggtitle("PCA Plot") +
  theme_minimal()

ggsave(
  filename = here("output", "pca_plot.pdf"),
  plot = pca_plot,
  width = 8,
  height = 6
)
# -------------------------
# PCA subplot without one normal sample
# -------------------------
# -------------------------
# Subset PCA plot (remove one normal sample)
# -------------------------
dds_subset <- dds[, colnames(dds) != "SRR36750771"]

vsd_subset <- vst(dds_subset, blind = FALSE)

pca_subset_plot <- plotPCA(vsd_subset, intgroup = "tumor_type") +
  ggtitle("PCA Plot (Subset: One Normal Sample Removed)") +
  theme_minimal()

ggsave(
  filename = here("output", "pca_subset.pdf"),
  plot = pca_subset_plot,
  width = 8,
  height = 6
)

# -------------------------
# Top 10 DEG table for report
# -------------------------
top10_deg <- head(deg_results[order(deg_results$adjusted_p_value), ], 10)
print(top10_deg)

write.csv(
  top10_deg,
  file = here("output", "top10_DEGs.csv"),
  row.names = FALSE
)
