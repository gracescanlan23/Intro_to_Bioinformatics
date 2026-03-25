library(DESeq2)
library(ggrepel)

source("scripts/normalize_counts.R")

dds <- DESeq(dds)
res <- results(dds)

# volcano plot

library(ggplot2)

# convert results to dataframe
# Explanation: res is a special object -> we convert it to a normal dataframe
# We remove rows with missing values for: adjusted p-values (padj) and log2 fold change. 
# Why: Some genes don't get valid statistics -> they would break plots or give weird points.
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# define significance
# What this does: Classifies each gene into: Upregulated (red), Downregulated (blue), or Not Significant (gray) based on adjusted p-value and fold change.
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# top 10 genes for labeling
# What this does: Sorts genes by smallest adjusted p-value (most significant first), selects top 10, and saves gene names from rownames into a column
#Why: these are the most significant genes that we want to label on the plots for better visualization and interpretation.
top_genes <- res_df[order(res_df$padj), ][1:10, ]
top_genes$gene <- rownames(top_genes)

# volcano plot with labels
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 100
  ) +
  scale_color_manual(values = c(
    "Down" = "blue",
    "Not Significant" = "gray",
    "Up" = "red"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
# Volcano Plot Code Explanation:
# Axes: x-axis is log2 fold change (how much a gene's expression changes), y-axis is -log10 of the adjusted p-value (significance of that change).
# Points: Each point is a gene, colored by significance (red for upregulated, blue for downregulated, gray for not significant).
# Labels: The most significant genes (top 10) are labeled with their names for easy identification.
# Threshold Lines: Dashed lines indicate thresholds for significance (p-value < 0.05) and fold change (log2 fold change > 1 or < -1).



# MA plot with labels
ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 100
  ) +
  scale_x_log10() +
  scale_color_manual(values = c(
    "Down" = "blue",
    "Not Significant" = "gray",
    "Up" = "red"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean Expression (log scale)",
    y = "Log2 Fold Change"
  )

# MA Plot Code Explanation:
# Axes: x-axis is the mean expression of each gene across all samples (on a log scale), y-axis is the log2 fold change (how much a gene's expression changes).
# Points: Each point is a gene, colored by significance (red for upregulated, blue for downregulated, gray for not significant).
# Labels: The most significant genes (top 10) are labeled with their names for easy identification.
# Threshold Line: A dashed line at y=0 indicates no change in expression (log2 fold change of 0 means the gene is not differentially expressed).

# Save plots
pdf("normalized_counts/volcano_plot.pdf")
print(volcano_plot)
dev.off()

pdf("normalized_counts/ma_plot.pdf")
print(ma_plot)
dev.off()