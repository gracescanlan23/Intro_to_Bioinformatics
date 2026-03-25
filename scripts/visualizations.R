library(DESeq2)
library(ggrepel)

source("scripts/normalize_counts.R")

dds <- DESeq(dds)
res <- results(dds)

# volcano plot

library(ggplot2)

# convert results to dataframe
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# define significance
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# top 10 genes for labeling
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

pdf("normalized_counts/volcano_plot.pdf")
print(volcano_plot)
dev.off()

pdf("normalized_counts/ma_plot.pdf")
print(ma_plot)
dev.off()