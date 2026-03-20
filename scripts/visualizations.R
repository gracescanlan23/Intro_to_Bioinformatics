library(DESeq2)

source("scripts/normalize_counts.R")

dds <- DESeq(dds)
res <- results(dds)

# volcano plot

library(ggplot2)

# Convert to dataframe
res_df <- as.data.frame(res)


# Remove NA values
res_df <- na.omit(res_df)

# Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )


  # MA plot

  res_df <- as.data.frame(res)
res_df <- na.omit(res_df)

ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean Expression (log scale)",
    y = "Log2 Fold Change"
  )

# save volcano and MA plots to file

res_df <- as.data.frame(res)

# remove rows with missing values needed for plotting
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# define significance groups
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

pdf("normalized_counts/volcano_plot.pdf")
print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ))
dev.off()

pdf("normalized_counts/ma_plot.pdf")
print(ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean Expression (log scale)",
    y = "Log2 Fold Change"
  ))
dev.off()