# Install packages if not already installed
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("tidyverse")) install.packages("tidyverse")

library(DESeq2)
library(tidyverse)
library(here)

# create normalized counts directory if it doesn't exist
if (!dir.exists(normalized_dir)) {
  dir.create(normalized_dir, recursive = TRUE)
}

# set paths to project root and counts directory
counts_path <- here("counts", "gene_counts.txt")
meta_path <- here("config", "metadata.tsv")
normalized_dir <- here("normalized_counts")

# Load gene counts
counts <- read.table(counts_path,
                     header = TRUE,
                     row.names = 1,
                     comment.char = "#")

# Keep only sample columns (drop featureCounts annotation cols)
counts <- counts[, 6:ncol(counts)]

# Clean counts column names: HISAT2.SRRxxxx.sorted(.bam) -> SRRxxxx
colnames(counts) <- gsub("^HISAT2\\.", "", colnames(counts))
colnames(counts) <- gsub("\\.sorted(\\.bam)?$", "", colnames(counts))
colnames(counts) <- trimws(colnames(counts))

# Load metadata 
metadata <- read.table(meta_path, header = TRUE, sep = "\t")

# Set metadata rownames from sample_id (and trim whitespace)
metadata$sample_id <- trimws(metadata$sample_id)
rownames(metadata) <- metadata$sample_id

# Check what's missing (THIS should be character(0)) this section created by AI to troubleshoot dataframe problem with column names.
missing <- setdiff(rownames(metadata), colnames(counts))
missing

# Stop early if something is missing
stopifnot(length(missing) == 0)

# Reorder/subset
counts <- counts[, rownames(metadata)]

# Load metadata
metadata <- read.table(meta_path,
                       header = TRUE,
                       row.names = 1,
                       sep = "\t")

metadata

###### Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ tumor_type
)

# Remove very low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

###### Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

head(normalized_counts)

# Save normalized counts to file
write.csv(normalized_counts,
          file = "normalized_counts/normalized_gene_counts.csv",
          row.names = TRUE)

######## Before vs after normalization plots

# Before boxplot
boxplot(log2(counts + 1),
        main = "Raw Counts",
        xlab = "Samples",
        ylab = "Log2(Counts + 1)",
        las = 2)

# After boxplot
boxplot(log2(normalized_counts + 1),
        main = "Normalized Counts",
        xlab = "Samples",
        ylab = "Log2(Normalized Counts + 1)",
        las = 2)

# Save boxplots to file
pdf("normalized_counts/boxplots.pdf")

boxplot(log2(counts + 1),
        main = "Raw Counts",
        las = 2)

boxplot(log2(normalized_counts + 1),
        main = "Normalized Counts",
        las = 2)

dev.off()

# Overlay Density Plot raw vs normalized (robust)

pdf(here("normalized_counts", "raw_vs_normalized_density.pdf"))

raw_vec  <- log2(as.vector(as.matrix(counts)) + 1)
norm_vec <- log2(as.vector(as.matrix(normalized_counts)) + 1)

# drop any weird values just in case
raw_vec  <- raw_vec[is.finite(raw_vec)]
norm_vec <- norm_vec[is.finite(norm_vec)]

plot(density(raw_vec),
     col = "red",
     main = "Raw vs Normalized Counts",
     xlab = "Log2(counts + 1)")

lines(density(norm_vec), col = "blue")

legend("topright",
       legend = c("Raw", "Normalized"),
       col = c("red", "blue"),
       lwd = 2)

dev.off()

# PCA plot
# subsetting the DESeq object
dds_subset <- dds[, colnames(dds) != "SRR36750771"]

vst_subset <- vst(dds_subset, blind = FALSE)
plotPCA(vst_subset, intgroup = "tumor_type") +
  ggtitle("PCA of Normalized Counts (Subset)") +
  theme_minimal()

pdf("normalized_counts/pca_subset.pdf")
print(plotPCA(vst_subset, intgroup = "tumor_type"))

dev.off()

vsd <- vst(dds, blind = FALSE)
dir.create(here("normalized_counts"), showWarnings = FALSE, recursive = TRUE)

saveRDS(vsd, file = here("normalized_counts", "vst_object.rds"))

plotPCA(vsd, intgroup = "tumor_type") +
  ggtitle("PCA of Normalized Counts") +
  theme_minimal()

# Save PCA plot to file
pdf("normalized_counts/pca.pdf")
print(plotPCA(vsd, intgroup = "tumor_type"))
dev.off()