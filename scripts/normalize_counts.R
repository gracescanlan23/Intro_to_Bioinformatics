# Install packages if not already installed
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("tidyverse")) install.packages("tidyverse")

library(DESeq2)
library(tidyverse)
library(here)

# set paths to project root and counts directory
counts_path <- here("counts", "gene_counts.txt")
meta_path <- here("config", "metadata.tsv")
normalized_dir <- here("normalized_counts")

###### Load gene counts
counts <- read.table(counts_path,
                     header = TRUE,
                     row.names = 1,
                     comment.char = "#")

# Remove annotation columns from Featurecounts
counts <- counts[, 6:ncol(counts)]

####### clean up column names #######

# Remove HISAT2. prefix
colnames(counts) <- sub("^HISAT2\\.", "", colnames(counts))

# Remove .sorted suffix
colnames(counts) <- sub("\\.sorted$", "", colnames(counts))

# Check result
head(colnames(counts))

head(counts)

# Load metadata
metadata <- read.table(meta_path,
                       header = TRUE,
                       row.names = 1,
                       sep = "\t")

metadata

# Make sure sample names in metadata match column names in counts
rownames(metadata) <- metadata$SampleID

counts <- counts[, rownames(metadata)]

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

# PCA plot

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "tumor_type") +
  ggtitle("PCA of Normalized Counts") +
  theme_minimal()

# AI used to troubleshoot dataframe problem with column names.