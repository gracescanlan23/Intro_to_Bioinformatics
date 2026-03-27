# Install packages if not already installed
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("apeglm")) BiocManager::install("apeglm")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if (!require("AnnotationDbi")) BiocManager::install("AnnotationDbi")

library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(here)

# set paths to project root and counts directory
counts_path <- here("counts", "gene_counts.txt")
meta_path <- here("config", "metadata.tsv")
normalized_dir <- here("normalized_counts")

# create normalized counts directory if it doesn't exist
if (!dir.exists(normalized_dir)) {
  dir.create(normalized_dir, recursive = TRUE)
}



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

# Sanity Checks
cat("Count matrix dimensions:\n")
print(dim(counts))
cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nPreview of count matrix:\n")
print(head(counts[, 1:min(5, ncol(counts))]))

cat("\nPreview of metadata:\n")
print(head(metadata))

missing <- setdiff(rownames(metadata), colnames(counts))
cat("\nSamples in metadata but not counts:\n")
print(missing)

stopifnot(length(missing) == 0)


# Reorder/subset
counts <- counts[, rownames(metadata)]

###### Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ tumor_type
)

# Make sure condition is a factor
dds$tumor_type <- factor(dds$tumor_type)

cat("\nTumor type levels:\n")
print(levels(dds$tumor_type))


# Remove very low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

cat("\nGenes remaining after filtering:\n")
print(nrow(dds))

####### Normalization Explanation
# DESeq2 was used for normalization because it provides a robust method for correcting differences
# in sequencing depth and library size across samples. RNA-seq datasets often contain variability due to techincal 
# factors, and DESeq2 addresses this by estimating size factors that scale counts to make samples comparable. 
# Unlike simple normalization methods, DESeq2 uses a median-of-ratios approach, which is less sensitive to extreme values and highly expressed genes. 
# This method assumes that most genes are not differentially expressed, allowing it to accurately estimate normalization
# factors across the dataset. Additionally, DESeq2 models count data using a negative bionomial distribution, which 
# is well-suited for RNA-seq data that exhibit overdispersion. This improves the reliability of downstream statistical
# testing for differential expression. Overall, DESeq2 normalization ensures that observed differences in gene expression are more
# likely to reflect true biological variation rather than technical bias. 


###### Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

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

