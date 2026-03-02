library(DESeq2)
library(pheatmap)
library(matrixStats)
library(here)

# Load normalized data
vsd <- readRDS(here("normalized_counts", "vst_object.rds"))
metadata <- read.table(here("config", "metadata.tsv"), header = TRUE, sep = "\t")

mat <- assay(vsd)

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

library(matrixStats)

# 1) make sure metadata rownames are the sample IDs
rownames(metadata) <- metadata$sample_id

# 2) build annotation df in the exact same order as the heatmap columns
anno <- metadata[colnames(mat), c("condition", "tumor_type"), drop = FALSE]

# 3) quick check
stopifnot(identical(rownames(anno), colnames(mat)))

# 4) heatmap with annotation
pheatmap(mat[topVarGenes, ],
         scale = "row",
         annotation_col = anno,
         main = "Top 50 Variable Genes Heatmap")

# k-means Clustering of Samples

set.seed(123)
k <- 3
km <- kmeans(t(mat), centers = k)

km$cluster

# Save plots to files

pdf(here("normalized_counts", "clustering_plots.pdf"))

plot(hc)
pheatmap(sampleDistMatrix)
pheatmap(mat[topVarGenes, ], scale="row")

dev.off()