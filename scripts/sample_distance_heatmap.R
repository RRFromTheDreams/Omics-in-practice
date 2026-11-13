
library(DESeq2)
library(pheatmap)

# Load previously transformed data
vsd <- readRDS("results/vsd_exploratory.rds")
vsd_mat <- assay(vsd)

# Calculate distances between samples
sample_dists <- dist(t(vsd_mat))
sample_dist_mat <- as.matrix(sample_dists)

# Use biological sample numbers as labels
mapping <- read.delim("metadata/sample_mapping.tsv")
mapping <- mapping[match(colnames(vsd_mat), mapping$barcode), ]

stopifnot(identical(colnames(vsd_mat), mapping$barcode))

labels <- paste0("Sample ", mapping$sample)
rownames(sample_dist_mat) <- labels
colnames(sample_dist_mat) <- labels

# Draw heatmap
png("figures/sample_distance_heatmap.png",
    width = 2200, height = 2000, res = 250)

pheatmap(
  sample_dist_mat,
  clustering_distance_rows = as.dist(sample_dist_mat),
  clustering_distance_cols = as.dist(sample_dist_mat),
  main = "BRB-seq sample-to-sample distances",
  fontsize = 8
)

dev.off()

cat("Sample-distance heatmap completed.\n")
