
library(DESeq2)
library(pheatmap)
library(data.table)

# Load variance-stabilized expression data
vsd <- readRDS("results/vsd_exploratory.rds")
vsd_mat <- assay(vsd)

# Load experimental metadata
metadata <- fread("metadata/sample_metadata.tsv")

metadata <- metadata[
  match(colnames(vsd_mat), metadata$sample),
]

stopifnot(identical(colnames(vsd_mat), metadata$sample))

# Genes included in the original heatmap
gene_ids <- c(
  MYC = "ENSG00000136997",
  UBTF = "ENSG00000108312",
  MAX = "ENSG00000125952",
  CDKN1A = "ENSG00000124762",
  MKI67 = "ENSG00000148773",
  BAX = "ENSG00000087088",
  BCL2 = "ENSG00000171791"
)

stable_ids <- sub("\\..*$", "", rownames(vsd_mat))
idx <- match(gene_ids, stable_ids)

if (anyNA(idx)) {
  stop(paste("Missing genes:", paste(names(gene_ids)[is.na(idx)],
                                   collapse = ", ")))
}

heat_mat <- vsd_mat[idx, , drop = FALSE]
rownames(heat_mat) <- names(gene_ids)

# Row z-score normalization
heat_z <- t(scale(t(heat_mat)))

# Annotate columns by treatment time
annotation <- data.frame(
  Treatment = factor(
    metadata$condition,
    levels = c("0h","6h","1d","2d","3d","4d","7d","8d")
  )
)

rownames(annotation) <- metadata$sample

# Name samples
colnames(heat_z) <- paste0("S", metadata$sample_number)
rownames(annotation) <- colnames(heat_z)

# Order by experimental time
ord <- order(annotation$Treatment, metadata$sample_number)
heat_z <- heat_z[, ord, drop = FALSE]
annotation <- annotation[ord, , drop = FALSE]

# Save heatmap
png("figures/MYC_UBTF_heatmap_conditions.png",
    width = 2500, height = 1700, res = 250)

pheatmap(
  heat_z,
  annotation_col = annotation,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  fontsize_row = 12,
  main = "MYC, UBTF and related genes across treatment time"
)

dev.off()

cat("Treatment-annotated heatmap completed.\n")
