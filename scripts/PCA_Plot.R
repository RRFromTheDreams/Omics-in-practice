
library(DESeq2)
library(data.table)
library(ggplot2)
library(ggrepel)

# Load existing transformed data
vsd <- readRDS("results/vsd_exploratory.rds")
metadata <- fread("metadata/sample_metadata.tsv")

# Match metadata to the transformed count matrix
metadata <- metadata[
  match(colnames(vsd), metadata$sample),
]

stopifnot(
  !anyNA(metadata$sample),
  identical(colnames(vsd), metadata$sample)
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("0h","6h","1d","2d","3d","4d","7d","8d")
)

# Calculate PCA from the VST data
pca <- plotPCA(
  vsd,
  intgroup = "sample_number",
  returnData = TRUE
)

# Add experimental conditions
pca$condition <- metadata$condition
pca$sample_number <- metadata$sample_number

percent_var <- round(100 * attr(pca, "percentVar"))

# Create improved figure
p <- ggplot(
  pca,
  aes(PC1, PC2, color = condition)
) +
  geom_point(size = 3.5) +
  geom_text_repel(
    aes(label = sample_number),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of K562 BRB-seq samples",
    subtitle = "Samples colored by drug-treatment time",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Treatment time"
  ) +
  coord_equal() +
  theme_bw(base_size = 13)

ggsave(
  "figures/PCA_with_conditions.png",
  p, width = 9, height = 7, dpi = 300
)

ggsave(
  "figures/PCA_with_conditions.pdf",
  p, width = 9, height = 7
)

cat("Updated PCA completed successfully.\n")
