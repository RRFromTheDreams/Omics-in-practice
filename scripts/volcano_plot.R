
library(data.table)
library(ggplot2)
library(ggrepel)

# Load DESeq2 results
res <- fread("results/DE_8d_vs_0h.tsv")

# Extract stable Ensembl IDs
res$stable_id <- sub("\\..*$", "", res$gene_id)

# Identify genes of interest
res$gene_label <- NA_character_

res$gene_label[
  res$stable_id == "ENSG00000136997"
] <- "MYC"

res$gene_label[
  res$stable_id == "ENSG00000108312"
] <- "UBTF"

# Define thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Classify genes
res$status <- "Not significant"

res$status[
  !is.na(res$padj) &
  res$padj < padj_cutoff &
  res$log2FoldChange >= lfc_cutoff
] <- "Upregulated"

res$status[
  !is.na(res$padj) &
  res$padj < padj_cutoff &
  res$log2FoldChange <= -lfc_cutoff
] <- "Downregulated"

# Convert adjusted p-values for plotting
res$minus_log10_padj <- -log10(
  pmax(res$padj, .Machine$double.xmin)
)

# Plot
p <- ggplot(
  res,
  aes(log2FoldChange, minus_log10_padj, color = status)
) +
  geom_point(alpha = 0.6, size = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  scale_color_manual(values = c(
    "Downregulated" = "blue",
    "Not significant" = "grey75",
    "Upregulated" = "red"
  )) +
  geom_text_repel(
    data = subset(res, !is.na(gene_label) &
                  is.finite(minus_log10_padj)),
    aes(label = gene_label),
    color = "black",
    size = 4,
    max.overlaps = Inf
  ) +
  labs(
    title = "K562 BRB-seq: 8 days vs 0 hours",
    subtitle = "DESeq2 differential expression",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Gene classification"
  ) +
  theme_bw(base_size = 13)

ggsave(
  "figures/volcano_8d_vs_0h.png",
  p, width = 9, height = 7, dpi = 300
)

ggsave(
  "figures/volcano_8d_vs_0h.pdf",
  p, width = 9, height = 7
)

# Display gene counts per category
print(table(res$status))

# Display MYC and UBTF results
print(
  res[!is.na(res$gene_label),
      c("gene_label", "log2FoldChange", "padj"),
      with = FALSE]
)

cat("Volcano plot completed successfully.\n")
