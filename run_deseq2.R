
library(DESeq2)
library(data.table)
library(apeglm)

# 1. Load raw counts and metadata
df <- fread("data/gene_counts.tsv",
            data.table = FALSE,
            check.names = FALSE)

counts <- as.matrix(df[, -1, drop = FALSE])
rownames(counts) <- df[[1]]
storage.mode(counts) <- "integer"

metadata <- fread("metadata/sample_metadata.tsv",
                  data.table = FALSE)

# 2. Match metadata to the count matrix
metadata <- metadata[
  match(colnames(counts), metadata$sample),
  ,
  drop = FALSE
]

stopifnot(
  !anyNA(metadata$sample),
  identical(colnames(counts), metadata$sample)
)

rownames(metadata) <- metadata$sample

# 3. Set experimental conditions
metadata$condition <- factor(
  metadata$condition,
  levels = c("0h", "6h", "1d", "2d",
             "3d", "4d", "7d", "8d")
)

# 4. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# 5. Filter genes with very low expression
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Genes retained:", nrow(dds), "\n")

# 6. Fit the DESeq2 model
dds <- DESeq(dds)

# 7. Save the fitted model
saveRDS(dds, "results/dds_fitted.rds")

# 8. Print available coefficients
print(resultsNames(dds))

# 9. Compare every time point against 0h
timepoints <- c("6h", "1d", "2d",
                "3d", "4d", "7d", "8d")

for (tp in timepoints) {

  contrast <- c("condition", tp, "0h")

  res <- results(
    dds,
    contrast = contrast,
    alpha = 0.05
  )

  # Shrink log2 fold changes
  coef_name <- paste0(
    "condition_", tp, "_vs_0h"
  )

  res_shrunk <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )

  # Export results
  res_df <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)

  res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
  ]

  filename <- paste0(
    "results/DE_", tp, "_vs_0h.tsv"
  )

  fwrite(res_df, filename, sep = "\t")

  cat("\nComparison:", tp, "vs 0h\n")
  print(summary(res))
}

cat("\nDESeq2 analysis completed successfully.\n")
