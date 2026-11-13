
library(data.table)

# Import the gene-count matrix
counts_df <- fread(
  "data/gene_counts.tsv",
  data.table = FALSE,
  check.names = FALSE
)

# Extract gene identifiers and counts
gene_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- gene_ids
storage.mode(counts) <- "numeric"

# Import barcode-to-sample mapping
mapping <- fread(
  "metadata/sample_mapping.tsv",
  data.table = FALSE
)

# Validate gene identifiers and count values
stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyNA(counts))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))

# Match barcode mapping to matrix columns
stopifnot(!anyDuplicated(mapping$barcode))
stopifnot(setequal(colnames(counts), mapping$barcode))

mapping <- mapping[
  match(colnames(counts), mapping$barcode),
  ,
  drop = FALSE
]

stopifnot(identical(colnames(counts), mapping$barcode))

# Display summary
cat("Genes:", nrow(counts), "\n")
cat("Samples:", ncol(counts), "\n")
cat("All checks passed!\n")

cat("\nLibrary sizes:\n")
print(summary(colSums(counts)))

cat("\nSample mapping:\n")
print(mapping)
