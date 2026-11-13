library(data.table)

mapping <- fread("metadata/sample_mapping.tsv")

# Confirm all 24 biological sample numbers are present
stopifnot(
  nrow(mapping) == 24,
  setequal(mapping$sample, 1:24),
  !anyDuplicated(mapping$barcode)
)

timepoints <- c(
  "0h", "6h", "1d", "2d",
  "3d", "4d", "7d", "8d"
)

# Samples 1-8, 9-16, 17-24 follow the same time sequence
metadata <- data.frame(
  sample = mapping$barcode,
  sample_number = mapping$sample,
  condition = timepoints[((mapping$sample - 1) %% 8) + 1],
  replicate_series = paste0(
    "Rep", ((mapping$sample - 1) %/% 8) + 1
  )
)

metadata$condition <- factor(
  metadata$condition,
  levels = timepoints
)

# Check three samples per time point
stopifnot(all(table(metadata$condition) == 3))

fwrite(
  metadata,
  "metadata/sample_metadata.tsv",
  sep = "\t"
)

print(metadata)
print(table(metadata$condition))

cat("Metadata created successfully.\n")
