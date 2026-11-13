
library(DESeq2)
library(apeglm)

# Load the fitted DESeq2 model
dds <- readRDS("results/dds_fitted.rds")

# Get the 8-day versus 0-hour result
res <- results(
  dds,
  contrast = c("condition", "8d", "0h"),
  alpha = 0.05
)

# Shrink fold changes, as in our volcano analysis
res_shrunk <- lfcShrink(
  dds,
  coef = "condition_8d_vs_0h",
  type = "apeglm"
)

# Create the MA plot
pdf("figures/MA_8d_vs_0h.pdf",
    width = 8, height = 6)

plotMA(
  res_shrunk,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "MA plot: 8 days vs 0 hours"
)

dev.off()

# Also save a PNG for PowerPoint
png("figures/MA_8d_vs_0h.png",
    width = 2400, height = 1800, res = 300)

plotMA(
  res_shrunk,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "MA plot: 8 days vs 0 hours"
)

dev.off()

cat("MA plots created successfully.\n")
