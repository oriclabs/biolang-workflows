#!/usr/bin/env Rscript

# Export Seurat's actual Annoy k.weight neighbours without repeating CCA or
# integration. Inputs are the already measured scored-anchor and weight-PCA
# artifacts.

suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: hbc_weight_neighbours_seurat.R CCA_DIR INTEGRATION_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
cca_dir <- args[[1]]
integration_dir <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

anchors <- read.csv(file.path(cca_dir, "anchors.csv"), check.names = FALSE)
pca_frame <- read.csv(file.path(integration_dir, "query-weight-pca.csv"), check.names = FALSE)
pc_columns <- paste0("PC_", seq_len(30))
pca <- as.matrix(pca_frame[, pc_columns, drop = FALSE])
rownames(pca) <- paste0("stim_", seq_len(nrow(pca)))
anchor_cells <- unique(anchors$right + 1L)
started <- proc.time()[["elapsed"]]
neighbours <- Seurat:::NNHelper(
  data = pca[anchor_cells, , drop = FALSE], query = pca, k = 100,
  method = "annoy", n.trees = 50, eps = 0
)
elapsed <- proc.time()[["elapsed"]] - started
sample_cells <- 1L + seq.int(0L, by = 28L, length.out = 512L)
write.csv(
  SeuratObject::Indices(neighbours)[sample_cells, , drop = FALSE],
  file.path(output_dir, "indices.csv"), row.names = FALSE, quote = FALSE
)
write.csv(
  SeuratObject::Distances(neighbours)[sample_cells, , drop = FALSE],
  file.path(output_dir, "distances.csv"), row.names = FALSE, quote = FALSE
)
write.csv(data.frame(
  query_cells = nrow(pca), unique_anchor_cells = length(anchor_cells),
  k_weight = 100, elapsed_seconds = elapsed
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
