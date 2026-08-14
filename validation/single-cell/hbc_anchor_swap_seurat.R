#!/usr/bin/env Rscript

# Seurat 5.5.1 oracle for the full neutral-matrix anchor-swap boundary. Both
# programs receive identical SCT residuals, scored anchors, weighting PCA, and
# cell order; only integration correction, PCA, graph, and Louvain differ.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: hbc_anchor_swap_seurat.R INPUT_DIR ANCHOR_DIR WEIGHT_DIR METADATA_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1L]]
anchor_dir <- args[[2L]]
weight_dir <- args[[3L]]
metadata_dir <- args[[4L]]
output_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

features <- read.csv(file.path(input_dir, "features.csv"), check.names = FALSE)
read_blmat_transposed <- function(path) {
  rows <- (file.info(path)$size - 24) / (8 * nrow(features))
  con <- file(path, "rb")
  on.exit(close(con))
  if (readChar(con, nchars = 8L, useBytes = TRUE) != "BLMATF64") {
    stop("invalid BLMATF64 magic: ", path)
  }
  seek(con, where = 24, origin = "start")
  matrix(
    readBin(con, what = "numeric", n = rows * nrow(features), size = 8L, endian = "little"),
    nrow = nrow(features), ncol = rows
  )
}
write_blmat <- function(path, matrix) {
  con <- file(path, "wb")
  on.exit(close(con))
  writeChar("BLMATF64", con, eos = NULL, useBytes = TRUE)
  writeBin(as.double(c(nrow(matrix), ncol(matrix))), con, size = 8L, endian = "little")
  # Header dimensions are uint64 in the protocol. These HBC dimensions are
  # exactly representable and R has no native uint64 writer, so emit the bytes.
  seek(con, where = 8L, origin = "start")
  dimensions <- c(nrow(matrix), ncol(matrix))
  for (value in dimensions) {
    bytes <- as.raw(c(
      value %% 256, (value %/% 256) %% 256, (value %/% 65536) %% 256,
      (value %/% 16777216) %% 256, 0, 0, 0, 0
    ))
    writeBin(bytes, con)
  }
  seek(con, where = 24L, origin = "start")
  writeBin(as.vector(t(matrix)), con, size = 8L, endian = "little")
}

ctrl <- read_blmat_transposed(file.path(input_dir, "ctrl.f64"))
stim <- read_blmat_transposed(file.path(input_dir, "stim.f64"))
rownames(ctrl) <- rownames(stim) <- features$gene
metadata <- read.csv(file.path(metadata_dir, "cells.csv"), check.names = FALSE)
if (nrow(metadata) != ncol(ctrl) + ncol(stim)) stop("cell metadata mismatch")
cell_names <- paste(metadata$sample, metadata$barcode, sep = "::")
colnames(ctrl) <- cell_names[seq_len(ncol(ctrl))]
colnames(stim) <- cell_names[ncol(ctrl) + seq_len(ncol(stim))]

anchors <- read.csv(file.path(anchor_dir, "anchors.csv"), check.names = FALSE)
filtered <- as.matrix(data.frame(
  cell1 = anchors$left + 1L, cell2 = anchors$right + 1L, score = anchors$score
))
weight_table <- read.csv(file.path(weight_dir, "query-weight-pca.csv"), check.names = FALSE)
weight_scores <- as.matrix(weight_table[paste0("PC_", seq_len(30L))])
rownames(weight_scores) <- colnames(stim)

started <- proc.time()[["elapsed"]]
anchor_cells <- unique(filtered[, "cell2"])
weight_neighbors <- Seurat:::NNHelper(
  data = weight_scores[anchor_cells, , drop = FALSE], query = weight_scores,
  k = 100L, method = "annoy", n.trees = 50L, eps = 0
)
distances <- 1 - SeuratObject::Distances(weight_neighbors) /
  SeuratObject::Distances(weight_neighbors)[, ncol(SeuratObject::Distances(weight_neighbors))]
integration_vectors <- t(stim[, filtered[, "cell2"], drop = FALSE]) -
  t(ctrl[, filtered[, "cell1"], drop = FALSE])
rownames(integration_vectors) <- colnames(stim)[filtered[, "cell2"]]
weights <- Seurat:::FindWeightsC(
  cells2 = 0:(ncol(stim) - 1L), distances = as.matrix(distances),
  anchor_cells2 = colnames(stim)[anchor_cells],
  integration_matrix_rownames = rownames(integration_vectors),
  cell_index = SeuratObject::Indices(weight_neighbors),
  anchor_score = filtered[, "score"], min_dist = 0, sd = 1,
  display_progress = FALSE
)
integrated_query <- Seurat:::IntegrateDataC(
  integration_matrix = as.sparse(integration_vectors),
  weights = as.sparse(weights), expression_cells2 = as.sparse(t(stim))
)
integrated <- as.matrix(cbind(ctrl, t(integrated_query)))
colnames(integrated) <- cell_names
correction_seconds <- proc.time()[["elapsed"]] - started
if (identical(tolower(Sys.getenv("HBC_WRITE_INTEGRATED_MATRIX")), "true")) {
  # Cell-major neutral matrix for independent PCA implementations. This is a
  # validation artifact only; normal runs avoid the additional 711 MB write.
  write_blmat(file.path(output_dir, "integrated.f64"), t(integrated))
}

pca_started <- proc.time()[["elapsed"]]
pca <- Seurat:::RunPCA.default(
  object = integrated, npcs = 50L, seed.use = 42L,
  approx = TRUE, verbose = FALSE
)
scores <- SeuratObject::Embeddings(pca)
loadings <- SeuratObject::Loadings(pca)
pca_seconds <- proc.time()[["elapsed"]] - pca_started
write_blmat(file.path(output_dir, "matrix.f64"), scores)
write_blmat(file.path(output_dir, "loadings.f64"), loadings)

graph_started <- proc.time()[["elapsed"]]
nn <- Seurat:::NNHelper(
  data = scores[, seq_len(40L), drop = FALSE], k = 20L,
  method = "annoy", n.trees = 50L, eps = 0
)
snn <- Seurat:::ComputeSNN(SeuratObject::Indices(nn), prune = 1 / 15)
clusters <- Seurat:::RunModularityClustering(
  SNN = snn, modularity = 1L, resolution = 0.8, algorithm = 1L,
  n.start = 10L, n.iter = 10L, random.seed = 0L, print.output = FALSE
)
graph_cluster_seconds <- proc.time()[["elapsed"]] - graph_started
write.csv(data.frame(
  sample = metadata$sample, barcode = metadata$barcode,
  cluster = as.integer(clusters)
), file.path(output_dir, "cells.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(
  cells = nrow(scores), features = nrow(integrated), anchors = nrow(filtered),
  pcs = ncol(scores), clusters = length(unique(clusters)),
  correction_seconds = correction_seconds, pca_seconds = pca_seconds,
  graph_cluster_seconds = graph_cluster_seconds,
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
