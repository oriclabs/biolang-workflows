#!/usr/bin/env Rscript

# Validation-only Seurat 5.5.1 kNN/SNN/Louvain oracle over fixed integrated PCs.

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: hbc_graph_seurat.R PC_INPUT_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1L]]
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_blmat <- function(path, rows) {
  con <- file(path, "rb")
  on.exit(close(con))
  if (readChar(con, nchars = 8L, useBytes = TRUE) != "BLMATF64") {
    stop("invalid BLMATF64 magic: ", path)
  }
  columns <- (file.info(path)$size - 24) / (8 * rows)
  if (columns != floor(columns)) stop("invalid BLMATF64 dimensions: ", path)
  seek(con, where = 24, origin = "start")
  matrix(
    readBin(con, what = "numeric", n = rows * columns, size = 8L, endian = "little"),
    nrow = rows, ncol = columns, byrow = TRUE
  )
}

metadata <- read.csv(file.path(input_dir, "cells.csv"), check.names = FALSE)
pcs <- read_blmat(file.path(input_dir, "matrix.f64"), nrow(metadata))
if (nrow(pcs) != nrow(metadata) || ncol(pcs) < 40L) {
  stop("PC matrix/metadata shape mismatch")
}
pcs <- pcs[, seq_len(40L), drop = FALSE]
rownames(pcs) <- paste(metadata$sample, metadata$barcode, sep = "::")

started <- proc.time()[["elapsed"]]
nn <- Seurat:::NNHelper(
  data = pcs, query = pcs, k = 20L, method = "annoy",
  n.trees = 50L, eps = 0
)
indices <- SeuratObject::Indices(nn)
distances <- SeuratObject::Distances(nn)
knn <- data.frame(
  source = rep(seq_len(nrow(indices)) - 1L, each = ncol(indices)),
  rank = rep(seq_len(ncol(indices)), times = nrow(indices)),
  target = as.vector(t(indices - 1L)),
  distance = as.vector(t(distances))
)
write.csv(knn, file.path(output_dir, "knn.csv"), row.names = FALSE, quote = FALSE)

snn <- Seurat:::ComputeSNN(indices, prune = 1 / 15)
snn_triplets <- summary(Matrix::triu(snn, k = 1L))
snn_edges <- data.frame(
  source = snn_triplets$i - 1L,
  target = snn_triplets$j - 1L,
  weight = snn_triplets$x
)
write.csv(snn_edges, file.path(output_dir, "snn.csv"), row.names = FALSE, quote = FALSE)

clusters <- Seurat:::RunModularityClustering(
  SNN = snn, modularity = 1L, resolution = 0.8,
  algorithm = 1L, n.start = 10L, n.iter = 10L,
  random.seed = 0L, print.output = FALSE
)
write.csv(data.frame(
  sample = metadata$sample,
  barcode = metadata$barcode,
  cluster = as.integer(clusters)
), file.path(output_dir, "cells.csv"), row.names = FALSE, quote = FALSE)
elapsed <- proc.time()[["elapsed"]] - started
write.csv(data.frame(
  cells = nrow(pcs), pcs = ncol(pcs), k = 20L,
  knn_rows = nrow(knn), snn_edges = nrow(snn_edges),
  clusters = length(unique(clusters)), elapsed_seconds = elapsed
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
