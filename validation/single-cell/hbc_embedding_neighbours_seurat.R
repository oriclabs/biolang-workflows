#!/usr/bin/env Rscript

# Validation-only export of the four neighbour-index matrices used by
# ScoreAnchors, from full-precision BLMATF64 embeddings.

suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: hbc_embedding_neighbours_seurat.R EMBEDDING_DIR OUTPUT_DIR")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_blmat <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  if (readChar(con, nchars = 8, useBytes = TRUE) != "BLMATF64") stop("invalid magic")
  header <- readBin(con, what = "raw", n = 16)
  low_u32 <- function(bytes) {
    sum(as.integer(bytes[1:4]) * c(1, 256, 65536, 16777216))
  }
  rows <- low_u32(header[1:8])
  columns <- low_u32(header[9:16])
  matrix(
    readBin(con, what = "numeric", n = rows * columns, size = 8, endian = "little"),
    nrow = rows, ncol = columns, byrow = TRUE
  )
}

left <- read_blmat(file.path(input_dir, "left-embedding.f64"))
right <- read_blmat(file.path(input_dir, "right-embedding.f64"))
rownames(left) <- paste0("left_", seq_len(nrow(left)))
rownames(right) <- paste0("right_", seq_len(nrow(right)))
started <- proc.time()[["elapsed"]]
indices <- function(data, query, k) {
  SeuratObject::Indices(Seurat:::NNHelper(
    data = data, query = query, k = k, method = "annoy", n.trees = 50
  ))
}
write.csv(indices(right, left, 30) - 1L,
          file.path(output_dir, "left-to-right-neighbours.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(indices(left, right, 30) - 1L,
          file.path(output_dir, "right-to-left-neighbours.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(indices(left, left, 31)[, seq_len(30), drop = FALSE] - 1L,
          file.path(output_dir, "left-within-neighbours.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(indices(right, right, 31)[, seq_len(30), drop = FALSE] - 1L,
          file.path(output_dir, "right-within-neighbours.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(data.frame(
  left_cells = nrow(left), right_cells = nrow(right), dimensions = ncol(left),
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "neighbour-summary.csv"), row.names = FALSE, quote = FALSE)
