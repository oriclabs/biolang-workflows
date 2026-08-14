#!/usr/bin/env Rscript

# Seurat 5.5.1 weight-PCA oracle over the same merged SCT residual matrix.

suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: hbc_weight_pca_seurat.R INPUT_DIR OUTPUT_DIR")
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
features <- read.csv(file.path(input_dir, "features.csv"), check.names = FALSE)
nfeatures <- nrow(features)

read_blmat_transposed <- function(path) {
  rows <- (file.info(path)$size - 24) / (8 * nfeatures)
  con <- file(path, "rb")
  on.exit(close(con))
  if (readChar(con, nchars = 8, useBytes = TRUE) != "BLMATF64") {
    stop("invalid BLMATF64 magic: ", path)
  }
  seek(con, where = 24, origin = "start")
  matrix(
    readBin(con, what = "numeric", n = rows * nfeatures, size = 8, endian = "little"),
    nrow = nfeatures, ncol = rows
  )
}

ctrl <- read_blmat_transposed(file.path(input_dir, "ctrl.f64"))
stim <- read_blmat_transposed(file.path(input_dir, "stim.f64"))
rownames(ctrl) <- rownames(stim) <- features$gene
colnames(ctrl) <- paste0("ctrl_", seq_len(ncol(ctrl)))
colnames(stim) <- paste0("stim_", seq_len(ncol(stim)))
merged <- cbind(ctrl, stim)
centered <- merged - rowMeans(merged)
rm(merged, ctrl)
invisible(gc())

started <- proc.time()[["elapsed"]]
pca <- Seurat:::RunPCA.default(
  object = centered, npcs = 30, seed.use = 42,
  approx = TRUE, verbose = FALSE
)
scores <- SeuratObject::Embeddings(pca)[colnames(stim), , drop = FALSE]
write.csv(
  data.frame(cell = seq_len(nrow(scores)) - 1L, scores),
  file.path(output_dir, "query-weight-pca.csv"), row.names = FALSE, quote = FALSE
)
write.csv(data.frame(
  query_cells = nrow(scores), dimensions = ncol(scores),
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "weight-summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
