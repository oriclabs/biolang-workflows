#!/usr/bin/env Rscript

# Validation-only export of Seurat 5.5.1's MIT Standardize stage.  The output
# is cells x features in BioLang's row-major BLMATF64 format.

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: hbc_standardize_oracle.R INPUT_DIR OUTPUT_DIR")
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
features <- read.csv(file.path(input_dir, "features.csv"), check.names = FALSE)
nfeatures <- nrow(features)

read_blmat_transposed <- function(path) {
  bytes <- file.info(path)$size
  cells <- (bytes - 24) / (8 * nfeatures)
  if (cells != floor(cells)) stop("invalid BLMATF64 dimensions: ", path)
  con <- file(path, "rb")
  on.exit(close(con))
  if (readChar(con, nchars = 8, useBytes = TRUE) != "BLMATF64") stop("invalid magic")
  seek(con, where = 24, origin = "start")
  matrix(
    readBin(con, what = "numeric", n = cells * nfeatures, size = 8, endian = "little"),
    nrow = nfeatures, ncol = cells
  )
}

write_blmat_cells <- function(path, matrix) {
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(charToRaw("BLMATF64"), con)
  for (value in c(ncol(matrix), nrow(matrix))) {
    writeBin(as.raw(c(
      value %% 256, (value %/% 256) %% 256,
      (value %/% 65536) %% 256, (value %/% 16777216) %% 256,
      0, 0, 0, 0
    )), con)
  }
  writeBin(as.numeric(matrix), con, size = 8, endian = "little")
}

started <- proc.time()[["elapsed"]]
ctrl <- Seurat:::Standardize(
  read_blmat_transposed(file.path(input_dir, "ctrl.f64")), display_progress = FALSE
)
stim <- Seurat:::Standardize(
  read_blmat_transposed(file.path(input_dir, "stim.f64")), display_progress = FALSE
)
write_blmat_cells(file.path(output_dir, "ctrl.f64"), ctrl)
write_blmat_cells(file.path(output_dir, "stim.f64"), stim)
write.csv(features, file.path(output_dir, "features.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(
  ctrl_cells = ncol(ctrl), stim_cells = ncol(stim), features = nrow(ctrl),
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
