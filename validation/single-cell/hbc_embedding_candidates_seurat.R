#!/usr/bin/env Rscript

# Validation-only crossover: run Seurat's Annoy mutual-neighbour search on an
# embedding produced by either Seurat or BioLang.  This isolates the neighbour
# engine from CCA construction.

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: hbc_embedding_candidates_seurat.R EMBEDDING_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_embedding <- function(path) {
  frame <- read.csv(path, check.names = FALSE)
  as.matrix(frame[, setdiff(names(frame), "cell"), drop = FALSE])
}

left <- read_embedding(file.path(input_dir, "left-embedding.csv"))
right <- read_embedding(file.path(input_dir, "right-embedding.csv"))
if (ncol(left) != ncol(right)) stop("embedding dimensions differ")
rownames(left) <- paste0("left_", seq_len(nrow(left)))
rownames(right) <- paste0("right_", seq_len(nrow(right)))

started <- proc.time()[["elapsed"]]
ab <- Seurat:::NNHelper(
  data = right, query = left, k = 30,
  method = "annoy", n.trees = 50
)
ba <- Seurat:::NNHelper(
  data = left, query = right, k = 30,
  method = "annoy", n.trees = 50
)
ab_index <- SeuratObject::Indices(ab)
ba_index <- SeuratObject::Indices(ba)
candidates <- do.call(rbind, lapply(seq_len(nrow(ab_index)), function(left_cell) {
  right_cell <- ab_index[left_cell, seq_len(5)]
  keep <- vapply(right_cell, function(value) {
    left_cell %in% ba_index[value, seq_len(5)]
  }, logical(1))
  if (!any(keep)) return(NULL)
  cbind(left = left_cell, right = right_cell[keep])
}))

write.csv(
  data.frame(left = candidates[, "left"] - 1L, right = candidates[, "right"] - 1L),
  file.path(output_dir, "candidate-anchors.csv"), row.names = FALSE, quote = FALSE
)
write.csv(data.frame(
  left_cells = nrow(left), right_cells = nrow(right), dimensions = ncol(left),
  candidates = nrow(candidates),
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
