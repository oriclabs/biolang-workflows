#!/usr/bin/env Rscript

# Regenerate the small MIT Seurat Louvain conformance fixture used by BioLang CI.
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: seurat_louvain_fixture.R SNN_EDGE_CSV")
edges <- read.csv(args[[1L]], check.names = FALSE)
if (all(c("i", "j", "w") %in% colnames(edges))) {
  edges <- data.frame(source = edges$i, target = edges$j, weight = edges$w)
}
n <- max(edges$source, edges$target) + 1L
snn <- Matrix::sparseMatrix(
  i = c(edges$source + 1L, edges$target + 1L),
  j = c(edges$target + 1L, edges$source + 1L),
  x = c(edges$weight, edges$weight), dims = c(n, n)
)
clusters <- Seurat:::RunModularityClustering(
  SNN = snn, modularity = 1L, resolution = 0.8,
  algorithm = 1L, n.start = 10L, n.iter = 10L,
  random.seed = 0L, print.output = FALSE
)
cat(paste(as.integer(clusters), collapse = ","), "\n")
