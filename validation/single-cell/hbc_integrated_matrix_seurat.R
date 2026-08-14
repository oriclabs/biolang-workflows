#!/usr/bin/env Rscript

# Validation-only Seurat 5.5.1 integration oracle over identical SCT matrices
# and the scored anchors exported by hbc_cca_seurat.R.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: hbc_integrated_matrix_seurat.R INPUT_DIR CCA_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1]]
cca_dir <- args[[2]]
output_dir <- args[[3]]
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
anchors <- read.csv(file.path(cca_dir, "anchors.csv"), check.names = FALSE)
required <- c("left", "right", "score")
if (!all(required %in% colnames(anchors))) {
  stop("CCA anchor artifact must include left, right, and score")
}
filtered_anchors <- as.matrix(data.frame(
  cell1 = anchors$left + 1L,
  cell2 = anchors$right + 1L,
  score = anchors$score
))

# Preserve the exact query weighting reduction as a second black-box boundary.
merged <- cbind(ctrl, stim)
centered <- merged - rowMeans(merged)
weight_pca <- Seurat:::RunPCA.default(
  object = centered, npcs = 30, seed.use = 42,
  approx = TRUE, verbose = FALSE
)
weight_scores <- SeuratObject::Embeddings(weight_pca)[colnames(stim), , drop = FALSE]
write.csv(
  data.frame(cell = seq_len(nrow(weight_scores)) - 1L, weight_scores),
  file.path(output_dir, "query-weight-pca.csv"), row.names = FALSE, quote = FALSE
)
rm(merged, centered, weight_pca)
invisible(gc())

started <- proc.time()[["elapsed"]]
anchor_cells <- unique(filtered_anchors[, "cell2"])
weight_neighbours <- Seurat:::NNHelper(
  data = weight_scores[anchor_cells, , drop = FALSE],
  query = weight_scores, k = 100,
  method = "annoy", n.trees = 50, eps = 0
)
distances <- SeuratObject::Distances(weight_neighbours)
distances <- 1 - distances / distances[, ncol(distances)]
integration_matrix <- t(stim[, filtered_anchors[, "cell2"], drop = FALSE]) -
  t(ctrl[, filtered_anchors[, "cell1"], drop = FALSE])
rownames(integration_matrix) <- colnames(stim)[filtered_anchors[, "cell2"]]
weights <- Seurat:::FindWeightsC(
  cells2 = 0:(ncol(stim) - 1),
  distances = as.matrix(distances),
  anchor_cells2 = colnames(stim)[anchor_cells],
  integration_matrix_rownames = rownames(integration_matrix),
  cell_index = SeuratObject::Indices(weight_neighbours),
  anchor_score = filtered_anchors[, "score"],
  min_dist = 0, sd = 1, display_progress = FALSE
)
diagnostic_cells <- 1L + seq.int(0L, by = 28L, length.out = 512L)
diagnostic_weights <- do.call(rbind, lapply(diagnostic_cells, function(cell) {
  column <- weights[, cell, drop = FALSE]
  entries <- summary(column)
  data.frame(
    query_cell = cell - 1L,
    anchor_index = entries$i - 1L,
    weight = entries$x
  )
}))
write.csv(
  diagnostic_weights,
  file.path(output_dir, "integration-weights.csv"),
  row.names = FALSE, quote = FALSE
)
integrated_query <- Seurat:::IntegrateDataC(
  integration_matrix = as.sparse(integration_matrix),
  weights = as.sparse(weights),
  expression_cells2 = as.sparse(t(stim))
)
elapsed <- proc.time()[["elapsed"]] - started

sample_cells <- diagnostic_cells
sample_features <- 1L + seq.int(0L, by = 6L, length.out = 500L)
sample_matrix <- as.matrix(integrated_query[sample_cells, sample_features, drop = FALSE])
write.csv(
  sample_matrix,
  file.path(output_dir, "integrated-query-sample.csv"),
  row.names = FALSE, quote = FALSE
)
write.csv(data.frame(
  query_cells = nrow(integrated_query), features = ncol(integrated_query),
  anchors = nrow(filtered_anchors), k_weight = 100,
  elapsed_seconds = elapsed
), file.path(output_dir, "integration-summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
