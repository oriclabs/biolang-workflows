#!/usr/bin/env Rscript

# Validation-only Seurat 5.5.1 CCA/anchor oracle over byte-identical matrices.

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: hbc_cca_seurat.R INPUT_DIR OUTPUT_DIR")
}
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("this oracle is pinned to Seurat 5.5.1")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_blmat_rows <- function(path, matrix) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("BLMATF64"), con)
  for (value in c(nrow(matrix), ncol(matrix))) {
    writeBin(as.raw(c(
      value %% 256, (value %/% 256) %% 256,
      (value %/% 65536) %% 256, (value %/% 16777216) %% 256,
      0, 0, 0, 0
    )), con)
  }
  # Transposing makes each original row contiguous in R's column-major layout.
  writeBin(as.double(t(matrix)), con, size = 8L, endian = "little")
}

features <- read.csv(file.path(input_dir, "features.csv"), check.names = FALSE)
nfeatures <- nrow(features)
read_blmat_transposed <- function(path) {
  bytes <- file.info(path)$size
  rows <- (bytes - 24) / (8 * nfeatures)
  if (rows != floor(rows)) stop("invalid BLMATF64 dimensions: ", path)
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readChar(con, nchars = 8, useBytes = TRUE)
  if (magic != "BLMATF64") stop("invalid BLMATF64 magic: ", path)
  seek(con, where = 24, origin = "start")
  values <- readBin(
    con, what = "numeric", n = rows * nfeatures,
    size = 8, endian = "little"
  )
  matrix(values, nrow = nfeatures, ncol = rows)
}

ctrl <- read_blmat_transposed(file.path(input_dir, "ctrl.f64"))
stim <- read_blmat_transposed(file.path(input_dir, "stim.f64"))
rownames(ctrl) <- rownames(stim) <- features$gene
colnames(ctrl) <- paste0("ctrl_", seq_len(ncol(ctrl)))
colnames(stim) <- paste0("stim_", seq_len(ncol(stim)))

started <- proc.time()[["elapsed"]]
cca <- Seurat:::RunCCA.default(
  object1 = ctrl, object2 = stim, standardize = TRUE,
  num.cc = 30, seed.use = 42, verbose = FALSE
)
raw <- cca$ccv
l2_rows <- function(matrix) {
  norms <- sqrt(rowSums(matrix * matrix))
  matrix / pmax(norms, .Machine$double.eps)
}
embedding <- l2_rows(raw)
left_embedding <- embedding[colnames(ctrl), , drop = FALSE]
right_embedding <- embedding[colnames(stim), , drop = FALSE]

ab <- Seurat:::NNHelper(
  data = right_embedding, query = left_embedding, k = 30,
  method = "annoy", n.trees = 50
)
ba <- Seurat:::NNHelper(
  data = left_embedding, query = right_embedding, k = 30,
  method = "annoy", n.trees = 50
)
ab_index <- SeuratObject::Indices(ab)
ba_index <- SeuratObject::Indices(ba)
candidates <- do.call(rbind, lapply(seq_len(nrow(ab_index)), function(left) {
  right <- ab_index[left, seq_len(5)]
  keep <- vapply(right, function(value) {
    left %in% ba_index[value, seq_len(5)]
  }, logical(1))
  if (!any(keep)) return(NULL)
  cbind(left = left, right = right[keep])
}))

left_raw <- raw[colnames(ctrl), , drop = FALSE]
right_raw <- raw[colnames(stim), , drop = FALSE]
loadings <- ctrl %*% left_raw + stim %*% right_raw
num_features <- vapply(seq_len(100), function(number) {
  length(unique(unlist(lapply(seq_len(30), function(dimension) {
    unlist(Seurat:::Top(
      data = loadings[, dimension, drop = FALSE],
      num = number, balanced = TRUE
    ))
  }))))
}, integer(1))
eligible <- num_features[num_features < 200]
max_per_dimension <- which.max(eligible)
top_features <- unique(unlist(lapply(seq_len(30), function(dimension) {
  unlist(Seurat:::Top(
    data = loadings[, dimension, drop = FALSE],
    num = max_per_dimension, balanced = TRUE
  ))
})))

left_filter <- Seurat:::L2Norm(t(ctrl[top_features, , drop = FALSE]), MARGIN = 1)
right_filter <- Seurat:::L2Norm(t(stim[top_features, , drop = FALSE]), MARGIN = 1)
filter_nn <- Seurat:::NNHelper(
  data = right_filter, query = left_filter, k = 200,
  method = "annoy", n.trees = 50
)
filter_index <- SeuratObject::Indices(filter_nn)
retained <- candidates[vapply(seq_len(nrow(candidates)), function(row) {
  candidates[row, "right"] %in% filter_index[candidates[row, "left"], ]
}, logical(1)), , drop = FALSE]

# ScoreAnchors uses the first k.score neighbours from the within- and
# cross-dataset searches. FindNeighbors asks for k.score + 1 within-dataset
# neighbours so that the self match is present in those first k.score rows.
aa <- Seurat:::NNHelper(
  data = left_embedding, query = left_embedding, k = 31,
  method = "annoy", n.trees = 50
)
bb <- Seurat:::NNHelper(
  data = right_embedding, query = right_embedding, k = 31,
  method = "annoy", n.trees = 50
)
aa_index <- SeuratObject::Indices(aa)
bb_index <- SeuratObject::Indices(bb)
offset <- nrow(left_embedding)
raw_scores <- vapply(seq_len(nrow(retained)), function(row) {
  left <- retained[row, "left"]
  right <- retained[row, "right"]
  neighbours_left <- c(aa_index[left, seq_len(30)], ab_index[left, seq_len(30)] + offset)
  neighbours_right <- c(ba_index[right, seq_len(30)], bb_index[right, seq_len(30)] + offset)
  length(intersect(neighbours_left, neighbours_right))
}, numeric(1))
low_score <- as.numeric(quantile(raw_scores, 0.01))
high_score <- as.numeric(quantile(raw_scores, 0.90))
scores <- pmin(1, pmax(0, (raw_scores - low_score) / (high_score - low_score)))

write.csv(
  data.frame(cell = seq_len(nrow(left_embedding)) - 1L, left_embedding),
  file.path(output_dir, "left-embedding.csv"), row.names = FALSE, quote = FALSE
)
write_blmat_rows(file.path(output_dir, "left-embedding.f64"), left_embedding)
write_blmat_rows(file.path(output_dir, "right-embedding.f64"), right_embedding)
write_blmat_rows(file.path(output_dir, "left-projection.f64"), left_raw)
write_blmat_rows(file.path(output_dir, "right-projection.f64"), right_raw)
write.csv(
  data.frame(cell = seq_len(nrow(right_embedding)) - 1L, right_embedding),
  file.path(output_dir, "right-embedding.csv"), row.names = FALSE, quote = FALSE
)
write.csv(
  data.frame(
    left = retained[, "left"] - 1L, right = retained[, "right"] - 1L,
    score = scores, raw_score = raw_scores
  ),
  file.path(output_dir, "anchors.csv"), row.names = FALSE, quote = FALSE
)
write.csv(
  data.frame(left = candidates[, "left"] - 1L, right = candidates[, "right"] - 1L),
  file.path(output_dir, "candidate-anchors.csv"), row.names = FALSE, quote = FALSE
)
write.csv(
  data.frame(
    rank = seq_along(top_features), gene = top_features,
    feature_index = match(top_features, rownames(ctrl)) - 1L
  ),
  file.path(output_dir, "filter-features.csv"), row.names = FALSE, quote = FALSE
)
write.csv(data.frame(
  ctrl_cells = ncol(ctrl), stim_cells = ncol(stim),
  features = nrow(ctrl), dimensions = ncol(embedding),
  candidates = nrow(candidates), retained = nrow(retained),
  elapsed_seconds = proc.time()[["elapsed"]] - started
), file.path(output_dir, "summary.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
