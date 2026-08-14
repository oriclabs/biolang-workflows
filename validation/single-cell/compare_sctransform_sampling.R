#!/usr/bin/env Rscript

# Validation-only diagnosis of the density-weighted step-one gene sample.
# It uses public R/sctransform functions in a separate process and exports only
# measurements; neither BioLang nor the GPL provider invokes this script.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: compare_sctransform_sampling.R INPUT_MEX PROVIDER_DIR ORACLE_DIR OUTPUT_JSON")
}
input_dir <- args[[1L]]
provider_dir <- args[[2L]]
oracle_dir <- args[[3L]]
output_path <- args[[4L]]

existing <- function(path) {
  if (file.exists(path)) path else if (file.exists(paste0(path, ".gz"))) paste0(path, ".gz") else stop("missing ", path)
}
umi <- as(readMM(existing(file.path(input_dir, "matrix.mtx"))), "CsparseMatrix")
features <- read.delim(existing(file.path(input_dir, "features.tsv")), header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(existing(file.path(input_dir, "barcodes.tsv")), header = FALSE, stringsAsFactors = FALSE)
rownames(umi) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
colnames(umi) <- as.character(barcodes[[1L]])

provider <- read.csv(file.path(provider_dir, "sampling.csv"), check.names = FALSE)
provider_fit <- read.csv(file.path(provider_dir, "fit-genes.csv"), check.names = FALSE)$gene
oracle_fit <- read.csv(file.path(oracle_dir, "fit-genes.csv"), check.names = FALSE)$gene

# These are exactly the public operations used by sctransform::vst after the
# candidate axis is established. Holding that axis fixed isolates density-grid
# and probability-sampling compatibility from upstream filtering.
log_gmean <- log10(sctransform:::row_gmean(umi[provider$gene, , drop = FALSE], eps = 1))
density_fit <- density(log_gmean, bw = "nrd", adjust = 1)
r_weights <- 1 / (approx(density_fit$x, density_fit$y, xout = log_gmean)$y + .Machine$double.eps)
manual_density_weights <- function(x, bw, accurate_tail) {
  n <- 512L
  from <- min(x) - 3 * bw
  to <- max(x) + 3 * bw
  lo <- from - 4 * bw
  up <- to + 4 * bw
  y <- .Call(stats:::C_BinDist, x, rep.int(1 / length(x), length(x)), lo, up, n)
  kords <- seq.int(0, ((2L * n - 1) / (n - 1)) * (up - lo), length.out = 2L * n)
  kords[(n + 2):(2 * n)] <- -kords[n:2]
  if (accurate_tail) {
    kernel <- dnorm(kords, sd = bw)
  } else {
    z <- kords / bw
    kernel <- 0.3989422804014327 * exp(-0.5 * z * z) / bw
  }
  convolved <- pmax.int(0, Re(fft(fft(y) * Conj(fft(kernel)), inverse = TRUE))[1L:n] / length(y))
  xords <- seq.int(lo, up, length.out = n)
  output_x <- seq.int(from, to, length.out = n)
  output_y <- approx(xords, convolved, output_x)$y
  1 / (approx(output_x, output_y, xout = x)$y + .Machine$double.eps)
}
naive_weights <- manual_density_weights(log_gmean, density_fit$bw, FALSE)
accurate_manual_weights <- manual_density_weights(log_gmean, density_fit$bw, TRUE)
from <- min(log_gmean) - 3 * density_fit$bw
to <- max(log_gmean) + 3 * density_fit$bw
r_sequence <- seq.int(from, to, length.out = 512L)
indices <- 0:511
step_sequence <- from + indices * ((to - from) / 511)
fraction_sequence <- from + (to - from) * (indices / 511)

set.seed(1448145L)
invisible(sample(colnames(umi), size = min(5000L, ncol(umi))))
r_fit <- sample(provider$gene, size = min(2000L, nrow(provider)), prob = r_weights)
set.seed(1448145L)
invisible(sample(colnames(umi), size = min(5000L, ncol(umi))))
r_fit_provider_weights <- sample(
  provider$gene,
  size = min(2000L, nrow(provider)),
  prob = provider$sampling_weight
)

signed_relative <- (provider$sampling_weight - r_weights) / pmax(abs(r_weights), .Machine$double.eps)
relative <- abs(signed_relative)
largest <- which.max(relative)
log_mean_error <- if ("log_geometric_mean" %in% colnames(provider)) {
  provider$log_geometric_mean - log_gmean
} else {
  rep(NA_real_, length(log_gmean))
}
largest_log_mean <- if (all(is.na(log_mean_error))) NA_integer_ else which.max(abs(log_mean_error))
metrics <- list(
  candidates = nrow(provider),
  r_bandwidth = unname(density_fit$bw),
  sequence_step_max_absolute_error = max(abs(r_sequence - step_sequence)),
  sequence_fraction_max_absolute_error = max(abs(r_sequence - fraction_sequence)),
  weight_pearson = unname(cor(provider$sampling_weight, r_weights)),
  provider_vs_naive_weight_max_relative_error = max(abs(provider$sampling_weight - naive_weights) / naive_weights),
  provider_vs_accurate_manual_weight_max_relative_error = max(abs(provider$sampling_weight - accurate_manual_weights) / accurate_manual_weights),
  accurate_manual_vs_density_max_relative_error = max(abs(accurate_manual_weights - r_weights) / r_weights),
  weight_relative_error_median = unname(median(relative)),
  weight_relative_error_p90 = unname(quantile(relative, 0.9)),
  weight_relative_error_max = unname(max(relative)),
  weight_signed_relative_error_median = unname(median(signed_relative)),
  weight_signed_relative_error_min = unname(min(signed_relative)),
  weight_signed_relative_error_max = unname(max(signed_relative)),
  largest_error_gene = provider$gene[[largest]],
  largest_error_provider_weight = provider$sampling_weight[[largest]],
  largest_error_r_weight = r_weights[[largest]],
  log_geometric_mean_max_absolute_error = max(abs(log_mean_error), na.rm = TRUE),
  largest_log_geometric_mean_error_gene = if (is.na(largest_log_mean)) NA_character_ else provider$gene[[largest_log_mean]],
  largest_log_geometric_mean_error_provider = if (is.na(largest_log_mean)) NA_real_ else provider$log_geometric_mean[[largest_log_mean]],
  largest_log_geometric_mean_error_r = if (is.na(largest_log_mean)) NA_real_ else log_gmean[[largest_log_mean]],
  largest_error_gene_log_mean_difference = log_mean_error[[largest]],
  r_sample_vs_oracle_intersection = length(intersect(r_fit, oracle_fit)),
  r_sample_vs_provider_intersection = length(intersect(r_fit, provider_fit)),
  r_provider_weight_sample_vs_oracle_intersection = length(intersect(r_fit_provider_weights, oracle_fit)),
  r_provider_weight_sample_vs_provider_intersection = length(intersect(r_fit_provider_weights, provider_fit)),
  provider_vs_oracle_intersection = length(intersect(provider_fit, oracle_fit)),
  sample_size = length(r_fit)
)
json <- jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE, digits = 16)
writeLines(json, output_path)
cat(json, "\n")
