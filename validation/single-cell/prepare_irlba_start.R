#!/usr/bin/env Rscript

# Validation-only export of the seeded normal start vector used by irlba.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: prepare_irlba_start.R LENGTH OUTPUT_CSV")
set.seed(42L)
write.csv(
  data.frame(value = rnorm(as.integer(args[[1L]]))), args[[2L]],
  row.names = FALSE, quote = FALSE
)
