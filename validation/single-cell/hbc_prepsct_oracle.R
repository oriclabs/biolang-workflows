#!/usr/bin/env Rscript

# Black-box Seurat boundary export for independent BioLang validation.
# This script is validation-only: BioLang and the GPL SCT provider never load,
# link, or invoke it. It records the exact matrices handed from
# PrepSCTIntegration to FindIntegrationAnchors.

validation_library <- Sys.getenv("BIOLANG_VALIDATION_R_LIB", unset = "")
if (nzchar(validation_library)) {
  .libPaths(c(normalizePath(validation_library, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: hbc_prepsct_oracle.R CTRL_RAW STIM_RAW OUTPUT_DIR")
}
ctrl_dir <- normalizePath(args[[1L]], mustWork = TRUE)
stim_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
if (dir.exists(output_dir) && length(list.files(output_dir)) > 0L) {
  stop("output directory must be absent or empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_blmat <- function(path, matrix) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw("BLMATF64"), connection)
  # Header dimensions are uint64. R has no native uint64 writer, and these
  # dimensions fit in the low 32 bits, so emit the little-endian bytes.
  for (value in c(ncol(matrix), nrow(matrix))) {
    writeBin(as.raw(c(
      value %% 256, (value %/% 256) %% 256,
      (value %/% 65536) %% 256, (value %/% 16777216) %% 256,
      0, 0, 0, 0
    )), connection)
  }
  # A genes-by-cells R matrix is already laid out as all genes for cell 1,
  # followed by all genes for cell 2: exactly row-major cells-by-genes.
  writeBin(as.double(matrix), connection, size = 8L, endian = "little")
}

started <- proc.time()[["elapsed"]]
set.seed(123456L)
ctrl <- CreateSeuratObject(Read10X(ctrl_dir, gene.column = 2L),
                           min.features = 100L, project = "ctrl")
stim <- CreateSeuratObject(Read10X(stim_dir, gene.column = 2L),
                           min.features = 100L, project = "stim")
merged <- JoinLayers(merge(ctrl, stim, add.cell.ids = c("ctrl", "stim")))
rm(ctrl, stim)
merged$nUMI <- merged$nCount_RNA
merged$nGene <- merged$nFeature_RNA
merged$log10GenesPerUMI <- log10(merged$nGene) / log10(merged$nUMI)
merged$mitoRatio <- PercentageFeatureSet(merged, pattern = "^MT-") / 100
filtered <- subset(
  merged,
  subset = nUMI >= 500 & nGene >= 250 &
    log10GenesPerUMI > 0.80 & mitoRatio < 0.20
)
counts <- GetAssayData(filtered, layer = "counts")
counts <- counts[Matrix::rowSums(counts > 0) >= 10L, , drop = FALSE]
metadata <- filtered[[]]
filtered <- CreateSeuratObject(counts, meta.data = metadata, project = "HBC")
rm(merged, counts, metadata)
if (nrow(filtered) != 14065L || ncol(filtered) != 29629L) {
  stop("filtered HBC checkpoint failed")
}

objects <- SplitObject(filtered, split.by = "orig.ident")
rm(filtered)
for (sample in names(objects)) {
  objects[[sample]] <- SCTransform(
    objects[[sample]], vars.to.regress = "mitoRatio", vst.flavor = "v2",
    conserve.memory = TRUE, verbose = FALSE
  )
  variables <- VariableFeatures(objects[[sample]])
  write.csv(
    data.frame(rank = seq_along(variables), gene = variables),
    file.path(output_dir, paste0(sample, "-variable-features.csv")),
    row.names = FALSE, quote = FALSE
  )
}

features <- SelectIntegrationFeatures(objects, nfeatures = 3000L, verbose = FALSE)
write.csv(
  data.frame(rank = seq_along(features), gene = features),
  file.path(output_dir, "features.csv"), row.names = FALSE, quote = FALSE
)
objects <- PrepSCTIntegration(
  objects, anchor.features = features, verbose = FALSE
)
for (sample in names(objects)) {
  residuals <- GetAssayData(objects[[sample]], assay = "SCT", layer = "scale.data")
  residuals <- residuals[features, , drop = FALSE]
  write_blmat(file.path(output_dir, paste0(sample, ".f64")), residuals)
}

write.csv(
  data.frame(
    seurat = as.character(packageVersion("Seurat")),
    sctransform = as.character(packageVersion("sctransform")),
    glm_gam_poi = as.character(packageVersion("glmGamPoi")),
    cells = sum(vapply(objects, ncol, integer(1L))),
    features = length(features),
    elapsed_seconds = proc.time()[["elapsed"]] - started
  ),
  file.path(output_dir, "manifest.csv"), row.names = FALSE, quote = FALSE
)
