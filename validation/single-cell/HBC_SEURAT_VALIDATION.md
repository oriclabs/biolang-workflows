# HBC Seurat validation boundary

`hbc_seurat_reference.R` is an independent, validation-only execution of the
public HBC control/stimulated PBMC course workflow. It is not an implementation
source for BioLang and is not called by BioLang, Cargo, package tests, notebook
execution, or book builds.

This separation is deliberate. Seurat and SeuratObject are MIT-licensed, while
parts of their dependency graph use GPL-family licences. The R packages may be
installed and executed in a separate environment to measure compatibility, but
they are not BioLang runtime dependencies. MIT-covered Seurat R/C++ files may
be inspected or ported with their copyright and licence notice retained;
copyleft dependency implementations must not be copied into BioLang.

From the BioLang repository root:

```powershell
$env:BIOLANG_VALIDATION_R_LIB = (Resolve-Path .validation-r-library)
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  validation/single-cell/hbc_seurat_reference.R `
  ctrl_raw stim_raw validation-results/hbc-seurat
```

The script follows the HBC lesson calls and records:

- SHA-256 hashes of all six 10x input files;
- exact cell and gene QC checkpoints;
- the 3,000 selected integration features;
- per-cell cluster, UMAP, and 40-PC manifests;
- the full cluster trajectory for resolutions 0.4, 0.6, 0.8, 1.0, and 1.4;
- R/Seurat session information, logs, artifact sizes, and hashes.

The HBC course leaves `dims` unspecified for `FindIntegrationAnchors`, which
means Seurat's 1:30 default is used there. The course then uses PCs 1:40 for
UMAP and `FindNeighbors`. Keeping those two choices distinct is necessary for
an honest course reproduction.

Generated results belong under `validation-results/`, which is ignored because
the PC and cell manifests can be large. Small, reviewed evidence snapshots may
be copied into the validated book only after a successful run; claims in the
book must name their originating artifact and input hashes.

## Measured built-in baseline (2026-08-10)

The independent Seurat 5.5.1 run and the final BioLang CPU run joined all
29,629 cells exactly. Seurat returned 19 clusters and BioLang returned 21. The
standalone SCTransform curves are highly correlated, but calibrated HBC
conformance currently fails on raw theta and residual-variance slope (see
`SCTRANSFORM_ALGORITHM_AUDIT.md`). The end-to-end partition also does not
match: ARI is 0.5750, adjusted mutual information is 0.7269, and optimal
one-to-one mapped-cell accuracy is 69.92%. Integration-feature overlap is
2,716/3,000 (90.53%) and integrated-PC 15-neighbour mean Jaccard is 0.2278.

Marker validation found 2,993 overlapping mapped cluster/gene pairs, 191
overlapping mapped top-50 pairs, and matching peak clusters for 7/15 canonical
genes. See `validation-results/hbc-comparison-sct99-final-cpu.json` in the
validation workspace for the complete generated record.

Using Windows process `PeakWorkingSet64`, the paired resource run measured:

| Runtime | Wall time | Peak host working set |
|---|---:|---:|
| BioLang release, CPU | 481.2 s | 16.06 GiB |
| R 4.5.2 / Seurat 5.5.1 | 1548.3 s | 12.50 GiB |

BioLang is 3.22 times faster in that run, but uses 1.28 times more peak host
memory. It therefore does not yet meet the combined target of end-to-end 95%
scientific agreement while being both faster and lower-memory than Seurat.

## Real-data external-SCT validation (2026-08-12)

The process-isolated GPL SCTransform provider was then run on both complete HBC
filtered matrices (14,847 control cells and 14,782 stimulated cells). It wrote
5,000-feature Pearson-residual matrices through the `BLMATF64` file interchange;
BioLang performed feature selection, CCA anchors, integration, PCA, SNN,
Louvain, UMAP, and RNA-scale marker testing. The Seurat oracle was not rerun or
linked into BioLang: its hashed 2026-08-09 artifacts were consumed only by the
independent Python comparison scripts.

This run establishes two different levels of agreement which should not be
collapsed into one percentage:

| Comparison | Measured result |
|---|---:|
| Cells joined by sample and barcode | 29,629 / 29,629 |
| Top-3,000 integration-feature overlap | 90.67% |
| Exact numeric partition ARI / AMI | 0.5154 / 0.6998 |
| Optimal one-to-one numeric-cluster cell accuracy | 62.46% |
| Broad PBMC identity ARI | 0.9111 |
| Broad PBMC identity exact cell agreement | **95.11%** |
| Median within-cell-type control-to-stimulated effect correlation | 0.9260 |
| Median top-100 stimulated-gene overlap | 84% |
| Interferon panel direction | 12 / 12 positive in both, in all 7 broad types |

Broad identities are assigned independently for each result from
cluster-average `log1p(CP10K)` expression. Each panel marker is standardized
across clusters before panel averaging. This avoids an earlier invalid method
that treated absent positive-marker rows as zero and silently resolved tied
scores to the first cell type. The condition comparison is descriptive
pseudobulk only: HBC has one library per condition, so it cannot support
replicate-aware differential-expression p-values.

An isolation probe held Seurat's integrated PCs 1:40 fixed and ran only
BioLang's 20-neighbour SNN and Louvain resolution 0.8. BioLang returned 18
clusters versus Seurat's 19, with ARI 0.8901, AMI 0.9033, and 94.57% optimally
mapped cell agreement. The corresponding full-pipeline PC-neighbour Jaccard is
only 0.2150. Therefore most remaining exact-partition drift is in
integration/PCA, not SCTransform or graph clustering.

Measured sequential resources were:

| Stage | Wall time | Peak host working set |
|---|---:|---:|
| GPL SCT control | 4.04 s | 0.695 GiB |
| GPL SCT stimulated | 3.64 s | 0.695 GiB |
| BioLang downstream, CPU | 433.87 s | 17.36 GiB |
| Combined external-SCT + BioLang wall time | 441.56 s | 17.36 GiB |
| Seurat comparison run | 1548.3 s | 12.50 GiB |

BioLang is 3.51 times faster in these measured runs, but its peak working set is
1.39 times Seurat's. The speed target is met; the lower-memory target is not.
The validated scientific claim is therefore “95.11% broad cell-type agreement,”
not “95% identical Seurat clusters” and not complete numerical equivalence.

The generated evidence is under the ignored directory
`validation-results/hbc-biological-2026-08-12/`. The key files are
`comparison-v3.json`, `biology-v4.json`, and
`seurat-pcs-cluster-comparison.json`.
