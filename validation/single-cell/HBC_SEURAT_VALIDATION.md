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

## Current-code stage audit (2026-08-12 to 2026-08-13)

The current comparison uses a hash-pinned Seurat 5.5.1 / sctransform 0.4.3 R
environment only as an external oracle. GPL-family packages are not linked,
imported, or referenced by BioLang. Each stage below either consumes neutral
CSV/`BLMATF64` artifacts or compares independently generated outputs.

### SCTransform wrapper

Both complete HBC samples were checked through the same VST, mitochondrial
regression, ranking clip, and returned-assay clip used by the Seurat HBC flow.

| SCTransform metric | Control | Stimulated |
|---|---:|---:|
| Modelled-gene set | 13,799 / 13,799 | 13,695 / 13,695 |
| Raw theta slope | 1.000000007 | 1.000248 |
| Raw theta median / p90 relative error | 0.0000079% / 0.0000099% | 0.0858% / 1.8370% |
| Residual-variance slope / rank correlation | 0.99999995 / 0.999999999995 | 1.002966 / 0.99999959 |
| Top-3,000 feature overlap | 3,000 / 3,000 | 2,996 / 3,000 |
| Post-mito residual RMSE / R residual SD | 0.00000154% | 0.1223% |
| Wrapper compute time, R / Rust | 97.28 s / 1.80 s | 96.00 s / 1.98 s |

### Identical-input CCA and anchors

The same 3,000-feature residual matrices were supplied byte-for-byte to both
engines. BioLang's matrix-free CCA uses a 32-vector guard subspace and 12
block-power passes; it never allocates the quadratic 29,629-cell cross-product.

| Metric | Seurat | BioLang / agreement |
|---|---:|---:|
| Candidate anchors | 29,974 | 29,970 |
| Retained anchors | 19,246 | 19,244 |
| Assigned component correlation, mean / median | - | 0.999999996 / 0.9999999998 |
| Principal-angle median / maximum | - | 0.00092 degrees / 0.02893 degrees |
| CCA Procrustes relative error | - | 0.00009049 |
| Candidate identity recall / Jaccard | - | 99.89% / 99.77% |
| Retained identity recall / Jaccard | - | 99.93% / 99.84% |
| Filter-feature overlap | - | 193 / 193 |
| Common-anchor normalized-score correlation | - | 0.99717 |
| Normalized-score median / p90 absolute error | - | 2.8e-16 / 0.0303 |

Correlation does not conceal the remaining score difference: 81.21% of common
anchors have an identical integer raw score and 18.79% differ, although the
median normalized absolute error is effectively zero. The common candidate
pairs have identical filter decisions. The remaining difference enters at the
approximate-neighbour boundary: 37 candidate pairs and 16/14 retained pairs
are unique to Seurat/BioLang, respectively, and those small changes alter the
shared-neighbour scores used by integration.

### Weighting PCA and compiled integration kernel

The fresh 30-PC reduction used for integration weights is effectively
identical: minimum component correlation 0.9999999984, maximum principal angle
0.00423 degrees, and Procrustes relative error 0.00000569.

Seurat's compiled `FindWeightsC` and `IntegrateDataC` were then invoked as an
external black box with identical residual matrices, scored anchors, weighting
PCA, and parameters. A deterministic 512-cell by 500-feature sample spans the
whole corrected query assay.

| Corrected-matrix metric | Result |
|---|---:|
| Observations | 256,000 |
| Pearson correlation | **0.9998949** |
| Regression slope / intercept | 0.9997747 / -0.0000209 |
| RMSE / Seurat residual SD | 1.4501% |
| Median / p90 absolute error | 0.0000000416 / 0.001874 |
| Median row cosine | 0.999999999999998 |
| Seurat compiled correction time | 23.25 s |
| BioLang correction time, GPU auto / GPU off | 8.99 s / 8.57 s |

The GPU-auto and GPU-off BioLang sample files have the same SHA-256 hash. The
4.5-fold BioLang speedup over its earlier 40.52-second correction loop comes
from parallelising independent complete cell rows; summation order and output
are unchanged. BioLang is 2.6 times faster than the isolated Seurat compiled
correction in this measurement.

### End-to-end HBC result

Near-exact numeric stages do not imply identical community boundaries. The
current eight-pass full run joined all 29,629 cells and measured:

| End-to-end metric | Result |
|---|---:|
| Seurat / BioLang clusters at resolution 0.8 | 19 / 23 |
| Exact partition ARI / AMI | 0.6556 / 0.7733 |
| Optimal one-to-one mapped-cell accuracy | 75.75% |
| Integration-feature overlap | 2,998 / 3,000 (**99.93%**) |
| Integrated-PC 15-neighbour mean Jaccard | 0.3160 |
| UMAP 15-neighbour mean Jaccard | 0.04069 |
| Mapped positive-marker recall / Jaccard | 91.21% / 72.68% |
| Mapped top-50 marker recall | 83.37% |
| Canonical peak-cluster matches | 13 / 15 |
| Broad PBMC identity ARI / exact cell agreement | 0.8871 / **94.85%** |
| Median within-type stimulation-effect correlation | 0.9584 |
| Median top-100 stimulated-gene overlap | 87% |
| Interferon direction | 12 / 12 positive in both, in all 7 broad types |

With Seurat's cluster labels held fixed, BioLang's marker statistics remain
near exact: marker-pair Jaccard 99.75%, top-50 overlap 99.47%, floating-point
log2-fold-change agreement, and detection-fraction correlations above
0.999999. Most end-to-end marker drift is inherited from the different
partition rather than marker-test arithmetic.

The final CPU resource run scopes full SCT matrices and the anchor set to the
integration stage, then releases them before graph, UMAP, and marker work. The
29,629-point SVG is opt-in (`--write-svg`); every coordinate is still written
to `cells.csv`, so presentation-string construction is not included in the
scientific resource measurement.

| Complete run | Wall time | Peak host working set |
|---|---:|---:|
| BioLang, CPU, current scoped/no-SVG flow | **445.97 s** | **11.072 GiB** |
| Seurat 5.5.1 comparison run | 1548.3 s | 12.50 GiB |

BioLang was 3.47 times faster and used 11.4% less peak host memory in these
complete measured runs. Relative to the prior completed BioLang q8 run, peak
memory fell from 17.592 to 11.072 GiB and wall time fell from 475.10 to 445.97
seconds. The new run has byte-identical features, cells (including clusters and
UMAP), and markers; its 40 PCs correlate component-for-component at least
0.9999999999999998 (the CSV column order differs).

Therefore the defensible conclusion is: low-level SCT, weighting PCA, marker
statistics, and integration correction closely match; broad biological labels
agree for 94.85% of cells; exact Seurat cluster equivalence is not achieved;
and this measured workflow now meets the faster/lower-memory resource target.

### Maximum-ARI investigation (2026-08-13)

The current code was re-run end to end on all 29,629 HBC cells using the
standalone GPL-3.0-only Rust SCT executable, targeted integration-feature
residual materialization, uncentered integrated PCA, Spotify Annoy 1.17.3, a
32-vector/12-pass CCA guard subspace, and the Seurat-compatible Louvain
implementation. Resolution remained fixed at Seurat's 0.8; no oracle-guided
resolution was used.

| Current independent run | Result |
|---|---:|
| Seurat / BioLang clusters | 19 / 22 |
| Exact partition ARI / AMI | **0.70562 / 0.79223** |
| Optimal one-to-one mapped-cell accuracy | 77.21% |
| Broad-type ARI / exact agreement | 0.88896 / **94.80%** |
| Integration-feature overlap | 2,994 / 2,996 (99.93%) |
| Wall time / peak host working set | **446.44 s / 8.683 GiB** |

This is the current independently measured result at Seurat's unchanged
resolution and seed, improving the original 0.6556 result. It remains 3.47
times faster and uses 30.5% less peak host memory than the measured Seurat run
(1548.3 seconds and 12.50 GiB). It is not exact Seurat partition parity.

This resource figure is the final cold repeat after exposing the CCA controls.
It reproduced the prior production partition exactly: 22 clusters, ARI
0.7056166113, AMI 0.7922256023, and 94.80% broad-type cell agreement. An earlier
223.52-second run remains a valid observation but is not used for the headline
speed ratio because the later repeat is the more conservative measurement.

Louvain itself is no longer an unexplained discrepancy. Given Seurat's exact
1,114,296-edge SNN, BioLang now reproduces all 19 clusters with ARI, AMI, and
mapped accuracy equal to 1.0. The new Apache-2.0 Spotify Annoy bridge also
reproduces the graph from Seurat's fixed PCs exactly: all 562,951 non-self kNN
pairs, all 1,114,296 SNN edges, their weights, and the 19-cluster partition
match (ARI 1.0). With Seurat's scored anchors and corrected matrix held fixed,
the current BioLang uncentered PCA, graph, and clustering path reproduces the
complete 20-cluster fixed-input partition: ARI, AMI, and mapped accuracy are
all 1.0. The earlier recorded ARI 0.91005 came from an older PCA implementation
and is superseded by this measured rerun. These isolation results locate the
dominant remaining full-run drift before the corrected matrix, at the
CCA/anchor boundary rather than PCA or Louvain.

On byte-identical residual inputs, the wider, longer CCA solve plus Annoy now
gives 99.93% retained-anchor recall and 99.84% Jaccard; CCA component
correlation averages 0.999999996 and the maximum principal angle is 0.029
degrees. BioLang retained 19,161 anchors in the final independent full run.

A controlled hybrid-anchor experiment quantifies the downstream sensitivity.
Using Seurat pairs with BioLang scores gives ARI 0.86257; BioLang pairs with
Seurat scores gives ARI 0.86016; both BioLang pairs and scores gives ARI
0.79556; and both Seurat pairs and scores gives ARI 1.0. Thus both the 30
retained-pair boundary differences and the changed scores are material. This
does not justify tuning against cluster labels: it identifies the remaining
algorithmic boundary that a future clean-room approximate-SVD/indexing study
would need to reproduce.

A disclosed resolution sweep over one fixed BioLang embedding and graph tested
0.40 through 1.20 in increments of 0.05. It is a diagnostic upper bound, not
an independent result: selecting a resolution using the oracle labels would be
test-set leakage. The final Annoy/provider embedding reached a maximum ARI of
0.70686 at resolution 0.95. An earlier targeted-feature/GPU-search diagnostic
reached 0.71556 at resolution 0.75. The production default remains 0.8 and no
acceptance claim uses either oracle-selected value.

Fresh R black-box SCT oracles were regenerated from the exact current Matrix
Market fixtures after discovering that the earlier stored oracle was stale.
Against those current artifacts, control top-3,000 overlap is 3,000/3,000 and
residual RMSE is 0.0000015% of oracle residual SD; stimulated overlap is
2,996/3,000 and residual RMSE is 0.122% of SD. The larger integration-feature
loss had a separate cause. Seurat selects the union
of per-object variable features and `PrepSCTIntegration` calculates residuals
for selected genes missing from an object's stored SCT assay. The current
file protocol previously materialized the top 5,000 columns per sample in
advance, forcing BioLang to replace 103 selected genes. Protocol v1 now accepts
an explicit feature request after the ranking pass and materializes only those
residual columns. The final run integrated 2,996 model-eligible genes, 2,994 of
which are in Seurat's 3,000-feature set, without allocating all-gene residuals.

The provider now also accepts a barcode-aligned mitochondrial fraction and a
per-feature regression mask. These preserve metadata calculated before gene
filtering and reproduce `PrepSCTIntegration`'s distinction between existing
regressed `scale.data` rows and newly recomputed union-feature residuals. A
diagnostic run using Seurat's exact 3,000-feature list reached ARI 0.71869, but
that is not an independent validation result and is not the number reported
above.

The new evidence records are
`hbc-review-20260813-maxari-prepsct-lexical-q12-annoy-comparison.json`,
`hbc-review-20260813-maxari-prepsct-lexical-q12-annoy/resources.json`,
`hbc-review-20260813-sct-current-ctrl-comparison.json`, and
`hbc-review-20260813-sct-current-stim-comparison.json`. Stage-isolation records
include `hbc-review-20260813-graph-annoy-comparison.json`,
`hbc-review-20260813-cca-annoy-q12o32-comparison.json`, and
`hbc-review-20260813-anchor-swap-biolang-uncentered-comparison.json`.
The final default-regression records are
`hbc-review-20260813-production-final-cca-controls/resources.json`,
`hbc-review-20260813-production-final-cca-controls-comparison.json`, and
`hbc-review-20260813-production-final-cca-controls-biology.json`.

### Weight, SCT-sampling, and PCA isolation (2026-08-13)

The integration correction is now excluded as a source of the remaining
partition difference. With the same residual matrices, scored anchors, and
query weighting PCA, all non-zero anchor weights agree to floating-point
precision (global Pearson correlation 1.0; RMSE `4.73e-17`) and the corrected
matrix has Pearson correlation 1.0 and RMSE `1.35e-15`. The apparent 0.989
support Jaccard comes only from 563 serialized zero-valued farthest entries;
it is not a numerical correction difference.

The stimulated SCT density-weighted sample differs by 64 of 2,000 genes even
though its density weights agree with R to roughly 15 significant digits.
Sequential sampling without replacement amplifies those last-bit differences.
An explicitly labelled validation-only run forced Seurat's ordered fit-gene
sample. It then reproduced all 3,000 integration features and reduced the SCT
residual RMSE to 0.0161% of the oracle residual SD, but downstream partition
agreement became worse: ARI fell from **0.70562** to **0.65347**. The external
fit list is therefore retained only as a diagnostic and is not a production
default or a justified accuracy tweak.

The PCA convergence ceiling was independently investigated rather than inferred
from an end-to-end duration. On the real 14,847 by 3,000 control residual
matrix, the previous ceiling converged after 132 sweeps while a bounded
50-sweep run took 28.17 seconds versus 80.48 seconds in the paired benchmark.
Against the 132-sweep result, the 50-component median correlation is 1.0, the
minimum is 0.999926, and the mean aligned per-cell cosine is 0.9999994. The
production full run retains exactly the independently reported ARI
(`0.7056166113`). BioLang now reports `sweeps` and `converged` in the PCA result;
`BIOLANG_PCA_MAX_SWEEPS` remains available for a deliberately stricter tail.

The integrated-PCA trajectory was also measured on Seurat's exact corrected
matrix rather than inferred from the residual-matrix benchmark. At the current
50-sweep ceiling, the minimum same-index correlation across 50 PCs is
0.99999930, the maximum principal angle is 0.0681 degrees, and the Procrustes
relative error is 0.0000767. Feeding those PCs through the current Annoy/SNN/
Louvain path reproduces the fixed Seurat partition exactly (20/20 clusters,
ARI 1.0). Runs at 1, 2, 3, 4, 5, 8, 12, 20, and 30 sweeps record the convergence
trajectory; 50 is retained as the bounded production default.

CCA controls were tested as an algorithmic sensitivity check, not selected by
cluster ARI. Raising the solve from 12 to 20 passes did not improve the
byte-identical-input comparison, while reducing the guard subspace from 32 to
7 worsened the maximum principal angle from 0.0289 to 0.729 degrees and reduced
anchor agreement. The production defaults therefore remain 12 passes and 32
oversampling vectors. `sc_find_anchors` now accepts `cca_sweeps` and
`cca_oversample` so this convergence can be reproduced explicitly.

The corresponding evidence is
`hbc-review-20260813-integration-weight-comparison.json`,
`hbc-review-20260813-integrated-matrix-comparison-weight-trace.json`,
`hbc-review-20260813-sct-stim-oracle-fit-comparison.json`,
`hbc-review-20260813-oracle-fit-full-comparison-v4.json`,
`hbc-review-20260813-production-pca50-comparison.json`, and
`hbc-review-20260813-pca-50-vs-300.json`. The corrected-matrix and anchor
boundary records are
`hbc-review-20260813-anchor-swap-current-pca-comparison.json`,
`hbc-review-20260813-anchor-swap-current-clusters.json`, and
`hbc-review-20260813-anchor-sensitivity.json`.

Generated evidence is in the ignored `validation-results/hbc-review-20260812-*`
and `validation-results/hbc-review-20260813-*` directories. Principal records
are `hbc-review-20260812-cca-comparison-q8-scored.json`,
`hbc-review-20260812-weight-pca-comparison.json`,
`hbc-review-20260813-integrated-matrix-comparison-parallel.json`,
`hbc-review-20260812-q8-comparison.json`,
`hbc-review-20260812-q8-biology.json`, and
`hbc-review-20260812-fixed-markers-comparison-v2.json`. The final resource
record is
`hbc-review-20260813-scoped-nosvg-cpu/resources.json`; deterministic q8
reproduction is recorded in `hbc-review-20260813-q8-repro-pca.json`.

## Fresh exact-cluster boundary result (2026-08-14)

The reference and BioLang stages were rerun from the current exact PrepSCT
matrices. The fresh reference partition contains 22 clusters; older 19- and
20-cluster stage fixtures must not be mixed with this run.

On full-precision reference CCA embeddings, BioLang reproduces all 29,927
candidate anchors, all 19,232 retained anchors, every cross and within-dataset
neighbour rank, every raw score, and every normalized score within `1e-12`.
This rerun exposed and fixed an Annoy search-budget error: `k.score=30` needs a
31-neighbour self search before the self hit is removed.

With those exact scored anchors, the exact weighting reduction, and seeded
direct-matrix restarted-Lanczos PCA, BioLang independently performs correction,
PCA, Annoy SNN construction, and Louvain clustering and returns 22 clusters with
**ARI 1.0000 for all 29,629 cells**. Its downstream run took 31.827 seconds,
versus 67.43 seconds for the corresponding R stages. This is an exact fixed-
boundary replay result, not an end-to-end native-parity claim.

The remaining native difference is specifically the iterative CCA numerical
boundary. Native block-subspace CCA has mean assigned-component correlation
0.9999999958 and 0.99797 retained-anchor Jaccard, but its 19 changed/missing and
20 added retained pairs propagate to 20 clusters and ARI 0.80802. The native
restarted-Lanczos CCA variant gives ARI 0.80541. Consequently, the default
MIT-native pipeline is still not bit-for-bit Seurat-equivalent on this fixture.
Strict equivalence currently requires an independently generated full-precision
CCA artifact. BioLang accepts that artifact through a generic numeric contract;
it does not load or link the producer.

That strict contract is now automatic when the separately installed
GPL-3.0-only `bl-seurat-provider` is selected. A fresh real-data run through
the actual process boundary reproduced all 29,927 candidate anchors, all
19,232 retained anchors and scores, and the complete 22-cluster partition
across 29,629 cells (ARI 1.0000). It took 782.443 seconds end to end. This is
an exact optional external-provider result; the default MIT-native CCA result
and its ARI 0.80802 remain unchanged and distinctly reported.

The concise measured record and artifact list are in
[`HBC_REVALIDATION_2026-08-13.md`](HBC_REVALIDATION_2026-08-13.md).
