# HBC revalidation after deterministic SCT sampling

This note records measured results only. R/Seurat ran as an independent
black-box validation process; neither BioLang nor its MIT libraries load or
link the R packages. The GPL SCTransform provider remains a separate process.

## SCTransform boundary

On both filtered HBC samples, ordinary density-weighted sampling selected the
same 2,000 fit genes as R (2,000/2,000). No oracle fit-gene list was supplied.
The platform-dependent one-ULP `exp` result that previously changed this
sequential sample is now replaced by a deterministic high-precision evaluation.

For both control and stimulated samples:

- all 3,000 top features and their order match;
- median raw-theta relative error is about `5e-8`;
- residual RMSE divided by the R residual standard deviation is between
  `1.54e-8` and `1.80e-8`;
- every scale-sensitive acceptance gate passes.

Artifacts:

- `hbc-review-20260813-sct-ctrl-repro-exp-comparison.json`
- `hbc-review-20260813-sct-stim-repro-exp-comparison.json`
- `hbc-review-20260813-sct-sampling-ctrl-repro-exp.json`
- `hbc-review-20260813-sct-sampling-stim-repro-exp.json`

## SelectIntegrationFeatures and PrepSCTIntegration

The original provider workflow incorrectly treated a capped `ranking.csv` as
the complete model axis. The provider now exports `modelled-genes.csv`, and
multi-object feature eligibility uses that complete axis. The independently
selected integration features then match Seurat exactly: 3,000/3,000, in the
same order.

A fresh Seurat 5.5.1 black-box run exported the matrices immediately after
`PrepSCTIntegration`. Comparing all values against the provider gives:

| Sample | relative matrix RMS error | median per-gene relative RMSE |
|---|---:|---:|
| control | `2.21e-8` | `1.79e-9` |
| stimulated | `1.81e-8` | `2.21e-9` |

Artifact: `hbc-review-20260813-prepsct-provider-vs-r.json`.

## CCA and anchors on the exact R PrepSCT matrices

The remaining difference begins after PrepSCT. On identical R-exported input:

- mean assigned CCA correlation: `0.99999999585`;
- median principal angle: `0.000932` degrees;
- candidate anchors: 29,882 common, `0.99850` recall of the smaller set;
- retained anchors: 19,213 common, `0.99901` recall of the smaller set;
- all 193 filter features match;
- common-anchor filter decisions agree exactly;
- R elapsed time: `729.46 s`; BioLang elapsed time: `215.083 s`
  (`3.39x` faster for this fixed-input boundary).

The residual mismatch is therefore a small CCA/nearest-neighbour boundary
difference, not SCT, feature selection, PrepSCT, filtering arithmetic, or the
fixed-input downstream graph/Louvain implementation.

Artifact: `hbc-review-20260813-cca-real-r-input-comparison.json`.

## End-to-end interpretation

Using the correct exact feature set produced 22 BioLang clusters versus 19 in
Seurat, ARI `0.65347`, broad-cell-type exact agreement `93.99%`, and mapped
marker-pair recall `90.56%`. An earlier ARI near `0.70` included a stale capped
feature axis and was partly accidental cancellation; it must not be presented
as the result of the corrected workflow.

Exact cluster equivalence is not yet achieved. The current defensible claim is
near-numerical SCT/PrepSCT equivalence and 99.9% retained-anchor recall, with
remaining instability at CCA neighbour boundaries amplified by community
detection.

## Exact-boundary isolation (2026-08-14)

A fresh Seurat 5.5.1 run was generated rather than reusing the older 19/20
cluster fixtures. With the current inputs and package versions, the independent
reference has **22 clusters** and 19,232 retained anchors. This supersedes the
stale cluster counts for this stage-isolation experiment.

Full-precision CCA embeddings were exchanged as neutral `BLMATF64` matrices.
On those fixed embeddings, BioLang now matches every discrete anchor operation:

| Boundary | Agreement |
|---|---:|
| Candidate anchors | 29,927 / 29,927 exact |
| Retained anchors | 19,232 / 19,232 exact |
| Cross-neighbour ranks | 888,870 / 888,870 exact |
| Within-dataset neighbour ranks | 888,870 / 888,870 exact |
| Raw anchor scores | 19,232 / 19,232 exact |
| Normalized scores within `1e-12` | 19,232 / 19,232 |

The scoring mismatch was a real implementation bug. Seurat requests 31 Annoy
neighbours for `k.score=30`, removes the self hit, and scores the first 30.
BioLang previously requested only 30 total candidates. Because Annoy derives
its default search budget from the requested count, the apparently discarded
extra neighbour changed some approximate-search results. BioLang now uses the
same 31-candidate budget.

Fresh fixed-boundary downstream results are also exact. With the reference
19,232 anchors and weighting reduction, BioLang's direct-matrix restarted-
Lanczos PCA followed by its own Annoy SNN and Louvain produces **22 clusters and
ARI 1.0000 across all 29,629 cells**. Feeding the exact 40 reference PCs to the
BioLang graph/Louvain stages independently gives the same exact partition.
The strict BioLang correction/PCA/graph/cluster run took 31.827 seconds; the
corresponding independent R stages took 67.43 seconds. The reference CCA export
itself took 711.9 seconds and is not included in either downstream comparison.

The ordinary MIT-native CCA remains the sole exact-partition blocker. On the
same byte-identical PrepSCT matrices it produced 29,930 candidates and 19,233
retained anchors. Candidate Jaccard is 0.99690, retained-anchor Jaccard is
0.99797, mean assigned component correlation is 0.9999999958, and the maximum
principal angle is 0.02923 degrees. Propagating those native anchors through
the same strict downstream path yields 20 clusters and ARI 0.80802. The
paper-derived restarted-Lanczos CCA variant yields 20 clusters and ARI 0.80541.
This demonstrates that correlation-level CCA agreement is not sufficient for
an exact community partition near unstable graph boundaries.

BioLang therefore exposes two honestly distinct modes:

- native mode is MIT-only, fast, and scientifically close, but is not an exact
  Seurat partition on this HBC fixture;
- strict replay accepts full-precision embeddings, filter indices, weighting
  reduction, and a seeded PCA configuration from an independently installed
  provider. BioLang then reproduces anchors, correction, graph, and clusters
  exactly without linking or embedding that provider.

Principal generated artifacts are
`hbc-review-20260814-fixed-r-binary-anchors-k31-comparison.json`,
`hbc-review-20260814-anchor-neighbours-k31-comparison.json`,
`hbc-review-20260814-anchor-swap-lanczos-labels.json`,
`hbc-review-20260814-internal-cca-lanczos-pca-labels.json`, and
`hbc-review-20260814-anchor-swap-lanczos-pca-comparison.json`.

## Automatic strict-provider run (2026-08-14)

The same boundary was then exercised through the user-facing automatic path,
not by manually loading oracle matrices. `sc_find_anchors(...,
{external_provider: "auto"})` launched the separately installed
`bl-seurat-provider` process for CCA and the weighting reduction; BioLang
performed anchor search/filter/score and correction, then invoked that provider
again for seeded PCA before performing SNN construction and Louvain itself.

Against the fresh Seurat 5.5.1 reference, this run reproduced all discrete
results:

| Boundary | Automatic BioLang result |
|---|---:|
| Candidate anchors | 29,927 / 29,927 exact |
| Retained anchors | 19,232 / 19,232 exact |
| Raw anchor scores | 19,232 / 19,232 exact |
| Normalized scores within `1e-12` | 19,232 / 19,232 |
| Cell labels | 29,629 / 29,629 in the same partition |
| Clusters / ARI | 22 / `1.0000` |

The complete automatic run took `782.443 s`: anchors and their external CCA/
weighting boundary `719.126 s`, correction `11.281 s`, external PCA `35.190 s`,
SNN graph `7.541 s`, and Louvain `7.767 s`. During the run, Windows reported
an observed peak working set of `4.806 GiB` for the R process and `3.594 GiB`
for the BioLang process. These are per-process observed peaks, not a sampled
whole-system peak and not necessarily simultaneous; they must not be added and
presented as a formal peak-memory benchmark.

This exact mode is deliberately not the native performance claim: it executes
the pinned R/Seurat numerical boundary in a separate GPL-3.0-only program.
`bl.exe`, its default integration path, the exchanged `BLMATF64` matrices, and
BioLang's anchor/correction/graph code remain MIT. The run artifact names the
external backend, and current provider manifests record protocol, provider,
Seurat, and irlba versions. WebAssembly continues to use the native path and
cannot launch the provider.

Automatic-run evidence is
`hbc-review-20260814-external-provider-automatic-run3/summary.csv`,
`hbc-review-20260814-external-provider-automatic-anchors.json`, and
`hbc-review-20260814-external-provider-automatic-labels.json`.
