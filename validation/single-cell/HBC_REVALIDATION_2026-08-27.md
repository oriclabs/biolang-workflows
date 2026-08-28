# HBC release revalidation - 2026-08-27

This is a fresh execution of the current BioLang 1.5.0 release build over all
29,629 filtered HBC control/stimulated PBMC cells. The independent Seurat 5.5.1
oracle was not rerun: its six input hashes still match the source Matrix Market
files, and neither the pinned R environment nor the oracle script changed.
BioLang CPU and GPU outputs, the current GPL SCTransform executable, and the
optional strict Seurat provider were rerun.

The generated evidence is under the ignored `validation-results/` directory.
Every BioLang run has a sibling `biolang.run/v1` record containing the script,
input, output, module, executable, seed, backend, elapsed-time, and peak-RSS
hashes or measurements. Generated artifacts are deliberately not committed.

## Measured result

| Boundary | Current measurement |
|---|---:|
| Cells joined by sample and barcode | 29,629 / 29,629 |
| Seurat / native BioLang clusters at resolution 0.8 | 19 / 21 |
| Native partition ARI / AMI | **0.69627 / 0.79914** |
| Optimally mapped numeric-cluster cell accuracy | 79.63% |
| Broad PBMC identity ARI / exact agreement | 0.89101 / **94.92%** |
| Integration-feature overlap | **3,000 / 3,000** |
| Integrated-PC 15-neighbour mean Jaccard | 0.32109 |
| Mapped top-50 marker recall | 80.63% |
| Canonical peak-cluster matches | 14 / 15 |
| Fixed-Seurat-label marker pair Jaccard | **99.75%** |
| Fixed-Seurat-label top-50 overlap | **99.47%** |

Marker arithmetic is not the cause of the end-to-end marker difference. With
Seurat's labels held fixed, shared effect sizes agree to floating-point
precision and the feature sets are nearly exact. The end-to-end marker drift
is inherited from the different native partition.

The current standalone GPL SCTransform executable was rebuilt from commit
`fe236f1` and rerun independently against the stored R 4.5.2 /
`sctransform` 0.4.3 oracles. Both real samples pass every scale-sensitive gate:

| SCTransform boundary | Control | Stimulated |
|---|---:|---:|
| Top-3,000 feature overlap | 100% | 100% |
| Raw-theta regression slope | 0.999999989 | 0.999999989 |
| Raw-theta median relative error | 0.0000051% | 0.0000050% |
| Residual-variance slope | 0.999999948 | 0.999999959 |
| Residual RMSE / oracle residual SD | 0.00000154% | 0.00000180% |

The control `genes.csv`, ranking, residual probe, fit genes, fit cells, and
sampling artifacts are byte-identical to the previously validated build. The
new commit changes naming and the exposed conserve-memory contract, not these
validated numeric results.

## CPU, GPU, and strict compatibility mode

| Run | Wall time | Peak host working set | Scientific result |
|---|---:|---:|---|
| BioLang native CPU | **368.545 s** | **8.386 GiB** | 21 clusters; ARI 0.69627 |
| BioLang native GPU (`RTX 3080`, Vulkan) | 457.731 s | 8.563 GiB plus unmeasured device memory | identical to CPU |
| Seurat 5.5.1 full oracle | 1548.3 s | 12.50 GiB | 19 clusters |
| Strict external Seurat provider | 576.639 s | BioLang 3.521 GiB; R CCA observed about 4.03 GiB | 22-cluster supplied boundary, ARI 1.0 |

Native CPU was 4.20 times faster than the measured Seurat run and used 32.9%
less peak host memory. GPU was 24.2% slower than CPU on this workflow. The GPU
and CPU runs nevertheless produced identical ordered features, numeric cluster
labels, marker pairs, and PCs (relative RMSE 0). CPU is therefore the
recommended backend for this HBC workflow; GPU remains optional and its device
memory is not hidden inside the host-RSS claim.

Strict external mode reproduced all 29,927 candidate anchors, all 19,232
retained anchors, and the complete supplied partition exactly. Its anchor CSVs
are byte-identical to the previous automatic-provider run. This is an optional
GPL process boundary, not native MIT parity and not the native performance
claim.

## What still does not match

Native MIT mode does **not** reproduce the exact Seurat partition. ARI 0.69627
is a useful regression baseline, not "95% cluster equivalence." Broad identity
agreement is 94.92%, narrowly below the aspirational 95% target. The first
material divergence remains native CCA/anchor selection: the current native
run retained 19,231 anchors from 29,934 candidates, whereas the strict supplied
boundary retained 19,232 from 29,927. Small neighbor-boundary differences are
amplified by integration and community boundaries.

No resolution was selected against oracle labels. Resolution stays at the HBC
value 0.8, and the validation gates distinguish regression floors from the
still-unmet 95% targets.

## Automated checks

- `compare_hbc_backends.py` requires CPU/GPU cluster ARI, ordered features,
  PCs, and marker pairs to agree.
- `validate_hbc_release.py` combines SCT, native biology, fixed-label markers,
  strict replay, elapsed time, and memory into explicit regression gates.
- `compare_hbc_results.py` now rejects marker oracles whose feature labels are
  absent from the retained assay, preventing the old R-row-name export from
  silently producing a false marker regression.

The aggregate record is
`validation-results/hbc-current-20260827-validation.json`; all regression gates
pass. The aspirational native ARI 0.95, broad agreement 0.95, and GPU-faster-
than-CPU targets remain explicitly unmet.
