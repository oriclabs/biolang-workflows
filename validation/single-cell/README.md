# Scanpy and Seurat validation

These checks generate a deterministic 10x MEX dataset locally. No count matrix
is committed or downloaded by the package.

Copy the installed package examples into a standalone working directory:

```text
bl examples singlecell --copy singlecell-examples
cd singlecell-examples
python make_demo_10x.py --output validation_data
bl run validation/biolang_reference.bl
python validation/compare_labels.py validation_data/truth.csv validation_biolang_labels.csv
python validation/scanpy_reference.py validation_data validation_scanpy_labels.csv
python validation/compare_labels.py validation_scanpy_labels.csv validation_biolang_labels.csv
Rscript validation/seurat_reference.R validation_data validation_seurat_labels.csv
python validation/compare_labels.py validation_seurat_labels.csv validation_biolang_labels.csv
```

The generated matrix and label CSV files are written under the copied working
directory. The full BioLang source repository is not required.

Python requires Scanpy and its Leiden dependencies. Use the pinned Python 3.12
environment below instead of a system Python installation:

```text
# Windows, from the copied example directory
py -3.12 -m venv .venv-scanpy
.venv-scanpy\Scripts\python -m pip install -r validation\requirements-scanpy.txt

# macOS/Linux, from the copied example directory
python3.12 -m venv .venv-scanpy
.venv-scanpy/bin/python -m pip install -r validation/requirements-scanpy.txt
```

Then substitute the environment's Python executable for `python` in the
commands above. The Numba constraint avoids a copied-function cache regression
seen with the Numba 0.66 Windows stack.

R requires `Seurat`, `Matrix`, and the Leiden backend used by Seurat. If
`Rscript` is not on `PATH`, invoke its absolute path. Keep both validation
environments separate; they are not runtime dependencies of the BioLang
package.

Cluster IDs themselves are arbitrary. The comparison joins cells by barcode and
uses adjusted Rand index, which compares partitions without requiring matching
numeric IDs. The default acceptance threshold is `0.70`; record the exact
package versions and inspect marker agreement before treating a result as
scientific equivalence.

## Verified reference run

The seeded fixture was validated on 2026-07-27:

| Runtime | Versions | Retained cells | Clusters | ARI vs BioLang |
|---|---|---:|---:|---:|
| BioLang | repository build after sparse single-cell changes | 220 | 4 | 1.0000 vs truth |
| Scanpy | Python 3.12.10, Scanpy 1.12.3, Numba 0.65.1, igraph 1.0.0, NumPy 2.4.6, pandas 3.0.5 | 220 | 4 | 1.0000 |
| Seurat | R 4.5.2, Seurat 5.5.1, Matrix 1.7.4, leidenbase 0.1.37 | 220 | 4 | 1.0000 |

This verifies the deterministic synthetic partition and workflow plumbing. It
does not establish equivalence for every tissue, protocol, QC threshold, or
real biological dataset.

## Advanced donor-aware validation

The advanced example uses four paired donors and deliberately induces an
interferon response in T cells:

```text
bl run singlecell/examples/advanced_analysis.bl
python singlecell/examples/validation/advanced_reference.py
Rscript singlecell/examples/validation/advanced_reference.R
```

Both reference scripts independently aggregate raw counts, calculate log2 CPM,
run a paired test, and compare the resulting effects and IFIT1 p-value with
`singlecell-results/paired-de.csv`. These scripts validate the transparent
paired exploratory model. Publication analyses should additionally run the
exported raw pseudobulk counts through DESeq2 or edgeR with an explicit
`~ donor + condition` design.

## Standalone SCTransform conformance

`run_sctransform_validation.py` runs the original R package and BioLang in two
separate processes against the same count matrix. The R package is a
validation-only installation and is not imported, linked, or needed by
BioLang. The driver records model parameters, residual variance, the ranked
feature list, a residual probe drawn from the top 3,000 features, transform and
whole-process time, and peak host working set.

Build BioLang in release mode, install `sctransform` separately in R, and run a
small deterministic check. The source option exposes the checkout's packages
without installing them globally:

```powershell
python validation/single-cell/run_sctransform_validation.py `
  --mode synthetic `
  --biolang-source ../biolang `
  --executable ../biolang/target/release/bl.exe `
  --r-library ../biolang/.validation-r-library `
  --output validation-results/sctransform-synthetic-current
```

For a real 10x MEX matrix:

```powershell
Rscript validation/single-cell/prepare_hbc_sctransform_fixture.R `
  ctrl_raw stim_raw validation-results/hbc-sctransform-input
python validation/single-cell/run_sctransform_validation.py `
  --mode tenx `
  --input validation-results/hbc-sctransform-input/ctrl `
  --biolang-source ../biolang `
  --executable ../biolang/target/release/bl.exe `
  --r-library ../biolang/.validation-r-library `
  --output validation-results/sctransform-hbc-ctrl-current
```

The output directory must be new so observations from different builds cannot
be mixed. `comparison.json` contains every numeric and resource metric plus
independent pass/fail gates. CPU is the reproducibility default; pass
`--gpu auto` or `--gpu on` only when intentionally validating that backend.
The calibrated oracle requires `glmGamPoi_offset` by default and stops if R
silently falls back to `nb_offset`. Use `--oracle-method nb_offset` only for a
separately named comparison of that backend.

To validate the residual assay from an HBC-style
`SCTransform(vars.to.regress = "mitoRatio")` call, pass a sixth `true` argument
to `sctransform_oracle.R`. The oracle then applies Seurat's post-VST
`ScaleData` regression before exporting the stratified residual probe. Compare
that fresh directory with a provider run made using `--regress-mito`; this
keeps core-VST and wrapper-level claims distinct. If mitochondrial metadata was
calculated before gene filtering, use the fixture's barcode-aligned
`cell-covariates.csv` with `--mito-fraction-file` instead of recomputing it.

Run the external provider with both `--regress-mito` and
`--seurat-conserve-memory` for that same HBC wrapper contract. The latter keeps
Seurat's distinct residual-variance ranking and returned-assay clips; omitting
it intentionally validates direct `sctransform::vst` behavior instead.
For a two-pass integration export, also pass each sample's original top-feature
list through `--regress-features-file`; this retains its regressed SCT rows
while newly requested union features follow `PrepSCTIntegration`'s model-only
residual path.

## Full HBC biological validation

`hbc_biological_external_sct.bl` exercises the process boundary used by the
optional GPL SCTransform executable, then runs the complete downstream analysis
in BioLang. `compare_hbc_biology.py` checks broad PBMC identities and the
control-to-stimulated response; `compare_hbc_results.py` retains the stricter
numeric cluster, feature, neighbour, UMAP, and marker comparisons.

To determine where partition drift enters, use
`prepare_hbc_seurat_pcs.py`, `hbc_cluster_seurat_pcs.bl`, and
`compare_hbc_cluster_probe.py`. That diagnostic fixes Seurat's integrated PCs
and tests only BioLang's SNN/Louvain stages. It is an isolation experiment, not
a mixed-engine analysis workflow.

`hbc_graph_seurat.R`, `hbc_graph_biolang.bl`, and `compare_hbc_graphs.py`
separate kNN identity, SNN identity/weights, and Louvain partitioning.
`hbc_louvain_seurat_graph.bl` tests BioLang clustering over the exact Seurat
SNN. `hbc_resolution_sweep.bl` and `compare_hbc_resolution_sweep.py` provide a
clearly labelled diagnostic sweep over one fixed embedding and graph; an
oracle-selected maximum from that sweep must not be presented as independent
validation or used to change the production default.

For numeric stage boundaries, `prepare_hbc_cca_inputs.py` creates byte-identical
3,000-feature residual matrices. `hbc_cca_seurat.R` and `hbc_cca_biolang.bl`
compare CCA coordinates, anchor identities, scores, and filter features;
`hbc_weight_pca_*` isolates the PCA used for anchor weights; and
`hbc_integrated_matrix_*` invokes the Seurat compiled weight/correction oracle
and BioLang over identical scored anchors. The Python comparators report both
correlation and scale-sensitive slope, RMSE, absolute-error, and neighbour
identity metrics. These scripts are validation tools and are not runtime
dependencies or examples of mixing engines in one scientific analysis.

`prepare_hbc_anchor_hybrids.py` creates validation-only pair/score hybrids and
`compare_hbc_anchor_sensitivity.py` reports how pair identity and score drift
separately affect the downstream partition. `hbc_anchor_swap_*` can optionally
export and consume corrected-matrix and PCA artifacts, allowing integrated
PCA, graph, and clustering to be tested without changing upstream inputs.
`hbc_cca_biolang.bl` records explicit `cca_sweeps` and `cca_oversample` values;
these are convergence diagnostics, not parameters to tune against Seurat
cluster labels.

For full-precision boundary replay, `hbc_cca_seurat.R` also writes `BLMATF64`
embeddings and projections. `hbc_fixed_embedding_anchors_biolang.bl` holds
those coordinates fixed while BioLang performs candidate search, filtering,
and scoring. `hbc_embedding_candidates_*` and
`hbc_embedding_neighbours_seurat.R` isolate candidate and `k.score` neighbour
ranks; `compare_hbc_embedding_candidates.py` and
`compare_hbc_anchor_neighbours.py` require exact discrete agreement. Use
`prepare_irlba_start.R` plus `BIOLANG_PCA_SOLVER=lanczos` only for a disclosed
strict-replay run. The native solver remains the untuned independent result.

The current fresh reference has 22 clusters. Exact supplied numeric boundaries
produce 22 BioLang clusters and ARI 1.0; native CCA produces 20 clusters and
ARI 0.80802 through the same downstream path. See
`HBC_REVALIDATION_2026-08-13.md` before quoting either number.

`strict_external_provider_smoke.bl` checks the public package call using
`compatibility: "external"`. `hbc_external_provider_biolang.bl` is the real-
data counterpart: it launches the separately installed provider automatically,
then writes anchors, candidate anchors, labels, stage timings, backend names,
and provider-version provenance. The 2026-08-14 run matched every anchor and
the complete 22-cluster reference partition (ARI 1.0 across 29,629 cells).
This is an optional GPL process boundary, not a claim that native CCA is exact.

`compare_sctransform_sampling.R` reconstructs the public density weights and
ordered fit-gene sample for diagnosis. `benchmark_pca.bl` times PCA separately
from integration and accepts `BIOLANG_PCA_INPUT`; optional
`BIOLANG_PCA_TRACE`, `BIOLANG_PCA_TRACE_SHIFT`, and
`BIOLANG_PCA_MAX_SWEEPS` settings expose convergence without changing the HBC
workflow. These are diagnostic controls, not parameters to tune against Seurat
cluster labels.

The latest measured real-data results, limitations, resource figures, and the
licensing boundary are recorded in [HBC_SEURAT_VALIDATION.md](HBC_SEURAT_VALIDATION.md).
