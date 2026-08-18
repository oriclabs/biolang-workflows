# BioLang Workflows

Applied, reproducible analyses written in BioLang. This repository contains
courses, notebooks, real-data workflows, benchmarks, and independent reference
validations. The BioLang compiler, runtime, scientific packages, and core
algorithms remain in the [BioLang repository](https://github.com/oriclabs/biolang).

## Repository boundary

This repository depends on BioLang; BioLang never depends on this repository.

| This repository | Main BioLang repository |
|---|---|
| HBC and other complete analyses | Compiler, runtime, and CLI |
| Practical learning material | `singlecell` and other package APIs |
| Real-data notebooks | SCTransform Rust implementation |
| External black-box validation | Unit tests and small conformance fixtures |
| Generated benchmark reports | Minimal package examples |

Generated sites, downloaded datasets, R/Python environments, and validation
results are deliberately not committed.

## Contents

- `courses/`: the HBC course adaptation and its measured validation companion.
- `examples/`: small runnable BioLang examples organized by domain.
- `workflows/single-cell/`: complete single-cell scripts and notebooks.
- `validation/single-cell/`: Seurat, Scanpy, and SCTransform black-box comparisons.
- `datasets/`: download manifests, checksums, and provenance records only.
- `benchmarks/`: benchmark definitions; generated measurements are ignored.

The language reference remains sourced and maintained in the main BioLang
repository. Generated HTML and PDFs remain ignored and never become
implementation source.

## Use with a released BioLang installation

Install BioLang 1.4.0 or newer and ensure `bl` is on `PATH`. Workflows importing
`singlecell` also require that package to be installed in BioLang's package
directory. Until a package registry is available, install it from a BioLang
source checkout:

```powershell
bl install ..\biolang\packages\singlecell
bl run workflows\single-cell\seurat_standard_workflow.bl
```

The install command copies the package to the user's BioLang package directory.
To avoid changing the user installation during development, point `BIOLANG_PATH`
at the source checkout instead:

```powershell
$env:BIOLANG_PATH = (Resolve-Path ..\biolang\packages)
& ..\biolang\target\release\bl.exe run workflows\single-cell\seurat_standard_workflow.bl
```

The wrappers under `scripts/` provide the second form without permanently
changing the calling shell. They default `BIOLANG_GPU` to `off` for
cross-machine reproducibility, but preserve an explicit user choice:

```powershell
$env:BIOLANG_GPU = "auto" # opt in when backend-specific results are acceptable
scripts/run-with-local-biolang.ps1 ..\biolang run examples\quickstart.bl
```

## Validate the repository

```powershell
python scripts/check_repository.py
```

Heavy R/Seurat and Scanpy comparisons are opt-in. They run as separate
processes and are not dependencies of BioLang or ordinary workflow execution.
See `validation/README.md` for the boundary and commands.

## Migration status

The first migration is intentionally non-destructive. Source material remains
in the BioLang checkout until this repository builds and validates independently.
See `MIGRATION.md` for the exact source-to-destination map and removal criteria.
