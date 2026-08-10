# External validation boundary

Validation programs execute BioLang, R/Seurat, Scanpy, and other references as
separate processes over common input files. They may compare numeric output but
must not link external reference implementations into BioLang or translate
their source into the MIT implementation.

External environments are optional and must never be required to build or use
BioLang. Preserve each tool's own licence notices and record exact package and
runtime versions in generated manifests.

Single-cell validation, including the standalone SCTransform black-box runner,
is under `single-cell/`. Generated output belongs under
`validation-results/`.

