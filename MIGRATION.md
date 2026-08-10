# Migration from the BioLang repository

The initial migration was copied from BioLang commit
`8d0db026f21ae178e700a1fbaeb8318efca5dbc3` on 2026-08-11. It is
non-destructive: no source content has yet been removed from BioLang.
The source checkout also contained pending validation-documentation,
comparison, and runner changes under `packages/singlecell/examples/validation`;
those workflow-facing changes were intentionally captured here and then
adapted to the external repository layout. No pending Rust/core source was
copied as an unrecorded algorithm change.

| BioLang source | Workflow destination |
|---|---|
| `books/practical-bioinformatics` | `books/practical-bioinformatics` |
| `books/biostatistics` | `books/biostatistics` |
| `books/msmb` | `books/msmb` |
| `books/single-cell-rna-seq` | `books/single-cell-rna-seq` |
| `books/hbc-scrnaseq` | `courses/hbc-scrnaseq` |
| `books/hbc-scrnaseq-validated` | `courses/hbc-scrnaseq-validated` |
| `examples` | `examples` |
| `packages/singlecell/examples/*` | `workflows/single-cell` |
| `packages/singlecell/examples/validation` | `validation/single-cell` |

The BioLang language reference remains in the main repository. Package source,
SCTransform Rust code, package tests, algorithm audits tied to core behavior,
and compact smoke examples also remain there.

## Criteria before deleting duplicated source

1. Every copied book builds from this repository.
2. Representative `.bl` and `.bln` workflows run against a released or
   explicitly selected BioLang checkout.
3. HBC and SCTransform black-box validation produces fresh reports here.
4. Links from the BioLang documentation point to stable locations here.
5. CI checks paths, large-file policy, and accidental generated artifacts.
6. The migration commit/tag is recorded in both repositories.

Only after these checks should a later BioLang change remove the duplicated
books and workflows and replace them with concise links.
