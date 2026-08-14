#!/usr/bin/env python3
"""Order Seurat's measured HBC labels on the filtered BioLang MEX cells."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_barcodes(path: Path) -> list[str]:
    return [line.strip().split("\t", 1)[0] for line in path.read_text().splitlines()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_cells", type=Path)
    parser.add_argument("ctrl_mex", type=Path)
    parser.add_argument("stim_mex", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    reference = pd.read_csv(args.seurat_cells)
    reference["key"] = reference["sample"].astype(str) + "::" + reference["barcode"].astype(str)
    if reference["key"].duplicated().any():
        raise RuntimeError("Seurat sample-plus-barcode keys are not unique")
    lookup = reference.set_index("key")["cluster"]

    rows: list[dict[str, object]] = []
    for sample, mex in (("ctrl", args.ctrl_mex), ("stim", args.stim_mex)):
        for barcode in read_barcodes(mex / "barcodes.tsv"):
            key = f"{sample}::{barcode}"
            if key not in lookup.index:
                raise RuntimeError(f"filtered MEX cell is absent from Seurat artifact: {key}")
            rows.append({"sample": sample, "barcode": barcode, "cluster": int(lookup[key])})

    output = pd.DataFrame(rows)
    if len(output) != len(reference) or set(output["sample"] + "::" + output["barcode"]) != set(reference["key"]):
        raise RuntimeError("Seurat and filtered-MEX cell sets differ")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_csv, index=False)
    print(f"wrote {len(output)} fixed Seurat cluster labels")


if __name__ == "__main__":
    main()
