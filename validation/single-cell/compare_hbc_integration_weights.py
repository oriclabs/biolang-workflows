#!/usr/bin/env python3
"""Compare sparse Seurat and BioLang integration-weight diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat", type=Path)
    parser.add_argument("biolang", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    seurat = pd.read_csv(args.seurat)
    biolang = pd.read_csv(args.biolang)
    keys = ["query_cell", "anchor_index"]
    joined = seurat.merge(biolang, on=keys, how="outer", suffixes=("_seurat", "_biolang"))
    joined[["weight_seurat", "weight_biolang"]] = joined[
        ["weight_seurat", "weight_biolang"]
    ].fillna(0.0)
    common = joined[(joined.weight_seurat > 0) & (joined.weight_biolang > 0)]
    support_intersection = int(len(common))
    support_union = int(len(joined))
    differences = joined.weight_biolang.to_numpy() - joined.weight_seurat.to_numpy()
    per_cell = joined.groupby("query_cell").apply(
        lambda frame: pd.Series({
            "support_intersection": int(np.sum((frame.weight_seurat > 0) & (frame.weight_biolang > 0))),
            "support_union": int(np.sum((frame.weight_seurat > 0) | (frame.weight_biolang > 0))),
            "l1": float(np.sum(np.abs(frame.weight_biolang - frame.weight_seurat))),
        }),
        include_groups=False,
    )
    result = {
        "sample_cells": int(joined.query_cell.nunique()),
        "seurat_nonzero_weights": int(len(seurat)),
        "biolang_nonzero_weights": int(len(biolang)),
        "support_intersection": support_intersection,
        "support_jaccard": support_intersection / support_union,
        "weight_pearson_with_zeros": float(np.corrcoef(joined.weight_seurat, joined.weight_biolang)[0, 1]),
        "weight_rmse_with_zeros": float(np.sqrt(np.mean(differences * differences))),
        "common_weight_pearson": float(np.corrcoef(common.weight_seurat, common.weight_biolang)[0, 1]),
        "common_weight_median_absolute_error": float(np.median(np.abs(common.weight_biolang - common.weight_seurat))),
        "per_cell_support_jaccard_median": float(np.median(per_cell.support_intersection / per_cell.support_union)),
        "per_cell_support_jaccard_p10": float(np.quantile(per_cell.support_intersection / per_cell.support_union, 0.1)),
        "per_cell_weight_l1_median": float(np.median(per_cell.l1)),
        "per_cell_weight_l1_p90": float(np.quantile(per_cell.l1, 0.9)),
        "per_cell_weight_l1_max": float(np.max(per_cell.l1)),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
