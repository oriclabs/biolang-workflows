#!/usr/bin/env python3
"""Compare clustering after feeding identical Seurat PCs to BioLang."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_cells", type=Path)
    parser.add_argument("biolang_cells", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    seurat = pd.read_csv(args.seurat_cells)
    biolang = pd.read_csv(args.biolang_cells)
    for frame in (seurat, biolang):
        frame["key"] = frame["sample"].astype(str) + "::" + frame["barcode"].astype(str)
        frame["cluster"] = frame["cluster"].astype(str)
    joined = seurat[["key", "cluster"]].merge(
        biolang[["key", "cluster"]],
        on="key",
        suffixes=("_seurat", "_biolang"),
        validate="one_to_one",
    )
    if len(joined) != len(seurat) or len(joined) != len(biolang):
        raise RuntimeError("cell identity mismatch in cluster isolation probe")

    contingency = pd.crosstab(joined["cluster_biolang"], joined["cluster_seurat"])
    rows, columns = linear_sum_assignment(-contingency.to_numpy())
    mapped_correct = int(contingency.to_numpy()[rows, columns].sum())
    result = {
        "joined_cells": len(joined),
        "fixed_input": "Seurat integrated PCs 1:40",
        "seurat_clusters": int(joined["cluster_seurat"].nunique()),
        "biolang_clusters": int(joined["cluster_biolang"].nunique()),
        "adjusted_rand_index": float(
            adjusted_rand_score(joined["cluster_seurat"], joined["cluster_biolang"])
        ),
        "adjusted_mutual_information": float(
            adjusted_mutual_info_score(
                joined["cluster_seurat"], joined["cluster_biolang"]
            )
        ),
        "one_to_one_mapped_accuracy": mapped_correct / len(joined),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
