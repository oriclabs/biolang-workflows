#!/usr/bin/env python3
"""Compare a BioLang resolution sweep with an independent Seurat partition."""

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
    parser.add_argument("sweep_dir", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    reference = pd.read_csv(args.seurat_cells)[["sample", "barcode", "cluster"]]
    observed = pd.read_csv(args.sweep_dir / "cells.csv")
    joined = reference.merge(observed, on=["sample", "barcode"], validate="one_to_one")
    if len(joined) != len(reference) or len(joined) != len(observed):
        raise RuntimeError("cell identity mismatch")

    rows = []
    for column in (name for name in observed.columns if name.startswith("r")):
        contingency = pd.crosstab(joined[column], joined.cluster)
        left, right = linear_sum_assignment(-contingency.to_numpy())
        matched = int(contingency.to_numpy()[left, right].sum())
        rows.append({
            "resolution": int(column[1:]) / 100,
            "clusters": int(joined[column].nunique()),
            "adjusted_rand_index": float(adjusted_rand_score(joined.cluster, joined[column])),
            "adjusted_mutual_information": float(adjusted_mutual_info_score(joined.cluster, joined[column])),
            "one_to_one_mapped_accuracy": matched / len(joined),
        })
    rows.sort(key=lambda row: row["resolution"])
    result = {
        "scope": "diagnostic sweep on one fixed BioLang PC embedding and SNN graph",
        "selection_warning": "Do not use the oracle-selected resolution as an independent validation result.",
        "results": rows,
        "maximum_ari": max(rows, key=lambda row: row["adjusted_rand_index"]),
    }
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
