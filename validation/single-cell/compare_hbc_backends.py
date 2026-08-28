#!/usr/bin/env python3
"""Compare two BioLang HBC runs without treating backend drift as harmless."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--min-cluster-ari", type=float, default=0.999999)
    parser.add_argument("--min-pc-correlation", type=float, default=0.999999)
    parser.add_argument("--max-pc-relative-rmse", type=float, default=1e-6)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def keyed(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["key"] = frame["sample"].astype(str) + "::" + frame["barcode"].astype(str)
    if frame["key"].duplicated().any():
        raise RuntimeError("sample-plus-barcode keys are not unique")
    return frame.set_index("key").sort_index()


def main() -> None:
    args = parse_args()
    left_cells = keyed(pd.read_csv(args.left / "cells.csv"))
    right_cells = keyed(pd.read_csv(args.right / "cells.csv"))
    if not left_cells.index.equals(right_cells.index):
        raise RuntimeError("backend runs do not contain the same ordered cell keys")

    left_labels = left_cells["cluster"].astype(str)
    right_labels = right_cells["cluster"].astype(str)
    cluster_ari = float(adjusted_rand_score(left_labels, right_labels))
    exact_cluster_labels = float(np.mean(left_labels == right_labels))

    left_features = pd.read_csv(args.left / "features.csv")["gene"].astype(str).tolist()
    right_features = pd.read_csv(args.right / "features.csv")["gene"].astype(str).tolist()

    left_pcs = keyed(pd.read_csv(args.left / "pcs.csv"))
    right_pcs = keyed(pd.read_csv(args.right / "pcs.csv"))
    if not left_pcs.index.equals(right_pcs.index):
        raise RuntimeError("backend PC manifests do not contain the same cell keys")
    pc_columns = sorted(
        (name for name in left_pcs.columns if name.startswith("pc_")),
        key=lambda name: int(name.removeprefix("pc_")),
    )
    right_pc_columns = sorted(
        (name for name in right_pcs.columns if name.startswith("pc_")),
        key=lambda name: int(name.removeprefix("pc_")),
    )
    if pc_columns != right_pc_columns:
        raise RuntimeError("backend PC columns differ")
    left_values = left_pcs[pc_columns].to_numpy(dtype=np.float64)
    right_values = right_pcs[pc_columns].to_numpy(dtype=np.float64)
    correlations = [
        float(np.corrcoef(left_values[:, index], right_values[:, index])[0, 1])
        for index in range(left_values.shape[1])
    ]
    pc_rmse = float(np.sqrt(np.mean((left_values - right_values) ** 2)))
    pc_scale = float(np.std(left_values))
    pc_relative_rmse = pc_rmse / pc_scale if pc_scale else pc_rmse

    left_markers = pd.read_csv(args.left / "markers.csv")
    right_markers = pd.read_csv(args.right / "markers.csv")
    left_pairs = set(zip(left_markers["cluster"].astype(str), left_markers["gene"].astype(str)))
    right_pairs = set(zip(right_markers["cluster"].astype(str), right_markers["gene"].astype(str)))

    results = {
        "left": str(args.left.resolve()),
        "right": str(args.right.resolve()),
        "cells": int(len(left_cells)),
        "cluster_ari": cluster_ari,
        "exact_numeric_cluster_label_fraction": exact_cluster_labels,
        "ordered_features_exact": left_features == right_features,
        "feature_count": len(left_features),
        "minimum_component_correlation": min(correlations),
        "pc_rmse": pc_rmse,
        "pc_relative_rmse": pc_relative_rmse,
        "marker_pairs_exact": left_pairs == right_pairs,
        "marker_pair_jaccard": len(left_pairs & right_pairs) / len(left_pairs | right_pairs),
        "file_sha256": {
            name: {
                "left": sha256(args.left / name),
                "right": sha256(args.right / name),
                "exact": sha256(args.left / name) == sha256(args.right / name),
            }
            for name in ("features.csv", "markers.csv", "summary.csv")
        },
    }
    gates = {
        "cluster_ari": cluster_ari >= args.min_cluster_ari,
        "ordered_features_exact": left_features == right_features,
        "pc_component_correlation": min(correlations) >= args.min_pc_correlation,
        "pc_relative_rmse": pc_relative_rmse <= args.max_pc_relative_rmse,
        "marker_pairs_exact": left_pairs == right_pairs,
    }
    results["gates"] = gates
    results["passed"] = all(gates.values())
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))
    if not results["passed"]:
        failed = ", ".join(name for name, passed in gates.items() if not passed)
        raise SystemExit(f"backend parity gates failed: {failed}")


if __name__ == "__main__":
    main()
