#!/usr/bin/env python3
"""Compare marker engines while holding the measured Seurat clusters fixed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_markers", type=Path)
    parser.add_argument("biolang_markers", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    seurat = pd.read_csv(args.seurat_markers)
    biolang = pd.read_csv(args.biolang_markers)
    for frame in (seurat, biolang):
        frame["cluster"] = frame["cluster"].astype(str)
        frame["gene"] = frame["gene"].astype(str)

    seurat_pairs = set(zip(seurat["cluster"], seurat["gene"]))
    biolang_pairs = set(zip(biolang["cluster"], biolang["gene"]))
    intersection = seurat_pairs & biolang_pairs
    union = seurat_pairs | biolang_pairs

    joined = seurat.merge(
        biolang,
        on=["cluster", "gene"],
        how="inner",
        suffixes=("_seurat", "_biolang"),
        validate="one_to_one",
    )
    top50_by_cluster: dict[str, int] = {}
    for cluster in sorted(set(seurat["cluster"]) | set(biolang["cluster"])):
        left = set(
            seurat[seurat["cluster"] == cluster]
            .sort_values(["avg_log2FC", "gene"], ascending=[False, True])
            .head(50)["gene"]
        )
        right = set(
            biolang[biolang["cluster"] == cluster]
            .sort_values(["avg_log2fc", "gene"], ascending=[False, True])
            .head(50)["gene"]
        )
        top50_by_cluster[cluster] = len(left & right)

    def correlation(left: str, right: str) -> float:
        return float(np.corrcoef(joined[left], joined[right])[0, 1])

    def rmse(left: str, right: str) -> float:
        return float(np.sqrt(np.mean((joined[left] - joined[right]) ** 2)))

    results = {
        "fixed_input": "Seurat HBC resolution-0.8 clusters",
        "seurat_positive_marker_rows": int(len(seurat)),
        "biolang_positive_marker_rows": int(len(biolang)),
        "marker_pair_intersection": len(intersection),
        "marker_pair_recall_min_set": len(intersection) / min(len(seurat_pairs), len(biolang_pairs)),
        "marker_pair_jaccard": len(intersection) / len(union),
        "top50_overlap_by_cluster": top50_by_cluster,
        "top50_overlap_total": sum(top50_by_cluster.values()),
        "top50_overlap_fraction": sum(top50_by_cluster.values()) / (50 * len(top50_by_cluster)),
        "joined_marker_pairs": int(len(joined)),
        "avg_log2fc_pearson": correlation("avg_log2FC", "avg_log2fc"),
        "avg_log2fc_rmse": rmse("avg_log2FC", "avg_log2fc"),
        "pct_1_pearson": correlation("pct.1", "pct_1"),
        "pct_1_rmse": rmse("pct.1", "pct_1"),
        "pct_2_pearson": correlation("pct.2", "pct_2"),
        "pct_2_rmse": rmse("pct.2", "pct_2"),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
