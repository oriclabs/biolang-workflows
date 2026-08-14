#!/usr/bin/env python3
"""Quantify HBC downstream sensitivity to anchor pairs and anchor scores."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score


def anchors(path: Path) -> dict[tuple[int, int], tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        result = {}
        for row in csv.DictReader(handle):
            pair = (int(row["left"]), int(row["right"]))
            result[pair] = (
                float(row.get("score", "nan")),
                float(row.get("raw_score", "nan")),
            )
        return result


def cluster_metrics(reference_path: Path, observed_path: Path) -> dict[str, float | int]:
    reference = pd.read_csv(reference_path)
    observed = pd.read_csv(observed_path)
    for frame in (reference, observed):
        frame["key"] = frame["sample"].astype(str) + "::" + frame["barcode"].astype(str)
        frame["cluster"] = frame["cluster"].astype(str)
    joined = reference[["key", "cluster"]].merge(
        observed[["key", "cluster"]], on="key", suffixes=("_reference", "_observed")
    )
    table = pd.crosstab(joined["cluster_observed"], joined["cluster_reference"])
    rows, columns = linear_sum_assignment(-table.to_numpy())
    return {
        "clusters": int(joined["cluster_observed"].nunique()),
        "ari": float(adjusted_rand_score(joined["cluster_reference"], joined["cluster_observed"])),
        "ami": float(adjusted_mutual_info_score(joined["cluster_reference"], joined["cluster_observed"])),
        "mapped_accuracy": float(table.to_numpy()[rows, columns].sum() / len(joined)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_anchor_dir", type=Path)
    parser.add_argument("biolang_anchor_dir", type=Path)
    parser.add_argument("seurat_cells", type=Path)
    parser.add_argument("exact_biolang_cells", type=Path)
    parser.add_argument("hybrid_results", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    seurat = anchors(args.seurat_anchor_dir / "anchors.csv")
    biolang = anchors(args.biolang_anchor_dir / "anchors.csv")
    seurat_candidates = set(anchors(args.seurat_anchor_dir / "candidate-anchors.csv"))
    biolang_candidates = set(anchors(args.biolang_anchor_dir / "candidate-anchors.csv"))
    common_candidates = seurat_candidates & biolang_candidates
    common = set(seurat) & set(biolang)
    score_errors = np.asarray(
        [abs(seurat[pair][0] - biolang[pair][0]) for pair in common]
    )
    score_changed = int(np.count_nonzero(score_errors > 1e-12))
    raw_changed = sum(seurat[pair][1] != biolang[pair][1] for pair in common)
    common_filter_disagreements = sum(
        ((pair in seurat) != (pair in biolang)) for pair in common_candidates
    )

    variants = {
        "both_seurat": args.exact_biolang_cells,
        "seurat_pairs_biolang_scores": args.hybrid_results
        / "hbc-review-20260813-anchor-hybrid-seurat-pairs-biolang-scores"
        / "cells.csv",
        "biolang_pairs_seurat_scores": args.hybrid_results
        / "hbc-review-20260813-anchor-hybrid-biolang-pairs-seurat-scores"
        / "cells.csv",
        "both_biolang": args.hybrid_results
        / "hbc-review-20260813-anchor-boundary-biolang-current"
        / "cells.csv",
    }
    result = {
        "candidate_pairs": {
            "seurat": len(seurat_candidates),
            "biolang": len(biolang_candidates),
            "common": len(common_candidates),
            "common_filter_decision_disagreements": common_filter_disagreements,
        },
        "retained_pairs": {
            "seurat": len(seurat),
            "biolang": len(biolang),
            "common": len(common),
            "seurat_only": len(set(seurat) - common),
            "biolang_only": len(set(biolang) - common),
        },
        "common_anchor_scores": {
            "materially_changed_normalized": score_changed,
            "materially_changed_normalized_fraction": score_changed / len(common),
            "normalized_median_absolute_error": float(np.median(score_errors)),
            "normalized_p90_absolute_error": float(np.quantile(score_errors, 0.9)),
            "changed_raw": raw_changed,
            "changed_raw_fraction": raw_changed / len(common),
        },
        "downstream_partition": {
            name: cluster_metrics(args.seurat_cells, path)
            for name, path in variants.items()
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
