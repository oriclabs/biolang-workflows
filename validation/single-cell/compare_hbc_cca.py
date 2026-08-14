#!/usr/bin/env python3
"""Compare Seurat and BioLang CCA embeddings, filter genes, and anchor identities."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import orthogonal_procrustes, subspace_angles
from scipy.optimize import linear_sum_assignment


def embedding(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    columns = [name for name in rows[0] if name != "cell"]
    return np.asarray([[float(row[name]) for name in columns] for row in rows])


def pairs(path: Path) -> set[tuple[int, int]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (int(row["left"]), int(row["right"]))
            for row in csv.DictReader(handle)
        }


def scored_pairs(path: Path) -> dict[tuple[int, int], tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "score" not in rows[0] or "raw_score" not in rows[0]:
        return {}
    return {
        (int(row["left"]), int(row["right"])): (float(row["score"]), float(row["raw_score"]))
        for row in rows
    }


def genes(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["gene"] for row in csv.DictReader(handle)}


def overlap(left: set, right: set) -> dict[str, float | int]:
    common = len(left & right)
    return {
        "left": len(left),
        "right": len(right),
        "intersection": common,
        "recall_min_set": common / min(len(left), len(right)),
        "jaccard": common / len(left | right),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat", type=Path)
    parser.add_argument("biolang", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    seurat = np.vstack(
        [embedding(args.seurat / "left-embedding.csv"), embedding(args.seurat / "right-embedding.csv")]
    )
    biolang = np.vstack(
        [embedding(args.biolang / "left-embedding.csv"), embedding(args.biolang / "right-embedding.csv")]
    )
    correlations = np.corrcoef(seurat.T, biolang.T)[:30, 30:]
    rows, cols = linear_sum_assignment(-np.abs(correlations))
    assigned = np.abs(correlations[rows, cols])
    rotation, _ = orthogonal_procrustes(biolang, seurat)
    aligned = biolang @ rotation
    row_cosines = np.sum(aligned * seurat, axis=1) / np.maximum(
        np.linalg.norm(aligned, axis=1) * np.linalg.norm(seurat, axis=1),
        np.finfo(float).eps,
    )
    angles = np.degrees(subspace_angles(seurat, biolang))

    seurat_candidates = pairs(args.seurat / "candidate-anchors.csv")
    biolang_candidates = pairs(args.biolang / "candidate-anchors.csv")
    seurat_retained = pairs(args.seurat / "anchors.csv")
    biolang_retained = pairs(args.biolang / "anchors.csv")
    common_candidates = seurat_candidates & biolang_candidates
    common_retention_agreement = sum(
        ((pair in seurat_retained) == (pair in biolang_retained))
        for pair in common_candidates
    ) / len(common_candidates)
    seurat_scores = scored_pairs(args.seurat / "anchors.csv")
    biolang_scores = scored_pairs(args.biolang / "anchors.csv")
    common_scored = sorted(seurat_scores.keys() & biolang_scores.keys())
    seurat_normalized = np.asarray([seurat_scores[pair][0] for pair in common_scored])
    biolang_normalized = np.asarray([biolang_scores[pair][0] for pair in common_scored])
    seurat_raw = np.asarray([seurat_scores[pair][1] for pair in common_scored])
    biolang_raw = np.asarray([biolang_scores[pair][1] for pair in common_scored])

    result = {
        "cells": int(seurat.shape[0]),
        "dimensions": int(seurat.shape[1]),
        "same_index_abs_correlation": np.abs(np.diag(correlations)).tolist(),
        "assigned_abs_correlation": assigned.tolist(),
        "assigned_abs_correlation_mean": float(np.mean(assigned)),
        "assigned_abs_correlation_median": float(np.median(assigned)),
        "principal_angles_degrees": angles.tolist(),
        "principal_angle_median_degrees": float(np.median(angles)),
        "procrustes_relative_frobenius_error": float(
            np.linalg.norm(aligned - seurat) / np.linalg.norm(seurat)
        ),
        "procrustes_row_cosine_mean": float(np.mean(row_cosines)),
        "procrustes_row_cosine_median": float(np.median(row_cosines)),
        "candidate_anchors": overlap(seurat_candidates, biolang_candidates),
        "retained_anchors": overlap(seurat_retained, biolang_retained),
        "common_candidate_filter_decision_agreement": common_retention_agreement,
        "common_anchor_scores": {
            "anchors": len(common_scored),
            "raw_exact_fraction": float(np.mean(seurat_raw == biolang_raw)),
            "raw_pearson": float(np.corrcoef(seurat_raw, biolang_raw)[0, 1]),
            "raw_median_absolute_error": float(np.median(np.abs(seurat_raw - biolang_raw))),
            "normalized_pearson": float(np.corrcoef(seurat_normalized, biolang_normalized)[0, 1]),
            "normalized_median_absolute_error": float(
                np.median(np.abs(seurat_normalized - biolang_normalized))
            ),
            "normalized_p90_absolute_error": float(
                np.quantile(np.abs(seurat_normalized - biolang_normalized), 0.9)
            ),
        },
        "filter_features": overlap(
            genes(args.seurat / "filter-features.csv"),
            genes(args.biolang / "filter-features.csv"),
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
