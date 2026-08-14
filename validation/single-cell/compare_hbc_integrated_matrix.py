#!/usr/bin/env python3
"""Compare matching samples of Seurat and BioLang integrated SCT matrices."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def read_matrix(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.asarray([[float(value) for value in row.values()] for row in rows])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat", type=Path)
    parser.add_argument("biolang", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    seurat = read_matrix(args.seurat)
    biolang = read_matrix(args.biolang)
    if seurat.shape != biolang.shape:
        raise ValueError(f"shape mismatch: {seurat.shape} != {biolang.shape}")
    expected = seurat.ravel()
    observed = biolang.ravel()
    slope, intercept = np.polyfit(expected, observed, 1)
    errors = observed - expected
    relative = np.abs(errors) / np.maximum(np.abs(expected), 1e-12)
    material = np.abs(expected) > 1
    row_cosines = np.sum(seurat * biolang, axis=1) / np.maximum(
        np.linalg.norm(seurat, axis=1) * np.linalg.norm(biolang, axis=1),
        np.finfo(float).eps,
    )
    result = {
        "sample_cells": int(seurat.shape[0]),
        "sample_features": int(seurat.shape[1]),
        "observations": int(expected.size),
        "pearson": float(np.corrcoef(expected, observed)[0, 1]),
        "regression_slope": float(slope),
        "regression_intercept": float(intercept),
        "rmse": float(np.sqrt(np.mean(errors * errors))),
        "rmse_fraction_of_seurat_sd": float(np.sqrt(np.mean(errors * errors)) / np.std(expected)),
        "median_absolute_error": float(np.median(np.abs(errors))),
        "p90_absolute_error": float(np.quantile(np.abs(errors), 0.9)),
        "median_relative_error": float(np.median(relative)),
        "material_observations": int(np.sum(material)),
        "material_median_relative_error": float(np.median(relative[material])),
        "material_p90_relative_error": float(np.quantile(relative[material], 0.9)),
        "row_cosine_mean": float(np.mean(row_cosines)),
        "row_cosine_median": float(np.median(row_cosines)),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
