#!/usr/bin/env python3
"""Rotation-invariant comparison of two cell x component PCA CSV files."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import re
import struct

import numpy as np
from scipy.linalg import orthogonal_procrustes, subspace_angles
from scipy.optimize import linear_sum_assignment


def matrix(path: Path) -> np.ndarray:
    if path.suffix.casefold() == ".f64":
        with path.open("rb") as handle:
            header = handle.read(24)
            if len(header) != 24 or header[:8] != b"BLMATF64":
                raise ValueError(f"invalid BLMATF64 matrix: {path}")
            rows, columns = struct.unpack("<QQ", header[8:])
            values = np.fromfile(handle, dtype="<f8", count=rows * columns)
        if values.size != rows * columns:
            raise ValueError(f"truncated BLMATF64 matrix: {path}")
        return values.reshape((rows, columns))

    opener = gzip.open if path.suffix.casefold() == ".gz" else Path.open
    with opener(path, mode="rt", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    columns = [
        name
        for name in rows[0]
        if re.search(r"(?:^|_)pc_?\d+$", name, re.I)
        or re.fullmatch(r"V\d+", name, re.I)
    ]
    if not columns:
        raise ValueError(f"no PC columns found in {path}")
    columns.sort(key=lambda name: int(re.search(r"(\d+)$", name).group(1)))
    return np.asarray([[float(row[name]) for name in columns] for row in rows])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat", type=Path)
    parser.add_argument("biolang", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    seurat = matrix(args.seurat)
    biolang = matrix(args.biolang)
    dimensions = min(seurat.shape[1], biolang.shape[1])
    seurat = seurat[:, :dimensions]
    biolang = biolang[:, :dimensions]
    correlations = np.corrcoef(seurat.T, biolang.T)[: seurat.shape[1], seurat.shape[1] :]
    rows, columns = linear_sum_assignment(-np.abs(correlations))
    assigned = np.abs(correlations[rows, columns])
    rotation, _ = orthogonal_procrustes(biolang, seurat)
    aligned = biolang @ rotation
    row_cosines = np.sum(aligned * seurat, axis=1) / np.maximum(
        np.linalg.norm(aligned, axis=1) * np.linalg.norm(seurat, axis=1),
        np.finfo(float).eps,
    )
    angles = np.degrees(subspace_angles(seurat, biolang))
    result = {
        "cells": int(seurat.shape[0]),
        "dimensions": int(seurat.shape[1]),
        "same_index_abs_correlation": np.abs(np.diag(correlations)).tolist(),
        "assigned_abs_correlation": assigned.tolist(),
        "assigned_abs_correlation_min": float(np.min(assigned)),
        "assigned_abs_correlation_median": float(np.median(assigned)),
        "principal_angles_degrees": angles.tolist(),
        "principal_angle_max_degrees": float(np.max(angles)),
        "principal_angle_median_degrees": float(np.median(angles)),
        "procrustes_relative_frobenius_error": float(
            np.linalg.norm(aligned - seurat) / np.linalg.norm(seurat)
        ),
        "procrustes_row_cosine_mean": float(np.mean(row_cosines)),
        "procrustes_row_cosine_median": float(np.median(row_cosines)),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
