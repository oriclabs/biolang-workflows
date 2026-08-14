#!/usr/bin/env python3
"""Replay Seurat's integration kernel from exported Annoy neighbours."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial.distance import cdist


def csv_matrix(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.asarray([[float(value) for value in row.values()] for row in rows])


def blmat(path: Path) -> np.memmap:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"BLMATF64":
        raise ValueError(f"invalid BLMATF64 matrix: {path}")
    rows, columns = struct.unpack("<QQ", header[8:])
    return np.memmap(path, dtype="<f8", mode="r", offset=24, shape=(rows, columns))


def anchors(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.asarray([int(row["left"]) for row in rows]),
        np.asarray([int(row["right"]) for row in rows]),
        np.asarray([float(row["score"]) for row in rows]),
    )


def pca(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    columns = [f"PC_{index}" for index in range(1, 31)]
    return np.asarray([[float(row[column]) for column in columns] for row in rows])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("cca_dir", type=Path)
    parser.add_argument("integration_dir", type=Path)
    parser.add_argument("neighbour_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    left = blmat(args.input_dir / "ctrl.f64")
    right = blmat(args.input_dir / "stim.f64")
    anchor_left, anchor_right, anchor_score = anchors(args.cca_dir / "anchors.csv")
    query_pca = pca(args.integration_dir / "query-weight-pca.csv")
    neighbour_indices = csv_matrix(args.neighbour_dir / "indices.csv").astype(int) - 1
    neighbour_distances = csv_matrix(args.neighbour_dir / "distances.csv")
    seurat = csv_matrix(args.integration_dir / "integrated-query-sample.csv")
    sample_cells = np.arange(512) * 28
    sample_features = np.arange(500) * 6

    unique_anchor_cells = np.asarray(list(dict.fromkeys(anchor_right.tolist())))
    anchors_by_cell: dict[int, list[int]] = {}
    for index, cell in enumerate(anchor_right):
        anchors_by_cell.setdefault(int(cell), []).append(index)

    replay = np.empty_like(seurat)
    for output_row, cell in enumerate(sample_cells):
        distances = neighbour_distances[output_row]
        similarities = 1 - distances / distances[-1] if distances[-1] > 1e-15 else np.ones(100)
        selected: list[int] = []
        weights: list[float] = []
        for neighbour_position, similarity in zip(neighbour_indices[output_row], similarities):
            anchor_cell = int(unique_anchor_cells[neighbour_position])
            for anchor_index in anchors_by_cell[anchor_cell]:
                if len(selected) >= 100:
                    break
                selected.append(anchor_index)
                weights.append(1 - np.exp(-similarity * anchor_score[anchor_index] / 4))
            if len(selected) >= 100:
                break
        selected_array = np.asarray(selected)
        weight_array = np.asarray(weights)
        differences = (
            left[np.ix_(anchor_left[selected_array], sample_features)]
            - right[np.ix_(anchor_right[selected_array], sample_features)]
        )
        replay[output_row] = right[cell, sample_features] + (
            weight_array[:, None] * differences
        ).sum(axis=0) / weight_array.sum()

    errors = replay - seurat
    exact_distances = cdist(query_pca[sample_cells], query_pca[unique_anchor_cells])
    exact_indices = np.argsort(exact_distances, axis=1)[:, :100]
    recalls = np.asarray([
        len(set(expected) & set(observed)) / 100
        for expected, observed in zip(neighbour_indices, exact_indices)
    ])
    result = {
        "sample_cells": 512,
        "sample_features": 500,
        "observations": int(seurat.size),
        "replay_pearson": float(np.corrcoef(seurat.ravel(), replay.ravel())[0, 1]),
        "replay_rmse": float(np.sqrt(np.mean(errors * errors))),
        "replay_max_absolute_error": float(np.max(np.abs(errors))),
        "seurat_annoy_vs_exact_neighbour_recall_mean": float(np.mean(recalls)),
        "seurat_annoy_vs_exact_neighbour_recall_median": float(np.median(recalls)),
        "seurat_annoy_vs_exact_neighbour_recall_min": float(np.min(recalls)),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
