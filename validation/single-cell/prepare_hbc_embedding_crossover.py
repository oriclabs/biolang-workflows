#!/usr/bin/env python3
"""Convert CCA embedding CSVs into BioLang's compact row-major matrix format."""

from __future__ import annotations

import argparse
import csv
import struct
from pathlib import Path

import numpy as np


def read_embedding(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        keep = [index for index, name in enumerate(header) if name != "cell"]
        values = [[float(row[index]) for index in keep] for row in reader]
    matrix = np.asarray(values, dtype="<f8")
    if matrix.ndim != 2 or matrix.size == 0:
        raise ValueError(f"empty or invalid embedding: {path}")
    return matrix


def write_blmat(path: Path, matrix: np.ndarray) -> None:
    matrix = np.ascontiguousarray(matrix, dtype="<f8")
    with path.open("wb") as handle:
        handle.write(b"BLMATF64")
        handle.write(struct.pack("<QQ", matrix.shape[0], matrix.shape[1]))
        handle.write(matrix.tobytes(order="C"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("embedding_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--features", type=Path)
    parser.add_argument("--filter-features", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    left = read_embedding(args.embedding_dir / "left-embedding.csv")
    right = read_embedding(args.embedding_dir / "right-embedding.csv")
    if left.shape[1] != right.shape[1]:
        raise ValueError("left and right embedding dimensions differ")
    write_blmat(args.output_dir / "left.f64", left)
    write_blmat(args.output_dir / "right.f64", right)
    if bool(args.features) != bool(args.filter_features):
        raise ValueError("--features and --filter-features must be supplied together")
    if args.features:
        with args.features.open(newline="", encoding="utf-8-sig") as handle:
            feature_rows = list(csv.DictReader(handle))
        by_gene = {row["gene"]: index for index, row in enumerate(feature_rows)}
        with args.filter_features.open(newline="", encoding="utf-8-sig") as handle:
            selected = list(csv.DictReader(handle))
        with (args.output_dir / "filter-features.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=["rank", "gene", "feature_index"])
            writer.writeheader()
            for rank, row in enumerate(selected, 1):
                writer.writerow(
                    {"rank": rank, "gene": row["gene"], "feature_index": by_gene[row["gene"]]}
                )
    print(f"left={left.shape} right={right.shape} output={args.output_dir}")


if __name__ == "__main__":
    main()
