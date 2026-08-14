#!/usr/bin/env python3
"""Select byte-identical HBC integration matrices from SCT provider artifacts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import struct

import numpy as np


def read_shape(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"BLMATF64":
        raise ValueError(f"invalid BLMATF64 matrix: {path}")
    return struct.unpack("<QQ", header[8:24])


def provider_columns(directory: Path) -> tuple[dict[str, int], list[dict[str, str]]]:
    with (directory / "genes.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row["gene"]: index for index, row in enumerate(rows)}, rows


def write_selected(
    provider: Path, features: list[str], output: Path
) -> list[dict[str, str]]:
    rows, columns = read_shape(provider / "matrix.f64")
    by_gene, metadata = provider_columns(provider)
    missing = [gene for gene in features if gene not in by_gene]
    if missing:
        raise ValueError(f"{provider} lacks {len(missing)} integration features: {missing[:5]}")
    selected = np.asarray([by_gene[gene] for gene in features], dtype=np.int64)
    matrix = np.memmap(
        provider / "matrix.f64",
        dtype="<f8",
        mode="r",
        offset=24,
        shape=(rows, columns),
        order="C",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        handle.write(b"BLMATF64")
        handle.write(struct.pack("<QQ", rows, len(features)))
        for start in range(0, rows, 256):
            block = np.asarray(matrix[start : start + 256, selected], dtype="<f8")
            handle.write(block.tobytes(order="C"))
    return [metadata[index] for index in selected]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ctrl_provider", type=Path)
    parser.add_argument("stim_provider", type=Path)
    parser.add_argument("integration_features", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.integration_features.open(newline="", encoding="utf-8") as handle:
        features = [row["gene"] for row in csv.DictReader(handle)]
    if len(features) != 3000 or len(set(features)) != len(features):
        raise ValueError("expected 3,000 unique integration features")

    ctrl_metadata = write_selected(args.ctrl_provider, features, args.output / "ctrl.f64")
    stim_metadata = write_selected(args.stim_provider, features, args.output / "stim.f64")
    with (args.output / "features.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", "gene", "ctrl_gene_index", "stim_gene_index"])
        writer.writeheader()
        for rank, (gene, ctrl, stim) in enumerate(zip(features, ctrl_metadata, stim_metadata), 1):
            writer.writerow(
                {
                    "rank": rank,
                    "gene": gene,
                    "ctrl_gene_index": ctrl["gene_index"],
                    "stim_gene_index": stim["gene_index"],
                }
            )
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
