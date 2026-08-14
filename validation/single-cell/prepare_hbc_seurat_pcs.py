#!/usr/bin/env python3
"""Convert Seurat's HBC PC artifact to BioLang's compact matrix interchange."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()

    frame = pd.read_csv(args.input_csv)
    seurat_columns = [f"PC_{index}" for index in range(1, 41)]
    biolang_columns = [f"pc_{index}" for index in range(1, 41)]
    matrix_columns = [f"V{index}" for index in range(1, 41)]
    pc_columns = next(
        columns
        for columns in (seurat_columns, biolang_columns, matrix_columns)
        if all(column in frame for column in columns)
    )
    missing = [column for column in pc_columns if column not in frame]
    if missing:
        raise RuntimeError(f"Seurat PC artifact is missing columns: {missing}")
    metadata = frame
    if not all(column in frame for column in ["sample", "barcode"]):
        if args.metadata is None:
            raise RuntimeError("PC artifact has no sample/barcode columns; pass --metadata")
        metadata = pd.read_csv(args.metadata)
        if len(metadata) != len(frame) or not all(
            column in metadata for column in ["sample", "barcode"]
        ):
            raise RuntimeError("metadata rows or columns do not match the PC artifact")
    values = frame[pc_columns].to_numpy(dtype="<f8", copy=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "matrix.f64").open("wb") as handle:
        handle.write(b"BLMATF64")
        handle.write(struct.pack("<QQ", values.shape[0], values.shape[1]))
        handle.write(values.tobytes(order="C"))
    metadata[["sample", "barcode"]].to_csv(args.output_dir / "cells.csv", index=False)
    print(f"wrote {values.shape[0]} x {values.shape[1]} PCs")


if __name__ == "__main__":
    main()
