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
    args = parser.parse_args()

    frame = pd.read_csv(args.input_csv)
    pc_columns = [f"PC_{index}" for index in range(1, 41)]
    missing = [column for column in ["sample", "barcode", *pc_columns] if column not in frame]
    if missing:
        raise RuntimeError(f"Seurat PC artifact is missing columns: {missing}")
    values = frame[pc_columns].to_numpy(dtype="<f8", copy=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "matrix.f64").open("wb") as handle:
        handle.write(b"BLMATF64")
        handle.write(struct.pack("<QQ", values.shape[0], values.shape[1]))
        handle.write(values.tobytes(order="C"))
    frame[["sample", "barcode"]].to_csv(args.output_dir / "cells.csv", index=False)
    print(f"wrote {values.shape[0]} x {values.shape[1]} PCs")


if __name__ == "__main__":
    main()
