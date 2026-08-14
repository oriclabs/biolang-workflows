#!/usr/bin/env python3
"""Compare every neighbour identity/rank used by anchor scoring."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = {}
    for stem in (
        "left-to-right-neighbours",
        "right-to-left-neighbours",
        "left-within-neighbours",
        "right-within-neighbours",
    ):
        reference = pd.read_csv(args.reference / f"{stem}.csv").to_numpy(dtype=np.int64)
        candidate = pd.read_csv(args.candidate / f"{stem}.csv").to_numpy(dtype=np.int64)
        if reference.shape != candidate.shape:
            raise ValueError(f"{stem}: shape mismatch {reference.shape} != {candidate.shape}")
        equal = reference == candidate
        result[stem] = {
            "shape": list(reference.shape),
            "exact_entries": int(equal.sum()),
            "total_entries": int(equal.size),
            "exact_fraction": float(equal.mean()),
            "rows_exact": int(np.all(equal, axis=1).sum()),
            "rows_total": int(equal.shape[0]),
        }
    rendered = json.dumps(result, indent=2)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
