"""Compare cluster partitions by barcode without third-party dependencies."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


def read_labels(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return {
            (f'{row["sample"]}::{row["barcode"]}' if "sample" in row else row["barcode"]):
            row["cluster"]
            for row in rows
        }


def choose2(value: int) -> int:
    return value * (value - 1) // 2


def adjusted_rand(left: list[str], right: list[str]) -> float:
    n = len(left)
    if n < 2:
        return 1.0
    contingency = Counter(zip(left, right))
    left_counts = Counter(left)
    right_counts = Counter(right)
    sum_pairs = sum(choose2(count) for count in contingency.values())
    left_pairs = sum(choose2(count) for count in left_counts.values())
    right_pairs = sum(choose2(count) for count in right_counts.values())
    total_pairs = choose2(n)
    expected = left_pairs * right_pairs / total_pairs
    maximum = 0.5 * (left_pairs + right_pairs)
    denominator = maximum - expected
    return 1.0 if math.isclose(denominator, 0.0) else (sum_pairs - expected) / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--min-ari", type=float, default=0.70)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    reference = read_labels(args.reference)
    candidate = read_labels(args.candidate)
    shared = sorted(reference.keys() & candidate.keys())
    if not shared:
        raise SystemExit("no shared barcodes")
    ari = adjusted_rand(
        [reference[barcode] for barcode in shared],
        [candidate[barcode] for barcode in shared],
    )
    result = {
        "shared_cells": len(shared),
        "reference_cells": len(reference),
        "candidate_cells": len(candidate),
        "reference_clusters": len({reference[cell] for cell in shared}),
        "candidate_clusters": len({candidate[cell] for cell in shared}),
        "adjusted_rand_index": ari,
    }
    print(
        f"shared={result['shared_cells']} reference={result['reference_cells']} "
        f"candidate={result['candidate_cells']} ARI={ari:.4f}"
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if ari < args.min_ari:
        raise SystemExit(
            f"ARI {ari:.4f} is below required threshold {args.min_ari:.4f}"
        )


if __name__ == "__main__":
    main()
