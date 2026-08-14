#!/usr/bin/env python3
"""Compare two fixed-embedding mutual-neighbour candidate sets."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def pairs(directory: Path) -> set[tuple[int, int]]:
    with (directory / "candidate-anchors.csv").open(newline="", encoding="utf-8-sig") as handle:
        return {(int(row["left"]), int(row["right"])) for row in csv.DictReader(handle)}


def scored(directory: Path) -> dict[tuple[int, int], tuple[float, float]]:
    path = directory / "anchors.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {
            (int(row["left"]), int(row["right"])): (
                float(row["score"]), float(row["raw_score"])
            )
            for row in csv.DictReader(handle)
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reference = pairs(args.reference)
    candidate = pairs(args.candidate)
    common = reference & candidate
    result = {
        "reference_count": len(reference),
        "candidate_count": len(candidate),
        "common_count": len(common),
        "reference_recall": len(common) / len(reference),
        "candidate_recall": len(common) / len(candidate),
        "jaccard": len(common) / len(reference | candidate),
        "reference_only": len(reference - candidate),
        "candidate_only": len(candidate - reference),
    }
    reference_scored = scored(args.reference)
    candidate_scored = scored(args.candidate)
    if reference_scored and candidate_scored:
        scored_common = set(reference_scored) & set(candidate_scored)
        result["retained"] = {
            "reference_count": len(reference_scored),
            "candidate_count": len(candidate_scored),
            "common_count": len(scored_common),
            "reference_only": len(set(reference_scored) - scored_common),
            "candidate_only": len(set(candidate_scored) - scored_common),
            "raw_score_exact": sum(
                reference_scored[pair][1] == candidate_scored[pair][1]
                for pair in scored_common
            ),
            "normalized_score_within_1e_12": sum(
                math.isclose(
                    reference_scored[pair][0], candidate_scored[pair][0],
                    rel_tol=0.0, abs_tol=1e-12
                )
                for pair in scored_common
            ),
        }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
