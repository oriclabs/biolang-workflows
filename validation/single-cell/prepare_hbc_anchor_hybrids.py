#!/usr/bin/env python3
"""Create validation-only anchor sets that isolate pair and score differences."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read(path: Path) -> tuple[list[str], dict[tuple[int, int], dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = {
            (int(row["left"]), int(row["right"])): row
            for row in reader
        }
        return list(reader.fieldnames or []), rows


def write(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.mkdir(parents=True, exist_ok=False)
    with (path / "anchors.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def with_scores(
    base: dict[tuple[int, int], dict[str, str]],
    scores: dict[tuple[int, int], dict[str, str]],
) -> list[dict[str, str]]:
    result = []
    for pair, original in base.items():
        row = dict(original)
        source = scores.get(pair)
        if source is not None:
            for field in ("score", "raw_score"):
                if field in source:
                    row[field] = source[field]
        result.append(row)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat", type=Path)
    parser.add_argument("biolang", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    seurat_fields, seurat = read(args.seurat / "anchors.csv")
    biolang_fields, biolang = read(args.biolang / "anchors.csv")
    fields = list(dict.fromkeys(seurat_fields + biolang_fields))
    common = set(seurat) & set(biolang)

    write(
        args.output / "seurat-pairs-biolang-scores",
        fields,
        with_scores(seurat, biolang),
    )
    write(
        args.output / "biolang-pairs-seurat-scores",
        fields,
        with_scores(biolang, seurat),
    )
    write(
        args.output / "common-pairs-seurat-scores",
        fields,
        [seurat[pair] for pair in sorted(common)],
    )
    write(
        args.output / "common-pairs-biolang-scores",
        fields,
        [biolang[pair] for pair in sorted(common)],
    )
    print(
        f"seurat={len(seurat)} biolang={len(biolang)} common={len(common)} "
        f"seurat_only={len(set(seurat) - common)} "
        f"biolang_only={len(set(biolang) - common)}"
    )


if __name__ == "__main__":
    main()
