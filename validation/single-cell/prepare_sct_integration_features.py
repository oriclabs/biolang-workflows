#!/usr/bin/env python3
"""Select Seurat-style integration features from provider ranking artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("control_provider", type=Path)
    parser.add_argument("stimulated_provider", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--n-features", type=int, default=3000)
    parser.add_argument(
        "--control-variable-output",
        type=Path,
        help="optional output for the control object's own top features",
    )
    parser.add_argument(
        "--stimulated-variable-output",
        type=Path,
        help="optional output for the stimulated object's own top features",
    )
    parser.add_argument(
        "--all-candidates",
        action="store_true",
        help="write the complete model-eligible union instead of only the selected prefix",
    )
    args = parser.parse_args()

    complete_rankings = [
        pd.read_csv(directory / "ranking.csv")
        for directory in (args.control_provider, args.stimulated_provider)
    ]
    rankings = [ranking.head(args.n_features) for ranking in complete_rankings]
    # `ranking.csv` may be intentionally capped to bound residual
    # materialization. Eligibility is the full min-cells model axis, not the
    # intersection of those caps; otherwise a top feature in one sample is
    # silently removed merely because it ranks below the other sample's cap.
    model_axes = []
    for directory, ranking in zip(
        (args.control_provider, args.stimulated_provider), complete_rankings
    ):
        manifest = directory / "modelled-genes.csv"
        axis = pd.read_csv(manifest) if manifest.exists() else ranking
        model_axes.append(set(axis.gene.astype(str)))
    eligible = model_axes[0] & model_axes[1]
    by_gene = []
    index_by_gene: dict[str, int] = {}
    for ranking in rankings:
        ranks: dict[str, int] = {}
        for position, row in enumerate(ranking.itertuples(index=False), start=1):
            gene = str(row.gene)
            ranks[gene] = position
            index_by_gene[gene] = int(row.gene_index)
        by_gene.append(ranks)

    candidates = set().union(*(set(ranks) for ranks in by_gene)) & eligible
    ordered = sorted(
        candidates,
        key=lambda gene: (
            -sum(gene in ranks for ranks in by_gene),
            sorted(ranks[gene] for ranks in by_gene if gene in ranks)[
                (sum(gene in ranks for ranks in by_gene) - 1) // 2
            ]
            if sum(gene in ranks for ranks in by_gene) % 2
            else sum(ranks[gene] for ranks in by_gene if gene in ranks)
            / sum(gene in ranks for ranks in by_gene),
            gene,
        ),
    )
    if not args.all_candidates:
        ordered = ordered[: args.n_features]
    output = pd.DataFrame(
        {
            "gene_index": [index_by_gene[gene] for gene in ordered],
            "gene": ordered,
            "rank": range(1, len(ordered) + 1),
        }
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_csv, index=False)
    for ranking, path in zip(
        rankings,
        (args.control_variable_output, args.stimulated_variable_output),
    ):
        if path is None:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        ranking.loc[:, ["gene_index", "gene"]].to_csv(path, index=False)
    print(f"wrote {len(output)} selected integration features")


if __name__ == "__main__":
    main()
