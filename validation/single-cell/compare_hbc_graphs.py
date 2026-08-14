#!/usr/bin/env python3
"""Compare fixed-PC kNN/SNN graphs and their Louvain partitions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score


def partition_metrics(reference: pd.DataFrame, observed: pd.DataFrame) -> dict[str, float | int]:
    keys = ["sample", "barcode"]
    joined = reference[keys + ["cluster"]].merge(
        observed[keys + ["cluster"]], on=keys, suffixes=("_reference", "_observed"),
        validate="one_to_one",
    )
    if len(joined) != len(reference) or len(joined) != len(observed):
        raise RuntimeError("cell identity mismatch")
    contingency = pd.crosstab(joined.cluster_observed, joined.cluster_reference)
    rows, columns = linear_sum_assignment(-contingency.to_numpy())
    matched = int(contingency.to_numpy()[rows, columns].sum())
    return {
        "reference_clusters": int(joined.cluster_reference.nunique()),
        "observed_clusters": int(joined.cluster_observed.nunique()),
        "adjusted_rand_index": float(adjusted_rand_score(
            joined.cluster_reference, joined.cluster_observed
        )),
        "adjusted_mutual_information": float(adjusted_mutual_info_score(
            joined.cluster_reference, joined.cluster_observed
        )),
        "one_to_one_mapped_accuracy": matched / len(joined),
    }


def edge_key(frame: pd.DataFrame) -> pd.Series:
    low = np.minimum(frame.source.to_numpy(), frame.target.to_numpy())
    high = np.maximum(frame.source.to_numpy(), frame.target.to_numpy())
    return (
        pd.Series(low, index=frame.index).astype(str)
        + ":"
        + pd.Series(high, index=frame.index).astype(str)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_dir", type=Path)
    parser.add_argument("biolang_dir", type=Path)
    parser.add_argument("biolang_same_graph_dir", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    seurat_knn = pd.read_csv(args.seurat_dir / "knn.csv")
    biolang_knn = pd.read_csv(args.biolang_dir / "knn.csv")
    seurat_knn = seurat_knn[seurat_knn.source != seurat_knn.target]
    reference_pairs = set(zip(seurat_knn.source, seurat_knn.target))
    observed_pairs = set(zip(biolang_knn.source, biolang_knn.target))
    common_pairs = reference_pairs & observed_pairs
    union_pairs = reference_pairs | observed_pairs

    seurat_snn = pd.read_csv(args.seurat_dir / "snn.csv")
    biolang_snn = pd.read_csv(args.biolang_dir / "snn.csv")
    seurat_snn["key"] = edge_key(seurat_snn)
    biolang_snn["key"] = edge_key(biolang_snn)
    seurat_snn = seurat_snn.drop_duplicates("key")
    biolang_snn = biolang_snn.drop_duplicates("key")
    joined_snn = seurat_snn[["key", "weight"]].merge(
        biolang_snn[["key", "weight"]], on="key",
        suffixes=("_seurat", "_biolang"), validate="one_to_one",
    )
    seurat_edges = set(seurat_snn.key)
    biolang_edges = set(biolang_snn.key)

    result = {
        "fixed_input": "Seurat integrated PCs 1:40",
        "knn": {
            "seurat_nonself_pairs": len(reference_pairs),
            "biolang_pairs": len(observed_pairs),
            "identity_recall": len(common_pairs) / len(reference_pairs),
            "identity_jaccard": len(common_pairs) / len(union_pairs),
        },
        "snn": {
            "seurat_edges": len(seurat_edges),
            "biolang_edges": len(biolang_edges),
            "identity_recall": len(seurat_edges & biolang_edges) / len(seurat_edges),
            "identity_jaccard": len(seurat_edges & biolang_edges) / len(seurat_edges | biolang_edges),
            "common_weight_correlation": float(joined_snn.weight_seurat.corr(
                joined_snn.weight_biolang
            )),
            "common_weight_median_absolute_error": float(np.median(np.abs(
                joined_snn.weight_seurat - joined_snn.weight_biolang
            ))),
        },
        "independent_graph_partition": partition_metrics(
            pd.read_csv(args.seurat_dir / "cells.csv"),
            pd.read_csv(args.biolang_dir / "cells.csv"),
        ),
        "same_seurat_graph_partition": partition_metrics(
            pd.read_csv(args.seurat_dir / "cells.csv"),
            pd.read_csv(args.biolang_same_graph_dir / "cells.csv"),
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
