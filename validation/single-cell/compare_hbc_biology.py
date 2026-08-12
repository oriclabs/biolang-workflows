#!/usr/bin/env python3
"""Compare broad PBMC identities and condition effects, not numeric clusters.

HBC contains one control and one stimulated library, so within-cell-type
results are descriptive pseudobulk log2-CPM effects, not replicate-aware
significance tests. This complements, rather than replaces, ARI and marker
comparison in compare_hbc_results.py.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score


PANELS = {
    "T": ["CD3D", "CD3E", "TRBC1", "TRBC2", "LTB"],
    "NK": ["GNLY", "NKG7", "KLRD1", "PRF1", "GZMB"],
    "B": ["MS4A1", "CD79A", "CD74", "CD37", "HLA-DRA"],
    "CD14_Mono": ["CD14", "S100A8", "S100A9", "CTSS", "LYZ"],
    "FCGR3A_Mono": ["FCGR3A", "MS4A7", "LST1", "IFITM3", "LGALS3BP"],
    "DC": ["FCER1A", "CD1C", "CST3", "CLEC10A", "HLA-DPA1"],
    "Platelet": ["PPBP", "PF4", "NRGN", "GNG11", "RGS18"],
}
IFN_GENES = [
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27",
    "MX1", "OAS1", "OAS2", "IRF7", "RSAD2", "XAF1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("seurat_dir", type=Path)
    parser.add_argument("biolang_dir", type=Path)
    parser.add_argument("ctrl_mex", type=Path)
    parser.add_argument("stim_mex", type=Path)
    parser.add_argument("output_json", type=Path)
    return parser.parse_args()


def load_cells(directory: Path) -> pd.DataFrame:
    path = directory / "cells.csv"
    if not path.exists():
        path = directory / "cells.csv.gz"
    cells = pd.read_csv(path)
    cells["sample"] = cells["sample"].astype(str)
    cells["barcode"] = cells["barcode"].astype(str)
    cells["cluster"] = cells["cluster"].astype(str)
    cells["key"] = cells["sample"] + "::" + cells["barcode"]
    return cells


def find_file(directory: Path, stem: str) -> Path:
    for name in (stem, f"{stem}.gz"):
        path = directory / name
        if path.exists():
            return path
    raise FileNotFoundError(f"{stem} not found under {directory}")


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle if line.strip()]


def read_mex(directory: Path) -> tuple[object, list[str], list[str]]:
    matrix_path = find_file(directory, "matrix.mtx")
    opener = gzip.open if matrix_path.suffix == ".gz" else open
    with opener(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr()
    feature_lines = read_lines(find_file(directory, "features.tsv"))
    genes = [
        (line.split("\t")[1] if len(line.split("\t")) > 1 else line.split("\t")[0])
        for line in feature_lines
    ]
    barcodes = [line.split("\t")[0] for line in read_lines(find_file(directory, "barcodes.tsv"))]
    if matrix.shape != (len(genes), len(barcodes)):
        raise RuntimeError(
            f"MEX shape mismatch: matrix={matrix.shape}, labels={(len(genes), len(barcodes))}"
        )
    return matrix, genes, barcodes


def labels_for_barcodes(
    cells: pd.DataFrame, sample: str, barcodes: list[str], column: str
) -> np.ndarray:
    subset = cells[cells["sample"] == sample]
    lookup = dict(zip(subset["barcode"], subset[column], strict=True))
    labels = np.array([lookup.get(barcode) for barcode in barcodes], dtype=object)
    if any(value is None for value in labels):
        raise RuntimeError(f"{sample} manifest does not cover the filtered MEX barcodes")
    return labels


def expression_type_map(
    ctrl_matrix,
    stim_matrix,
    genes: list[str],
    ctrl_barcodes: list[str],
    stim_barcodes: list[str],
    cells: pd.DataFrame,
) -> tuple[dict[str, str], dict[str, dict[str, float]]]:
    """Assign broad identities from cluster-average normalized expression.

    Positive-marker tables can omit every panel gene for a cluster. Scoring those
    tables silently assigns the first panel on ties, so use expression itself.
    Each marker is standardized across clusters before panel averaging; this
    measures marker specificity rather than rewarding generally abundant genes.
    """
    gene_index = {gene: index for index, gene in enumerate(genes)}
    selected_genes = list(dict.fromkeys(gene for panel in PANELS.values() for gene in panel))
    missing = [gene for gene in selected_genes if gene not in gene_index]
    if missing:
        raise RuntimeError(f"cell-type panel genes absent from MEX: {missing}")
    selected_indices = [gene_index[gene] for gene in selected_genes]
    clusters = sorted(cells["cluster"].unique())
    sums = {cluster: np.zeros(len(selected_genes), dtype=np.float64) for cluster in clusters}
    counts = {cluster: 0 for cluster in clusters}

    for sample, matrix, barcodes in (
        ("ctrl", ctrl_matrix, ctrl_barcodes),
        ("stim", stim_matrix, stim_barcodes),
    ):
        cluster_labels = labels_for_barcodes(cells, sample, barcodes, "cluster")
        library_sizes = np.asarray(matrix.sum(axis=0)).ravel()
        scales = np.divide(
            10_000.0,
            library_sizes,
            out=np.zeros_like(library_sizes, dtype=np.float64),
            where=library_sizes > 0,
        )
        normalized = matrix[selected_indices, :].tocsc().multiply(scales).tocsc()
        normalized.data = np.log1p(normalized.data)
        for cluster in clusters:
            indices = np.flatnonzero(cluster_labels == cluster)
            if len(indices):
                sums[cluster] += np.asarray(normalized[:, indices].sum(axis=1)).ravel()
                counts[cluster] += len(indices)

    means = np.vstack([sums[cluster] / counts[cluster] for cluster in clusters])
    centre = means.mean(axis=0)
    spread = means.std(axis=0)
    standardized = (means - centre) / np.where(spread > 0, spread, 1.0)
    selected_lookup = {gene: index for index, gene in enumerate(selected_genes)}

    by_cluster: dict[str, str] = {}
    score_artifact: dict[str, dict[str, float]] = {}
    for row_index, cluster in enumerate(clusters):
        scores = {
            label: float(
                np.mean([standardized[row_index, selected_lookup[gene]] for gene in panel])
            )
            for label, panel in PANELS.items()
        }
        by_cluster[str(cluster)] = max(scores, key=scores.get)
        score_artifact[str(cluster)] = scores
    return by_cluster, score_artifact


def aggregate_by_type(matrix, labels: np.ndarray, minimum_cells: int = 50) -> dict[str, np.ndarray]:
    result = {}
    for label in sorted(set(labels)):
        indices = np.flatnonzero(labels == label)
        if len(indices) >= minimum_cells:
            result[str(label)] = np.asarray(matrix[:, indices].sum(axis=1)).ravel()
    return result


def log2_cpm(counts: np.ndarray) -> np.ndarray:
    return np.log2((counts + 0.5) / (counts.sum() + 1.0) * 1_000_000.0)


def condition_effects(
    ctrl_matrix,
    stim_matrix,
    ctrl_labels: np.ndarray,
    stim_labels: np.ndarray,
) -> dict[str, np.ndarray]:
    ctrl = aggregate_by_type(ctrl_matrix, ctrl_labels)
    stim = aggregate_by_type(stim_matrix, stim_labels)
    return {
        label: log2_cpm(stim[label]) - log2_cpm(ctrl[label])
        for label in sorted(ctrl.keys() & stim.keys())
    }


def effect_summary(
    genes: list[str],
    seurat_effects: dict[str, np.ndarray],
    biolang_effects: dict[str, np.ndarray],
) -> dict:
    gene_index = {gene: index for index, gene in enumerate(genes)}
    common_types = sorted(seurat_effects.keys() & biolang_effects.keys())
    per_type = {}
    for label in common_types:
        left = seurat_effects[label]
        right = biolang_effects[label]
        correlation = float(np.corrcoef(left, right)[0, 1])
        left_top = set(np.argsort(left)[-100:])
        right_top = set(np.argsort(right)[-100:])
        ifn_indices = [gene_index[gene] for gene in IFN_GENES if gene in gene_index]
        per_type[label] = {
            "cells_effect_correlation": correlation,
            "top100_stimulated_gene_overlap": len(left_top & right_top) / 100.0,
            "seurat_ifn_median_log2fc": float(np.median(left[ifn_indices])),
            "biolang_ifn_median_log2fc": float(np.median(right[ifn_indices])),
            "seurat_ifn_positive": int(np.sum(left[ifn_indices] > 0)),
            "biolang_ifn_positive": int(np.sum(right[ifn_indices] > 0)),
            "ifn_genes_tested": len(ifn_indices),
        }
    return {
        "common_cell_types": common_types,
        "per_cell_type": per_type,
        "median_effect_correlation": float(
            np.median([row["cells_effect_correlation"] for row in per_type.values()])
        ),
        "median_top100_stimulated_gene_overlap": float(
            np.median([row["top100_stimulated_gene_overlap"] for row in per_type.values()])
        ),
    }


def main() -> None:
    args = parse_args()
    seurat_cells = load_cells(args.seurat_dir)
    biolang_cells = load_cells(args.biolang_dir)

    ctrl_matrix, ctrl_genes, ctrl_barcodes = read_mex(args.ctrl_mex)
    stim_matrix, stim_genes, stim_barcodes = read_mex(args.stim_mex)
    if ctrl_genes != stim_genes:
        raise RuntimeError("control and stimulated MEX gene axes differ")

    seurat_map, seurat_scores = expression_type_map(
        ctrl_matrix, stim_matrix, ctrl_genes, ctrl_barcodes, stim_barcodes, seurat_cells
    )
    biolang_map, biolang_scores = expression_type_map(
        ctrl_matrix, stim_matrix, ctrl_genes, ctrl_barcodes, stim_barcodes, biolang_cells
    )
    seurat_cells["cell_type"] = seurat_cells["cluster"].map(seurat_map)
    biolang_cells["cell_type"] = biolang_cells["cluster"].map(biolang_map)
    joined = seurat_cells[["key", "cell_type"]].merge(
        biolang_cells[["key", "cell_type"]],
        on="key",
        suffixes=("_seurat", "_biolang"),
        validate="one_to_one",
    )

    seurat_effects = condition_effects(
        ctrl_matrix,
        stim_matrix,
        labels_for_barcodes(seurat_cells, "ctrl", ctrl_barcodes, "cell_type"),
        labels_for_barcodes(seurat_cells, "stim", stim_barcodes, "cell_type"),
    )
    biolang_effects = condition_effects(
        ctrl_matrix,
        stim_matrix,
        labels_for_barcodes(biolang_cells, "ctrl", ctrl_barcodes, "cell_type"),
        labels_for_barcodes(biolang_cells, "stim", stim_barcodes, "cell_type"),
    )

    results = {
        "joined_cells": len(joined),
        "broad_cell_type_ari": float(
            adjusted_rand_score(joined["cell_type_seurat"], joined["cell_type_biolang"])
        ),
        "broad_cell_type_exact_agreement": float(
            np.mean(joined["cell_type_seurat"] == joined["cell_type_biolang"])
        ),
        "seurat_cluster_to_type": seurat_map,
        "biolang_cluster_to_type": biolang_map,
        "panel_scoring_method": (
            "cluster-average log1p(CP10K) expression; each marker standardized "
            "across clusters before averaging within its broad-cell-type panel"
        ),
        "seurat_expression_panel_scores": seurat_scores,
        "biolang_expression_panel_scores": biolang_scores,
        "condition_effects": effect_summary(ctrl_genes, seurat_effects, biolang_effects),
        "interpretation": (
            "One library per condition: log2-CPM effects are descriptive and "
            "must not be reported as replicate-aware differential-expression p-values."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
