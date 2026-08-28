#!/usr/bin/env python3
"""Aggregate HBC evidence into regression gates and a compact measured table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-dir", type=Path, required=True)
    parser.add_argument("--gpu-dir", type=Path, required=True)
    parser.add_argument("--cpu-comparison", type=Path, required=True)
    parser.add_argument("--cpu-biology", type=Path, required=True)
    parser.add_argument("--backend-comparison", type=Path, required=True)
    parser.add_argument("--fixed-marker-comparison", type=Path, required=True)
    parser.add_argument("--sct-control", type=Path, required=True)
    parser.add_argument("--sct-stimulated", type=Path, required=True)
    parser.add_argument("--exact-dir", type=Path, required=True)
    parser.add_argument("--exact-labels", type=Path, required=True)
    parser.add_argument("--exact-bl-peak-gib", type=float, required=True)
    parser.add_argument("--exact-r-observed-peak-gib", type=float, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-markdown", type=Path, required=True)
    parser.add_argument("--seurat-wall-seconds", type=float, default=1548.3)
    parser.add_argument("--seurat-peak-gib", type=float, default=12.50)
    parser.add_argument("--baseline-wall-seconds", type=float, default=446.439)
    parser.add_argument("--baseline-peak-gib", type=float, default=8.683)
    args = parser.parse_args()

    cpu_resources = load_json(args.cpu_dir / "resources.json")
    gpu_resources = load_json(args.gpu_dir / "resources.json")
    cpu_summary = load_csv_row(args.cpu_dir / "summary.csv")
    gpu_summary = load_csv_row(args.gpu_dir / "summary.csv")
    comparison = load_json(args.cpu_comparison)
    biology = load_json(args.cpu_biology)
    backends = load_json(args.backend_comparison)
    fixed_markers = load_json(args.fixed_marker_comparison)
    sct_control = load_json(args.sct_control)
    sct_stimulated = load_json(args.sct_stimulated)
    exact_summary = load_csv_row(args.exact_dir / "summary.csv")
    exact_labels = load_json(args.exact_labels)

    cpu_wall = float(cpu_resources["wall_seconds"])
    cpu_peak = float(cpu_resources["peak_working_set_gib"])
    gpu_wall = float(gpu_resources["wall_seconds"])
    gpu_peak = float(gpu_resources["peak_working_set_gib"])
    broad_agreement = float(biology["broad_cell_type_exact_agreement"])
    native_ari = float(comparison["adjusted_rand_index"])

    gates = {
        "cpu_run_succeeded": int(cpu_resources["exit_code"]) == 0,
        "gpu_run_succeeded": int(gpu_resources["exit_code"]) == 0,
        "all_cells_joined": int(comparison["joined_cells"]) == 29629,
        "integration_feature_overlap_at_least_0_99": float(comparison["feature_overlap_min_set"]) >= 0.99,
        "native_cluster_ari_regression_floor_0_65": native_ari >= 0.65,
        "broad_identity_agreement_regression_floor_0_94": broad_agreement >= 0.94,
        "condition_effect_correlation_at_least_0_90": float(biology["condition_effects"]["median_effect_correlation"]) >= 0.90,
        "fixed_marker_pair_jaccard_at_least_0_99": float(fixed_markers["marker_pair_jaccard"]) >= 0.99,
        "fixed_marker_top50_overlap_at_least_0_99": float(fixed_markers["top50_overlap_fraction"]) >= 0.99,
        "sct_control_all_scale_sensitive_gates": bool(sct_control["passed"]),
        "sct_stimulated_all_scale_sensitive_gates": bool(sct_stimulated["passed"]),
        "cpu_gpu_scientific_parity": bool(backends["passed"]),
        "strict_external_partition_exact": float(exact_labels["adjusted_rand_index"]) >= 0.999999,
        "strict_external_anchor_counts_exact": int(exact_summary["candidates"]) == 29927 and int(exact_summary["retained"]) == 19232,
        "cpu_faster_than_seurat": cpu_wall < args.seurat_wall_seconds,
        "cpu_lower_host_memory_than_seurat": cpu_peak < args.seurat_peak_gib,
        "cpu_wall_no_more_than_10pct_over_baseline": cpu_wall <= args.baseline_wall_seconds * 1.10,
        "cpu_peak_no_more_than_10pct_over_baseline": cpu_peak <= args.baseline_peak_gib * 1.10,
    }
    targets = {
        "native_exact_cluster_ari_at_least_0_95": native_ari >= 0.95,
        "broad_identity_agreement_at_least_0_95": broad_agreement >= 0.95,
        "gpu_faster_than_cpu": gpu_wall < cpu_wall,
    }
    result = {
        "schema": "biolang.hbc.validation/v1",
        "passed": all(gates.values()),
        "gates": gates,
        "targets": targets,
        "measured": {
            "cells": int(comparison["joined_cells"]),
            "seurat_clusters": int(comparison["seurat_clusters"]),
            "biolang_native_clusters": int(cpu_summary["clusters"]),
            "native_ari": native_ari,
            "native_ami": float(comparison["adjusted_mutual_information"]),
            "mapped_cell_accuracy": float(comparison["one_to_one_mapped_accuracy"]),
            "broad_identity_agreement": broad_agreement,
            "broad_identity_ari": float(biology["broad_cell_type_ari"]),
            "integration_feature_overlap": float(comparison["feature_overlap_min_set"]),
            "fixed_marker_pair_jaccard": float(fixed_markers["marker_pair_jaccard"]),
            "fixed_marker_top50_overlap": float(fixed_markers["top50_overlap_fraction"]),
            "cpu_wall_seconds": cpu_wall,
            "cpu_peak_host_gib": cpu_peak,
            "gpu_wall_seconds": gpu_wall,
            "gpu_peak_host_gib": gpu_peak,
            "gpu_device_memory": "not measured",
            "seurat_wall_seconds": args.seurat_wall_seconds,
            "seurat_peak_host_gib": args.seurat_peak_gib,
            "cpu_speedup_over_seurat": args.seurat_wall_seconds / cpu_wall,
            "cpu_host_memory_reduction_vs_seurat": 1.0 - cpu_peak / args.seurat_peak_gib,
            "cpu_gpu_cluster_ari": float(backends["cluster_ari"]),
            "cpu_gpu_pc_relative_rmse": float(backends["pc_relative_rmse"]),
            "strict_external_clusters": int(exact_summary["clusters"]),
            "strict_external_ari": float(exact_labels["adjusted_rand_index"]),
            "strict_external_wall_seconds": float(exact_summary["elapsed_seconds"]),
            "strict_external_bl_peak_host_gib": args.exact_bl_peak_gib,
            "strict_external_r_cca_observed_peak_host_gib": args.exact_r_observed_peak_gib,
            "sct_control_top3000_overlap": float(sct_control["metrics"]["top_feature_overlap"]),
            "sct_stimulated_top3000_overlap": float(sct_stimulated["metrics"]["top_feature_overlap"]),
            "sct_control_residual_rmse_over_sd": float(sct_control["metrics"]["residual_rmse_over_oracle_sd"]),
            "sct_stimulated_residual_rmse_over_sd": float(sct_stimulated["metrics"]["residual_rmse_over_oracle_sd"]),
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    status = lambda value: "pass" if value else "not met"
    measured = result["measured"]
    lines = [
        "# HBC validation summary",
        "",
        "| Boundary | Measured result |",
        "|---|---:|",
        f"| Cells joined | {measured['cells']:,} / 29,629 |",
        f"| Seurat / native BioLang clusters | {measured['seurat_clusters']} / {measured['biolang_native_clusters']} |",
        f"| Native ARI / AMI | {measured['native_ari']:.5f} / {measured['native_ami']:.5f} |",
        f"| Broad identity exact agreement | {measured['broad_identity_agreement'] * 100:.2f}% |",
        f"| Integration-feature overlap | {measured['integration_feature_overlap'] * 100:.2f}% |",
        f"| Fixed-label marker pair Jaccard | {measured['fixed_marker_pair_jaccard'] * 100:.2f}% |",
        f"| Fixed-label top-50 marker overlap | {measured['fixed_marker_top50_overlap'] * 100:.2f}% |",
        f"| CPU wall / peak host memory | {cpu_wall:.3f} s / {cpu_peak:.3f} GiB |",
        f"| GPU wall / peak host memory | {gpu_wall:.3f} s / {gpu_peak:.3f} GiB + unmeasured device memory |",
        f"| Seurat wall / peak host memory | {args.seurat_wall_seconds:.1f} s / {args.seurat_peak_gib:.2f} GiB |",
        f"| Strict external clusters / ARI | {measured['strict_external_clusters']} / {measured['strict_external_ari']:.4f} |",
        "",
        "## Regression gates",
        "",
    ]
    lines.extend(f"- {status(passed)}: `{name}`" for name, passed in gates.items())
    lines.extend(["", "## Aspirational targets", ""])
    lines.extend(f"- {status(passed)}: `{name}`" for name, passed in targets.items())
    lines.append("")
    args.output_markdown.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        failed = ", ".join(name for name, passed in gates.items() if not passed)
        raise SystemExit(f"HBC regression gates failed: {failed}")


if __name__ == "__main__":
    main()
