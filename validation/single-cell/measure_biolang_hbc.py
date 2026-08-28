#!/usr/bin/env python3
"""Run a BioLang HBC workflow and record wall time and peak host memory."""

from __future__ import annotations

import argparse
import ctypes
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


class _ProcessMemoryCounters(ctypes.Structure):
    _fields_ = [
        ("cb", ctypes.c_ulong),
        ("PageFaultCount", ctypes.c_ulong),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


def peak_working_set(process: subprocess.Popen[bytes]) -> int:
    """Return the Windows peak working set, or zero on unsupported systems."""
    if os.name != "nt":
        return 0
    counters = _ProcessMemoryCounters()
    counters.cb = ctypes.sizeof(counters)
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    if not query(int(process._handle), ctypes.byref(counters), counters.cb):
        return 0
    return int(counters.PeakWorkingSetSize)


def normalized_environment() -> dict[str, str]:
    """Collapse case-duplicate Windows keys such as Path/PATH."""
    if os.name != "nt":
        return dict(os.environ)
    result: dict[str, tuple[str, str]] = {}
    for key, value in os.environ.items():
        result[key.casefold()] = (key, value)
    return {key: value for key, value in result.values()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--executable",
        default=os.environ.get("BIOLANG_EXE", "bl.exe" if os.name == "nt" else "bl"),
    )
    parser.add_argument("--biolang-source")
    parser.add_argument("--ctrl-input")
    parser.add_argument("--stim-input")
    parser.add_argument("--ctrl-sct")
    parser.add_argument("--stim-sct")
    parser.add_argument(
        "--notebook",
        default="workflows/single-cell/hbc_course_validation.bln",
    )
    parser.add_argument(
        "--output", default="validation-results/hbc-biolang-current"
    )
    parser.add_argument("--gpu", choices=("auto", "off", "on"), default="off")
    parser.add_argument("--write-svg", action="store_true")
    parser.add_argument(
        "--seed",
        type=int,
        default=123456,
        help="seed BioLang's runtime random stream (default: 123456)",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    requested_executable = Path(args.executable)
    candidate = (
        requested_executable
        if requested_executable.is_absolute()
        else repo / requested_executable
    )
    executable_text = str(candidate.resolve()) if candidate.is_file() else shutil.which(args.executable)
    if not executable_text:
        raise SystemExit(
            f"BioLang executable does not exist and is not on PATH: {args.executable}"
        )
    executable = Path(executable_text).resolve()
    notebook = (repo / args.notebook).resolve()
    output = (repo / args.output).resolve()
    if output.exists():
        raise SystemExit(f"output directory already exists: {output}")
    output.mkdir(parents=True)

    environment = normalized_environment()
    if args.biolang_source:
        source = Path(args.biolang_source).resolve()
        packages = source / "packages"
        if not (packages / "singlecell" / "biolang.toml").is_file():
            raise SystemExit(f"BioLang singlecell package not found under: {packages}")
        existing = environment.get("BIOLANG_PATH", "")
        environment["BIOLANG_PATH"] = (
            str(packages) if not existing else str(packages) + os.pathsep + existing
        )
    environment["BIOLANG_GPU"] = args.gpu
    environment["BIOLANG_HBC_OUTPUT"] = str(output)
    environment["BIOLANG_HBC_WRITE_SVG"] = "true" if args.write_svg else "false"
    required_inputs = {
        "BIOLANG_HBC_CTRL_INPUT": args.ctrl_input,
        "BIOLANG_HBC_STIM_INPUT": args.stim_input,
        "BIOLANG_HBC_CTRL_SCT": args.ctrl_sct,
        "BIOLANG_HBC_STIM_SCT": args.stim_sct,
    }
    resolved_inputs: list[Path] = []
    for name, value in required_inputs.items():
        if value:
            resolved = Path(value).resolve()
            if not resolved.exists():
                raise SystemExit(f"{name} path does not exist: {resolved}")
            environment[name] = str(resolved)
            resolved_inputs.append(resolved)
        elif name not in environment:
            raise SystemExit(f"pass --{name.removeprefix('BIOLANG_HBC_').lower().replace('_', '-')} or set {name}")
        else:
            resolved = Path(environment[name]).resolve()
            if not resolved.exists():
                raise SystemExit(f"{name} path does not exist: {resolved}")
            environment[name] = str(resolved)
            resolved_inputs.append(resolved)
    # A .bln file is an executable notebook, whereas `bl notebook file.bl`
    # renders plain source as notebook content and exits successfully without
    # evaluating it.  Select the command from the artifact type so a measured
    # "success" always means that the requested analysis actually ran.
    subcommand = "notebook" if notebook.suffix.casefold() == ".bln" else "run"
    command = [str(executable), subcommand, str(notebook)]
    record = output.parent / f"{output.name}-run.json"
    if record.exists():
        raise SystemExit(f"run record already exists: {record}")
    if subcommand == "run":
        # Keep logs and rolling resource checkpoints out of the declared
        # artifacts. The CLI hashes only stable scientific outputs after the
        # script succeeds, while this wrapper independently samples peak RSS.
        command.extend(["--record", str(record), "--seed", str(args.seed)])
        for path in resolved_inputs:
            command.extend(["--input", str(path)])
        expected_outputs = [
            output / "features.csv",
            output / "markers.csv",
            output / "cells.csv",
            output / "summary.csv",
            output / "timings.csv",
            output / "pcs.csv",
        ]
        if args.write_svg:
            expected_outputs.append(output / "umap.svg")
        for path in expected_outputs:
            command.extend(["--output", str(path)])
    creation_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0

    started = time.perf_counter()
    with (output / "stdout.log").open("wb") as stdout, (
        output / "stderr.log"
    ).open("wb") as stderr:
        process = subprocess.Popen(
            command,
            cwd=repo,
            env=environment,
            stdout=stdout,
            stderr=stderr,
            creationflags=creation_flags,
        )
        peak_bytes = 0
        last_checkpoint = 0.0
        while process.poll() is None:
            peak_bytes = max(peak_bytes, peak_working_set(process))
            # Persist a rolling sample so an external timeout or interrupted
            # terminal still leaves the observed peak and elapsed time.
            elapsed = time.perf_counter() - started
            if elapsed - last_checkpoint >= 5.0:
                partial = {
                    "executable": str(executable),
                    "notebook": str(notebook),
                    "subcommand": subcommand,
                    "gpu": args.gpu,
                    "status": "running",
                    "elapsed_seconds": round(elapsed, 3),
                    "peak_working_set_bytes_so_far": peak_bytes,
                    "peak_working_set_gib_so_far": round(peak_bytes / 1024**3, 3),
                }
                (output / "resources.partial.json").write_text(
                    json.dumps(partial, indent=2) + "\n", encoding="utf-8"
                )
                last_checkpoint = elapsed
            time.sleep(0.25)
        peak_bytes = max(peak_bytes, peak_working_set(process))

    measurement = {
        "executable": str(executable),
        "notebook": str(notebook),
        "subcommand": subcommand,
        "gpu": args.gpu,
        "write_svg": args.write_svg,
        "seed": args.seed,
        "run_record": str(record) if subcommand == "run" else None,
        "exit_code": process.returncode,
        "wall_seconds": round(time.perf_counter() - started, 3),
        "peak_working_set_bytes": peak_bytes,
        "peak_working_set_gib": round(peak_bytes / 1024**3, 3),
    }
    (output / "resources.json").write_text(
        json.dumps(measurement, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(measurement, indent=2))
    if process.returncode:
        print(f"validation failed; see {output / 'stderr.log'}", file=sys.stderr)
    return int(process.returncode or 0)


if __name__ == "__main__":
    raise SystemExit(main())
