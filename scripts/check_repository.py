#!/usr/bin/env python3
"""Dependency-free structural and large-artifact checks."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


REQUIRED = (
    "LICENSE",
    "workflow.toml",
    "books/practical-bioinformatics",
    "books/biostatistics",
    "courses/hbc-scrnaseq",
    "courses/hbc-scrnaseq-validated",
    "examples",
    "workflows/single-cell/hbc_course_validation.bln",
    "validation/single-cell/run_sctransform_validation.py",
)

FORBIDDEN_DIRECTORIES = {
    "validation-results",
    "benchmark-results",
    "__pycache__",
    ".validation-r-library",
}

MAX_COMMITTED_FILE_BYTES = 10 * 1024 * 1024


def candidate_files(root: Path) -> list[Path]:
    if (root / ".git").exists():
        completed = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=root,
            check=True,
            stdout=subprocess.PIPE,
        )
        return [
            root / value.decode("utf-8")
            for value in completed.stdout.split(b"\0")
            if value
        ]
    return [path for path in root.rglob("*") if path.is_file()]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    errors: list[str] = []
    for relative in REQUIRED:
        if not (root / relative).exists():
            errors.append(f"missing required path: {relative}")

    for path in candidate_files(root):
        relative = path.relative_to(root)
        if any(part in FORBIDDEN_DIRECTORIES for part in relative.parts):
            errors.append(f"generated directory must not be committed: {relative}")
        elif path.stat().st_size > MAX_COMMITTED_FILE_BYTES:
            errors.append(
                f"file exceeds 10 MiB policy: {relative} ({path.stat().st_size} bytes)"
            )

    if errors:
        print("repository check failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print("repository check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
