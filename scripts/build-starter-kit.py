#!/usr/bin/env python3
"""Build the HBC starter archive from workflow and BioLang source checkouts."""

import argparse
import os
from pathlib import Path
import shutil
import tempfile
import zipfile


REPOSITORY = Path(__file__).resolve().parents[1]
BOOK = REPOSITORY / "courses" / "hbc-scrnaseq" / "src"
OUTPUT = BOOK / "downloads" / "singlecell-starter.zip"
KIT = "singlecell-starter"


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--biolang-source",
        type=Path,
        default=Path(os.environ.get("BIOLANG_SOURCE", REPOSITORY.parent / "biolang")),
        help="BioLang source checkout containing packages/singlecell",
    )
    return parser.parse_args()


def main():
    args = arguments()
    package = args.biolang_source.resolve() / "packages" / "singlecell"
    downloads = BOOK / "downloads"

    for required in (package, downloads / "get-data.py"):
        if not required.exists():
            raise SystemExit(f"missing: {required}")

    with tempfile.TemporaryDirectory() as temporary:
        temporary_path = Path(temporary)
        stage = temporary_path / KIT
        stage.mkdir()

        shutil.copytree(
            package,
            stage / "singlecell",
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".git"),
        )
        shutil.copy2(downloads / "get-data.py", stage)

        for source in sorted(downloads.iterdir()):
            if source.suffix in {".bl", ".bln"}:
                shutil.copy2(source, stage)

        readme = downloads / "starter-kit-README.md"
        if readme.exists():
            shutil.copy2(readme, stage / "README.md")

        count = 0
        with zipfile.ZipFile(
            OUTPUT, "w", zipfile.ZIP_DEFLATED, compresslevel=9
        ) as archive:
            for folder, dirs, files in os.walk(stage):
                dirs[:] = [name for name in dirs if name != "__pycache__"]
                for name in sorted(files):
                    source = Path(folder) / name
                    relative = source.relative_to(temporary_path).as_posix()
                    archive.write(source, relative)
                    count += 1

    with zipfile.ZipFile(OUTPUT) as archive:
        invalid = [name for name in archive.namelist() if "\\" in name]
        if invalid:
            raise SystemExit(f"archive has backslash paths: {invalid[:3]}")

    print(f"{OUTPUT.name}: {count} files, {OUTPUT.stat().st_size // 1024} KB")


if __name__ == "__main__":
    main()
