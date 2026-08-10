#!/usr/bin/env sh
set -eu

if [ "$#" -lt 2 ]; then
    echo "usage: $0 BIOLANG_REPO BL_ARGUMENT..." >&2
    exit 2
fi

biolang_repo=$1
shift
executable="$biolang_repo/target/release/bl"
packages="$biolang_repo/packages"

if [ ! -x "$executable" ]; then
    echo "release executable not found: $executable; run cargo build --release -p bl-cli" >&2
    exit 2
fi
if [ ! -f "$packages/singlecell/biolang.toml" ]; then
    echo "BioLang singlecell package not found under: $packages" >&2
    exit 2
fi

BIOLANG_PATH=$packages
export BIOLANG_PATH
if [ -z "${BIOLANG_GPU+x}" ]; then
    BIOLANG_GPU=off
    export BIOLANG_GPU
fi
exec "$executable" "$@"
