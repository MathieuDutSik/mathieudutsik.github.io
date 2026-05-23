#!/usr/bin/env bash
# Copy the built WASM artifacts from the sibling wasm_polyhedral repo into
# this site's wasm/ directory. After running, commit the updated files.
#
# Override the source path with WASM_POLYHEDRAL=/path/to/wasm_polyhedral

set -euo pipefail

cd "$(dirname "$0")"

SRC_DEFAULT="$(cd ../../GITwasm_polyhedral/wasm_polyhedral && pwd)"
SRC="${WASM_POLYHEDRAL:-$SRC_DEFAULT}"

if [ ! -d "$SRC/dist" ]; then
    echo "wasm_polyhedral/dist/ not found at $SRC" >&2
    echo "set WASM_POLYHEDRAL env var, or run build.sh in that repo first" >&2
    exit 1
fi

mkdir -p wasm

# Copy each binding's .js + .wasm pair. The current module bundles every
# WASM entry point into one binary; add new module basenames here if/when
# you split the build.
for binding in polyhedral; do
    for ext in js wasm; do
        src="$SRC/dist/${binding}.${ext}"
        dst="wasm/${binding}.${ext}"
        if [ ! -f "$src" ]; then
            echo "missing artifact: $src" >&2
            exit 1
        fi
        cp "$src" "$dst"
        echo "  $src  ->  $dst"
    done
done

echo ""
echo "Don't forget to: git add wasm/ && git commit"
