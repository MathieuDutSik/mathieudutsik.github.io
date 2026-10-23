#!/usr/bin/env bash
# Copy the zone diagram demo (page + WASM module) from the sibling
# zone_computation repo into Demo/zone_diagram/. After running, commit the
# updated files.
#
# Build it first in that repo with `make wasm`. Override the source path with
# ZONE_COMPUTATION=/path/to/zone_computation

set -euo pipefail

cd "$(dirname "$0")"

SRC_DEFAULT="$(cd ../zone_computation && pwd)"
SRC="${ZONE_COMPUTATION:-$SRC_DEFAULT}"

if [ ! -f "$SRC/web/zone2d.wasm" ]; then
    echo "zone_computation/web/zone2d.wasm not found at $SRC" >&2
    echo "set ZONE_COMPUTATION env var, or run 'make wasm' in that repo first" >&2
    exit 1
fi

mkdir -p Demo/zone_diagram

for f in index.html style.css app.js zone2d.js zone2d.wasm; do
    cp "$SRC/web/$f" "Demo/zone_diagram/$f"
    echo "  $SRC/web/$f  ->  Demo/zone_diagram/$f"
done

echo ""
echo "Don't forget to: git add Demo/zone_diagram && git commit"
