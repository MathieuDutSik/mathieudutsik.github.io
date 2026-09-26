#!/usr/bin/env bash
# Build the two-dimensional k-covering simulation in Demo/two_dim_simulation/
# and write the static bundle into Demo/two_dim_simulation/dist/, which is
# what the Demo page embeds. After running, commit the updated files.
#
# Requires node and pnpm. Dependencies are installed into node_modules/ for
# the build and removed afterwards to keep the checkout small; pnpm keeps
# them in its content-addressable store, so the next install is fast.

set -euo pipefail

cd "$(dirname "$0")/Demo/two_dim_simulation"

echo "[build-demo.sh] Installing dependencies..."
pnpm install --frozen-lockfile

# The site serves the app from a subdirectory, so the bundle must use
# relative asset URLs (--base ./) instead of Vite's default absolute /assets/.
./node_modules/.bin/vite build --base ./ --outDir dist --emptyOutDir

rm -rf node_modules

echo ""
echo "Installed into Demo/two_dim_simulation/dist:"
find dist -type f | sort
echo ""
echo "Don't forget to: git add Demo/ && git commit"
