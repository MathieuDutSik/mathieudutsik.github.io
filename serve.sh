#!/usr/bin/env bash
# Local development helper for the Jekyll site.
#
# Usage:
#   ./serve.sh serve     start dev server with live reload (default)
#   ./serve.sh build     one-shot build into _site/
#   ./serve.sh clean     remove _site/ and .jekyll-cache/
#   ./serve.sh install   (re)install Ruby gem dependencies
#
# The site at https://mathieudutsik.github.io/ is published from the
# default branch by GitHub Pages on push. This script is only for local
# preview; it never deploys anywhere.

set -euo pipefail

cd "$(dirname "$0")"

# Homebrew Ruby (Apple Silicon). The user shell may add this via ~/.zshrc,
# but a script invoked from a clean environment does not pick that up.
export PATH="/opt/homebrew/opt/ruby/bin:$PATH"

HOST="127.0.0.1"
PORT="4000"

usage() {
    cat <<'EOF'
Local development helper for the Jekyll site.

Usage:
  ./serve.sh serve     start dev server with live reload (default)
  ./serve.sh build     one-shot build into _site/
  ./serve.sh clean     remove _site/ and .jekyll-cache/
  ./serve.sh install   (re)install Ruby gem dependencies

The site at https://mathieudutsik.github.io/ is published from the
default branch by GitHub Pages on push. This script is only for local
preview; it never deploys anywhere.
EOF
}

ensure_deps() {
    if [ ! -d vendor/bundle ]; then
        echo "[serve.sh] Installing Ruby dependencies into vendor/bundle/ (one-time)..."
        bundle config set --local path 'vendor/bundle' >/dev/null
        bundle install
    fi
}

stop_existing_server() {
    local pids
    pids=$(lsof -ti tcp:"$PORT" 2>/dev/null || true)
    if [ -n "$pids" ]; then
        echo "[serve.sh] Port $PORT busy — stopping existing process(es): $pids"
        kill $pids 2>/dev/null || true
        sleep 1
    fi
}

cmd="${1:-serve}"
case "$cmd" in
    serve)
        ensure_deps
        stop_existing_server
        echo "[serve.sh] Serving on http://$HOST:$PORT/  (Ctrl-C to stop)"
        exec bundle exec jekyll serve \
            --host "$HOST" \
            --port "$PORT" \
            --livereload
        ;;
    build)
        ensure_deps
        exec bundle exec jekyll build
        ;;
    clean)
        rm -rf _site/ .jekyll-cache/
        echo "[serve.sh] Removed _site/ and .jekyll-cache/"
        ;;
    install)
        bundle config set --local path 'vendor/bundle' >/dev/null
        bundle install
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        echo "Unknown command: $cmd" >&2
        usage >&2
        exit 1
        ;;
esac
