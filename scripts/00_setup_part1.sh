#!/usr/bin/env bash
# Recreate the locked Part 1 snapshot and the Python env. Part 1 folder is read-only; we only clone it.
set -euo pipefail
cd "$(dirname "$0")/.."
P1="${PART1_DIR:-../coessentiality_protein complexes}"
LOCK=c7f95ee007576be9311c0a5f08c6ac63ca23a009
if [ ! -d external/part1_locked ]; then
  git clone -q --no-hardlinks "$P1" external/part1_locked
fi
git -C external/part1_locked checkout -q "$LOCK"
test "$(git -C external/part1_locked rev-parse HEAD)" = "$LOCK"
# datasets: Part 1 code reads datasets/ relative to its root -> symlink to ours (identical checksums, DATA.md)
[ -e external/part1_locked/datasets ] || ln -s "$(pwd)/datasets" external/part1_locked/datasets
[ -d .venv ] || "$P1/.env/bin/python" -m venv --system-site-packages .venv
echo "part1 locked at $LOCK"
