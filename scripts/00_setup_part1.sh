#!/usr/bin/env bash
# Recreate the locked Part 1 snapshot and the Python env. Part 1 folder is read-only; we only clone it.
#
# Two routes to the same locked code:
#   1. PART1_DIR is the original Part 1 repository: clone it and check out the locked commit c7f95ee.
#   2. Otherwise (public reproduction): clone the public release of the method,
#      https://github.com/cancergenetics/ai-hackathon-Coessentiality, at 3c684c9. Its history is fresh, so c7f95ee is
#      not in it; every Part 1 path this project reads must have the same git object id as at c7f95ee
#      (LOCKED_OBJECTS in src/testeval/paths.py and src/part2/paths.py), and this script refuses otherwise.
set -euo pipefail
cd "$(dirname "$0")/.."
P1="${PART1_DIR:-../coessentiality_protein complexes}"
LOCK=c7f95ee007576be9311c0a5f08c6ac63ca23a009
MIRROR_URL=https://github.com/cancergenetics/ai-hackathon-Coessentiality.git
MIRROR=3c684c972f1e5c3ead722b4b3e2c508fa3cb5ddf
if [ ! -d external/part1_locked ]; then
  if git -C "$P1" cat-file -e "$LOCK^{commit}" 2>/dev/null; then
    git clone -q --no-hardlinks "$P1" external/part1_locked
  else
    git clone -q "$MIRROR_URL" external/part1_locked
  fi
fi
if git -C external/part1_locked cat-file -e "$LOCK^{commit}" 2>/dev/null; then
  git -C external/part1_locked checkout -q "$LOCK"
  test "$(git -C external/part1_locked rev-parse HEAD)" = "$LOCK"
else
  git -C external/part1_locked checkout -q "$MIRROR"
  python3 - <<'PY'
import subprocess, sys
sys.path.insert(0, "src")
from testeval.paths import LOCKED_OBJECTS, PART1
bad = [p for p, oid in LOCKED_OBJECTS.items()
       if subprocess.run(["git", "-C", str(PART1), "rev-parse", f"HEAD:{p}"], capture_output=True,
                         text=True).stdout.strip() != oid]
sys.exit(f"public mirror differs from the locked commit at: {bad}" if bad else 0)
PY
fi
# datasets: Part 1 code reads datasets/ relative to its root -> symlink to ours (identical checksums, DATA.md)
[ -e external/part1_locked/datasets ] || ln -s "$(pwd)/datasets" external/part1_locked/datasets
if [ ! -d .venv ]; then
  if [ -x "$P1/.env/bin/python" ]; then
    "$P1/.env/bin/python" -m venv --system-site-packages .venv
  else
    python3 -m venv .venv && .venv/bin/pip install -q -r requirements.txt
  fi
fi
echo "part1 locked at $(git -C external/part1_locked rev-parse HEAD)"
