#!/usr/bin/env bash
# Verify dataset checksums against DATA.md
set -euo pipefail
cd "$(dirname "$0")/.."
grep -E '^\| [A-Za-z]+/[^ ]+ \|' DATA.md | awk -F'|' '{gsub(/ /,"",$2); gsub(/ /,"",$5); print $5"  datasets/"$2}' | shasum -a 256 -c -
