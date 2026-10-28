#!/usr/bin/env bash
# FZ-C score runs on the construction set A (DESIGN v3.1 §2.3): zavg_BEF then GLS on A, 20 PatientID-blocked 80%
# subsamples (seeds 1000-1019) and R=4 split halves (seeds 2000-2003). Serial (one ARPACK job at a time, round-2
# algorithms SHOULD 10); one process per run so each sidecar's peak RSS is that run's.
#   PART1_DIR=... bash scripts/34_run_A_blocks.sh
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONHASHSEED=0
LOG=results/part2/construction/runs
mkdir -p "$LOG"
SETS="A $(for s in $(seq 1000 1019); do printf 'boot%s ' $s; done)$(for s in 2000 2001 2002 2003; do printf 'half%sa half%sb ' $s $s; done)"
for method in zavg_BEF gls; do
  for set in $SETS; do
    echo "$(date '+%F %T') start $method $set"
    .venv/bin/python scripts/31_score_blocks.py --set "$set" --method "$method" > "$LOG/run_${set}_${method}.log" 2>&1
    echo "$(date '+%F %T') done  $method $set"
  done
done
echo ALL_DONE
