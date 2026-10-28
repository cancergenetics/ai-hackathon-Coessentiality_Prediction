#!/usr/bin/env bash
# Score one pair set with the LOCKED Part 1 predictor (external/part1_locked, commit c7f95ee) for every method of
# docs/TEST_PROTOCOL.md §2. Label-free: the predictor only receives an ID-only pair file (scripts/19_make_pair_files.py).
#
#   scripts/20_score_test_pairs.sh train     # all hu.MAP train pairs (fit set for thresholds / calibration)
#   scripts/20_score_test_pairs.sh dryrun    # Part 1 dev and IFT pairs (dry run on development data)
#   scripts/20_score_test_pairs.sh test      # hu.MAP TEST pairs; refuses unless docs/TEST_PROTOCOL.md is committed
#
# Outputs: data_cache/test_eval/scores/<set>/<method>.tsv (gitignored); run records (git commit of the clone, input
# SHA-256s, null-moment checks, mapping counts, runtimes) -> results/test_eval{,_dryrun}/score_records/<set>/.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$(pwd)"
PY="$ROOT/.venv/bin/python"
P1="$ROOT/external/part1_locked"
LOCK=c7f95ee007576be9311c0a5f08c6ac63ca23a009
MODE="${1:?usage: $0 train|dryrun|test}"

test "$(git -C "$P1" rev-parse HEAD)" = "$LOCK" || { echo "Part 1 clone is not at $LOCK" >&2; exit 1; }
if [ -n "$(git -C "$P1" status --porcelain -- src scripts results)" ]; then
  echo "Part 1 clone has local modifications" >&2; exit 1
fi

case "$MODE" in
  train)  SETS="train";               REC="$ROOT/results/test_eval/score_records" ;;
  dryrun) SETS="dryrun_dev dryrun_ift"; REC="$ROOT/results/test_eval_dryrun/score_records" ;;
  test)   SETS="test";                REC="$ROOT/results/test_eval/score_records"
          "$PY" -c "import sys; sys.path.insert(0, 'src'); from testeval.protocol import require_frozen; require_frozen()" ;;
  *) echo "unknown mode $MODE" >&2; exit 2 ;;
esac

# ID files (labels never reach the predictor)
"$PY" scripts/19_make_pair_files.py --set "$( [ "$MODE" = dryrun ] && echo dryrun || echo "$MODE" )" > /dev/null

# method name used in results  ->  locked predictor --method
METHODS="zavg_BEF:zavg:B,E,F B:single gls_signed:gls spearman_signed:spearman pearson_signed:pearson"
for S in $SETS; do
  IDS="$ROOT/data_cache/test_eval/ids/$S.tsv"
  OUT="$ROOT/data_cache/test_eval/scores/$S"
  mkdir -p "$OUT" "$REC/$S"
  for M in $METHODS; do
    NAME="${M%%:*}"; ARG="${M#*:}"
    echo "[$(date +%H:%M:%S)] $S $NAME (--method $ARG)"
    ( cd "$P1" && /usr/bin/time -p "$PY" scripts/predict_pairs.py --pairs "$IDS" --method "$ARG" \
        --out "$OUT/$NAME.tsv" --record "$REC/$S/$NAME.run.json" ) 2>&1 | grep -E "^real|n_mapped|runtime_s" || true
  done
  echo "[$(date +%H:%M:%S)] $S spec_k20_tinf_std (post hoc reference)"
  /usr/bin/time -p "$PY" scripts/20b_score_posthoc_reference.py --pairs "$IDS" --out "$OUT/spec_k20_tinf_std.tsv" \
      --record "$REC/$S/spec_k20_tinf_std.run.json" 2>&1 | grep -E "^real|runtime_s" || true
done
echo "[$(date +%H:%M:%S)] done: $MODE"
