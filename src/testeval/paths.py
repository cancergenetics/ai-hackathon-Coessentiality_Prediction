"""All file locations used by the test-evaluation harness (one place, no hard-coded paths elsewhere)."""
from __future__ import annotations

import os
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
PART1 = PROJECT / "external" / "part1_locked"            # locked clone (read-only), commit LOCKED_COMMIT
LOCKED_COMMIT = "c7f95ee007576be9311c0a5f08c6ac63ca23a009"
# original Part 1 folder (read-only): development labels + its data_cache (read, never written)
PART1_ORIG = Path(os.environ.get("PART1_DIR", PROJECT.parent / "coessentiality_protein complexes"))
TRAIN_POS = PART1_ORIG / "datasets" / "HuMAP" / "Train_positive_pairs.txt"
TRAIN_NEG = PART1_ORIG / "datasets" / "HuMAP" / "Train_negative_pairs.txt"
P1_ASSIGNMENT = PART1_ORIG / "data_cache" / "splits" / "protein_assignment.csv"
P1_ASSIGNMENT_SHA = PART1 / "results" / "splits" / "protein_assignment.sha256"
P1_LOOK01 = PART1 / "results" / "ift" / "look01_2026-10-27.json"
P1_POSTHOC = PART1 / "results" / "revision" / "ift_posthoc_spec_k20_tinf_std.json"

PROTOCOL = PROJECT / "docs" / "TEST_PROTOCOL.md"
CACHE = PROJECT / "data_cache" / "test_eval"           # gitignored: pair-ID files and score TSVs
SCORES = CACHE / "scores"                               # SCORES/<set>/<method>.tsv
RESULTS = PROJECT / "results" / "test_eval"             # committed results of the single test look
DRYRUN = PROJECT / "results" / "test_eval_dryrun"       # committed dry-run results (development data only)
FITTED_RULES = RESULTS / "fitted_rules.json"           # decision rules + calibrators fitted on train pairs only
LOOKS = RESULTS / "looks.jsonl"


def test_files() -> tuple[Path, Path]:
    """hu.MAP test pair files. Only scripts/20_score_test_pairs.sh (mode test) and scripts/21_evaluate_test.py
    (mode test) may call this (tests/test_testeval.py::test_only_entry_points_touch_test_files)."""
    d = PROJECT / "datasets" / "huMAP"
    return d / "Test_positive_pairs.txt", d / "Test_negative_pairs.txt"
