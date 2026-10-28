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
# pinned (Part 1 DATA.md, locked commit): PART1_DIR may relocate the folder but can never substitute the files
TRAIN_SHA256 = {"Train_positive_pairs.txt": "77276dbac8ca108cac4a9c196b0d26bfac3912b616d9647b44adacdc88092fc9",
                "Train_negative_pairs.txt": "c37e27dfe12ddfcecc34e084719a4f8b2aaf47e106455e3da46d7199aed836f3"}
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


def verify_train_files(pos=None, neg=None) -> dict:
    """Raise unless the train files hash to the pinned Part 1 values (called by every entry point that reads them)."""
    import hashlib
    out = {}
    for f in (pos or TRAIN_POS, neg or TRAIN_NEG):
        h = hashlib.sha256(Path(f).read_bytes()).hexdigest()
        if TRAIN_SHA256.get(Path(f).name) != h:
            raise SystemExit(f"{f}: SHA-256 {h} is not the pinned Part 1 train file hash: refused")
        out[Path(f).name] = h
    return out
