"""Pair files -> one canonical labelled pair table (docs/TEST_PROTOCOL.md §3), and score joining.

Rules (fixed a priori, applied identically to every method):
  R1 canonical unordered pair: (a, b) with a < b lexicographically on the UniProt accession string (Part 1
     ``data.load_pairs`` rule); duplicates (same or reversed orientation, within or across lines) collapse to one.
  R2 self pairs (a == b after stripping): excluded from every analysis (not a pair of two proteins; Part 1's loader
     rejects them); counted.
  R3 label conflict (canonical pair in both the positive and the negative file): excluded from every analysis
     (label undefined); counted.
  R4 pairs also present in a train file: kept in the primary analysis; flagged (in_train_pos / in_train_neg); the
     sensitivity set S2 that drops them is analysed only if any exist (docs/part1_audit.md found 0).
  Strata C1 / C2 / C3: both / one / neither protein among the proteins of the train POSITIVE file (audit definition).
  R5 unmapped pairs (a protein without a unique DepMap gene, both proteins on the same gene, or an undefined score):
     kept in the primary analysis at the locked predictor's floor score (min of mapped scores in the scored file - 1;
     intention-to-score); the mapped-only analysis is secondary.
Pair files are two tab-separated UniProt columns, no header; blank lines skipped, extra columns ignored (as Part 1
``predict.read_pairs``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_pair_file(path) -> pd.DataFrame:
    """Raw lines -> DataFrame(x, y) of stripped accessions (Part 1 predict.read_pairs semantics)."""
    xs, ys = [], []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 2 or not f[0].strip() or not f[1].strip():
                raise ValueError(f"{path}:{n}: expected two tab-separated UniProt IDs")
            xs.append(f[0].strip())
            ys.append(f[1].strip())
    return pd.DataFrame(dict(x=xs, y=ys), dtype=object)


def canonical(df: pd.DataFrame) -> pd.DataFrame:
    x, y = df.x.values.astype(str), df.y.values.astype(str)
    a = np.where(x < y, x, y)
    b = np.where(x < y, y, x)
    return pd.DataFrame(dict(a=a, b=b), dtype=object)


def key(a, b) -> np.ndarray:
    return np.char.add(np.char.add(np.asarray(a, dtype=str), "\t"), np.asarray(b, dtype=str))


def labelled_table(pos_path, neg_path, train_paths: tuple | None = None) -> tuple[pd.DataFrame, dict]:
    """Canonical pair table over pos ∪ neg (R1-R4 flags) and an audit of counts (no pair lists).

    Columns: a, b, label (1/0; -1 for conflicts), self_pair, conflict, n_lines_pos, n_lines_neg, in_train_pos,
    in_train_neg, analysed (= not self_pair and not conflict). Sorted by (a, b)."""
    parts, audit = {}, {}
    for name, path in (("pos", pos_path), ("neg", neg_path)):
        raw = read_pair_file(path)
        c = canonical(raw)
        audit[f"{name}_lines"] = int(len(raw))
        audit[f"{name}_self_lines"] = int((c.a == c.b).sum())
        cnt = c.groupby(["a", "b"]).size().rename(f"n_lines_{name}")
        audit[f"{name}_unique_pairs"] = int(len(cnt))
        audit[f"{name}_duplicate_lines"] = int(len(raw) - len(cnt))
        audit[f"{name}_reversed_duplicates"] = int(_reversed_dups(raw))
        parts[name] = cnt
    t = pd.concat([parts["pos"], parts["neg"]], axis=1).fillna(0).astype(np.int64).reset_index()
    t = t.sort_values(["a", "b"], kind="mergesort").reset_index(drop=True)
    t["self_pair"] = (t.a == t.b).values
    t["conflict"] = ((t.n_lines_pos > 0) & (t.n_lines_neg > 0)).values
    t["label"] = np.where(t.conflict, -1, (t.n_lines_pos > 0).astype(np.int8)).astype(np.int8)
    t["analysed"] = ~(t.self_pair | t.conflict)
    if train_paths is not None:
        tp, tn = (set(key(*canonical(read_pair_file(p))[["a", "b"]].values.T)) for p in train_paths)
        k = key(t.a.values, t.b.values)
        t["in_train_pos"] = np.fromiter((x in tp for x in k), bool, len(k))
        t["in_train_neg"] = np.fromiter((x in tn for x in k), bool, len(k))
        an = t.analysed.values
        tr_prot, tr_pos_prot = set(), set()
        for i, p in enumerate(train_paths):
            r = read_pair_file(p)
            tr_prot |= set(r.x) | set(r.y)
            if i == 0:                                          # train_paths[0] = the train POSITIVE file
                tr_pos_prot = set(r.x) | set(r.y)
        prots = pd.unique(t.loc[an, ["a", "b"]].values.ravel())
        t["a_in_train"] = t.a.isin(tr_prot).values
        t["b_in_train"] = t.b.isin(tr_prot).values
        # audit strata (docs/part1_audit.md §4): number of the pair's proteins among train-POSITIVE proteins
        n_tp = t.a.isin(tr_pos_prot).values.astype(int) + t.b.isin(tr_pos_prot).values.astype(int)
        t["c_stratum"] = np.array(["C3", "C2", "C1"], dtype=object)[n_tp]
        y = t.label.values == 1
        audit.update(
            train_overlap=dict(
                pos_in_train_pos=int((an & y & t.in_train_pos).sum()), pos_in_train_neg=int((an & y & t.in_train_neg).sum()),
                neg_in_train_pos=int((an & ~y & t.in_train_pos).sum()), neg_in_train_neg=int((an & ~y & t.in_train_neg).sum())),
            proteins_analysed=int(len(prots)),
            proteins_also_in_train=int(len(set(prots) & tr_prot)),
            proteins_also_in_train_positives=int(len(set(prots) & tr_pos_prot)),
            c_strata={c: dict(pos=int((an & y & (t.c_stratum == c)).sum()), neg=int((an & ~y & (t.c_stratum == c)).sum()))
                      for c in ("C1", "C2", "C3")},
            pairs_both_proteins_unseen_in_train=dict(
                pos=int((an & y & ~t.a_in_train & ~t.b_in_train).sum()),
                neg=int((an & ~y & ~t.a_in_train & ~t.b_in_train).sum())))
    audit.update(unique_canonical_pairs=int(len(t)), self_pairs=int(t.self_pair.sum()), conflicts=int(t.conflict.sum()),
                 analysed_pos=int((t.analysed & (t.label == 1)).sum()),
                 analysed_neg=int((t.analysed & (t.label == 0)).sum()))
    audit["analysed_prevalence"] = audit["analysed_pos"] / max(audit["analysed_pos"] + audit["analysed_neg"], 1)
    return t, audit


def _reversed_dups(raw: pd.DataFrame) -> int:
    """Lines whose reversed orientation also occurs as a line (counted once per such line)."""
    fw = set(key(raw.x.values, raw.y.values))
    rv = key(raw.y.values, raw.x.values)
    return int(sum(1 for k, (x, y) in zip(rv, zip(raw.x.values, raw.y.values)) if x != y and k in fw))


def write_id_file(tables, out_path) -> dict:
    """Unique canonical non-self pairs of one or more tables (columns a, b), IDs ONLY (no labels), sorted, written
    as the predictor's input (tab-separated, no header). One file per evaluated set, so that the predictor's
    unmapped floor (min of mapped scores in the file - 1) is shared by positives and negatives."""
    t = pd.concat([x[["a", "b"]] for x in tables], ignore_index=True).drop_duplicates()
    t = t[t.a != t.b].sort_values(["a", "b"], kind="mergesort")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out_path, sep="\t", header=False, index=False, lineterminator="\n")
    return dict(path=str(out_path), n_pairs=int(len(t)))


def read_scores(path) -> pd.DataFrame:
    """Predictor TSV (uniprot_a, uniprot_b, score, mapped) -> canonical a, b, score, mapped (one row per pair;
    duplicate canonical rows must agree exactly)."""
    s = pd.read_csv(path, sep="\t", dtype={"uniprot_a": str, "uniprot_b": str})
    c = canonical(s.rename(columns={"uniprot_a": "x", "uniprot_b": "y"}))
    c["score"] = s.score.values.astype(np.float64)
    c["mapped"] = s.mapped.values.astype(bool)
    g = c.groupby(["a", "b"]).agg(score_min=("score", "min"), score_max=("score", "max"), mapped=("mapped", "min"))
    if (g.score_min != g.score_max).any():
        raise ValueError(f"{path}: duplicated pairs with different scores")
    out = g.rename(columns={"score_min": "score"}).drop(columns="score_max").reset_index()
    if np.isnan(out.score.values).any():
        raise ValueError(f"{path}: NaN scores")
    return out


def attach_scores(table: pd.DataFrame, scores: pd.DataFrame, name: str) -> tuple[np.ndarray, np.ndarray]:
    """Scores and mapped flags aligned to ``table`` rows (every table pair must be scored)."""
    m = table[["a", "b"]].merge(scores, on=["a", "b"], how="left", validate="one_to_one")
    if m.score.isna().any():
        raise ValueError(f"{name}: {int(m.score.isna().sum())} pairs without a score")
    return m.score.values.astype(np.float64), m.mapped.values.astype(bool)
