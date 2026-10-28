#!/usr/bin/env python
"""Audit: were hu.MAP TEST pairs / labels / test-derived information accessed during Part 1?

Read-only with respect to Part 1 (working tree and .git are only read; git history operations run in the read-only
clone external/part1_locked, whose history is checked to equal Part 1's). The hu.MAP test pair lists are read ONLY for
set-membership / overlap checks. No score, metric, AP, AUROC or prediction is computed on test pairs.

Usage:  .venv/bin/python scripts/verify/audit_part1.py
Output: results/audit/part1_audit.json (+ JSON summary on stdout)
"""
from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[2]
P1 = Path("/Users/metinyazar/Desktop/autoresearch/coessentiality_protein complexes")
CLONE = ROOT / "external" / "part1_locked"
TEST_POS = ROOT / "datasets/huMAP/Test_positive_pairs.txt"
TEST_NEG = ROOT / "datasets/huMAP/Test_negative_pairs.txt"
SHA_POS = "c4e3e68edcfc088516b3334dbfae248fb8f276c6475c13dc621abc8b4bba2d96"
SHA_NEG = "c4a3e5468da0790690ee910af1822dcaaa301dea1afaad45816e34f5e40ae9de"
OUT = ROOT / "results/audit/part1_audit.json"
ENV = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
ACC = re.compile(r"(?<![A-Za-z0-9])(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?![A-Za-z0-9])")


def iso(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path, bs=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def git(repo, *args, text=True):
    r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, env=ENV)
    if r.returncode not in (0, 1):
        raise RuntimeError(f"git {args}: {r.stderr.decode()[:300]}")
    return r.stdout.decode(errors="replace") if text else r.stdout


def read_pairs(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != "":
                rows.append((p[0].strip(), p[1].strip()))
    return rows


def canon(a, b):
    return (a, b) if a <= b else (b, a)


# --------------------------------------------------------------------------------------------- 1. timeline
def check_timeline():
    st = {}
    for name, p in (("test_positive", TEST_POS), ("test_negative", TEST_NEG)):
        s = os.stat(p)
        st[name] = dict(path=str(p), birth_utc=iso(s.st_birthtime), mtime_utc=iso(s.st_mtime), birth_epoch=s.st_birthtime,
                        size=s.st_size, sha256=sha256(p))
    st["sha256_match"] = (st["test_positive"]["sha256"] == SHA_POS, st["test_negative"]["sha256"] == SHA_NEG)
    birth = min(st["test_positive"]["birth_epoch"], st["test_negative"]["birth_epoch"])

    fmt = "%H|%cI|%aI|%s"
    p1_log = git(P1, "log", "--all", f"--format={fmt}").strip().splitlines()
    cl_log = git(CLONE, "log", "--all", f"--format={fmt}").strip().splitlines()
    commits = []
    for ln in p1_log:
        h, c, a, s = ln.split("|", 3)
        commits.append(dict(hash=h[:7], committer_date=c, author_date=a, subject=s))
    cts = [datetime.fromisoformat(c["committer_date"]).timestamp() for c in commits]
    ats = [datetime.fromisoformat(c["author_date"]).timestamp() for c in commits]
    branches = [b.strip(" *") for b in git(P1, "branch", "-a").splitlines()]
    head = git(CLONE, "rev-parse", "HEAD").strip()
    # Part 1 working-tree / data mtimes
    after_birth, after_last_commit, n_files, newest = [], [], 0, (0, None)
    last_commit = max(cts)
    for dp, dn, fn in os.walk(P1, followlinks=False):
        dn[:] = [d for d in dn if d != ".git"]
        for f in fn:
            fp = os.path.join(dp, f)
            try:
                s = os.lstat(fp)
            except OSError:
                continue
            n_files += 1
            rel = os.path.relpath(fp, P1)
            if s.st_mtime > newest[0] and f != ".DS_Store":
                newest = (s.st_mtime, rel)
            if s.st_mtime > last_commit + 1:
                after_last_commit.append((rel, iso(s.st_mtime)))
            if s.st_mtime > birth:
                after_birth.append((rel, iso(s.st_mtime)))
    nonfinder = lambda L: [x for x in L if not x[0].endswith(".DS_Store")]
    # .git internals
    git_after = []
    for dp, dn, fn in os.walk(P1 / ".git"):
        for f in fn:
            fp = os.path.join(dp, f)
            m = os.lstat(fp).st_mtime
            if m > last_commit + 120:
                git_after.append((os.path.relpath(fp, P1), iso(m)))
    # Part 1 train files
    tr = {}
    for f in sorted((P1 / "datasets/HuMAP").iterdir()):
        tr[f.name] = dict(birth_utc=iso(os.stat(f).st_birthtime), mtime_utc=iso(os.stat(f).st_mtime), sha256=sha256(f))
    return dict(
        test_files={k: v for k, v in st.items() if k != "sha256_match"} | dict(sha256_match_expected=list(st["sha256_match"])),
        part1_n_commits=len(commits), part1_branches=branches,
        part1_first_commit=min(commits, key=lambda c: c["committer_date"])["committer_date"],
        part1_last_commit=max(commits, key=lambda c: datetime.fromisoformat(c["committer_date"]))["committer_date"],
        clone_head=head, clone_history_equals_part1=(sorted(cl_log) == sorted(p1_log)),
        commits_after_test_birth=sum(1 for t in cts + ats if t > birth),
        gap_last_commit_to_test_birth_minutes=round((birth - last_commit) / 60, 1),
        part1_files_walked=n_files,
        part1_files_mtime_after_last_commit_excl_DS_Store=nonfinder(after_last_commit),
        part1_files_mtime_after_test_birth_excl_DS_Store=nonfinder(after_birth),
        part1_ds_store_after_test_birth=[x for x in after_birth if x[0].endswith(".DS_Store")],
        part1_newest_non_finder_file=dict(path=newest[1], mtime_utc=iso(newest[0])),
        part1_git_internals_modified_after_last_commit=git_after,
        part1_datasets_HuMAP_files=tr,
        commits=commits,
    )


# --------------------------------------------------------------------------------------------- 2. references
STRICT = [r"Test_positive", r"Test_negative", r"Test_pairs", r"(huMAP|HuMAP|humap)/Test",
          SHA_POS[:12], SHA_NEG[:12], r"Test_positive_pairs\.txt", r"Test_negative_pairs\.txt"]
BROAD = r"[Tt]est[_ ](pos|neg)"


def check_references():
    revs = git(CLONE, "rev-list", "--all").split()
    strict_pat = "|".join(STRICT)
    out = git(CLONE, "grep", "-I", "-n", "-i", "-E", strict_pat, *revs, "--", ".")
    strict_hits = Counter()
    for ln in out.splitlines():
        # rev:path:lineno:text
        parts = ln.split(":", 3)
        if len(parts) == 4:
            strict_hits[(parts[1], parts[3].strip()[:160])] += 1
    strict_classified = []
    for (path, text), n in strict_hits.items():
        innocent = bool(re.search(r"test_negative_control", text, re.I)) or bool(re.search(r"def test_", text))
        strict_classified.append(dict(path=path, text=text, in_n_commits=n, innocent_test_function_name=innocent))

    out2 = git(CLONE, "grep", "-I", "-c", "-E", BROAD, *revs, "--", ".")
    broad_files = Counter()
    for ln in out2.splitlines():
        parts = ln.rsplit(":", 2)
        if len(parts) == 3:
            broad_files[parts[1]] += int(parts[2])
    cat = Counter()
    for p in broad_files:
        if p.startswith("results/") or p.startswith("manuscript/") or p.startswith("data_cache"):
            cat["results_json: internal CV fold counts (test_pos/test_neg keys)"] += 1
        elif p.startswith("tests/") or p.startswith("scripts/") or p.startswith("src/"):
            cat["code: fold-count variables / test function names"] += 1
        else:
            cat["other"] += 1
    other = [p for p in broad_files if not (p.startswith(("results/", "manuscript/", "data_cache", "tests/", "scripts/", "src/")))]

    # working tree (excluding .env/.git; binary skipped); strict pattern, case-insensitive
    rx = re.compile(strict_pat, re.I)
    wt = []
    n_scanned = 0
    for dp, dn, fn in os.walk(P1):
        dn[:] = [d for d in dn if d not in (".env", ".git", ".pytest_cache", "__pycache__")]
        for f in fn:
            fp = Path(dp) / f
            if f.endswith((".npy", ".npz", ".png", ".tif", ".pdf", ".docx", ".parquet", ".pkl")) or f == ".DS_Store":
                continue
            try:
                if fp.stat().st_size > 400_000_000:
                    continue
                opener = gzip.open if f.endswith(".gz") else open
                with opener(fp, "rt", errors="replace") as fh:
                    n_scanned += 1
                    for i, ln in enumerate(fh, 1):
                        if rx.search(ln):
                            wt.append(dict(path=str(fp.relative_to(P1)), line=i, text=ln.strip()[:160]))
            except (OSError, EOFError):
                continue
    for h in wt:
        h["innocent_test_function_name"] = bool(re.search(r"test_negative_control|def test_", h["text"], re.I))
    # copies of the test files by size + hash
    sizes = {os.stat(TEST_POS).st_size: SHA_POS, os.stat(TEST_NEG).st_size: SHA_NEG}
    copies = []
    for dp, dn, fn in os.walk(P1):
        dn[:] = [d for d in dn if d not in (".env", ".git")]
        for f in fn:
            fp = os.path.join(dp, f)
            try:
                if os.lstat(fp).st_size in sizes and sha256(fp) in sizes.values():
                    copies.append(fp)
            except OSError:
                pass
    return dict(
        n_revisions_searched=len(revs), strict_pattern=strict_pat,
        strict_hits_history=strict_classified,
        strict_hits_history_real=[h for h in strict_classified if not h["innocent_test_function_name"]],
        broad_pattern=BROAD, broad_hit_files_history=len(broad_files), broad_categories=dict(cat),
        broad_files_outside_expected_dirs=other,
        working_tree_text_files_scanned=n_scanned,
        working_tree_strict_hits=wt,
        working_tree_strict_hits_real=[h for h in wt if not h["innocent_test_function_name"]],
        files_in_part1_identical_to_test_files=copies,
    )


# --------------------------------------------------------------------------------------------- mapping helper
def uniprot_map(entrez_cols, hgnc_path):
    h = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=["symbol", "status", "entrez_id", "uniprot_ids"])
    h = h[h.status == "Approved"][["entrez_id", "uniprot_ids"]].dropna()
    h["u"] = h.uniprot_ids.str.split("|")
    h = h.explode("u")
    cand = h.groupby("u").entrez_id.apply(lambda s: sorted(set(s)))
    dep = {e: i for i, e in enumerate(entrez_cols)}
    status, gidx = {}, {}
    for u, es in cand.items():
        inn = [e for e in es if e in dep]
        if len(inn) == 1:
            status[u], gidx[u] = "mapped", dep[inn[0]]
        elif len(inn) == 0:
            status[u] = "not_in_depmap"
        else:
            status[u] = "ambiguous"
    return status, gidx


# --------------------------------------------------------------------------------------------- 3/4. pairs
def load_all():
    tp_raw, tn_raw = read_pairs(TEST_POS), read_pairs(TEST_NEG)
    rp_raw = read_pairs(P1 / "datasets/HuMAP/Train_positive_pairs.txt")
    rn_raw = read_pairs(P1 / "datasets/HuMAP/Train_negative_pairs.txt")
    return tp_raw, tn_raw, rp_raw, rn_raw


def descriptives(tp_raw, tn_raw, rp_raw, rn_raw):
    def stats(raw):
        c = [canon(a, b) for a, b in raw]
        cs = Counter(c)
        return dict(lines=len(raw), unique_unordered=len(cs), duplicate_extra_lines=len(raw) - len(cs),
                    self_pairs=sum(1 for a, b in raw if a == b)), set(k for k in cs if k[0] != k[1])
    stp, TP = stats(tp_raw)
    stn, TN = stats(tn_raw)
    str_, RP = stats(rp_raw)
    strn, RN = stats(rn_raw)
    both = TP & TN
    trprot_pos = {x for p in RP for x in p}
    trprot = trprot_pos | {x for p in RN for x in p}
    teprot_pos = {x for p in TP for x in p}
    teprot = teprot_pos | {x for p in TN for x in p}
    d = dict(
        test_pos=stp, test_neg=stn, part1_train_pos=str_, part1_train_neg=strn,
        test_pos_and_neg_both=len(both),
        test_pos_in_train_pos=len(TP & RP), test_pos_in_train_neg=len(TP & RN),
        test_neg_in_train_pos=len(TN & RP), test_neg_in_train_neg=len(TN & RN),
        train_pos_and_neg_both=len(RP & RN),
        test_pos_only=len(TP - RP - RN), test_pos_not_in_train_pos=len(TP - RP),
        test_neg_only=len(TN - RP - RN),
        labels_agree_test_pos_in_train_pos=len(TP & RP), labels_disagree_test_pos_in_train_neg=len(TP & RN),
        labels_disagree_test_neg_in_train_pos=len(TN & RP), labels_agree_test_neg_in_train_neg=len(TN & RN),
        test_prevalence=round(len(TP) / (len(TP) + len(TN)), 5),
        train_prevalence=round(len(RP) / (len(RP) + len(RN)), 5),
        test_proteins=len(teprot), test_pos_proteins=len(teprot_pos),
        test_proteins_in_train_proteins=len(teprot & trprot), test_proteins_not_in_train=len(teprot - trprot),
        test_proteins_in_train_pos_proteins=len(teprot & trprot_pos),
        train_proteins=len(trprot), train_pos_proteins=len(trprot_pos),
        test_pos_pairs_C1_both_proteins_in_train_pos=sum(1 for a, b in TP if a in trprot_pos and b in trprot_pos),
        test_pos_pairs_C2_one=sum(1 for a, b in TP if (a in trprot_pos) != (b in trprot_pos)),
        test_pos_pairs_C3_none=sum(1 for a, b in TP if a not in trprot_pos and b not in trprot_pos),
    )
    return d, TP, TN, RP, RN, teprot, trprot


def depmap_mapping(TP, TN, RP, RN, teprot, trprot):
    hdr = next(csv.reader(open(ROOT / "datasets/DepMAP/CRISPRGeneEffect.csv")))[1:]
    ent = [re.match(r"^.*\((\d+)\)$", c).group(1) for c in hdr]
    status, gidx = uniprot_map(ent, ROOT / "datasets/HGNC/hgnc_complete_set.txt")
    st = lambda u: status.get(u, "not_in_hgnc")
    mapped = lambda u: st(u) == "mapped"
    def scorable(p):
        a, b = p
        return mapped(a) and mapped(b) and gidx[a] != gidx[b]
    res = dict(
        test_proteins_mapped=sum(mapped(u) for u in teprot), test_proteins=len(teprot),
        test_proteins_status=dict(Counter(st(u) for u in teprot)),
        train_proteins_mapped=sum(mapped(u) for u in trprot), train_proteins=len(trprot),
        test_pos_scorable=sum(scorable(p) for p in TP), test_pos=len(TP),
        test_neg_scorable=sum(scorable(p) for p in TN), test_neg=len(TN),
        train_pos_scorable=sum(scorable(p) for p in RP), train_neg_scorable=sum(scorable(p) for p in RN),
    )
    # self-consistency of this mapping with Part 1's own uniprot_map.csv (train proteins)
    m1 = pd.read_csv(P1 / "data_cache/uniprot_map.csv", dtype=str).set_index("uniprot")
    comm = [u for u in m1.index if u in trprot]
    agree = sum((st(u) == m1.loc[u, "status"]) and (not mapped(u) or str(gidx[u]) == str(m1.loc[u, "gene_idx"])) for u in comm)
    res["mapping_agrees_with_part1_uniprot_map_csv"] = dict(n=len(comm), agree=int(agree))
    return res, status, gidx, ent


def pairs_universe_check(TP, TN, RP, RN, gidx, teprot, trprot):
    """Part 1 label/split artefacts in UniProt space and gene-index space."""
    res = {}
    tonly_p = TP - RP - RN
    tonly_n = TN - RP - RN
    # data_cache/pairs.pkl : full labelled pair table of Part 1
    pk = pd.read_pickle(P1 / "data_cache/pairs.pkl")
    pk_pairs = {canon(a, b): int(l) for a, b, l in zip(pk.a, pk.b, pk.label)}
    res["pairs_pkl"] = dict(rows=len(pk), equals_train_pos_union_neg=(set(pk_pairs) == (RP | RN)),
                            pos_rows=int((pk.label == 1).sum()),
                            test_only_pos_in_pkl=len(tonly_p & set(pk_pairs)), test_only_neg_in_pkl=len(tonly_n & set(pk_pairs)),
                            test_pos_in_pkl_as_pos=sum(1 for p in TP if pk_pairs.get(p) == 1),
                            test_pos_in_pkl_as_neg=sum(1 for p in TP if pk_pairs.get(p) == 0),
                            test_neg_in_pkl_as_pos=sum(1 for p in TN if pk_pairs.get(p) == 1),
                            test_neg_in_pkl_as_neg=sum(1 for p in TN if pk_pairs.get(p) == 0))
    # protein-level files
    prot_files = {}
    for rel in ("data_cache/splits/protein_assignment.csv", "results/splits/protein_assignment.csv", "data_cache/uniprot_map.csv"):
        col = "protein" if "assignment" in rel else "uniprot"
        s = set(pd.read_csv(P1 / rel, dtype=str)[col])
        prot_files[rel] = dict(n=len(s), test_only_proteins_present=len(s & (teprot - trprot)), train_proteins_covered=len(s & trprot))
    res["protein_level_files"] = prot_files
    # gene-index pair universe caches (ia, ib)
    N = 20000
    def enc(a, b):
        a, b = np.minimum(a, b), np.maximum(a, b)
        return a.astype(np.int64) * N + b
    def gene_keys(S):
        k = [enc(np.array([gidx[a]]), np.array([gidx[b]]))[0] for a, b in S if a in gidx and b in gidx and gidx[a] != gidx[b]]
        return np.unique(np.array(k, dtype=np.int64))
    kp, kn = gene_keys(tonly_p), gene_keys(tonly_n)
    cache = {}
    for rel in sorted(Path(P1 / "data_cache").rglob("*.pairs.npz")):
        z = np.load(rel, allow_pickle=False)
        keys = np.unique(enc(z["ia"], z["ib"]))
        key_pk = np.unique(enc(pk.ia.values, pk.ib.values))
        cache[str(rel.relative_to(P1))] = dict(
            n_pairs=len(keys), equals_pairs_pkl_gene_pairs=bool(np.array_equal(keys, np.unique(enc(pk[pk.scorable].ia.values, pk[pk.scorable].ib.values)))),
            test_only_pos_gene_pair_present=int(np.isin(kp, keys).sum()), test_only_neg_gene_pair_present=int(np.isin(kn, keys).sum()))
    res["pair_universe_npz_gene_space"] = cache
    res["test_only_pos_with_scorable_gene_pair"] = len(kp)
    res["test_only_neg_with_scorable_gene_pair"] = len(kn)
    return res


def text_scan(tonly_p, tonly_n, TP, TN):
    """Scan every text-like file in Part 1 (working tree + all historic blobs) for test-only pairs co-occurring on one
    line (or as a 2-element JSON list). Test pairs that ARE in the train files are excluded by construction."""
    lab = {}
    adj = defaultdict(set)
    for p in tonly_p:
        lab[p] = "pos"
        adj[p[0]].add(p[1]); adj[p[1]].add(p[0])
    for p in tonly_n:
        lab.setdefault(p, "neg")
        adj[p[0]].add(p[1]); adj[p[1]].add(p[0])
    hits = defaultdict(lambda: Counter())
    n_files = 0
    files_with_acc = 0

    def scan_lines(name, lines):
        nonlocal n_files, files_with_acc
        n_files += 1
        seen_acc = False
        for ln in lines:
            accs = ACC.findall(ln)
            if not accs:
                continue
            seen_acc = True
            if 2 <= len(accs) <= 12:
                s = set(accs)
                for a in s:
                    pr = adj.get(a)
                    if pr:
                        for b in pr & s:
                            if a < b:
                                hits[name][lab[(a, b)]] += 1
        files_with_acc += seen_acc

    SKIP = (".npy", ".npz", ".png", ".tif", ".pdf", ".docx", ".parquet", ".pkl", ".DS_Store", ".pyc", ".gif", ".jpg", ".xlsx")
    for dp, dn, fn in os.walk(P1):
        dn[:] = [d for d in dn if d not in (".env", ".git", ".pytest_cache", "__pycache__")]
        for f in fn:
            if f.endswith(SKIP):
                continue
            fp = Path(dp) / f
            rel = "WT:" + str(fp.relative_to(P1))
            try:
                if fp.stat().st_size > 400_000_000:
                    continue
                opener = gzip.open if f.endswith(".gz") else open
                with opener(fp, "rt", errors="replace") as fh:
                    scan_lines(rel, fh)
            except (OSError, EOFError):
                continue
    # parquet
    pq = []
    for fp in P1.rglob("*.parquet"):
        if ".env" in fp.parts:
            continue
        df = pd.read_parquet(fp)
        pq.append(dict(path=str(fp.relative_to(P1)), columns=list(df.columns), rows=len(df)))
    # historic blobs
    objs = git(CLONE, "rev-list", "--objects", "--all").splitlines()
    path_of = {}
    for ln in objs:
        p = ln.split(" ", 1)
        if len(p) == 2:
            path_of[p[0]] = p[1]
    inp = "\n".join(path_of).encode()
    r = subprocess.run(["git", "-C", str(CLONE), "cat-file", "--batch-check"], input=inp, capture_output=True, env=ENV)
    blobs = []
    for ln in r.stdout.decode().splitlines():
        sha, typ, size = ln.split()
        if typ == "blob" and int(size) < 60_000_000 and not path_of[sha].endswith(SKIP):
            blobs.append(sha)
    proc = subprocess.Popen(["git", "-C", str(CLONE), "cat-file", "--batch"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, env=ENV)
    nb = 0
    for sha in blobs:
        proc.stdin.write((sha + "\n").encode()); proc.stdin.flush()
        hdr = proc.stdout.readline().decode().split()
        size = int(hdr[2])
        data = proc.stdout.read(size); proc.stdout.read(1)
        nb += 1
        if path_of[sha].endswith(".gz"):
            try:
                data = gzip.decompress(data)
            except Exception:
                continue
        scan_lines("GIT:" + path_of[sha], io.StringIO(data.decode(errors="replace")))
    proc.stdin.close(); proc.wait()
    return dict(n_test_only_pos_pairs=len(tonly_p), n_test_only_neg_pairs=len(tonly_n),
                working_tree_and_history_text_units_scanned=n_files, historic_blobs_scanned=nb,
                units_containing_any_uniprot_accession=files_with_acc,
                # presence only (booleans): per-file counts of test pairs inside label-free score listings would amount to a
                # recall/precision-at-K on test pairs, which the hard rule forbids
                units_with_test_only_pair_on_one_line={k: {lab_: True for lab_ in v} for k, v in hits.items()},
                parquet_files_seen=pq)


# --------------------------------------------------------------------------------------------- code (label-free)
def check_scoring_code():
    files = ["mvp/src/coessential/method.py", "mvp/scripts/predict_all_pairs.py", "scripts/predict_pairs.py", "src/predict.py",
             "mvp/src/coessential/all_pairs.py", "mvp/src/coessential/data.py"]
    res = {}
    label_rx = re.compile(r"HuMAP|humap|Train_|POS_FILE|NEG_FILE|pairs\.pkl|read_pairs|load_pairs|\bsplits?\b", re.I)
    for f in files:
        txt = (CLONE / f).read_text()
        tree = ast.parse(txt)
        imps = sorted({(n.module or "") if isinstance(n, ast.ImportFrom) else a.name
                       for n in ast.walk(tree) for a in (n.names if isinstance(n, (ast.Import, ast.ImportFrom)) else [None]) if a is not None or isinstance(n, ast.ImportFrom)})
        # code lines (not docstring/comment) that touch label files
        code_hits = []
        for i, ln in enumerate(txt.splitlines(), 1):
            s = ln.strip()
            if s.startswith("#") or not label_rx.search(ln):
                continue
            code_hits.append((i, s[:140]))
        data_calls = sorted(set(re.findall(r"\bdata\.([a-zA-Z_]+)", txt))) if f == "src/predict.py" else None
        res[f] = dict(imports=imps, label_file_token_lines=code_hits, src_data_functions_used=data_calls,
                      docstring_declares_label_free=("label-free" in txt.lower() or "no pair label" in txt.lower() or "no pair labels" in txt.lower()))
    # src/data.py functions used by predict.py must not touch labels
    dtxt = (CLONE / "src/data.py").read_text()
    funcs = {n.name: ast.get_source_segment(dtxt, n) for n in ast.parse(dtxt).body if isinstance(n, ast.FunctionDef)}
    used = {}
    for fn in ("load_matrix", "build_matrix_cache", "build_uniprot_map"):
        seg = funcs.get(fn, "")
        used[fn] = dict(refs_label_files=bool(re.search(r"POS_FILE|NEG_FILE|HuMAP|Train_", seg)))
    res["src.data functions used by predict.py"] = used
    cfg = (CLONE / "src/config.py").read_text()
    res["src_config_label_file_constants"] = [ln.strip()[:140] for ln in cfg.splitlines() if "HuMAP" in ln]
    # run records of label-free scoring runs
    rr = {}
    for f in sorted((CLONE / "results/predict").glob("*.run.json")):
        d = json.loads(f.read_text())
        rr[f.name] = dict(label_free=d.get("label_free"), input=d.get("input") or d.get("pairs") or None)
    res["results_predict_run_records"] = rr
    # AST guard: only final_eval imports IFT loader (test exists); IFT looks log
    looks = [json.loads(l) for l in open(CLONE / "results/ift_looks.jsonl")]
    res["ift_looks"] = [dict(look=l["look"], what=l["what"][:110], status=l.get("status")) for l in looks]
    res["n_trial_records"] = sum(1 for _ in open(CLONE / "results/trials.jsonl"))
    return res


def check_hparams():
    tri = [json.loads(l) for l in open(CLONE / "results/trials.jsonl")]
    c = Counter((t.get("universe"), t.get("scheme")) for t in tri if t.get("event") == "start")
    ift_dir = sorted(p.name for p in (CLONE / "results/ift").glob("*"))
    return dict(trials_started_by_universe_scheme={f"{k[0]}|{k[1]}": v for k, v in c.items()},
                ift_result_files=ift_dir,
                design_docs={f: (CLONE / f).exists() for f in ("DESIGN.md", "RESEARCH_TREE.md", "DATA.md", "PROMPT.md")},
                data_md_withheld_statement=[l.strip()[:120] for l in (CLONE / "DATA.md").read_text().splitlines() if "withheld" in l.lower()],
                prompt_statement=[s.strip()[:260] for s in re.findall(r"The huMAP Test[^.]*\.[^.]*\.", (CLONE / "PROMPT.md").read_text())])


def check_other_channels():
    rx = re.compile(r"withheld|test set|hu\.?MAP test|test pairs|Test pairs", re.I)
    hits = []
    for f in ("LITERATURE.md", "DESIGN.md", "REPORT.md", "DATA.md", "docs/AGENTS.md"):
        for i, ln in enumerate((CLONE / f).read_text().splitlines(), 1):
            if rx.search(ln):
                hits.append(dict(file=f, line=i, text=ln.strip()[:200]))
    pub = [h for h in hits if h["file"] == "LITERATURE.md"]
    return dict(doc_mentions=hits, literature_mentions=pub)


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    R = {"generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
         "hard_rule": "no score/metric/prediction computed on test pairs; test lists used for set membership only"}
    R["check1_timeline"] = check_timeline()
    R["check2_references"] = check_references()
    tp_raw, tn_raw, rp_raw, rn_raw = load_all()
    desc, TP, TN, RP, RN, teprot, trprot = descriptives(tp_raw, tn_raw, rp_raw, rn_raw)
    mapres, status, gidx, ent = depmap_mapping(TP, TN, RP, RN, teprot, trprot)
    R["check4_overlap_descriptives"] = desc | dict(depmap_mapping=mapres)
    # integrity of Part 1 train files vs recorded hashes in DATA.md
    dm = (P1 / "DATA.md").read_text()
    R["part1_train_file_hashes_match_DATA_md"] = {
        f: (sha256(P1 / "datasets/HuMAP" / f) in dm) for f in ("Train_positive_pairs.txt", "Train_negative_pairs.txt")}
    R["new_project_inputs_match_part1_DATA_md"] = {
        "CRISPRGeneEffect.csv": sha256(ROOT / "datasets/DepMAP/CRISPRGeneEffect.csv") in dm,
        "hgnc_complete_set.txt": sha256(ROOT / "datasets/HGNC/hgnc_complete_set.txt") in dm}
    tonly_p, tonly_n = TP - RP - RN, TN - RP - RN
    R["check3_content_leakage"] = dict(
        label_artefact_checks=pairs_universe_check(TP, TN, RP, RN, gidx, teprot, trprot),
        text_scan=text_scan(tonly_p, tonly_n, TP, TN),
        scoring_code=check_scoring_code())
    R["check5_hyperparameter_selection"] = check_hparams()
    R["check6_other_channels"] = check_other_channels()
    OUT.write_text(json.dumps(R, indent=1, default=str))
    # concise stdout summary
    t = R["check1_timeline"]
    summ = dict(
        timeline=dict(test_birth=t["test_files"]["test_positive"]["birth_utc"], part1_last_commit=t["part1_last_commit"],
                      gap_min=t["gap_last_commit_to_test_birth_minutes"], commits_after=t["commits_after_test_birth"],
                      p1_files_modified_after_last_commit=len(t["part1_files_mtime_after_last_commit_excl_DS_Store"])),
        references_real_hits_history=len(R["check2_references"]["strict_hits_history_real"]),
        references_real_hits_worktree=len(R["check2_references"]["working_tree_strict_hits_real"]),
        copies_of_test_files=R["check2_references"]["files_in_part1_identical_to_test_files"],
        overlap=dict(test_pos=desc["test_pos"]["unique_unordered"], test_neg=desc["test_neg"]["unique_unordered"],
                     pos_in_train_pos=desc["test_pos_in_train_pos"], pos_in_train_neg=desc["test_pos_in_train_neg"],
                     neg_in_train_pos=desc["test_neg_in_train_pos"], neg_in_train_neg=desc["test_neg_in_train_neg"],
                     test_proteins=desc["test_proteins"], in_train=desc["test_proteins_in_train_proteins"],
                     mapped=mapres["test_proteins_mapped"]),
        text_scan_hits=R["check3_content_leakage"]["text_scan"]["units_with_test_only_pair_on_one_line"])
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
