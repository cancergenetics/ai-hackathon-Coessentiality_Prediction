"""Verifier A: independent re-implementation of the TEST_PROTOCOL v1.2 evaluation (no import of src/testeval or
scripts/21_*). Inputs: the two test label files, the label-free score TSVs, fitted_rules.json, the train positive
file (for C1-C3 and S2), and Part 1's build_clusters (allowed by §8.3). Part 1's metric functions and decide() are
imported ONLY as a cross-check of this file's own implementation.

Run:  PYTHONHASHSEED=0 .venv/bin/python results/test_eval/verify/verify_independent.py
Writes results/test_eval/verify/verifier_results.json
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

if os.environ.get("PYTHONHASHSEED") != "0":
    raise SystemExit("run with PYTHONHASHSEED=0 (Louvain)")

ROOT = Path(__file__).resolve().parents[3]
P1 = ROOT / "external" / "part1_locked"
sys.path.insert(0, str(P1))
from src.clusters import build_clusters  # noqa: E402  (Part 1, imported unchanged; allowed by §8.3)

OUT = Path(os.environ.get("VERIFY_OUT", ROOT / "results" / "test_eval" / "verify" / "verifier_results.json"))
POS = ROOT / "datasets" / "huMAP" / "Test_positive_pairs.txt"
NEG = ROOT / "datasets" / "huMAP" / "Test_negative_pairs.txt"
TRAIN_POS = Path(os.environ.get("PART1_DIR", ROOT.parent / "coessentiality_protein complexes")) / "datasets" / "HuMAP" / "Train_positive_pairs.txt"
TRAIN_NEG = TRAIN_POS.with_name("Train_negative_pairs.txt")
SCORES = ROOT / "data_cache" / "test_eval" / "scores" / "test"
RULES = json.loads((ROOT / "results" / "test_eval" / "fitted_rules.json").read_text())
METHODS = ["zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed", "spec_k20_tinf_std"]
CANDS = ["zavg_BEF", "B"]
BASES = ["gls_signed", "spearman_signed", "pearson_signed"]
REF = "spec_k20_tinf_std"
B_BOOT, SEED = int(os.environ.get("VERIFY_B", 2000)), 0
PI_TRAIN = 3314 / 634779
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ------------------------------------------------------------------------------------------ R1-R5
def read_pairs(path):
    d = pd.read_csv(path, sep="\t", header=None, names=["x", "y"], dtype=str)
    a = np.where(d.x.values < d.y.values, d.x.values, d.y.values)
    b = np.where(d.x.values < d.y.values, d.y.values, d.x.values)
    return pd.DataFrame(dict(a=a, b=b))


pos_raw, neg_raw = read_pairs(POS), read_pairs(NEG)
audit = dict(pos_lines=len(pos_raw), neg_lines=len(neg_raw),
             pos_self=int((pos_raw.a == pos_raw.b).sum()), neg_self=int((neg_raw.a == neg_raw.b).sum()))
pos = pos_raw[pos_raw.a != pos_raw.b].drop_duplicates()
neg = neg_raw[neg_raw.a != neg_raw.b].drop_duplicates()
audit["pos_duplicates"] = int(len(pos_raw) - audit["pos_self"] - len(pos))
audit["neg_duplicates"] = int(len(neg_raw) - audit["neg_self"] - len(neg))
pk = set(zip(pos.a, pos.b))
conf = [k for k in zip(neg.a, neg.b) if k in pk]
audit["conflicts"] = len(conf)
conf = set(conf)
pos = pos[[k not in conf for k in zip(pos.a, pos.b)]]
neg = neg[[k not in conf for k in zip(neg.a, neg.b)]]
pairs = pd.concat([pos.assign(label=1), neg.assign(label=0)]).sort_values(["a", "b"]).reset_index(drop=True)
audit.update(itt_pos=int(pairs.label.sum()), itt_neg=int((pairs.label == 0).sum()),
             itt_prevalence=float(pairs.label.mean()))

# R4 / S2: train overlap (flag only)
tr = pd.concat([read_pairs(TRAIN_POS), read_pairs(TRAIN_NEG)])
trk = set(zip(tr.a, tr.b))
audit["test_pairs_in_train"] = int(sum(k in trk for k in zip(pairs.a, pairs.b)))
trp = read_pairs(TRAIN_POS)
train_pos_prot = set(trp.a) | set(trp.b)
del tr, trk
log("pairs", audit)

# ------------------------------------------------------------------------------- scores + ITT floor
y_all = pairs.label.values.astype(bool)
S, MAPPED, floor_info = {}, {}, {}
for m in METHODS:
    t = pd.read_csv(SCORES / f"{m}.tsv", sep="\t", dtype={"uniprot_a": str, "uniprot_b": str})
    t = t.rename(columns={"uniprot_a": "a", "uniprot_b": "b"})
    j = pairs[["a", "b"]].merge(t, on=["a", "b"], how="left", validate="one_to_one")
    assert len(t) == len(pairs) and j.score.notna().sum() + j.mapped.isna().sum() >= 0
    if j.mapped.isna().any():
        raise RuntimeError(f"{m}: {j.mapped.isna().sum()} analysed pairs missing from score file")
    sc = j.score.values.astype(np.float64)
    mp = (j.mapped.values.astype(int) == 1) & np.isfinite(sc)
    fl = float(sc[mp].min() - 1.0)           # floor over the scored file (= all analysed pairs here, no conflicts)
    s = np.where(mp, sc, fl)
    floor_info[m] = dict(floor=fl, n_mapped=int(mp.sum()), n_unmapped=int((~mp).sum()),
                         n_mapped_pos=int((mp & y_all).sum()), n_mapped_neg=int((mp & ~y_all).sum()),
                         max_abs_tsv_minus_floor_on_unmapped=float(np.nanmax(np.abs(sc[~mp] - fl))) if (~mp).any() else 0.0,
                         n_nan_in_tsv=int(np.isnan(sc).sum()))
    S[m], MAPPED[m] = s, mp
log("floors", {k: (v["floor"], v["n_mapped"]) for k, v in floor_info.items()})
s1_mask = np.logical_and.reduce([MAPPED[m] for m in METHODS])
audit["mapped_flags_identical_across_methods"] = bool(all((MAPPED[m] == MAPPED[METHODS[0]]).all() for m in METHODS))
audit["s1_pos"], audit["s1_neg"] = int((s1_mask & y_all).sum()), int((s1_mask & ~y_all).sum())


# ---------------------------------------------------------------------------------- own metric code
class Sorted:
    """Descending sort with tie blocks; AP / AUROC / P@R with multiplicity weights (sklearn tie semantics)."""

    def __init__(self, y, s):
        self.o = np.argsort(-s, kind="mergesort")
        self.ys = y[self.o]
        ss = s[self.o]
        self.be = np.flatnonzero(np.r_[ss[1:] != ss[:-1], True])

    def _cum(self, w):
        ws = np.ones(len(self.ys)) if w is None else w[self.o]
        tp = np.cumsum(np.where(self.ys, ws, 0.0))[self.be]
        fp = np.cumsum(np.where(self.ys, 0.0, ws))[self.be]
        return tp, fp

    def ap(self, w=None):
        tp, fp = self._cum(w)
        if tp[-1] <= 0:
            return np.nan
        prec = tp / np.maximum(tp + fp, 1e-300)
        return float(np.sum(np.diff(np.r_[0.0, tp / tp[-1]]) * prec))

    def auroc(self, w=None):
        tp, fp = self._cum(w)
        P, N = tp[-1], fp[-1]
        if P <= 0 or N <= 0:
            return np.nan
        tpb, fpb = np.diff(np.r_[0.0, tp]), np.diff(np.r_[0.0, fp])
        return float(np.sum(fpb * (tp - tpb + 0.5 * tpb)) / (P * N))

    def p_at_r(self, levels=(0.1, 0.25, 0.5)):
        tp, fp = self._cum(None)
        rec, prec = tp / tp[-1], tp / (tp + fp)
        return {f"P@R{lv:g}": float(prec[min(np.searchsorted(rec, lv - 1e-12), len(prec) - 1)]) for lv in levels}


def ap_subset(y, s, keep):
    if keep.sum() == 0 or y[keep].sum() == 0:
        return np.nan
    return Sorted(y[keep], s[keep]).ap()


class Macro:
    """Macro per-cluster AP: clusters whose touching pairs (>=1 endpoint in the cluster) contain >=1 pos and >=1 neg
    (at unit weights). Replicate version: cluster c counts m_c times; its AP uses the replicate pair weights;
    clusters with m_c = 0 or no positive weight drop out. Segmented vectorised AP (integer weights: exact cumsums)."""

    def __init__(self, y, s, ca, cb):
        df = pd.DataFrame(dict(c=np.r_[ca, cb], i=np.r_[np.arange(len(ca)), np.arange(len(cb))])).drop_duplicates()
        segs, cids = [], []
        for c, g in df.groupby("c", sort=True):
            ix = g.i.values
            if y[ix].any() and (~y[ix]).any():
                ix = np.sort(ix)
                ix = ix[np.argsort(-s[ix], kind="mergesort")]
                segs.append(ix)
                cids.append(c)
        self.cids = np.array(cids)
        self.idx = np.concatenate(segs)
        lens = np.array([len(x) for x in segs])
        self.seg = np.repeat(np.arange(len(segs)), lens)
        self.start = np.r_[0, np.cumsum(lens)[:-1]]
        ss = s[self.idx]
        last_in_seg = np.zeros(len(self.idx), bool)
        last_in_seg[np.cumsum(lens) - 1] = True
        self.be = np.flatnonzero(np.r_[ss[1:] != ss[:-1], True] | last_in_seg)
        self.ys = y[self.idx]
        self.nseg = len(segs)

    def value(self, w=None, mult=None):
        ws = np.ones(len(self.idx)) if w is None else w[self.idx]
        ctp = np.cumsum(np.where(self.ys, ws, 0.0))
        cfp = np.cumsum(np.where(self.ys, 0.0, ws))
        off_tp = np.r_[0.0, ctp][self.start][self.seg]
        off_fp = np.r_[0.0, cfp][self.start][self.seg]
        be = self.be
        tp, fp, sg = ctp[be] - off_tp[be], cfp[be] - off_fp[be], self.seg[be]
        P = np.zeros(self.nseg)
        P[sg] = tp                                   # last block of each segment wins (ascending order)
        prev_tp = np.r_[0.0, tp[:-1]]
        first = np.r_[True, sg[1:] != sg[:-1]]
        prev_tp[first] = 0.0
        prec = tp / np.maximum(tp + fp, 1e-300)
        Pb = P[sg]
        contrib = np.where(Pb > 0, (tp - prev_tp) / np.where(Pb > 0, Pb, 1.0) * prec, 0.0)
        aps = np.bincount(sg, weights=contrib, minlength=self.nseg)
        ok = P > 0
        mc = np.ones(self.nseg) if mult is None else mult
        ok &= mc > 0
        if not ok.any():
            return np.nan, 0
        return float(np.sum(mc[ok] * aps[ok]) / np.sum(mc[ok])), int(ok.sum())


def prf(y, s, t, w=None):
    w = np.ones(len(y)) if w is None else w
    pp = s >= t
    tp = float(np.sum(w[pp & y]))
    fp = float(np.sum(w[pp & ~y]))
    P_ = float(np.sum(w[y]))
    p = tp / (tp + fp) if tp + fp > 0 else np.nan
    r = tp / P_ if P_ > 0 else np.nan
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
    return dict(threshold=float(t), tp=tp, fp=fp, n_pred_pos=tp + fp, precision=p, recall=r, f1=f1)


def iso_prob(m, s_eval, mapped):
    r = RULES["methods"][m]
    x, yk = np.asarray(r["isotonic"]["x"], float), np.asarray(r["isotonic"]["y"], float)
    s_cal = np.where(mapped, s_eval, r["fit_floor"])
    return np.interp(s_cal, x, yk)


def platt_prob(m, s_eval, mapped):
    r = RULES["methods"][m]
    pl = r["platt"]
    s_cal = np.where(mapped, s_eval, r["fit_floor"])
    z = pl["a"] + pl["b"] * (s_cal - pl["mu"]) / pl["sd"]
    return 1.0 / (1.0 + np.exp(-z))


def logit_fit(y, x, iters=100):
    """MLE logistic regression y ~ a + b x (Newton)."""
    X = np.c_[np.ones(len(x)), x]
    beta = np.zeros(2)
    for _ in range(iters):
        p = 1 / (1 + np.exp(-(X @ beta)))
        g = X.T @ (y - p)
        H = (X * (p * (1 - p))[:, None]).T @ X
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-12:
            break
    return float(beta[0]), float(beta[1])


def calib(y, p):
    yf = y.astype(float)
    n = len(y)
    brier = float(np.mean((p - yf) ** 2))
    prev = yf.mean()
    o = np.argsort(p, kind="stable")
    ece_m = sum(abs(p[b].sum() - yf[b].sum()) for b in np.array_split(o, 20)) / n
    bins = np.minimum((p * 15).astype(int), 14)
    ece_w = sum(abs(p[bins == k].sum() - yf[bins == k].sum()) for k in range(15)) / n
    pc = np.clip(p, 1e-6, 1 - 1e-6)
    a, b = logit_fit(yf, np.log(pc / (1 - pc)))
    return dict(brier=brier, bss_train_prev=1 - brier / float(np.mean((PI_TRAIN - yf) ** 2)),
                bss_test_prev=1 - brier / float(np.mean((prev - yf) ** 2)),
                citl=float(p.sum() / yf.sum()), cal_intercept=a, cal_slope=b,
                ece_equal_mass20=float(ece_m), ece_equal_width15=float(ece_w))


def holm(pv):
    items = sorted(pv.items(), key=lambda kv: kv[1])
    out, run = {}, 0.0
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (len(items) - i) * p))
        out[k] = run
    return out


def ci(x):
    lo, hi = np.nanquantile(x, [0.025, 0.975])
    return float(lo), float(hi)


def pval(d):
    return float(min(1.0, 2 * min(np.mean(d <= 0), np.mean(d >= 0))))


# ------------------------------------------------------------------------------------ one analysed set
def analyse(name, mask, do_boot=True):
    log(f"== {name}")
    idx = np.flatnonzero(mask)
    y = y_all[idx]
    a, b = pairs.a.values[idx], pairs.b.values[idx]
    cl = build_clusters(pd.DataFrame(dict(a=a, b=b, label=y.astype(int))))
    cmap, comp = dict(zip(cl.protein, cl.cluster)), dict(zip(cl.protein, cl.component))
    ca, cb = pd.Series(a).map(cmap).values, pd.Series(b).map(cmap).values
    pa, pb = pd.Series(a).map(comp).values, pd.Series(b).map(comp).values
    csz = cl.groupby("component").size()
    largest = int(csz[csz == csz.max()].index.min())
    keep_nl = (pa != largest) & (pb != largest)
    cu, cn = np.unique(ca[y & (ca == cb)], return_counts=True)
    top10 = cu[np.argsort(-cn, kind="stable")][:10]
    keep_nt = ~(np.isin(ca, top10) | np.isin(cb, top10))
    hard = y | (pa == pb)
    sizes = cl.groupby("cluster").size()
    info = dict(n_pos=int(y.sum()), n_neg=int((~y).sum()), n_proteins=len(cl), n_components=int(cl.component.nunique()),
                n_clusters=int(cl.cluster.nunique()), n_singleton_clusters=int((sizes == 1).sum()),
                largest_component=largest, largest_component_size=int(csz.max()),
                components_split=int(cl.loc[cl.louvain_split, "component"].nunique()),
                n_pos_no_largest=int(y[keep_nl].sum()), top10_clusters=[int(x) for x in top10],
                n_pos_between_clusters=int((y & (ca != cb)).sum()))
    log(info)
    per, srt, mac, probs = {}, {}, {}, {}
    for m in METHODS:
        s = S[m][idx]
        mp = MAPPED[m][idx]
        so = Sorted(y, s)
        srt[m] = so
        mac[m] = Macro(y, s, ca, cb)
        mv, nmac = mac[m].value()
        r = dict(ap=so.ap(), auroc=so.auroc(), macro_cluster_ap=mv, n_clusters_macro=nmac,
                 ap_no_largest_component=ap_subset(y, s, keep_nl), ap_no_top10_clusters=ap_subset(y, s, keep_nt),
                 ap_hard_neg=ap_subset(y, s, hard), n_hard_neg=int((hard & ~y).sum()))
        r.update(so.p_at_r())
        r["rule_b"] = prf(y, s, RULES["methods"][m]["f1max"]["threshold"])
        K = int(round(PI_TRAIN * len(y)))
        r["rule_c"] = dict(prf(y, s, np.sort(s)[::-1][K - 1]), K=K)
        if m == "zavg_BEF":
            r["rule_a_nonblind_D1"] = prf(y, s, 3.0)
        p_iso = iso_prob(m, s, mp)
        probs[m] = p_iso
        r["isotonic"] = calib(y, p_iso)
        r["platt"] = calib(y, platt_prob(m, s, mp))
        # strata C1-C3 (point AP)
        nt = np.array([x in train_pos_prot for x in a], int) + np.array([x in train_pos_prot for x in b], int)
        for k, lab in ((2, "C1"), (1, "C2"), (0, "C3")):
            kk = nt == k
            r[f"ap_{lab}"] = ap_subset(y, s, kk)
            r[f"n_pos_{lab}"], r[f"n_neg_{lab}"] = int(y[kk].sum()), int((~y[kk]).sum())
        per[m] = r
        log(m, {k: r[k] for k in ("ap", "auroc", "macro_cluster_ap", "ap_no_largest_component")}, "F1b", r["rule_b"]["f1"])
    # Part 1 cross-check of point metrics (w = 1)
    from src import metrics as p1m
    xchk = {}
    for m in METHODS:
        s = S[m][idx]
        su = p1m.metric_suite(y, s, np.ones(len(y)), ca, cb, pa, pb, largest_component=largest)
        xchk[m] = {k: float(abs(su[k2] - per[m][k])) for k, k2 in (("ap", "ap_ht"), ("auroc", "auroc_ht"),
                   ("macro_cluster_ap", "macro_cluster_ap"), ("ap_no_largest_component", "ap_ht_no_largest_component"),
                   ("ap_hard_neg", "ap_ht_hard_neg"))}
        xchk[m]["P@R"] = max(abs(su[f"ht_P@R{lv}"] - per[m][f"P@R{lv}"]) for lv in ("0.1", "0.25", "0.5"))
    log("part1 cross-check max", max(max(v.values()) for v in xchk.values()))
    out = dict(info=info, metrics=per, part1_metric_suite_abs_diff=xchk)
    if not do_boot:
        return out, None

    # ----------------------------------------------------------------- paired cluster bootstrap, seed 0
    uniq, inv = np.unique(np.r_[ca, cb], return_inverse=True)
    ia, ib = inv[:len(ca)], inv[len(ca):]
    same = ia == ib
    C = len(uniq)
    # macro: position of each Macro segment cluster in uniq
    rng = np.random.default_rng(SEED)
    stats = ("ap", "auroc", "ap_no_largest", "macro", "f1_b", "prec_b", "rec_b", "brier_iso")
    reps = {m: {k: np.empty(B_BOOT) for k in stats} for m in METHODS}
    srt_nl = {m: Sorted(y[keep_nl], S[m][idx][keep_nl]) for m in METHODS}
    mac_pos = {m: np.searchsorted(uniq, mac[m].cids) for m in METHODS}
    yf = y.astype(float)
    tb = time.time()
    for r_ in range(B_BOOT):
        mm = np.bincount(rng.integers(0, C, C), minlength=C).astype(np.float64)
        w = np.where(same, mm[ia], mm[ia] * mm[ib])
        wnl = w[keep_nl]
        for m in METHODS:
            so = srt[m]
            reps[m]["ap"][r_] = so.ap(w)
            reps[m]["auroc"][r_] = so.auroc(w)
            reps[m]["ap_no_largest"][r_] = srt_nl[m].ap(wnl)
            reps[m]["macro"][r_] = mac[m].value(w, mm[mac_pos[m]])[0]
            q = prf(y, S[m][idx], RULES["methods"][m]["f1max"]["threshold"], w)
            reps[m]["f1_b"][r_], reps[m]["prec_b"][r_], reps[m]["rec_b"][r_] = q["f1"], q["precision"], q["recall"]
            reps[m]["brier_iso"][r_] = float(np.sum(w * (probs[m] - yf) ** 2) / np.sum(w))
        if r_ % 200 == 0:
            log(f"boot {r_}/{B_BOOT} ({time.time() - tb:.0f}s)")
    point = {m: dict(ap=per[m]["ap"], auroc=per[m]["auroc"], ap_no_largest=per[m]["ap_no_largest_component"],
                     macro=per[m]["macro_cluster_ap"], f1_b=per[m]["rule_b"]["f1"], prec_b=per[m]["rule_b"]["precision"],
                     rec_b=per[m]["rule_b"]["recall"], brier_iso=per[m]["isotonic"]["brier"]) for m in METHODS}
    boot = dict(n_boot=B_BOOT, seed=SEED, n_clusters=int(C), per_method={}, diff={}, ratio={})
    for m in METHODS:
        boot["per_method"][m] = {k: dict(point=point[m][k], ci=ci(reps[m][k]), boot_sd=float(np.nanstd(reps[m][k], ddof=1)))
                                 for k in stats}

    def dif(a_, b_, k):
        d = reps[a_][k] - reps[b_][k]
        rel = reps[a_][k] / reps[b_][k] - 1
        return dict(diff=point[a_][k] - point[b_][k], ci=ci(d), p_two_sided=pval(d),
                    rel_gain=point[a_][k] / point[b_][k] - 1, rel_ci=ci(rel))

    for c in CANDS:
        for bb in BASES:
            boot["diff"][f"ap:{c}-{bb}"] = dif(c, bb, "ap")
            for k in ("macro", "ap_no_largest"):
                boot["diff"][f"{k}:{c}-{bb}"] = dif(c, bb, k)
                rr = reps[c][k] / reps[bb][k]
                boot["ratio"][f"{k}:{c}/{bb}"] = dict(ratio=point[c][k] / point[bb][k], ci=ci(rr))
        for bb in BASES:
            if c == "zavg_BEF":
                boot["diff"][f"auroc:{c}-{bb}"] = dif(c, bb, "auroc")
    boot["diff"]["ap:B-zavg_BEF"] = dif("B", "zavg_BEF", "ap")
    boot["diff"][f"ap:{REF}-zavg_BEF"] = dif(REF, "zavg_BEF", "ap")

    # ------------------------------------------------------------- §7.1 locked verdict (own) + §7.2
    verdict = dict(pairs={}, candidates={})
    primary_ok = None
    for i, c in enumerate(CANDS):
        for bb in BASES:
            rel = per[c]["ap"] / per[bb]["ap"] - 1
            gate = bool(np.isfinite(rel) and rel >= 0.10)
            ni = {k: bool(np.isfinite(per[c][k]) and np.isfinite(per[bb][k]) and per[c][k] >= 0.95 * per[bb][k])
                  for k in ("macro_cluster_ap", "ap_no_largest_component")}
            d = boot["diff"][f"ap:{c}-{bb}"]
            verdict["pairs"][f"{c}-{bb}"] = dict(rel_gain=rel, gate_ii=gate, ni=ni, passes=bool(gate and all(ni.values())),
                                                 ci_excludes_0=bool(d["ci"][0] > 0), borderline=bool(abs(d["ci"][0]) <= 0.011))
        beats = all(verdict["pairs"][f"{c}-{bb}"]["passes"] for bb in BASES)
        if i == 0:
            primary_ok, status = beats, "confirmatory"
        else:
            status = "confirmatory" if primary_ok else "exploratory"
        verdict["candidates"][c] = dict(clearly_beats=beats, status=status, confirmed=bool(beats and status == "confirmatory"),
                                        ci_excludes_0_all=all(verdict["pairs"][f"{c}-{bb}"]["ci_excludes_0"] for bb in BASES))
    verdict["holm_zavg_dAP"] = holm({bb: boot["diff"][f"ap:zavg_BEF-{bb}"]["p_two_sided"] for bb in BASES})
    # cross-check against Part 1's decide(), called unchanged
    sys.path.insert(0, str(P1 / "scripts"))
    import final_eval as p1fe
    per1 = {m: dict(ap_ht=per[m]["ap"], macro_cluster_ap=per[m]["macro_cluster_ap"],
                    ap_ht_no_largest_component=per[m]["ap_no_largest_component"]) for m in METHODS}
    boot1 = dict(diff={})
    for c in CANDS:
        for bb in BASES:
            d = boot["diff"][f"ap:{c}-{bb}"]
            boot1["diff"][f"{c}-{bb}"] = dict(diff=d["diff"], rel_gain=d["rel_gain"], ci_low=d["ci"][0], ci_high=d["ci"][1],
                                              rel_ci_low=d["rel_ci"][0], rel_ci_high=d["rel_ci"][1], p_two_sided=d["p_two_sided"])
    dec = p1fe.decide(per1, boot1, CANDS, BASES)
    verdict["part1_decide"] = {k: v["passes"] for k, v in dec["pairs"].items()}
    verdict["part1_decide_candidates"] = {k: dict(clearly_beats=v["clearly_beats_ift"], status=v["status"],
                                                  confirmed=v["confirmed"]) for k, v in dec["candidates"].items()}
    verdict["own_equals_part1_decide"] = all(verdict["part1_decide"][k] == v["passes"] for k, v in verdict["pairs"].items())
    log(name, "verdict", {c: v for c, v in verdict["candidates"].items()})
    return out, dict(bootstrap=boot, verdict=verdict)


res = dict(verifier="A (independent; no testeval import; test_results.json not opened)",
           protocol="docs/TEST_PROTOCOL.md FROZEN v1.2", audit=audit, floors=floor_info, sets={})
for nm, mask in (("ITT", np.ones(len(pairs), bool)), ("S1", s1_mask)):
    o, bt = analyse(nm, mask)
    res["sets"][nm] = dict(**o, **(bt or {}))
res["runtime_s"] = round(time.time() - T0, 1)


def _j(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    raise TypeError(type(x))


OUT.write_text(json.dumps(res, indent=1, default=_j))
log("wrote", OUT)
