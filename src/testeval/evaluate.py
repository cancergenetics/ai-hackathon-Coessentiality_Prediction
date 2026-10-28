"""Core of the test evaluation (docs/TEST_PROTOCOL.md). The same functions run the dry run on development data and
the single test look, so the code path of the look is the one that was verified.

  load_set      labelled pair files + the locked predictor's score TSVs -> EvalSet (analysed pairs only)
  fit_rules     decision thresholds + calibrators from the FIT set (train pairs) only
  evaluate      every pre-specified metric, decision rule, calibration summary, bootstrap and the verdict
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from testeval import bootstrap as B
from testeval import calibration as C
from testeval import pairs as P
from testeval.metrics import Ranked, auroc, top_cluster_ids
from testeval import part1
from testeval.part1 import SEEDS, clusters as p1c, data as p1d, metrics as p1m, universe as p1u

PRIMARY = "zavg_BEF"
SECONDARY = "B"
BASELINES = ("gls_signed", "spearman_signed", "pearson_signed")
REFERENCE = "spec_k20_tinf_std"
METHODS = (PRIMARY, SECONDARY) + BASELINES + (REFERENCE,)
PREDICTOR_ARG = {"zavg_BEF": "zavg:B,E,F", "B": "single", "gls_signed": "gls", "spearman_signed": "spearman",
                 "pearson_signed": "pearson", "spec_k20_tinf_std": "scripts/20b_score_posthoc_reference.py"}
Z_RULE = {PRIMARY: 3.0}            # Part 1 release cutoff (results/predictions/README.md: top_pairs = zavg_BEF >= 3)
MARGIN, NI_RATIO = 0.10, 0.95      # Part 1 DESIGN §3.7 (ii) relative gain >= 10 %, (iii) non-inferiority >= 0.95 x
NI_METRICS = ("macro_cluster_ap", "ap_no_largest_component")
MC_SPREAD = 0.011                  # dry run B: max seed-to-seed change of the dAP CI endpoints -> 'borderline' band
N_BOOT = 2000
BOOT_SEED = SEEDS[0]               # primary bootstrap seed; SEEDS[1:] = Monte-Carlo stability check
N_TOP = 10
RECALL_LEVELS = (0.1, 0.25, 0.5)


@dataclass
class EvalSet:
    name: str
    table: pd.DataFrame                  # analysed pairs: a, b, label, flags
    scores: dict                         # method -> float64 array aligned to table
    mapped: dict                         # method -> bool array (predictor's 'mapped')
    w: np.ndarray                        # design weights (1 = natural prevalence; HT weights in the IFT dry run)
    audit: dict = field(default_factory=dict)

    @property
    def y(self):
        return (self.table.label.values == 1)

    @property
    def mapped_all(self):
        return np.logical_and.reduce([self.mapped[m] for m in self.scores])

    def subset(self, keep, name) -> "EvalSet":
        keep = np.asarray(keep, bool)
        return EvalSet(name, self.table[keep].reset_index(drop=True), {k: v[keep] for k, v in self.scores.items()},
                       {k: v[keep] for k, v in self.mapped.items()}, self.w[keep], dict(self.audit))


def load_set(name, pos, neg, score_dir, methods=METHODS, train_paths=None, weights_fn=None) -> EvalSet:
    t, audit = P.labelled_table(pos, neg, train_paths)
    t = t[t.analysed].reset_index(drop=True)
    sc, mp = {}, {}
    for m in methods:
        s, mm = P.attach_scores(t, P.read_scores(Path(score_dir) / f"{m}.tsv"), f"{name}/{m}")
        sc[m], mp[m] = s, mm
    w = np.ones(len(t)) if weights_fn is None else np.asarray(weights_fn(t, mp), np.float64)
    audit["mapped_counts"] = {m: dict(pos=int((mp[m] & (t.label.values == 1)).sum()),
                                      neg=int((mp[m] & (t.label.values == 0)).sum())) for m in methods}
    return EvalSet(name, t, sc, mp, w, audit)


# ------------------------------------------------------------------------------------------ annotation
def gene_index(ids, genes) -> pd.Series:
    """UniProt -> DepMap column index (-1 unmapped), Part 1 rule (data.build_uniprot_map)."""
    return p1d.build_uniprot_map(ids, genes).gene_idx


def annotate_clusters(es: EvalSet, assignment: pd.DataFrame | None = None) -> dict:
    """Per-pair cluster / component ids. Protocol (assignment None): Part 1 ``clusters.build_clusters`` on this set's
    analysed pairs (all proteins, incl. unmapped; positive graph; components > 40 proteins split by networkx Louvain,
    seed SEEDS[0], PYTHONHASHSEED=0). Reproduction: Part 1's committed assignment."""
    t = es.table
    if assignment is None:
        if os.environ.get("PYTHONHASHSEED") != "0":
            raise RuntimeError("cluster construction must run with PYTHONHASHSEED=0 (Part 1 Louvain tie-breaking)")
        cl = p1c.build_clusters(pd.DataFrame(dict(a=t.a.values, b=t.b.values, label=es.y.astype(int))))
        stats = p1c.cluster_edge_stats(pd.DataFrame(dict(a=t.a, b=t.b, label=es.y.astype(int))), cl)
    else:
        cl = assignment
        stats = dict(source="Part 1 committed assignment")
    c = cl.set_index("protein")
    ca, cb = c.cluster.reindex(t.a).values, c.cluster.reindex(t.b).values
    oa, ob = c.component.reindex(t.a).values, c.component.reindex(t.b).values
    largest = int(cl.loc[cl.component_size.idxmax(), "component"])
    wp = pd.Series(ca[es.y & (ca == cb)]).value_counts()
    stats.update(n_clusters_with_within_pos=int(len(wp)),
                 kish_neff_clusters=p1c.kish_neff(wp.values) if len(wp) else 0.0,
                 largest_component=largest,
                 largest_component_size=int(cl.component_size.max()),
                 largest_component_pos_pairs=int((es.y & (oa == largest) & (ob == largest)).sum()))
    return dict(ca=np.asarray(ca), cb=np.asarray(cb), comp_a=np.asarray(oa), comp_b=np.asarray(ob),
                largest=largest, table=cl, stats=stats)


def essential_flags(es: EvalSet, X, genes) -> tuple[np.ndarray, np.ndarray]:
    ess = p1u.common_essential(X)
    gi = gene_index(pd.unique(es.table[["a", "b"]].values.ravel()), genes)
    ia = gi.reindex(es.table.a).values.astype(np.int64)
    ib = gi.reindex(es.table.b).values.astype(np.int64)
    ea = np.where(ia >= 0, ess[np.maximum(ia, 0)], False)
    eb = np.where(ib >= 0, ess[np.maximum(ib, 0)], False)
    return ea, eb


# ------------------------------------------------------------------------------------------ fitting
def weighted_prevalence(y, w) -> float:
    return float(np.sum(w * y) / np.sum(w))


def fit_rules(fit: EvalSet) -> dict:
    """Per method, from the fit set only: unmapped floor of the fit file, F1-maximising threshold (rule b), isotonic
    and Platt calibrators; plus the fit-set prevalence used by the top-k rule (rule c) and Brier skill."""
    y, w = fit.y, fit.w
    out = dict(fit_set=fit.name, n_pos=int(y.sum()), n_neg=int((~y).sum()),
               fit_prevalence=weighted_prevalence(y, w), weights="design weights of the fit set", methods={})
    for m, s in fit.scores.items():
        mp = fit.mapped[m]
        floor = float(s[~mp].max()) if (~mp).any() else float(s.min() - 1.0)
        if (~mp).any() and not np.all(s[~mp] == floor):
            raise RuntimeError(f"{m}: unmapped pairs do not share one floor score")
        r = Ranked(y, s)
        f1 = r.f1_max(w)
        out["methods"][m] = dict(fit_floor=floor, f1max=f1, isotonic=C.fit_isotonic(s, y, w),
                                 platt=C.fit_platt(s, y, w), fit_ap=r.ap(w))
        if m in Z_RULE:
            out["methods"][m]["z_rule"] = dict(threshold=Z_RULE[m], fit=r.confusion(Z_RULE[m], w))
    return out


# ------------------------------------------------------------------------------------------ evaluation
def topk_threshold(r: Ranked, w, prevalence: float) -> tuple[float, float]:
    """Rule c: K = round(prevalence_fit x total design weight); threshold = score of the first tie block whose
    cumulative weight reaches K (w = 1: the K-th ranked score)."""
    ws = np.asarray(w, np.float64)[r.order]
    cw = np.cumsum(ws)[r.block_end]
    K = float(np.round(prevalence * ws.sum()))
    k = int(np.searchsorted(cw, K - 1e-9))
    return float(r.block_score[min(k, len(cw) - 1)]), K


def point_metrics(es: EvalSet, cl: dict, ess=None, rules: dict | None = None) -> dict:
    y, w = es.y, es.w
    out = {}
    top = top_cluster_ids(y, cl["ca"], cl["cb"], N_TOP)
    keep_top = ~(np.isin(cl["ca"], top) | np.isin(cl["cb"], top))
    for m, s in es.scores.items():
        suite = p1m.metric_suite(y, s, w, cl["ca"], cl["cb"], cl["comp_a"], cl["comp_b"],
                                 None if ess is None else ess[0], None if ess is None else ess[1],
                                 largest_component=cl["largest"])
        suite.pop("ap_reweighted", None)          # Part 1 train-universe reweighting: not defined for this set
        r = Ranked(y, s)
        d = dict(ap=suite["ap_ht"], auroc=auroc(y, s, w), ap_unweighted=suite["ap_raw"],
                 macro_cluster_ap=suite["macro_cluster_ap"], n_clusters_macro=suite["n_clusters_macro"],
                 ap_no_largest_component=suite.get("ap_ht_no_largest_component"),
                 ap_no_top10_clusters=p1m.weighted_ap(y[keep_top], s[keep_top], w[keep_top]),
                 n_pos_no_top10_clusters=int(y[keep_top].sum()),
                 ap_hard_neg=suite.get("ap_ht_hard_neg"),
                 **{f"P@R{lv:g}": suite[f"ht_P@R{lv:g}"] for lv in RECALL_LEVELS},
                 part1_metric_suite=suite)
        if ess is not None:
            d.update({k.replace("ap_ht_", "ap_"): v for k, v in suite.items() if k.startswith("ap_ht_both")
                      or k == "ap_ht_mixed"})
        if rules is not None:
            rm = rules["methods"][m]
            dec = dict(rule_b_f1max_dev=dict(threshold=rm["f1max"]["threshold"], **r.confusion(rm["f1max"]["threshold"], w)))
            t, K = topk_threshold(r, w, rules["fit_prevalence"])
            dec["rule_c_topk"] = dict(threshold=t, K=K, **r.confusion(t, w))
            if "z_rule" in rm:
                dec["rule_a_z3"] = dict(threshold=rm["z_rule"]["threshold"], **r.confusion(rm["z_rule"]["threshold"], w))
            d["decision_rules"] = dec
            cal = {}
            for kind in ("isotonic", "platt"):
                p = C.apply(rm[kind], C.floor_to_fit(s, es.mapped[m], rm["fit_floor"]))
                cal[kind] = dict(**C.summary(p, y, w, rules["fit_prevalence"]),
                                 reliability_mass20=C.reliability(p, y, w, 20, "mass"),
                                 reliability_log16=C.reliability(p, y, w, 16, "log"))
            d["calibration"] = cal
        out[m] = d
    return out


def _method_fn(y, s, w_unused, thr, keep_nl, cal_p):
    r = Ranked(y, s)
    r_nl = Ranked(y[keep_nl], s[keep_nl])
    yf = y.astype(np.float64)

    def f(w):
        sm = r.summary(w, thr)
        d = dict(ap=sm["ap"], auroc=sm["auroc"], ap_no_largest_component=r_nl.ap(w[keep_nl]))
        if thr is not None:
            c = sm["confusion"]
            d.update(f1_rule_b=c["f1"], precision_rule_b=c["precision"], recall_rule_b=c["recall"])
        if cal_p is not None:
            d["brier_isotonic"] = float(np.sum(w * (cal_p - yf) ** 2) / np.sum(w))
        return d
    return f


def _ap_fn(y, s):
    r = Ranked(y, s)
    return lambda w: dict(ap=r.ap(w))


def comparisons(primary_only=False):
    d = [(f"{PRIMARY}|ap", f"{b}|ap") for b in BASELINES]
    if primary_only:
        return d
    d += [(f"{SECONDARY}|ap", f"{b}|ap") for b in BASELINES]
    d += [(f"{SECONDARY}|ap", f"{PRIMARY}|ap"), (f"{REFERENCE}|ap", f"{PRIMARY}|ap")]
    d += [(f"{PRIMARY}|auroc", f"{b}|auroc") for b in BASELINES]
    d += [(f"{c}|{k}", f"{b}|{k}") for c in (PRIMARY, SECONDARY) for b in BASELINES
          for k in ("macro_cluster_ap", "ap_no_largest_component")]
    return d


def _macro_fn(y, scores: dict, ca, cb, res: B.ClusterResampler):
    """Bootstrap version of Part 1 ``macro_cluster_ap`` (unweighted mean over clusters whose touching pairs hold >= 1
    positive and >= 1 negative of per-cluster AP with the pair weights). A cluster drawn m_c times counts m_c times;
    with all multiplicities 1 it equals Part 1's value (tested)."""
    groups = p1m._groups_touching(np.asarray(ca), np.asarray(cb))
    cl = []
    for c, ix in groups.items():
        if y[ix].any() and (~y[ix]).any():
            k = int(np.searchsorted(res.uniq, c))
            cl.append((k, ix, {m: Ranked(y[ix], s[ix]) for m, s in scores.items()}))

    def f(w):
        m = res.m if res.m is not None else np.ones(len(res.uniq))
        num = {k: 0.0 for k in scores}
        den = {k: 0.0 for k in scores}
        for k, ix, rk in cl:
            mc = m[k]
            if mc == 0:
                continue
            wi = w[ix]
            for name, r in rk.items():
                a = r.ap(wi)
                if np.isfinite(a):
                    num[name] += mc * a
                    den[name] += mc
        return {f"{n}|macro_cluster_ap": (num[n] / den[n] if den[n] > 0 else np.nan) for n in scores}
    return f


def bootstraps(es: EvalSet, cl: dict, rules: dict | None, n_boot=N_BOOT, seed=BOOT_SEED, level="all",
               mc_seeds=()) -> dict:
    """level 'all' (primary ITT set): paired cluster bootstrap of every statistic (AP, AUROC, AP without the largest
    component, macro per-cluster AP, F1/P/R at rule b, isotonic Brier) + protein-level and pair-level bootstraps (AP)
    + Monte-Carlo seeds. 'cluster' (mapped-only, the locked-verdict set): the cluster bootstrap of every statistic.
    'ap' (strata): cluster bootstrap of AP and paired dAP only."""
    y, w = es.y, es.w
    keep_nl = (cl["comp_a"] != cl["largest"]) & (cl["comp_b"] != cl["largest"])
    present = [m for m in METHODS if m in es.scores]
    out, rt = {}, {}
    res = B.ClusterResampler(cl["ca"], cl["cb"], n_boot, seed)
    diffs = [d for d in comparisons() if d[0].split("|")[0] in present and d[1].split("|")[0] in present]
    if level == "ap":
        fns = {m: _ap_fn(y, es.scores[m]) for m in present}
        diffs = [d for d in diffs if d[0].endswith("|ap") and d[1].endswith("|ap")]
    else:
        fns = {}
        for m in present:
            thr = rules["methods"][m]["f1max"]["threshold"] if rules else None
            cal_p = None
            if rules:
                rm = rules["methods"][m]
                cal_p = C.apply(rm["isotonic"], C.floor_to_fit(es.scores[m], es.mapped[m], rm["fit_floor"]))
            fns[m] = _method_fn(y, es.scores[m], w, thr, keep_nl, cal_p)
        fns["*macro"] = _macro_fn(y, {m: es.scores[m] for m in present}, cl["ca"], cl["cb"], res)
    t0 = time.perf_counter()
    out["cluster"] = B.run(res, w, fns, diffs)
    out["cluster"].update(scheme="paired cluster bootstrap (Part 1 src/metrics.cluster_bootstrap scheme)", seed=seed,
                          n_units=int(len(res.uniq)))
    rt["cluster"] = round(time.perf_counter() - t0, 1)
    if level == "all":
        ap_fns = {m: _ap_fn(y, es.scores[m]) for m in present}
        apd = [d for d in diffs if d[0].endswith("|ap") and d[1].endswith("|ap")]
        t0 = time.perf_counter()
        out["protein"] = B.run(B.protein_weights(es.table.a.values, es.table.b.values, n_boot, seed), w, ap_fns, apd)
        out["protein"].update(scheme="protein (node) bootstrap, pair weight m_a * m_b", seed=seed)
        rt["protein"] = round(time.perf_counter() - t0, 1)
        t0 = time.perf_counter()
        out["pair"] = B.run(B.pair_weights(y, n_boot, seed), w, ap_fns, apd)
        out["pair"].update(scheme="stratified pair bootstrap (Part 1 src/metrics._boot_weights)", seed=seed)
        rt["pair"] = round(time.perf_counter() - t0, 1)
    if mc_seeds:
        sub = [PRIMARY, *BASELINES]
        ap_fns = {m: _ap_fn(y, es.scores[m]) for m in sub if m in es.scores}
        t0 = time.perf_counter()
        out["cluster_mc_seeds"] = {
            int(s_): B.run(B.cluster_weights(cl["ca"], cl["cb"], n_boot, int(s_)), w, ap_fns, comparisons(True))["diff"]
            for s_ in mc_seeds}
        rt["cluster_mc_seeds"] = round(time.perf_counter() - t0, 1)
    out["runtime_s"] = rt
    return out


def locked_verdict(pm: dict, boot: dict) -> dict:
    """Part 1's locked IFT decision rule, computed by Part 1's own ``scripts/final_eval.decide`` (DESIGN v3.1 §3.7, IFT
    parts): per candidate c (fixed sequence zavg_BEF, then B) and baseline b, relative AP gain >= 10% (point
    estimate) and non-inferiority (macro per-cluster AP and AP without the largest component >= 0.95 x baseline);
    intersection-union over the three baselines; bootstrap CIs are information only. Inputs: Part 1 metric_suite
    values (ap_ht = AP with this set's design weights) and the paired cluster bootstrap in Part 1's format."""
    fe = part1.final_eval()
    cands = [c for c in (PRIMARY, SECONDARY) if c in pm]
    per = {m: pm[m]["part1_metric_suite"] for m in (*cands, *BASELINES)}
    diff = {}
    for c in cands:
        for b in BASELINES:
            d = boot["cluster"]["diff"][f"{c}|ap - {b}|ap"]
            diff[f"{c}-{b}"] = {k: d[k] for k in ("diff", "ci_low", "ci_high", "rel_gain", "rel_ci_low", "rel_ci_high",
                                                  "p_two_sided")}
    out = fe.decide(per, dict(diff=diff), cands, list(BASELINES))
    out["source"] = "external/part1_locked/scripts/final_eval.py::decide (called unchanged)"
    return out


def ci_statistic(pm: dict, boot: dict) -> dict:
    """Additional pre-specified statistic (TEST_PROTOCOL §7.2), reported regardless of the locked verdict: for each
    candidate vs baseline, whether the 95% paired cluster-bootstrap CI of dAP excludes 0, 'borderline' if the lower
    bound lies within MC_SPREAD of 0, and the descriptive bootstrap CIs of the non-degradation ratios."""
    out = {}
    cd = boot["cluster"]["diff"]
    for c in (PRIMARY, SECONDARY):
        if c not in pm:
            continue
        for b in BASELINES:
            d = cd[f"{c}|ap - {b}|ap"]
            ratios = {}
            for k in NI_METRICS:
                r = cd.get(f"{c}|{k} - {b}|{k}")
                if r:
                    ratios[k] = dict(ratio=1 + r["rel_gain"], ci_low=1 + r["rel_ci_low"], ci_high=1 + r["rel_ci_high"])
            out[f"{c} vs {b}"] = dict(diff=d["diff"], ci_low=d["ci_low"], ci_high=d["ci_high"],
                                     p_two_sided=d["p_two_sided"], ci_excludes_0=bool(d["ci_low"] > 0),
                                     borderline=bool(abs(d["ci_low"]) <= MC_SPREAD), non_degradation_ratio_ci=ratios)
    pv = {k: v["p_two_sided"] for k, v in out.items() if k.startswith(PRIMARY)}
    joint = {c: all(out[f"{c} vs {b}"]["ci_excludes_0"] for b in BASELINES) for c in (PRIMARY, SECONDARY) if c in pm}
    return dict(pairs=out, ci_excludes_0_vs_all_baselines=joint, holm_adjusted_p_primary=B.holm(pv) if pv else {})


def evaluate(es: EvalSet, rules: dict | None, cl: dict, ess=None, n_boot=N_BOOT, level="all", mc_seeds=()) -> dict:
    t0 = time.perf_counter()
    pm = point_metrics(es, cl, ess, rules)
    t1 = time.perf_counter()
    bs = bootstraps(es, cl, rules, n_boot=n_boot, level=level, mc_seeds=mc_seeds)
    t2 = time.perf_counter()
    res = dict(set=es.name, n_pos=int(es.y.sum()), n_neg=int((~es.y).sum()),
               prevalence_design=weighted_prevalence(es.y, es.w), prevalence_unweighted=float(es.y.mean()),
               n_unmapped_any=int((~es.mapped_all).sum()), cluster_stats=cl["stats"], metrics=pm, bootstrap=bs)
    if all(m in pm for m in (PRIMARY, *BASELINES)):
        res["locked_verdict"] = locked_verdict(pm, bs)
        res["ci_statistic"] = ci_statistic(pm, bs)
    res["runtime_s"] = dict(point_metrics=round(t1 - t0, 1), bootstraps=round(t2 - t1, 1))
    return res
