#!/usr/bin/env python
"""R-1 (referee M1): matched-set reference. POST HOC, DESCRIPTIVE, A only (results/part2/revision/INTENT.md).

  (a) matched random sets drawn from S genes outside every frozen unit (generators nbr / unif / ess, the acceptance
      rule of calib.matched_random_sets), scored and tested as the registered calibration; p < 0.05 rates against the
      real units and the complex-member matched sets.
  (b) selection-aware reference for "beats all matched sets": 10 random atlases (unit -> its r-th matched member set)
      run through the discovery selection (BB -> permutation re-test of the q < 0.1 screen -> BB); the symmetric flag
      rate and the |t| / |beta| distributions conditional on discovery.

No run_association call, no ledger or trials line, no lockbox data. Outputs: results/part2/revision/r1_*.
Usage (PART1_DIR set, PYTHONHASHSEED=0):  .venv/bin/python scripts/58_rev_matched_sets.py [--workers 7]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import multiprocessing as mp
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, calib, covariates as covmod, cscore, registry, revision as rv  # noqa: E402

OUT = ROOT / "results" / "part2" / "revision"
DIRS = ROOT / "results" / "part2" / "directions"
MAP = ROOT / "results" / "part2" / "fzx" / "zavg_engine_membership.tsv"
AIDS = ("H3a_lineage_atlas", "H3b_driver_family")
R_SETS, N_PERM, SCREEN_Q = 10, 10_000, 0.1


def runner():
    spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
    R = importlib.util.module_from_spec(spec)
    sys.modules["run_direction"] = R
    spec.loader.exec_module(R)
    return R


def alt_of(entry, c):
    return {"+": "less", "-": "greater", "two-sided": "two-sided"}[registry.contexts_of(entry)[c]["sign"]]


class Dir:
    """Per-direction context: entry, contexts, covariates, alternatives, score cache over named sets."""

    def __init__(self, env, R, aid):
        self.env, self.aid = env, aid
        self.entry, _ = R.entry_of(env, aid)
        self.q = float(self.entry["fdr_family"]["q"])
        self.seed = int(self.entry["calibration"]["seed"])
        self.perm_seed = int(self.entry["seeds"]["permutation"])
        self.ctxs = R.primary_contexts(self.entry)
        self.covs, _ = covmod.per_context(self.entry, self.ctxs, env.cov)
        self.alts = {c: alt_of(self.entry, c) for c in self.ctxs}
        items = {c["id"]: c for c in self.entry["context_family"]["contexts"]}
        self.lineage = {c for c in self.ctxs if cscore.is_lineage_context(c, items.get(c))}
        self.rem = {c: cscore.removal_for_context(c, env.homdel.get(c)) for c in self.ctxs}

    def test_sets(self, sets: dict) -> tuple[pd.DataFrame, dict]:
        """Score `sets` per context (removals) and test (covariates, alternative). Returns the long test table and the
        score frames per removal."""
        cache, rows = {}, []
        for c in self.ctxs:
            r = self.rem[c]
            if r not in cache:
                cache[r] = calib.score_sets(self.env.E_A, sets, r, self.env.annot)
            f = calib.test_columns(cache[r], self.env.ctx[c], self.covs[c], self.alts[c])
            rows.append(f.assign(context=c, set=f.index, unit=[k.split("#")[0] for k in f.index]))
        return pd.concat(rows, ignore_index=True), cache

    def perm_payload(self, cache, c, setname, n_perm, seedseq):
        env = self.env
        Y = cache[self.rem[c]]
        X, _, ok = assoc._design(env.ctx[c].reindex(Y.index), self.covs[c].reindex(Y.index))
        y = Y[setname].values[ok]
        m = ~np.isnan(y)
        ids = Y.index[ok][m]
        st = None if c in self.lineage else env.strata.reindex(ids).fillna("NA").astype(str).values
        return (y[m], X[m], list(ids), self.alts[c], n_perm, seedseq, st)


def run_perms(pool, payloads: dict) -> dict:
    if not payloads:
        return {}
    args = [(k,) + v for k, v in payloads.items()]
    return dict(pool.imap_unordered(rv.perm_task, args, chunksize=4))


# ------------------------------------------------------------------------------------------------ R-1(a)
def nonmember_inputs(env, R):
    """Pool N (S genes that are A columns and in no unit of the zavg membership table), scored unit members and gene
    essentiality, as calib.matched_random_sets uses them."""
    ps = R.load_pair_scores(env.proj.root)
    map_genes = set(pd.read_csv(MAP, sep="\t").gene)
    cols = set(env.E_A.columns)
    pool = sorted(g for g in ps.index if g in cols and g not in map_genes)
    units = env.fmap.coarse()
    allg = {g for u in units for g in env.fmap.members[u] if g in cols}
    unit_members = {u: [g for g in env.fmap.members[u] if g in allg] for u in units}
    info = dict(n_S=int(len(ps.index)), n_map_genes=len(map_genes), n_pool=len(pool))
    return ps, pool, units, unit_members, calib.clip_essentiality(env.E_A), info


def _gen_task(a):
    aid, gen, pool, ess, ps, unit_members, units, seed = a
    t0 = time.time()
    out = rv.pool_matched_sets(pool, ess, ps, unit_members, units, R_SETS, seed, gen)
    return (aid, gen), out + (round(time.time() - t0, 1),)


def part_a(d: Dir, R, real, member_tests, out_tests, quality_rows, rate_rows, js, gen_sets):
    env = d.env
    units = env.fmap.coarse()
    real_ev = real[real.evaluable.astype(bool)]
    for gen in rv.GENERATORS:
        t0 = time.time()
        sets, st, qual, gen_s = gen_sets[(d.aid, gen)]
        ft, _ = d.test_sets(sets) if sets else (pd.DataFrame(columns=["context", "set", "unit", "beta", "t", "p"]), None)
        out_tests.append(ft[["context", "set", "unit", "beta", "t", "p"]].assign(direction=d.aid, generator=gen))
        quality_rows.append(dict(
            direction=d.aid, generator=gen, n_units=len(units), **{f"status_{k}": int(v) for k, v in
                                                                  st.status.value_counts().items()},
            n_sets=len(sets), med_abs_d_ess=float(qual.d_ess.abs().median()) if len(qual) else None,
            p95_abs_d_ess=float(qual.d_ess.abs().quantile(0.95)) if len(qual) else None,
            med_abs_d_coh=float(qual.d_coh.abs().median()) if len(qual) else None,
            p95_abs_d_coh=float(qual.d_coh.abs().quantile(0.95)) if len(qual) else None,
            median_target_coh_available=float(st[st.n_sets > 0].target_coh.median()) if (st.n_sets > 0).any() else None,
            median_target_coh_unavailable=float(st[st.n_sets == 0].target_coh.median()) if (st.n_sets == 0).any()
            else None, generation_s=gen_s, test_s=round(time.time() - t0, 1)))
        have = set(st.unit[st.n_sets > 0])
        bb = calib.bb_rejections(ft.p.values.astype(float), ft.context.values, d.q) if len(ft) else {}
        js.setdefault("bb", {})[gen] = bb
        for scope, keep in (("all_units", None), ("common_units", have)):
            for typ, tab in (("real", real_ev), ("member_sets", member_tests), (f"nonmember_{gen}", ft)):
                if typ in ("real", "member_sets") and scope == "all_units" and gen != rv.GENERATORS[0]:
                    continue
                tb = tab if keep is None else tab[tab.unit.isin(keep)]
                scope_lab = scope if keep is None else f"common_units_{gen}"
                rate_rows.append(dict(direction=d.aid, scope=scope_lab, set_type=typ, context="pooled",
                                      n_units=int(tb.unit.nunique()), **rv.frac_lt(tb.p)))
                for c, g in tb.groupby("context"):
                    rate_rows.append(dict(direction=d.aid, scope=scope_lab, set_type=typ, context=c,
                                          n_units=int(g.unit.nunique()), **rv.frac_lt(g.p)))


# ------------------------------------------------------------------------------------------------ R-1(b)
def part_b(d: Dir, real, sets, member_tests, cache, pool, js, max_hours: float = 2.0, t_start: float = 0.0,
           approx: bool = False):
    """Random atlases r = 0..9 through the discovery selection. Returns the random test table with discovery flags."""
    mt = member_tests.copy()
    mt["rep"] = mt.set.str.split("#r").str[1].astype(int)
    tabs = []
    for r in range(R_SETS):
        t0 = time.time()
        res = mt[mt.rep == r].copy()
        res["context"] = pd.Categorical(res.context, d.ctxs, ordered=True)
        res = res.sort_values(["context", "unit"]).reset_index(drop=True)
        res["context"] = res.context.astype(str)
        res["evaluable"] = np.isfinite(res.p.astype(float))
        if approx:
            res["discovery"] = res.p < 0.05
            res["p_final"] = res.p
        else:
            def perm_fn(idx, res=res, r=r):
                pay = {i: d.perm_payload(cache, res["context"][i], res["set"][i], N_PERM, (d.perm_seed, r + 1, int(i)))
                       for i in idx}
                return run_perms(pool, pay)
            res = rv.discover(res, d.q, perm_fn, SCREEN_Q)
        res["rep"] = r
        tabs.append(res)
        el = time.time() - t0
        js.setdefault("replicates", []).append(dict(rep=r, n_tests=int(res.evaluable.sum()),
                                                    n_screened=int(res.p_perm.notna().sum()) if "p_perm" in res else 0,
                                                    n_discovery=int(res.discovery.sum()), runtime_s=round(el, 1)))
        print(f"  {d.aid} rep {r}: {int(res.discovery.sum())} discoveries, {el:.0f} s", flush=True)
        if r == 0 and not approx and d.aid == AIDS[0] and (time.time() - t_start) + el * R_SETS * 1.5 > max_hours * 3600:
            raise SystemExit(f"cost rule: projected > {max_hours} h after replicate 0 ({el:.0f} s); rerun with --approx")
    return pd.concat(tabs, ignore_index=True)


def summarise_b(d: Dir, real, rand, js) -> list:
    rd = real[real.evaluable.astype(bool)].copy()
    rd["abs_t"], rd["abs_beta"] = rd.t.abs(), rd.beta.abs()
    rand = rand.copy()
    rand["abs_t"], rand["abs_beta"] = rand.t.abs(), rand.beta.abs()
    sym = rv.symmetric_flag(rd[["context", "unit", "abs_t"]], rand[["context", "unit", "rep", "abs_t"]])
    rand["exceeds_all_others"] = sym.exceeds_all_others.values
    rand["max_abs_t_others"] = sym.max_abs_t_others.values
    mx = rand.groupby(["context", "unit"]).abs_t.max()
    rdisc = rd[rd.discovery.astype(bool)].copy()
    rdisc["max_abs_t_matched"] = [mx.get((c, u), np.nan) for c, u in zip(rdisc.context, rdisc.unit)]
    rdisc["beats_all"] = rdisc.abs_t > rdisc.max_abs_t_matched
    with_sets = rdisc[rdisc.max_abs_t_matched.notna()]
    xdisc = rand[rand.discovery.astype(bool)]
    rdisc["p_beats_rand_disc_t_ctx"] = rv.conditional_auc(rdisc, xdisc, "abs_t", ("context",))
    rdisc["p_beats_rand_disc_beta_ctx"] = rv.conditional_auc(rdisc, xdisc, "abs_beta", ("context",))
    rdisc["p_beats_rand_disc_t_unit"] = rv.conditional_auc(rdisc, xdisc, "abs_t", ("context", "unit"))
    rdisc["p_beats_rand_disc_beta_unit"] = rv.conditional_auc(rdisc, xdisc, "abs_beta", ("context", "unit"))
    # the same, unconditional (all random tests of the context): the reference without selection
    rdisc["p_beats_rand_all_t_ctx"] = rv.conditional_auc(rdisc, rand[np.isfinite(rand.abs_t)], "abs_t", ("context",))
    s = dict(
        n_real_discoveries=int(len(rdisc)), n_real_with_sets=int(len(with_sets)),
        real_beats_all=int(with_sets.beats_all.sum()),
        real_beats_all_rate=float(with_sets.beats_all.mean()),
        n_random_discoveries_total=int(len(xdisc)),
        random_discoveries_per_rep_mean=float(len(xdisc) / R_SETS),
        random_symmetric_flag=int(xdisc.exceeds_all_others.sum()),
        random_symmetric_flag_rate=float(xdisc.exceeds_all_others.mean()) if len(xdisc) else None,
        random_symmetric_flag_rate_all_tests=float(rand[np.isfinite(rand.abs_t)].exceeds_all_others.mean()),
        null_rate_one_of_11=1 / 11,
        mean_p_beats_rand_disc_t_ctx=float(rdisc.p_beats_rand_disc_t_ctx.mean()),
        mean_p_beats_rand_disc_beta_ctx=float(rdisc.p_beats_rand_disc_beta_ctx.mean()),
        mean_p_beats_rand_disc_t_unit=float(rdisc.p_beats_rand_disc_t_unit.mean()),
        n_real_with_unit_rand_disc=int(rdisc.p_beats_rand_disc_t_unit.notna().sum()),
        mean_p_beats_rand_disc_beta_unit=float(rdisc.p_beats_rand_disc_beta_unit.mean()),
        mean_p_beats_rand_all_t_ctx=float(rdisc.p_beats_rand_all_t_ctx.mean()),
        frac_real_above_median_rand_disc_t_ctx=float((rdisc.p_beats_rand_disc_t_ctx > 0.5).mean()),
        median_abs_t_real_disc=float(rdisc.abs_t.median()), median_abs_t_rand_disc=float(xdisc.abs_t.median()),
        median_abs_beta_real_disc=float(rdisc.abs_beta.median()),
        median_abs_beta_rand_disc=float(xdisc.abs_beta.median()))
    js["summary"] = s
    per = []
    for c, g in rdisc.groupby("context"):
        x = xdisc[xdisc.context == c]
        per.append(dict(direction=d.aid, context=c, n_real_disc=len(g), n_rand_disc_10reps=len(x),
                        rand_disc_per_rep=len(x) / R_SETS, real_beats_all_rate=float(g.beats_all[g.max_abs_t_matched
                                                                                              .notna()].mean()),
                        rand_symmetric_flag_rate=float(x.exceeds_all_others.mean()) if len(x) else np.nan,
                        mean_p_beats_rand_disc_t=float(g.p_beats_rand_disc_t_ctx.mean()),
                        mean_p_beats_rand_disc_beta=float(g.p_beats_rand_disc_beta_ctx.mean()),
                        median_abs_t_real=float(g.abs_t.median()),
                        median_abs_t_rand=float(x.abs_t.median()) if len(x) else np.nan))
    return per, rdisc, rand


def check_member_reproduction(d: Dir, mt: pd.DataFrame) -> dict:
    ct = pd.read_csv(DIRS / d.aid / "calibration_D8_tests.tsv.gz", sep="\t")
    ct = ct[ct.calibration == "matched_random_set"]
    m = ct.merge(mt[["context", "set", "t", "p"]], on=["context", "set"], how="outer", suffixes=("_c", "_r"),
                 indicator=True)
    both = m[m._merge == "both"]
    dt = (both.t_c - both.t_r).abs()
    nan_mismatch = int((both.t_c.isna() != both.t_r.isna()).sum())
    # the committed table stores t with 6 significant digits (%.6g): agreement is checked to that precision
    rel = dt / np.maximum(1.0, both.t_c.abs())
    out = dict(n_committed=int(len(ct)), n_recomputed=int(len(mt)), n_matched=int(len(both)),
               max_abs_dt=float(dt.max()), max_rel_dt=float(rel.max()), n_nan_mismatch=nan_mismatch,
               tolerance="|dt| <= 1e-5 * max(1, |t|) (6 significant digits stored)")
    if len(both) != len(ct) or len(ct) != len(mt) or not (rel.fillna(0) <= 1e-5).all() or nan_mismatch:
        raise SystemExit(f"member matched sets do not reproduce the committed calibration: {out}")
    return out


def check_discover_on_real(d: Dir, real: pd.DataFrame) -> dict:
    """The selection code reproduces the committed discovery flags when fed the committed permutation p."""
    res = real[["context", "unit", "p", "evaluable"]].copy()
    pp = real.p_perm
    out = rv.discover(res, d.q, lambda idx: {i: pp[i] for i in idx if np.isfinite(pp[i])}, SCREEN_Q)
    n_scr = int((out.q_adj.fillna(1) < SCREEN_Q).sum())
    ok = bool((out.discovery.values == real.discovery.astype(bool).values).all())
    r = dict(n_discovery=int(out.discovery.sum()), committed=int(real.discovery.sum()), identical=ok,
             n_screened=n_scr, committed_perm_retested=int(real.p_perm.notna().sum()))
    if not ok or n_scr != r["committed_perm_retested"]:
        raise SystemExit(f"discovery selection does not reproduce: {r}")
    return r


def check_perm_on_real(d: Dir, real: pd.DataFrame, pool, n: int = 4) -> list:
    """A few committed permutation p re-computed with the run's seed (perm_seed + row) on the real units."""
    env = d.env
    units = env.fmap.coarse()
    rows = real.index[real.p_perm.notna()][:n]
    cache = {}
    pay = {}
    for i in rows:
        c = real.context[i]
        rm = d.rem[c]
        if rm not in cache:
            cache[rm] = cscore.unit_scores(env.E_A, env.fmap, units, rm, env.annot)[0]
        Y = cache[rm]
        X, _, ok = assoc._design(env.ctx[c].reindex(Y.index), d.covs[c].reindex(Y.index))
        y = Y[real.unit[i]].values[ok]
        m = ~np.isnan(y)
        ids = Y.index[ok][m]
        st = None if c in d.lineage else env.strata.reindex(ids).fillna("NA").astype(str).values
        pay[int(i)] = (y[m], X[m], list(ids), real.alternative[i], N_PERM, (d.perm_seed + int(i),), st)
    got = run_perms(pool, pay)
    out = [dict(row=int(i), committed=float(real.p_perm[i]), recomputed=float(got[int(i)])) for i in rows]
    if not all(abs(o["committed"] - o["recomputed"]) < 1e-12 for o in out):
        raise SystemExit(f"permutation p does not reproduce: {out}")
    return out


def runner_entry_seed(env, R, aid):
    return R.entry_of(env, aid)[0]["calibration"]["seed"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--approx", action="store_true", help="documented approximation (cost rule)")
    a = ap.parse_args()
    t_start = time.time()
    R = runner()
    env = R.real_env()
    OUT.mkdir(parents=True, exist_ok=True)
    js = dict(label="POST HOC, DESCRIPTIVE (revision R-1; results/part2/revision/INTENT.md); A only, no lockbox data",
              A_models=len(env.E_A), directions={})
    out_tests, quality, rates, per_ctx_b, rdisc_all, rand_all = [], [], [], [], [], []
    ctx = mp.get_context("spawn")
    with ctx.Pool(a.workers, initializer=rv.init_perm_worker, initargs=(env.blocks,)) as pool:
        ps, npool, units, unit_members, gene_ess, pinfo = nonmember_inputs(env, R)
        js["nonmember_pool"] = pinfo
        tg = time.time()
        seeds = {aid: int(runner_entry_seed(env, R, aid)) for aid in AIDS}
        tasks = [(aid, gen, npool, gene_ess, ps.reindex(index=sorted(set(npool) | {g for m in unit_members.values()
                                                                                for g in m}),
                                                       columns=sorted(set(npool) | {g for m in unit_members.values()
                                                                                  for g in m})),
                  unit_members, units, seeds[aid]) for aid in AIDS for gen in rv.GENERATORS]
        del ps
        gen_sets = dict(pool.imap_unordered(_gen_task, tasks))
        del tasks
        js["nonmember_generation_wall_s"] = round(time.time() - tg, 1)
        print(f"non-member sets generated in {js['nonmember_generation_wall_s']} s", flush=True)
        for aid in AIDS:
            print(aid, flush=True)
            d = Dir(env, R, aid)
            dj = js["directions"][aid] = {}
            real = R.committed_discovery(env, aid).reset_index(drop=True)
            dj["check_discovery_selection"] = check_discover_on_real(d, real)
            dj["check_permutation_p"] = check_perm_on_real(d, real, pool)
            ess = calib.clip_essentiality(env.E_A)
            ps = R.load_pair_scores(env.proj.root)
            sets, st = calib.matched_random_sets(env.E_A, env.fmap, env.fmap.coarse(), ess, ps, R_SETS, d.seed)
            del ps
            mt, cache = d.test_sets(sets)
            dj["check_member_sets"] = check_member_reproduction(d, mt)
            dj["member_sets_status"] = st.status.value_counts().to_dict()
            # ---- (a)
            ta = time.time()
            part_a(d, R, real, mt, out_tests, quality, rates, dj.setdefault("a", {}), gen_sets)
            dj["a"]["runtime_s"] = round(time.time() - ta, 1)
            # ---- (b)
            tb = time.time()
            bj = dj.setdefault("b", {})
            rand = part_b(d, real, sets, mt, cache, pool, bj, t_start=t_start, approx=a.approx)
            bj["mode"] = "approximation: raw p < 0.05 screen" if a.approx else "full discovery selection"
            per, rdisc, rand = summarise_b(d, real, rand, bj)
            bj["runtime_s"] = round(time.time() - tb, 1)
            per_ctx_b += per
            rdisc_all.append(rdisc.assign(direction=aid))
            rand_all.append(rand.assign(direction=aid))
            print(json.dumps(bj["summary"], indent=1), flush=True)
    pd.DataFrame(quality).to_csv(OUT / "r1a_matching_quality.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(rates).to_csv(OUT / "r1a_rates.tsv", sep="\t", index=False, float_format="%.4f")
    pd.concat(out_tests).to_csv(OUT / "r1a_nonmember_tests.tsv.gz", sep="\t", index=False, float_format="%.6g",
                                compression=dict(method="gzip", mtime=0))
    pd.DataFrame(per_ctx_b).to_csv(OUT / "r1b_per_context.tsv", sep="\t", index=False, float_format="%.4f")
    cols = ["direction", "context", "unit", "beta", "t", "p_final", "abs_t", "max_abs_t_matched", "beats_all",
            "p_beats_rand_disc_t_ctx", "p_beats_rand_disc_beta_ctx", "p_beats_rand_disc_t_unit",
            "p_beats_rand_disc_beta_unit", "p_beats_rand_all_t_ctx"]
    pd.concat(rdisc_all)[cols].to_csv(OUT / "r1b_real_discoveries.tsv", sep="\t", index=False, float_format="%.6g")
    rc = ["direction", "rep", "context", "set", "unit", "beta", "t", "p", "q_adj", "p_perm", "p_final", "q_adj_perm",
          "discovery", "exceeds_all_others"]
    ra = pd.concat(rand_all)
    ra[ra.discovery.astype(bool)][rc].to_csv(OUT / "r1b_random_discoveries.tsv.gz", sep="\t", index=False,
                                             float_format="%.6g", compression=dict(method="gzip", mtime=0))
    js["runtime_s"] = round(time.time() - t_start, 1)
    js["peak_rss_bytes_main"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    js["peak_rss_bytes_workers_max"] = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (OUT / "r1_summary.json").write_text(json.dumps(js, indent=1, default=float) + "\n")
    print(f"done {js['runtime_s']} s")


if __name__ == "__main__":
    main()
