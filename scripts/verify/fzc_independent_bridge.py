#!/usr/bin/env python
"""Independent re-implementation (reviewer, FZ-C) of the tiers and the bridge endpoint from the COMMITTED maps.

Does not import src/part2/validation.py or scripts/37. Re-uses only producer code for inputs that are not under test:
the gene universe S + A-essentiality fraction (runsets.S_mask / genefilter on A) and the cached CCLE proteomics
parse (contexts.read_proteomics). Everything downstream (BioPlex parsing, union network, strata, nulls, BH, tier
rule, D5 coverage matching, bin weights, size-standardised precision, protein-node bootstrap) is written here.

Bootstrap: (i) "same-stream" = rng.multinomial per replicate as in the DESIGN wording (seed 5000/5001), and
(ii) "alt-stream" = index resampling + bincount (seed 5000), to show the CI is not an artefact of one RNG stream.
Tier nulls use DIFFERENT seeds (7100+level, 7200+level) and a simpler stratum-merging rule, so tier counts are a
Monte-Carlo robustness check, not a bit-for-bit reproduction.

  PART1_DIR=... PYTHONHASHSEED=0 .venv/bin/python scripts/verify/fzc_independent_bridge.py --out <dir>
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
CON = ROOT / "results" / "part2" / "construction"
BINS = [(3, 4), (5, 8), (9, 16), (17, 60)]


def sbin(n):
    return next((i for i, (lo, hi) in enumerate(BINS) if lo <= n <= hi), -1)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    m = len(ok)
    if m == 0:
        return q
    o = ok[np.argsort(p[ok])]
    adj = p[o] * m / np.arange(1, m + 1)
    adj = np.minimum(1, np.minimum.accumulate(adj[::-1])[::-1])
    q[o] = adj
    return q


def load_universe():
    from part2 import data, genefilter, runsets
    S = runsets.S_mask()
    A = list(runsets.A())
    _, st = genefilter.dependency_bearing_for_models(A, return_stats=True)
    g = data.protein_genes()
    idx = np.flatnonzero(S)
    fb = st["frac_below"][idx]
    ess = np.select([fb >= 0.9, fb >= 0.05], ["pan", "selective"], "other")
    return pd.DataFrame(dict(gene_index=idx, symbol=g.symbol.values[idx], entrez=g.entrez.values[idx].astype(np.int64),
                             ess=ess)), A


def load_bioplex():
    e = []
    for f in ("BioPlex_293T_Network_10K_Dec_2019.tsv", "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"):
        d = pd.read_csv(ROOT / "datasets" / "BioPlex" / f, sep="\t", dtype=str)
        a = pd.to_numeric(d["GeneA"], errors="coerce")
        b = pd.to_numeric(d["GeneB"], errors="coerce")
        k = a.notna() & b.notna() & (a != b)
        e.append(np.c_[np.minimum(a[k], b[k]), np.maximum(a[k], b[k])].astype(np.int64))
    e = np.unique(np.vstack(e), axis=0)
    deg = pd.Series(e.ravel()).value_counts()
    return e, deg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    U, A = load_universe()
    ent = U.entrez.values
    # ---- BioPlex network restricted to S (observed = Entrez node of the union network)
    E, deg = load_bioplex()
    obs = np.isin(ent, deg.index.values)
    opos = np.full(len(U), -1)
    opos[obs] = np.arange(obs.sum())
    e2o = pd.Series(opos[obs], index=ent[obs])
    e2o = e2o[~e2o.index.duplicated()]
    keep = np.isin(E[:, 0], e2o.index) & np.isin(E[:, 1], e2o.index)
    no = int(obs.sum())
    ADJ = np.zeros((no, no), bool)
    i, j = e2o.reindex(E[keep, 0]).values, e2o.reindex(E[keep, 1]).values
    ADJ[i, j] = ADJ[j, i] = True
    odeg = deg.reindex(ent[obs]).values.astype(int)
    oess = U.ess.values[obs]
    info = dict(bioplex_union_edges=int(len(E)), bioplex_union_nodes=int(len(deg)), S=int(len(U)),
                S_observed=no, edges_within_S=int(ADJ.sum() // 2))
    # ---- proteomics (A models only), symbol-mapped, >= 20 measured models
    from part2.contexts import read_proteomics
    P = read_proteomics()
    P = P.loc[P.index.isin(A)]
    s2l = pd.Series(np.arange(len(U)), index=U.symbol.values)
    s2l = s2l[~s2l.index.duplicated()]
    cols = [c for c in P.columns if c in s2l.index and P[c].notna().sum() >= 20]
    mloc = s2l[cols].values
    o = np.argsort(mloc)
    cols, mloc = [cols[k] for k in o], mloc[o]
    R = P[cols].astype(float).corr(method="pearson", min_periods=3).values      # pandas pairwise-complete
    ppos = np.full(len(U), -1)
    ppos[mloc] = np.arange(len(mloc))
    info.update(prot_A_models=int(len(P)), prot_S_measured=int(len(mloc)))
    print(json.dumps(info), flush=True)

    # ---- clusters from the committed maps
    def clusters(d):
        tab = pd.read_csv(CON / d / "clusters.tsv", sep="\t", keep_default_na=False)
        mem = pd.read_csv(CON / d / "membership.tsv", sep="\t")
        g2l = pd.Series(np.arange(len(U)), index=U.gene_index.values)
        mem["l"] = g2l.reindex(mem.gene_index.values).values
        assert mem.l.notna().all()
        ms = {c: np.sort(v.astype(int).values) for c, v in mem.groupby("cluster_id").l}
        assert all(len(ms[c]) == s for c, s in zip(tab.cluster_id, tab["size"]))
        return tab, ms

    maps = {"zavg": clusters("zavg_map"), "gls": clusters("gls_map")}
    for m, (tab, ms) in maps.items():
        l1 = np.concatenate([ms[c] for c in tab.cluster_id[tab.level == 1]])
        info[f"{m}_level1_disjoint"] = bool(len(l1) == len(np.unique(l1)))

    # ---- tiers (own null: degree-bin x ess-class strata, bins < 10 genes pooled down; different seeds)
    lb = np.floor(np.log2(np.maximum(odeg, 1))).astype(int)
    lab = np.array([f"{c}|{b}" for c, b in zip(oess, lb)], object)
    for _ in range(100):
        u, cnt = np.unique(lab, return_counts=True)
        small = u[cnt < 10]
        if not len(small):
            break
        s = small[0]
        c, b = s.split("|")
        same = sorted(int(x.split("|")[1]) for x in u if x.startswith(c + "|") and x != s)
        if same:
            tgt = max([x for x in same if x < int(b)], default=min(same))
            lab[lab == s] = f"{c}|{tgt}"
        else:
            lab[lab == s] = f"pooled|{b}"
    pools = {s: np.flatnonzero(lab == s) for s in np.unique(lab)}
    mi = np.flatnonzero(ppos >= 0)

    def t1(members, rng, B=10_000):
        p = np.sort(opos[members][opos[members] >= 0])
        n = len(p)
        if n < 3:
            return n, np.nan, np.nan
        x = int(np.triu(ADJ[np.ix_(p, p)], 1).sum())
        dr = np.stack([pools[lab[g]][rng.integers(0, len(pools[lab[g]]), B)] for g in p], 1)
        for _ in range(30):
            ss = np.sort(dr, 1)
            dup = (ss[:, 1:] == ss[:, :-1]).any(1)
            if not dup.any():
                break
            dr[dup] = np.stack([pools[lab[g]][rng.integers(0, len(pools[lab[g]]), dup.sum())] for g in p], 1)
        ia, ib = np.triu_indices(n, 1)
        null = ADJ[dr[:, ia], dr[:, ib]].sum(1)
        return n, x, (1 + (null >= x).sum()) / (1 + B)

    def t2(members, rng, B=10_000):
        p = np.sort(ppos[members][ppos[members] >= 0])
        n = len(p)
        if n < 3:
            return n, np.nan
        ia, ib = np.triu_indices(n, 1)
        x = np.nanmean(R[p[ia], p[ib]])
        dr = rng.integers(0, len(mi), (B, n))
        for _ in range(60):
            ss = np.sort(dr, 1)
            dup = (ss[:, 1:] == ss[:, :-1]).any(1)
            if not dup.any():
                break
            dr[dup] = rng.integers(0, len(mi), (int(dup.sum()), n))
        with np.errstate(all="ignore"):
            null = np.nanmean(R[dr[:, ia], dr[:, ib]], 1)
        return n, (1 + (null >= x).sum()) / (1 + B)

    tiers, tcount = {}, {}
    for m, (tab, ms) in maps.items():
        rows = []
        for lev in (1, 2):
            sub = tab[tab.level == lev]
            r1, r2 = np.random.default_rng(7100 + lev), np.random.default_rng(7200 + lev)
            res = [(c, *t1(ms[c], r1), *t2(ms[c], r2)) for c in sub.cluster_id]
            d = pd.DataFrame(res, columns=["cluster_id", "n_obs", "edges", "p1", "n_meas", "p2"])
            d["q1"], d["q2"] = bh(d.p1.values), bh(d.p2.values)
            ok1 = (d.q1 < 0.05) & (d.edges >= 2)
            ok2 = ~ok1 & (d.q2 < 0.05)
            tested = (d.n_obs >= 3) | (d.n_meas >= 3)
            d["tier"] = np.select([ok1, ok2, tested], ["T1", "T2", "T3-neg"], "T3-nc")
            d["level"] = lev
            rows.append(d)
        tiers[m] = pd.concat(rows, ignore_index=True).merge(tab[["cluster_id", "size", "locus"]], on="cluster_id")
        tiers[m].to_csv(out / f"{m}_tiers_reviewer.tsv", sep="\t", index=False, float_format="%.5g")
        tcount[m] = {f"level{lev}": tiers[m][tiers[m].level == lev].tier.value_counts().to_dict() for lev in (1, 2)}
        prod = pd.read_csv(CON / "tiers" / f"{m}_tiers.tsv", sep="\t")
        mm = tiers[m].merge(prod[["cluster_id", "tier", "edges", "n_obs", "n_meas"]], on="cluster_id", suffixes=("", "_prod"))
        tcount[m]["agreement"] = dict(
            tier_same=int((mm.tier == mm.tier_prod).sum()), n=int(len(mm)),
            edges_identical=bool(np.array_equal(mm.edges.fillna(-1).values, mm.edges_prod.fillna(-1).values)),
            n_obs_identical=bool((mm.n_obs == mm.n_obs_prod).all()), n_meas_identical=bool((mm.n_meas == mm.n_meas_prod).all()),
            discordant=mm.loc[mm.tier != mm.tier_prod, ["cluster_id", "level", "tier", "tier_prod"]].values.tolist())
    print(json.dumps(tcount, default=str), flush=True)

    # ---- bridge (D5)
    def match(ta, tb):
        ca, cb = ta["size"].sum(), tb["size"].sum()

        def trim(t, target):
            t = t.sort_values(["stability", "cohesion", "cluster_id"], ascending=[True, True, False])
            cs = t["size"].sum() - t["size"].cumsum()                     # coverage after dropping rows 0..i
            n_drop = int(np.argmax(cs.values <= target)) + 1 if t["size"].sum() > target else 0
            return t.iloc[n_drop:], t.iloc[:n_drop]
        if ca > cb:
            ka, da = trim(ta, cb)
            return ka, tb, da
        if cb > ca:
            kb, db = trim(tb, ca)
            return ta, kb, db
        return ta, tb, ta.iloc[:0]

    def pairs(tab, ms, posmap, value):
        rows = []
        for c, sz in zip(tab.cluster_id, tab["size"]):
            p = np.sort(posmap[ms[c]][posmap[ms[c]] >= 0])
            if len(p) < 2:
                continue
            ia, ib = np.triu_indices(len(p), 1)
            v = ADJ[p[ia], p[ib]].astype(float) if value == "edge" else R[p[ia], p[ib]]
            rows.append(np.c_[np.full(len(ia), sbin(sz)), p[ia], p[ib], v])
        x = np.vstack(rows)
        return x[np.isfinite(x[:, 3])]

    def stdz(x, w, mult=None):
        ww = np.ones(len(x)) if mult is None else mult[x[:, 1].astype(int)] * mult[x[:, 2].astype(int)]
        per = np.array([(ww[x[:, 0] == b] * x[x[:, 0] == b, 3]).sum() / ww[x[:, 0] == b].sum()
                        if ww[x[:, 0] == b].sum() > 0 else np.nan for b in range(4)])
        ok = np.isfinite(per) & (w > 0)
        return (per[ok] * w[ok]).sum() / w[ok].sum(), per

    def bridge(tz, tg, value="edge", seed=5000, B=2000):
        msz, msg = maps["zavg"][1], maps["gls"][1]
        posmap = opos if value == "edge" else ppos
        nn = no if value == "edge" else len(mi)
        xz, xg = pairs(tz, msz, posmap, value), pairs(tg, msg, posmap, value)
        sb = np.array([sbin(s) for s in list(tz["size"]) + list(tg["size"])])
        w = np.array([(sb == b).mean() for b in range(4)])
        ez, pz = stdz(xz, w)
        eg, pg = stdz(xg, w)
        used = np.unique(np.r_[xz[:, 1], xz[:, 2], xg[:, 1], xg[:, 2]]).astype(int)
        res = dict(zavg=ez, gls=eg, diff=ez - eg, per_bin_zavg=pz.tolist(), per_bin_gls=pg.tolist(), w=w.tolist(),
                   n_clusters=[len(tz), len(tg)], genes=[int(tz["size"].sum()), int(tg["size"].sum())],
                   n_pairs=[len(xz), len(xg)], n_nodes=len(used))
        for stream in ("same", "alt"):
            rng = np.random.default_rng(seed)
            d = np.empty(B)
            for t in range(B):
                mult = np.zeros(nn)
                if stream == "same":
                    mult[used] = rng.multinomial(len(used), np.full(len(used), 1 / len(used)))
                else:
                    mult[used] = np.bincount(rng.integers(0, len(used), len(used)), minlength=len(used))
                d[t] = stdz(xz, w, mult)[0] - stdz(xg, w, mult)[0]
            res[f"ci95_{stream}"] = np.nanpercentile(d, [2.5, 97.5]).tolist()
            res[f"boot_sd_{stream}"] = float(np.nanstd(d))
        return res

    tz, tg = [maps[m][0][maps[m][0].level == 1] for m in ("zavg", "gls")]
    mz, mg, dropped = match(tz, tg)
    br = dict(primary=bridge(mz, mg))
    print(json.dumps(dict(primary=br["primary"]), default=float), flush=True)
    br["no_matching"] = bridge(tz, tg)
    nl = lambda t: t[~t.locus.astype(str).str.lower().eq("true")]
    mz2, mg2, _ = match(nl(tz), nl(tg))
    br["locus_excluded"] = bridge(mz2, mg2)
    br["protein_coabundance"] = bridge(mz, mg, value="r", seed=5001)
    # ---- coverage-matching bias probe: dropped vs kept zavg clusters
    pr = lambda t: stdz(pairs(t, maps["zavg"][1], opos, "edge"), np.array(br["primary"]["w"]))
    dsz = dropped["size"].map(sbin).value_counts().sort_index().to_dict()
    ksz = mz["size"].map(sbin).value_counts().sort_index().to_dict()
    zt = tiers["zavg"][tiers["zavg"].level == 1].set_index("cluster_id")
    br["matching_probe"] = dict(
        n_dropped=int(len(dropped)), genes_dropped=int(dropped["size"].sum()),
        dropped_size_bins=dsz, kept_size_bins=ksz,
        dropped_precision_std=pr(dropped)[0], dropped_per_bin=pr(dropped)[1].tolist(),
        kept_precision_std=br["primary"]["zavg"],
        dropped_T1_frac=float((zt.loc[dropped.cluster_id].tier == "T1").mean()),
        kept_T1_frac=float((zt.loc[mz.cluster_id].tier == "T1").mean()),
        dropped_locus_frac=float(dropped.locus.astype(str).str.lower().eq("true").mean()),
        kept_locus_frac=float(mz.locus.astype(str).str.lower().eq("true").mean()),
        dropped_median_stability=float(dropped.stability.median()), kept_median_stability=float(mz.stability.median()))
    # random-drop alternative: drop random zavg clusters to matched coverage (50 draws) -> distribution of the diff
    rs = np.random.default_rng(9100)
    rd = []
    for _ in range(50):
        perm = tz.sample(frac=1, random_state=int(rs.integers(1 << 31)))
        cs = perm["size"].cumsum().values
        keepn = int(np.searchsorted(cs, tg["size"].sum(), side="right"))
        sub = perm.iloc[:keepn]
        x = pairs(sub, maps["zavg"][1], opos, "edge")
        sb = np.array([sbin(s) for s in list(sub["size"]) + list(tg["size"])])
        w = np.array([(sb == b).mean() for b in range(4)])
        rd.append(stdz(x, w)[0] - stdz(pairs(tg, maps["gls"][1], opos, "edge"), w)[0])
    br["random_drop_matching"] = dict(n=50, mean_diff=float(np.mean(rd)), q025=float(np.quantile(rd, .025)),
                                      q975=float(np.quantile(rd, .975)))
    rep = dict(info=info, tiers=tcount, bridge=br, runtime_s=round(time.perf_counter() - t0, 1))
    (out / "independent_bridge.json").write_text(json.dumps(rep, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    print(json.dumps({k: v for k, v in br.items() if k != "primary"}, default=float, indent=1), flush=True)


if __name__ == "__main__":
    main()
