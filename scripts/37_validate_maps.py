#!/usr/bin/env python
"""FZ-C validation (after the map commit; never changes the maps): tiers and the pre-registered bridge endpoint.

Refuses to run unless results/part2/construction/{zavg,gls}_map/MANIFEST.sha256 are committed at HEAD and the map
files match them (validation.require_frozen), so validation inputs are first read after the FZ-C commit.

Tiers (DESIGN v3.1 §2.3): T1 BioPlex (293T U HCT116, Entrez) within-cluster edges vs node-restricted degree x
essentiality-preserving null (10^4, seed 4100), >= 3 observed members, >= 2 edges, BH per level q < 0.05; T2 (not
T1) CCLE protein co-abundance on A models vs size-matched null (10^4, seed 4200), BH per level q < 0.05; T3-neg /
T3-nc. Computed for the zavg map (FZ-C) and, for comparison, the GLS map.
Bridge (primary): size-standardised within-cluster BioPlex precision, level-1 clusters, gene coverage matched
(decision D5), paired protein-node bootstrap B = 2000, seed 5000; "better" iff the 95% CI excludes 0.
Secondary: protein co-abundance (same standardisation, seed 5001); stability by size bin; S-wide pair-level AP
against BioPlex for zavg vs GLS vs Pearson on A. Sensitivities: no coverage matching; locus clusters excluded.

  PYTHONHASHSEED=0 .venv/bin/python scripts/37_validate_maps.py
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from part2 import construct, data, genefilter, paths, runsets, scores  # noqa: E402
from part2 import validation as V  # noqa: E402
from part2.consensus import ensure_hashseed  # noqa: E402

CON = paths.RESULTS / "construction"
MAPS = {"zavg_BEF": CON / "zavg_map", "gls": CON / "gls_map"}
MIN_PROT_MODELS = 20


def gene_universe():
    S = runsets.S_mask()
    _, st = genefilter.dependency_bearing_for_models(list(runsets.A()), return_stats=True)
    genes = data.protein_genes()
    sidx = np.flatnonzero(S)
    fb = st["frac_below"][sidx]
    ess = np.where(fb >= 0.9, "pan", np.where(fb >= 0.05, "selective", "other"))
    return S, sidx, genes.entrez.values[sidx].astype(np.int64), genes.symbol.values[sidx], ess


def clusters_of(d: Path, sidx):
    tab = pd.read_csv(d / "clusters.tsv", sep="\t", keep_default_na=False)
    memb = pd.read_csv(d / "membership.tsv", sep="\t")
    memb["loc"] = np.searchsorted(sidx, memb.gene_index.values)
    assert (sidx[memb["loc"].values] == memb.gene_index.values).all()
    mem = {k: np.sort(v.values) for k, v in memb.groupby("cluster_id")["loc"]}
    return tab, mem


def prot_pairs(R, mpos, clusters):
    out = []
    for ci, m in enumerate(clusters):
        p = np.sort(mpos[m][mpos[m] >= 0])
        if len(p) < 2:
            continue
        ia, ib = np.triu_indices(len(p), 1)
        r = R[p[ia], p[ib]]
        ok = np.isfinite(r)
        out.append(pd.DataFrame(dict(cluster=ci, bin=V.size_bin(len(m)), a=p[ia][ok], b=p[ib][ok], r=r[ok])))
    return pd.concat(out, ignore_index=True)


def main():
    ensure_hashseed()
    t0 = time.perf_counter()
    frozen = V.require_frozen(MAPS.values())
    S, sidx, entrez, symbol, ess = gene_universe()
    timings = {}
    # ---------------------------------------------------------------- validation inputs (first read here)
    t = time.perf_counter()
    keys, deg = V.load_bioplex()
    net = V.Network(entrez, keys, deg)
    from part2.contexts import read_proteomics
    P = read_proteomics()
    A = list(runsets.A())
    P = P.loc[P.index.isin(A)]
    sym2loc = pd.Series(np.arange(len(symbol)), index=symbol)
    sym2loc = sym2loc[~sym2loc.index.duplicated()]
    cols = [c for c in P.columns if c in sym2loc.index]
    Pm = P[cols]
    nobs = Pm.notna().sum(0)
    cols = [c for c in cols if nobs[c] >= MIN_PROT_MODELS]
    measured = np.zeros(len(sidx), bool)
    mloc = sym2loc[cols].values
    measured[mloc] = True
    order = np.argsort(mloc)
    R = scores.pairwise_pearson(Pm[cols].values[:, order])          # local order = np.flatnonzero(measured)
    mpos = np.full(len(sidx), -1)
    mpos[np.flatnonzero(measured)] = np.arange(measured.sum())
    timings["inputs_s"] = round(time.perf_counter() - t, 1)
    inputs = dict(bioplex=dict(files=[str(f.relative_to(paths.PROJECT)) for f in V.BIOPLEX_FILES],
                               sha256=[data.sha256_file(f) for f in V.BIOPLEX_FILES], n_edges_union=len(keys),
                               n_nodes_union=int(len(deg)), S_observed=int(net.observed.sum()),
                               edges_within_S=net.n_edges_universe),
                  proteomics=dict(n_A_models=int(len(P)), n_S_measured=int(measured.sum()),
                                  min_models_per_protein=MIN_PROT_MODELS, gene_mapping="Gene_Symbol == genes.tsv symbol"))
    out_t = CON / "tiers"
    out_b = CON / "bridge"
    out_t.mkdir(parents=True, exist_ok=True)
    out_b.mkdir(parents=True, exist_ok=True)
    # ---------------------------------------------------------------- tiers
    t = time.perf_counter()
    tabs, mems, tier_summ = {}, {}, {}
    for method, d in MAPS.items():
        tab, mem = clusters_of(d, sidx)
        tabs[method], mems[method] = tab, mem
        rows = []
        for lev in (1, 2):
            sub = tab[tab.level == lev]
            cl = [mem[c] for c in sub.cluster_id]
            t1 = V.t1_test(net, cl, ess, B=10_000, seed=4100 + lev)
            t2 = V.t2_test(R, measured, cl, B=10_000, seed=4200 + lev)
            tt = V.assign_tiers(t1, t2)
            tt.insert(0, "cluster_id", sub.cluster_id.values[tt.i.values])
            tt.insert(1, "level", lev)
            rows.append(tt.drop(columns="i"))
        tiers = pd.concat(rows, ignore_index=True)
        tiers = tiers.merge(tab[["cluster_id", "size", "locus", "members"]], on="cluster_id")
        name = "zavg" if method == "zavg_BEF" else "gls"
        tiers.to_csv(out_t / f"{name}_tiers.tsv", sep="\t", index=False, float_format="%.5g")
        s = {}
        for lev in (1, 2):
            x = tiers[tiers.level == lev]
            s[f"level{lev}"] = dict(n=int(len(x)), **{k: int(v) for k, v in x.tier.value_counts().items()},
                                    T1_nonlocus=int(((x.tier == "T1") & ~x.locus.astype(bool)).sum()),
                                    bioplex_testable=int((x.n_obs >= 3).sum()), prot_testable=int((x.n_meas >= 3).sum()),
                                    expected_false_T1=round(0.05 * int((x.tier == "T1").sum()), 1))
        tier_summ[name] = s
    timings["tiers_s"] = round(time.perf_counter() - t, 1)
    (out_t / "summary.json").write_text(json.dumps(dict(map_commits=frozen, inputs=inputs, tiers=tier_summ,
                                                         rules=V.__doc__), indent=1))
    print(json.dumps(tier_summ), flush=True)
    # ---------------------------------------------------------------- bridge
    t = time.perf_counter()
    l1 = {m: tabs[m][tabs[m].level == 1] for m in MAPS}
    mz, mg, cov = V.match_coverage(l1["zavg_BEF"], l1["gls"])

    def run(tz, tg, value="edge", seed=5000):
        cz = [mems["zavg_BEF"][c] for c in tz.cluster_id]
        cg = [mems["gls"][c] for c in tg.cluster_id]
        w = V.bin_weights(tz["size"], tg["size"])
        if value == "edge":
            pz, pg = V.pair_table(net, cz), V.pair_table(net, cg)
            n_nodes = len(net.obs_idx)
        else:
            pz, pg = prot_pairs(R, mpos, cz), prot_pairs(R, mpos, cg)
            n_nodes = int(measured.sum())
        res = V.bridge_bootstrap(pz, pg, n_nodes, w, B=2000, seed=seed, value=value)
        res.update(n_clusters_zavg=len(cz), n_clusters_gls=len(cg), genes_zavg=int(tz["size"].sum()),
                   genes_gls=int(tg["size"].sum()))
        return res

    primary = run(mz, mg)
    primary["coverage_matching"] = cov
    nl = lambda t_: t_[~t_.locus.astype(bool)]
    mz2, mg2, cov2 = V.match_coverage(nl(l1["zavg_BEF"]), nl(l1["gls"]))
    sens = dict(no_coverage_matching=run(l1["zavg_BEF"], l1["gls"]),
                locus_excluded=dict(run(mz2, mg2), coverage_matching=cov2))
    prot = run(mz, mg, value="r", seed=5001)
    # stability by size bin (level 1, coverage-matched maps and full maps)
    stab = []
    for m, tt in (("zavg_BEF", l1["zavg_BEF"]), ("gls", l1["gls"])):
        b = tt["size"].map(V.size_bin)
        for k, (lo, hi) in enumerate(V.SIZE_BINS):
            x = tt.stability[b == k]
            stab.append(dict(method=m, bin=f"{lo}-{hi}", n=int(len(x)), mean_stability=float(x.mean()) if len(x) else np.nan,
                             median_stability=float(x.median()) if len(x) else np.nan))
    stab = pd.DataFrame(stab)
    stab.to_csv(out_b / "stability_by_size.tsv", sep="\t", index=False, float_format="%.4f")
    timings["bridge_s"] = round(time.perf_counter() - t, 1)
    # pair-level AP over S pairs with both genes BioPlex-observed
    t = time.perf_counter()
    from sklearn.metrics import average_precision_score, roc_auc_score
    obs = net.obs_idx
    iu = np.triu_indices(len(obs), 1)
    y = net.adj[iu]
    ap = dict(n_pairs=int(len(y)), n_pos=int(y.sum()), prevalence=float(y.mean()))
    for m in ("zavg_BEF", "gls", "pearson"):
        blk = construct.load("A", m, S)
        v = blk.vals[np.ix_(obs, obs)][iu].astype(np.float64)
        del blk
        ok = np.isfinite(v)
        ap[m] = dict(AP=float(average_precision_score(y[ok], v[ok])), AUROC=float(roc_auc_score(y[ok], v[ok])),
                     n_scored=int(ok.sum()), AP_over_prevalence=float(average_precision_score(y[ok], v[ok]) / y[ok].mean()))
    timings["pair_ap_s"] = round(time.perf_counter() - t, 1)
    timings["total_s"] = round(time.perf_counter() - t0, 1)
    timings["peak_rss_gb"] = round(scores.peak_rss_gb(), 2)
    res = dict(map_commits=frozen, primary=primary, decision=("zavg better" if primary["excludes_0"] and primary["diff"] > 0
                                                              else "GLS better" if primary["excludes_0"] else "no difference claimed"),
               sensitivity=sens, secondary=dict(protein_coabundance=prot, pair_level_AP_bioplex=ap,
                                                stability_by_size=stab.to_dict("records")),
               inputs=inputs, timings=timings)
    (out_b / "bridge.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    print(json.dumps(dict(primary={k: primary[k] for k in ("zavg", "gls", "diff", "ci95", "excludes_0")},
                          decision=res["decision"], ap={m: ap[m]["AP"] for m in ("zavg_BEF", "gls", "pearson")},
                          timings=timings), default=float), flush=True)


if __name__ == "__main__":
    main()
