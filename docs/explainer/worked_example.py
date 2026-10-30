"""Worked example for docs/explainer/METHOD.md: one complex x one context traced through every Part 2 step.

Chosen example: complex Z0303 (SNAPC1-SNAPC5, the snRNA-activating protein complex) x context lineage:Bone
(one-vs-rest, H3a lineage atlas). See METHOD.md "Worked example" for why.

What this script does (construction set A only; it recomputes, it never re-decides):
  step 1  counts of A / L lines and of Bone lines in each, from the committed FZ-LB files (model labels only);
  step 2  member pair scores (zavg_BEF recomputed on A, cached S x S block), mutual top-k / z0 edge test on the A graph
          and on the 20 PatientID-blocked 80% subsample graphs that fed the consensus, co-clustering frequencies from
          the cached stage-1 Leiden runs, near-miss genes, locus flag, attached genes (committed table);
  step 3  the complex score per A line (cscore.unit_scores, removal rule of the context) and its distribution;
  step 4  the HC3 / Satterthwaite fit (engine pieces assoc._design + stats.HC3Design), the Freedman-Lane permutation
          re-run with the engine's seed, and Benjamini-Bogomolov recomputed from the committed p-values; each compared
          with the committed discovery table;
  step 5  the 10 matched random sets regenerated with the registered rule and seed and re-tested (compared with the
          committed calibration t), plus the committed never-expressed and pseudo-context results;
  step 6  lockbox result, robustness, strongest-findings status and post hoc labels, copied ONLY from committed result
          files; LOSO / leave-top-2-out / breadth re-run on A.
It does NOT call run_association (which would append to the trials ledger), read_lockbox or resume_lockbox_read, and
never opens data_cache/part2/lockbox_reads, any lockbox dependency value, or any hu.MAP file.

Usage (project root):
    PART1_DIR=".../coessentiality_protein complexes" PYTHONHASHSEED=0 .venv/bin/python docs/explainer/worked_example.py
Writes docs/explainer/data/worked_example.json, data/worked_example_scores_A.tsv and the worked-example figures
(via make_figures.py).
"""
from __future__ import annotations

import gzip
import io
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

from part2 import assoc, calib, covariates as covmod, cscore, data, locus, registry, robust, runsets, scores, stats  # noqa: E402
from part2.cmap import FrozenMap, symbol  # noqa: E402
from part2.consensus import gamma_grid  # noqa: E402
from part2.scores import select_topk  # noqa: E402

UNIT, CTX, AID = "Z0303", "lineage:Bone", "H3a_lineage_atlas"
Z0, K = 6.0, 10
N_PERM, PERM_SEED, CAL_SEED, R_SETS = 10_000, 4000, 4600, 10
A_BLOCK = ROOT / "data_cache/part2/blocks/zavg_BEF_bd28f78d612bef9f_S8c264ca513f8.npz"
LEVELS = ROOT / "data_cache/part2/construction/zavg_BEF_levels.npz"
OUT = HERE / "data"
FORBIDDEN = ("lockbox_reads", "huMAP", "test_eval")


def rd(path, **kw):
    p = ROOT / path
    assert not any(f in str(p) for f in FORBIDDEN), p
    return pd.read_csv(p, sep="\t", **kw)


def tsv_gz(path):
    with gzip.open(ROOT / path, "rt") as fh:
        return pd.read_csv(io.StringIO(fh.read()), sep="\t")


def f(x, nd=6):
    """JSON-friendly float (None for NaN)."""
    if x is None:
        return None
    x = float(x)
    return None if not np.isfinite(x) else float(f"{x:.{nd}g}")


def main():
    t0 = time.time()
    J = dict(unit=UNIT, context=CTX, analysis_id=AID)
    proj = assoc.Project()

    # ---------------------------------------------------------------------------------------- step 1: A / L split
    A = assoc.a_models(proj)                                       # manifest-verified A list
    L = [l for l in (ROOT / "results/part2/lockbox/lockbox_models.txt").read_text().splitlines() if l]
    man = json.loads((ROOT / "results/part2/lockbox/manifest.json").read_text())
    models = data.load_models()
    lin = models.OncotreeLineage
    strata = rd("results/part2/lockbox/strata.tsv")
    bone_row = strata[strata.stratum == CTX].iloc[0]
    pat = models.PatientID.reindex(A + L)
    J["step1"] = dict(
        n_U=man["n_U"], n_A=len(A), n_L=len(L), seed=man["seed"], frac_L=man["frac_L"],
        n_patients_straddling=man["n_patients_straddling"],
        n_patients_A=int(models.PatientID.reindex(A).nunique()), n_patients_L=int(models.PatientID.reindex(L).nunique()),
        n_bone_A=int((lin.reindex(A) == "Bone").sum()), n_bone_L=int((lin.reindex(L) == "Bone").sum()),
        n_bone_U_strata=int(bone_row.n_U), target_bone_L=int(bone_row.target_L),
        bone_diseases_A=lin.reindex(A).eq("Bone").pipe(lambda s: models.OncotreePrimaryDisease.reindex(A)[s.values])
        .value_counts().to_dict(),
        bone_diseases_L=lin.reindex(L).eq("Bone").pipe(lambda s: models.OncotreePrimaryDisease.reindex(L)[s.values])
        .value_counts().to_dict(),
        note="L counts are model labels only (Model.csv); no L dependency value is read.")
    assert J["step1"]["n_bone_A"] == 44 and J["step1"]["n_bone_L"] == 15
    del pat

    # ---------------------------------------------------------------------------------------- step 2: complex map
    genes = data.protein_genes()
    blk = scores.load_block(A_BLOCK)
    assert blk.meta["model_set_hash"] == "bd28f78d612bef9f"
    keys = genes.depmap_column.values[blk.s_index]
    syms = np.array([symbol(k) for k in keys])
    V = blk.vals                                                    # S x S float32, NaN diagonal
    clusters = rd("results/part2/construction/zavg_map/clusters.tsv", dtype={"parent": str})
    crow = clusters[clusters.cluster_id == UNIT].iloc[0]
    members = crow.members.split(";")
    loc = {s: i for i, s in enumerate(syms)}
    mi = np.array([loc[s] for s in members])
    sub = V[np.ix_(mi, mi)].astype(float)
    I_A, _ = select_topk(V[mi], K)                                  # top-10 within S of each member row (A)
    Ifull, Vfull = select_topk(V[mi], 60)

    def rank_in_row(i_loc, j_loc):
        row = V[i_loc].astype(float)
        v = row[j_loc]
        return int(np.sum(np.nan_to_num(row, nan=-np.inf) > v)) + 1

    pairs = []
    for a in range(len(mi)):
        for b in range(a + 1, len(mi)):
            s = float(sub[a, b])
            ra, rb = rank_in_row(mi[a], mi[b]), rank_in_row(mi[b], mi[a])
            pairs.append(dict(a=members[a], b=members[b], score=f(s, 5), rank_of_b_in_a=ra, rank_of_a_in_b=rb,
                              mutual_topk=bool(ra <= K and rb <= K), above_z0=bool(s >= Z0),
                              edge_A=bool(ra <= K and rb <= K and s >= Z0)))
    # the 20 subsample graphs that fed the consensus
    S_mask = np.zeros(len(genes), bool)
    S_mask[blk.s_index] = True
    boot_edges = {(p["a"], p["b"]): 0 for p in pairs}
    boot_names = [f"boot{s}" for s in runsets.BOOT_SEEDS]
    for nm in boot_names:
        p = scores.cache_path("zavg_BEF", runsets.models(nm), S_mask)
        b = scores.load_block(p)
        assert np.array_equal(b.s_index, blk.s_index)
        Vb = b.vals
        Ib, _ = select_topk(Vb[mi], K)
        for a in range(len(mi)):
            for c in range(a + 1, len(mi)):
                s = Vb[mi[a], mi[c]]
                if s >= Z0 and mi[c] in Ib[a] and mi[a] in Ib[c]:
                    boot_edges[(members[a], members[c])] += 1
        del b, Vb
    for p_ in pairs:
        p_["edge_in_subsample_graphs"] = boot_edges[(p_["a"], p_["b"])]
    # co-clustering frequencies from the cached stage-1 Leiden runs (20 graphs x 5 seeds per gamma)
    lv = np.load(LEVELS)
    gam = list(lv["gammas"])
    cocl = []
    for j, g in enumerate(gam):
        R = lv[f"runs_{j}"][:, mi]
        all5 = float(np.mean([len(set(r)) == 1 for r in R]))
        pf = [float(np.mean(R[:, a] == R[:, c])) for a in range(len(mi)) for c in range(a + 1, len(mi))]
        cocl.append(dict(j=j, gamma=f(g, 4), frac_runs_all_five_together=f(all5, 4), min_pair_freq=f(min(pf), 4),
                         mean_pair_freq=f(np.mean(pf), 4)))
    # near-miss genes: highest mean A score to the members, outside the complex
    memb = rd("results/part2/construction/zavg_map/membership.tsv")
    cl_of = dict(zip(memb.symbol, memb.cluster_id))
    meanto = np.nanmean(V[:, mi].astype(float), axis=1)
    meanto[mi] = -np.inf
    cand = np.argsort(-np.nan_to_num(meanto, nan=-np.inf))[:6]
    near = []
    R2 = lv["runs_2"]                                               # gamma index of the frozen cluster (j_start = 2)
    for g in cand:
        Ig, _ = select_topk(V[[g]], K)
        sc = V[g, mi].astype(float)
        mut = [bool(g in I_A[a] and mi[a] in Ig[0]) for a in range(len(mi))]
        edges = [bool(m_ and s >= Z0) for m_, s in zip(mut, sc)]
        own_top = [syms[x] for x in Ig[0][:5]]
        near.append(dict(gene=syms[g], cluster=cl_of.get(syms[g], ""), mean_score_to_members=f(meanto[g], 4),
                         max_score_to_members=f(np.nanmax(sc), 4),
                         n_members_scoring_ge_z0=int(np.sum(sc >= Z0)),
                         n_mutual_topk_with_members=int(sum(mut)), n_edges_to_members_A=int(sum(edges)),
                         rank_of_best_member_in_its_row=int(min(rank_in_row(g, m) for m in mi)),
                         its_own_top5=own_top,
                         coclustering_with_members_gamma0p1=f(np.mean([np.mean(R2[:, g] == R2[:, m]) for m in mi]), 3)))
    # each member's top-10 within S on A (who the members are near)
    top10 = {members[a]: [syms[x] for x in I_A[a]] for a in range(len(mi))}
    # locus flag inputs
    gl = locus.gene_locations(genes)
    gl_S = gl.iloc[blk.s_index].reset_index(drop=True)
    lf = locus.locus_flag(mi, gl_S.arm.values, gl_S.band.values, chrom=gl_S.chrom.values, pos=gl_S.pos.values)
    att = rd("results/part2/construction/zavg_map/attachments.tsv")
    att = att[att.cluster == UNIT]
    tiers = rd("results/part2/construction/tiers/zavg_tiers.tsv")
    trow = tiers[tiers.cluster_id == UNIT].iloc[0]
    var = rd("results/part2/construction/variability/zavg_BEF_per_complex_jaccard.tsv")
    vrow = var[var.cluster_id == UNIT].iloc[0]
    sel = json.loads((ROOT / "results/part2/construction/selection/selection.json").read_text())
    J["step2"] = dict(
        members=members, z0=Z0, k=K, n_S=int(len(syms)), n_level1=int((clusters.level == 1).sum()),
        pair_scores=pairs, member_score_matrix=[[f(x, 4) for x in r] for r in sub],
        members_top10_in_S_on_A=top10,
        n_edges_A_graph=int(sum(p["edge_A"] for p in pairs)), n_pairs=len(pairs),
        n_subsample_graphs=len(boot_names), coclustering_by_gamma=cocl, gamma_grid=[f(g, 4) for g in gamma_grid()],
        cluster_row=dict(gamma=f(crow.gamma), j_start=int(crow.j_start), j_end=int(crow.j_end),
                         n_adjacent=int(crow.n_adjacent), support=f(crow.support), stability=f(crow.stability),
                         cohesion=f(crow.cohesion), n_pan=int(crow.n_pan), n_selective=int(crow.n_selective),
                         locus=bool(crow.locus), top_arm=crow.top_arm, arm_share=f(crow.arm_share),
                         window_share=f(crow.window_share)),
        cohesion_recomputed=f(np.nanmean(sub[np.triu_indices(len(mi), 1)]), 6),
        locus_recomputed=dict(lf, member_locations=dict(zip(members, gl_S.location.values[mi].tolist())),
                              member_chrom=dict(zip(members, gl_S.chrom.values[mi].tolist()))),
        near_miss=near,
        attached=att[["symbol", "x", "size", "q", "locus_attach"]].to_dict("records"),
        tier=dict(tier=trow.tier, bioplex_edges=int(trow.edges), bioplex_pairs=int(trow.n_pairs),
                  null_mean_edges=f(trow.null_mean_edges), q_bioplex=f(trow.q_bioplex), mean_prot_r=f(trow.mean_r),
                  q_prot=f(trow.q_prot), note="validation input, read after the map was frozen"),
        variability=dict(jaccard_median=f(vrow.jaccard_median), jaccard_q25=f(vrow.jaccard_q25),
                         jaccard_q75=f(vrow.jaccard_q75), frac_ge_0_5=f(vrow.frac_ge_0_5),
                         random_p95=f(vrow.random_S_p95)),
        selection=dict(selected=sel["zavg_BEF"].get("selected"),
                       grid=[dict(z0=g_["z0"], k=g_["k"], stat=f(g_["stat"], 4), sd=f(g_["stat_sd"], 3))
                             for g_ in sel["zavg_BEF"]["grid"]]))
    del V, blk, lv

    # ---------------------------------------------------------------------------------------- step 3: score
    E_A = assoc.load_A_effects(proj=proj)
    assert set(E_A.index) <= set(A) and len(E_A) == 1007
    fmap = FrozenMap.from_membership(rd("results/part2/fzx/zavg_engine_membership.tsv"))
    ctx_full = pd.read_parquet(ROOT / "data_cache/part2/contexts/context_table.parquet")
    ctx = ctx_full.reindex(E_A.index)
    annot = cscore.load_annotation(E_A.columns)
    rem = cscore.removal_for_context(CTX)
    Ysc, info, _ = cscore.unit_scores(E_A, fmap, [UNIT], rem, annot)
    y = Ysc[UNIT]
    c = ctx[CTX]
    cov = covmod.source_table(ctx_full, E_A.index)
    mem_keys = fmap.members[UNIT]
    raw = E_A[mem_keys].astype(float)
    ex = pd.DataFrame(dict(score=y, bone=c, name=models.StrippedCellLineName.reindex(E_A.index),
                           disease=models.OncotreePrimaryDisease.reindex(E_A.index),
                           lineage=models.OncotreeLineage.reindex(E_A.index)))
    for k_ in mem_keys:
        ex[symbol(k_)] = raw[k_]
    ex.index.name = "ModelID"
    OUT.mkdir(exist_ok=True)
    ex.to_csv(OUT / "worked_example_scores_A.tsv", sep="\t", float_format="%.5g")
    inb, outb = ex[ex.bone == 1], ex[ex.bone == 0]
    pick = pd.concat([inb.nsmallest(3, "score"), inb.nlargest(2, "score"),
                      outb.iloc[(outb.score - outb.score.median()).abs().argsort()[:2]],
                      outb.nsmallest(1, "score"), outb.nlargest(1, "score")])
    show = pick[["name", "lineage", "disease", "score"] + [symbol(k_) for k_ in mem_keys]].reset_index()
    desc = lambda s: dict(n=int(s.notna().sum()), mean=f(s.mean(), 4), median=f(s.median(), 4), sd=f(s.std(), 4),  # noqa: E731
                          q10=f(s.quantile(.1), 4), q90=f(s.quantile(.9), 4), min=f(s.min(), 4), max=f(s.max(), 4))
    by_dis = {d: desc(g.score) for d, g in inb.groupby("disease")}
    J["step3"] = dict(
        context_definition="OncotreeLineage == Bone vs all other A lines (one-vs-rest); lineage contexts remove no member",
        context_counts=dict(n_context=int((c == 1).sum()), n_ref=int((c == 0).sum()), n_nan=int(c.isna().sum())),
        removal=dict(rule=str(rem), removed=info.removed.iloc[0], n_used=int(info.n_used.iloc[0])),
        score_rule="mean of members' CRISPRGeneEffect, each clipped to [-3, 3]; NaN if < 50% of members observed",
        n_score_nan=int(y.isna().sum()),
        example_lines=[{k_: (f(v, 4) if isinstance(v, (float, np.floating)) else v) for k_, v in r.items()}
                       for r in show.to_dict("records")],
        dist_in=desc(inb.score), dist_out=desc(outb.score), by_bone_disease=by_dis,
        raw_mean_difference=f(inb.score.mean() - outb.score.mean(), 4),
        member_mean_effect_A={symbol(k_): f(raw[k_].clip(-3, 3).mean(), 4) for k_ in mem_keys},
        member_frac_below_m05_A={symbol(k_): f((raw[k_] < -0.5).mean(), 4) for k_ in mem_keys})

    # ---------------------------------------------------------------------------------------- step 4: association
    entry, _ = registry.load_entry(proj.registry_dir, AID, root=proj.root)
    covs, dropped = covmod.per_context(entry, [CTX], cov)
    X, names, ok, ddrop = assoc._design(c.reindex(y.index), covs[CTX].reindex(y.index), return_dropped=True)
    yy = y.values[ok]
    m = ~np.isnan(yy)
    d = stats.HC3Design(X[m], 1)
    b, se, t = d.fit(yy[m])
    p = float(stats.pvalue_t(t, d.df, "two-sided"))
    # committed discovery row (and its row index = the engine's permutation seed offset)
    tests = tsv_gz("results/part2/directions/H3a_lineage_atlas/tests.tsv.gz")
    pos = int(np.flatnonzero((tests.context == CTX) & (tests.unit == UNIT))[0])
    cr = tests.iloc[pos]
    ids = y.index[ok][m]
    tp = time.time()
    p_perm = assoc.permutation_p(yy[m], X[m], ids, "two-sided", N_PERM, np.random.default_rng(PERM_SEED + pos),
                                 cov.PatientID, None)
    perm_s = time.time() - tp
    # exceedance count behind the GPD tail (same permutations, re-drawn with the same seed)
    bl = cov.PatientID.reindex(ids)
    bl = bl.where(bl.notna(), pd.Series(["__" + str(i) for i in ids], index=ids)).astype(str).values
    perms = stats.blocked_permutations(N_PERM, np.random.default_rng(PERM_SEED + pos), blocks=bl, strata=None,
                                       n=len(ids))
    Z = np.delete(X[m], 1, axis=1)
    Qz, _ = np.linalg.qr(Z)
    r_ = yy[m] - Qz @ (Qz.T @ yy[m])
    Rp = r_[perms].T
    tnull = (d.g @ Rp) / np.sqrt((d.w[:, None] ** 2 * d.resid(Rp) ** 2).sum(0))
    n_exceed = int(np.sum(np.abs(tnull) >= abs(t)))
    np.save(OUT / "worked_example_perm_tnull.npy", tnull.astype(np.float32))
    # Benjamini-Bogomolov from the committed per-test p-values (p_final = permutation p where re-tested)
    ev = tests[tests.evaluable.astype(str) == "True"].copy()
    bb = stats.benjamini_bogomolov(ev.p_final.values, ev.context.values, 0.05)
    ev["q_re"], ev["rej_re"] = bb["adj"], bb["reject"]
    me = ev[(ev.context == CTX) & (ev.unit == UNIT)].iloc[0]
    fam = bb["families"]
    bf = fam[fam.family == CTX].iloc[0]
    bb_an = stats.benjamini_bogomolov(ev.p.values, ev.context.values, 0.05)
    q_an = float(bb_an["adj"][np.flatnonzero((ev.context == CTX).values & (ev.unit == UNIT).values)[0]])
    from scipy import stats as sst
    lo, hi = b - sst.t.ppf(0.975, d.df) * se, b + sst.t.ppf(0.975, d.df) * se
    J["step4"] = dict(
        model="score ~ 1 + context + covariates (OLS); HC3 SE; Satterthwaite df; two-sided",
        covariates_registered=entry["covariates"], covariates_dropped_by_rule=dropped[CTX],
        design_columns=names, n_design_columns=len(names), design_dropped=ddrop,
        n=int(m.sum()), n_context=int(X[m, 1].sum()),
        refit=dict(beta=f(b), se_hc3=f(se), t=f(t), df=f(d.df, 5), p=f(p), ci95=[f(lo), f(hi)]),
        committed=dict(beta=f(cr.beta), se=f(cr.se), t=f(cr.t), df=f(cr.df, 5), p=f(cr.p), p_perm=f(cr.p_perm),
                       q_adj_analytic=f(cr.q_adj), q_final=f(cr.q_final), discovery=str(cr.discovery),
                       ci95=[f(cr.ci95_lo), f(cr.ci95_hi)], row_index=pos),
        reproduces=dict(beta=bool(abs(b - cr.beta) < 1e-5), t=bool(abs(t - cr.t) < 1e-4),
                        df=bool(abs(d.df - cr.df) < 1e-3), p_perm=bool(np.isclose(p_perm, cr.p_perm, rtol=1e-4)),
                        q=bool(np.isclose(me.q_re, cr.q_final, rtol=1e-4))),
        permutation=dict(n_perm=N_PERM, seed=PERM_SEED + pos, scheme="Freedman-Lane; whole PatientID blocks exchanged; "
                         "unrestricted across lineages (lineage context)", p_perm_refit=f(p_perm),
                         n_null_abs_t_ge_obs=n_exceed, p_empirical=f((n_exceed + 1) / (N_PERM + 1)),
                         max_abs_null_t=f(np.max(np.abs(tnull)), 4), gpd_used=bool(p_perm < (n_exceed + 1) / (N_PERM + 1)),
                         runtime_s=round(perm_s, 2)),
        fdr=dict(method="Benjamini-Bogomolov, families = 23 contexts, q = 0.05", M=int(bb["M"]), R=int(bb["R"]),
                 family_simes_p=f(bf.simes_p), family_q=f(bf.family_q), family_n_tests=int(bf.n_tests),
                 family_n_rejected=int(bf.n_rejected), within_family_level=f(0.05 * bb["R"] / bb["M"]),
                 q_recomputed=f(me.q_re), q_analytic_recomputed=f(q_an), n_discoveries_total=int(bb["reject"].sum()),
                 m_evaluable=int(len(ev)), rank_in_bone_family_by_p=int((ev[ev.context == CTX].p_final < me.p_final).sum()) + 1))

    # ---------------------------------------------------------------------------------------- step 5: calibration
    coarse = fmap.coarse()
    ui = coarse.index(UNIT)
    ess = calib.clip_essentiality(E_A)
    blk = scores.load_block(A_BLOCK)
    ps = pd.DataFrame(blk.vals, index=keys, columns=keys)
    del blk
    sets, status = matched_sets_one_unit(E_A, fmap, coarse, ui, ess, ps, R_SETS, CAL_SEED)
    Ys = calib.score_sets(E_A, sets, rem, annot)
    ft = calib.test_columns(Ys, c, covs[CTX], "two-sided")
    cal = tsv_gz("results/part2/directions/H3a_lineage_atlas/calibration_D8_tests.tsv.gz")
    cb = cal[(cal.context == CTX)]
    ms = cb[cb.unit == UNIT].set_index("set")
    nev = cb[cb.calibration == "never_expressed"]
    cj = json.loads((ROOT / "results/part2/directions/H3a_lineage_atlas/calibration_D8.json").read_text())
    pj = json.loads((ROOT / "results/part2/directions/H3a_lineage_atlas/calibration.json").read_text())
    pc = next(r for r in pj["per_context"] if r["context"] == CTX)
    allm = cb[cb.calibration == "matched_random_set"]
    J["step5"] = dict(
        matched_rule="10 sets drawn from other frozen complexes' members, matched on size (5), mean clipped A effect "
                     "(+/-0.1) and mean pairwise zavg_BEF on A (+/-1.0); seed default_rng([4600, unit index])",
        unit_index=ui, target=status,
        matched_sets=[dict(set=k_, genes=[symbol(g) for g in v], t_refit=f(ft.loc[k_, "t"], 5),
                           t_committed=f(ms.loc[k_, "t"], 5), p_committed=f(ms.loc[k_, "p"], 4))
                      for k_, v in sets.items()],
        matched_reproduces=bool(np.allclose(ft.t.values, ms.loc[list(sets), "t"].values, rtol=1e-4)),
        t_observed=f(t, 5),
        n_matched_with_abs_t_ge_obs=int((ms.t.abs() >= abs(t)).sum()), rank_of_observed=int((ms.t.abs() >= abs(t)).sum()) + 1,
        matched_bone_all_units_frac_p_lt_005=f((allm.p < 0.05).mean(), 4),
        real_bone_frac_p_lt_005=f(cj["matched_random_sets"]["per_context"][CTX]["real"]["frac_p_lt_005"], 4),
        random_bone_frac_p_lt_005=f(cj["matched_random_sets"]["per_context"][CTX]["random"]["frac_p_lt_005"], 4),
        never_expressed=dict(n_sets=int(len(nev)), frac_p_lt_005=f((nev.p < 0.05).mean(), 4),
                             t_values=[f(x, 4) for x in nev.t.values], max_abs_t=f(nev.t.abs().max(), 4),
                             pooled_all_contexts=f(cj["never_expressed_sets"]["pooled"]["frac_p_lt_005"], 4)),
        pseudo_contexts=dict(n_pseudo=pc["n_pseudo"], n_tests=pc["n_tests"], frac_p_lt_005=f(pc["frac_p_lt_005"], 4),
                             pooled_all_contexts=f(pj["frac_p_lt_005"], 4)))
    pd.DataFrame(dict(t=allm.t.values)).to_csv(OUT / "worked_example_matched_t_bone_all_units.tsv", sep="\t",
                                                index=False, float_format="%.5g")

    # ---------------------------------------------------------------------------------------- step 6: lockbox etc.
    conf = rd("results/part2/lockbox_confirm/H3a_lineage_atlas/confirmation.tsv")
    lr = conf[(conf.context == CTX) & (conf.unit == UNIT)].iloc[0]
    with gzip.open(ROOT / "results/part2/lockbox_confirm/primary_relabelled_R7.tsv.gz", "rt") as fh:
        r7 = pd.read_csv(fh, sep="\t")
    r7r = r7[(r7.entry == AID) & (r7.context == CTX) & (r7.unit == UNIT)].iloc[0]
    rob = rd("results/part2/robustness/robustness_per_hit.tsv")
    rr = rob[(rob.context == CTX) & (rob.unit == UNIT)].iloc[0]
    sf = rd("results/part2/robustness/strongest_findings.tsv")
    sr = sf[(sf.context == CTX) & (sf.unit == UNIT)].iloc[0]
    lead = rd("results/part2/revision/r2_leads_corrected.tsv")
    lrow = lead[(lead.context == CTX) & (lead.unit == UNIT)].iloc[0]
    d2 = rd("results/part2/revision/r3_demeter2_per_hit.tsv")
    drow = d2[(d2.context == CTX) & (d2.unit == UNIT)].iloc[0]
    bench = rd("results/part2/lockbox_confirm/BENCH_complex_vs_gene/bench_pairs.tsv")
    brow = bench[(bench.context == CTX) & (bench.unit == UNIT)].iloc[0]
    rk = rd("results/part2/revision/r5_ranking_alternatives.tsv")
    rkr = rk[(rk.context == CTX) & (rk.unit == UNIT)]
    # A-side robustness re-run with the engine functions (A only)
    groups = cov.OncotreePrimaryDisease
    lo_ = robust.loso(y, c, covs[CTX], groups, "two-sided")
    t2 = robust.leave_top2_out(y, c, covs[CTX], "two-sided")
    br = robust.breadth(E_A, mem_keys, c, covs[CTX], np.sign(b))
    per_member = {}
    for k_ in mem_keys:
        fo = robust.fit_one(E_A[k_].astype(float).clip(-3, 3), c, covs[CTX], "less" if b < 0 else "greater")
        per_member[symbol(k_)] = dict(beta=f(fo["beta"], 4), t=f(fo["t"], 4), p_one_sided=f(fo["p"], 3))
    J["step6"] = dict(
        lockbox=dict(source="results/part2/lockbox_confirm/H3a_lineage_atlas/confirmation.tsv",
                     beta_L=f(lr.beta_L), se_L=f(lr.se_L), df_L=f(lr.df_L, 5), n_L=int(lr.n_L),
                     n_context_L=int(lr.n_context_L), p_L_hc3=f(lr.p_L_hc3), p_L_perm=f(lr.p_L_perm), p_L=f(lr.p_L),
                     ci95_L=[f(lr.ci95_L_lo), f(lr.ci95_L_hi)], same_sign=str(lr.same_sign), status=lr.status,
                     primary_call_R7=r7r.primary_call, corrected_call_R7=r7r.corrected_call,
                     defect_sensitive=str(r7r.defect_sensitive)),
        robustness_committed=dict(loso_unit=rr.loso_unit, loso_folds=int(rr.loso_folds), loso_pass=int(rr.loso_pass),
                                  loso_skipped=int(rr.loso_skipped), top2_beta=f(rr.top2_beta), top2_p=f(rr.top2_p),
                                  breadth=f(rr.breadth), breadth_n=int(rr.breadth_n)),
        robustness_refit_A=dict(loso=lo_, top2={k_: (f(v) if not isinstance(v, (bool, np.bool_)) else bool(v))
                                                for k_, v in t2.items()}, breadth=br, per_member=per_member),
        strongest=dict(strongest=str(rr.strongest), rank=int(sr["rank"]), tier=sr.tier,
                       rank_alt_abs_t=int(rkr.iloc[0].rank_alt1_abs_t) if len(rkr) else None,
                       beats_all_matched_sets=str(sr.beats_all_matched_sets), n_strongest_total=int(len(sf))),
        scale_artefact=dict(essentiality=f(lrow.essentiality, 4), fitted_beta_nl=f(lrow.fitted_nl, 4),
                            residual_beta_nl=f(lrow.residual_beta_nl, 4), residual_z_nl=f(lrow.residual_z_nl, 4),
                            concordant_nl=str(lrow.concordant_nl), scale_artefact_possible=str(lrow.scale_artefact_possible),
                            lead_corrected=str(lrow.lead_corrected)),
        demeter2=dict(beta_D2=f(drow.beta_D2), se_D2=f(drow.se_D2), t_D2=f(drow.t_D2), df_D2=f(drow.df_D2, 4),
                      p_dir_D2=f(drow.p_dir_D2), n_D2=int(drow.n_D2), n_context_D2=int(drow.n_context_D2),
                      same_sign=str(drow.same_sign)),
        bench=dict(gstar=brow.gstar, delta_tL_complex_minus_gstar=f(brow.delta), tercile=brow.tercile,
                   breadth=f(brow.breadth)))
    if len(rkr):
        J["step6"]["strongest"]["ranks_alternatives"] = rkr.iloc[0].to_dict()

    J["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "worked_example.json").write_text(json.dumps(J, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    print(json.dumps({k_: J["step4"][k_] for k_ in ("refit", "committed", "reproduces", "permutation")}, indent=1, default=str))
    print("matched sets reproduce:", J["step5"]["matched_reproduces"], "runtime", J["runtime_s"])
    import make_figures
    make_figures.worked_figures(J, ex)


def matched_sets_one_unit(E_A, fmap, coarse, ui, essentiality, pair_scores, R, seed, tol_ess=0.1, tol_coh=1.0,
                          max_attempts=20_000):
    """calib.matched_random_sets restricted to one unit (same pool, same rng = default_rng([seed, ui]), same rule).
    The registered function loops over all 367 units; restricting it keeps the seed of this unit unchanged."""
    from part2.h1 import coherence
    u = coarse[ui]
    allg = sorted({g for v in coarse for g in fmap.members[v] if g in E_A.columns})
    gi = {g: i for i, g in enumerate(allg)}
    ess_v = essentiality.reindex(allg).values.astype(float)
    PS = pair_scores.reindex(index=allg, columns=allg).values.astype(float)
    pools = {v: np.array([gi[g] for g in fmap.members[v] if g in gi]) for v in coarse}
    mem = [g for g in fmap.members[u] if g in gi]
    k = len(mem)
    idx = np.array([gi[g] for g in mem])
    t_ess = float(np.nanmean(ess_v[idx]))
    t_coh = coherence(mem, pair_scores)
    iu = np.triu_indices(k, 1)
    others = [pools[v] for v in coarse if v != u and len(pools[v])]
    big = [p for p in others if len(p) >= k]
    rng = np.random.default_rng([seed, ui])
    status, sets = "unavailable", []
    for mult in (1.0, 2.0):
        seen, sets = set(), []
        for _ in range(max_attempts):
            if big:
                cand = big[rng.integers(len(big))]
            else:
                break
            s = np.sort(rng.choice(cand, k, replace=False))
            key = s.tobytes()
            if key in seen:
                continue
            if abs(float(np.nanmean(ess_v[s])) - t_ess) > tol_ess * mult:
                continue
            if np.isfinite(t_coh) and not abs(float(np.nanmean(PS[np.ix_(s, s)][iu])) - t_coh) <= tol_coh * mult:
                continue
            seen.add(key)
            sets.append(s)
            if len(sets) >= R:
                break
        if len(sets) >= R:
            status = "ok" if mult == 1.0 else "tolerances doubled"
            break
    out = {f"{u}#r{r}": [allg[i] for i in s] for r, s in enumerate(sets)}
    return out, dict(unit=u, size=k, target_ess=float(t_ess), target_coh=float(t_coh), status=status)


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    main()
