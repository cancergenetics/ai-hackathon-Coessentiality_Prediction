"""End-to-end synthetic pipeline for the Part 2 association engine (helper for tests/test_part2_e2e.py; also runnable:
``.venv/bin/python tests/part2_synth_pipeline.py [outdir]`` prints metrics and stage runtimes as JSON).

Builds a throw-away git project from part2.synth.simulate (freeze manifests, context table, frozen map, registry
entries), then runs exactly the real sequence: controls (one run, kill tally) -> H3a, H3b discovery (BB, permutation
re-entry, robustness) -> committed trials -> H2 (stage-wise inside discovered parents) -> H1 (events, pooled BH,
classes with TOST and matched null) -> per entry: committed target list, lockbox read, one-sided L test, close ->
BENCH read (gated: only after every confirmation is closed). No DepMap data is touched.
"""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, bench, confirm, controls, covariates as covmod, cscore, h1, h2, registry, robust, synth  # noqa
from part2.cmap import symbol  # noqa: E402

REL_CTX = "data_cache/part2/contexts/context_table.parquet"


# ------------------------------------------------------------------------------------------------ repo
def _git(root, *a):
    r = registry.git(root, *a, check=False)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return r


def commit_all(root, msg):
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", msg)


def build_repo(root: Path, sim: synth.Synth, events=None) -> dict:
    """Throw-away git project with FZ-LB / context / FZ-C manifests and the synthetic registry entries."""
    from part2 import h1 as h1m, h2 as h2m
    root.mkdir(parents=True, exist_ok=True)
    _git(root, "init", "-q")
    _git(root, "config", "user.email", "synthetic@example.org")
    _git(root, "config", "user.name", "synthetic")
    lb = root / "results" / "part2" / "lockbox"
    lb.mkdir(parents=True)
    A = sorted(sim.fold.index[sim.fold == "A"])
    L = sorted(sim.fold.index[sim.fold == "L"])
    tA, tL = "\n".join(A) + "\n", "\n".join(L) + "\n"
    (lb / "A_models.txt").write_text(tA)
    (lb / "lockbox_models.txt").write_text(tL)
    (lb / "manifest.json").write_text(json.dumps(dict(files={"A_models.txt": assoc._sha_text(tA),
                                                             "lockbox_models.txt": assoc._sha_text(tL)})))
    cache = root / "data_cache" / "part2" / "contexts"
    cache.mkdir(parents=True)
    tab = pd.concat([sim.contexts, sim.meta.drop(columns=["TP53_altered", "WGD"]),
                     sim.lof.add_prefix("H1loss:")], axis=1)
    tab.to_parquet(cache / "context_table.parquet")
    cdir = root / "results" / "part2" / "contexts"
    cdir.mkdir(parents=True)
    rel_ctx = "data_cache/part2/contexts/context_table.parquet"
    (cdir / "checksums.json").write_text(json.dumps(dict(files={rel_ctx: assoc._sha_file(root / rel_ctx)})))
    mdir = root / "results" / "part2" / "map"
    mdir.mkdir(parents=True)
    sim.fmap.to_membership().to_csv(mdir / "membership.tsv", sep="\t", index=False)
    (mdir / "manifest.json").write_text(json.dumps(dict(files={"membership.tsv": assoc._sha_file(mdir / "membership.tsv")})))
    fin = ["results/part2/contexts/checksums.json", "results/part2/lockbox/manifest.json",
           "results/part2/map/manifest.json"]
    reg = root / "docs" / "registry"
    reg.mkdir(parents=True)
    coarse = sim.fmap.coarse()
    h3a_ctx = [dict(id=f"lineage:{l}", type="lineage") for l in synth.LINEAGES]
    (reg / "H3a_lineage_atlas.yaml").write_text(synth.entry_yaml(
        "H3a_lineage_atlas", "H3a", h3a_ctx, coarse, "S_C_mean_winsor3", synth.H3A_COV,
        dict(method="benjamini_bogomolov", families="contexts", q=0.05), fin))
    (reg / "H3b_driver_family.yaml").write_text(synth.entry_yaml(
        "H3b_driver_family", "H3b", [dict(id=c, covariates=synth.H3B_COV[c]) for c in synth.H3B_CONTEXTS], coarse,
        "S_C_mean_winsor3", synth.BASE_COV, dict(method="benjamini_bogomolov", families="contexts", q=0.05), fin))
    subs = h2m.h2_units(sim.fmap)
    h2_ctx = [dict(id=f"lineage:{l}", type="lineage", covariates=synth.H3A_COV + ["NNMD", "CasActivity"]) for l in synth.LINEAGES]
    h2_ctx += [dict(id=c, covariates=synth.H3B_COV[c] + ["NNMD", "CasActivity"]) for c in synth.H3B_CONTEXTS]
    (reg / "H2_submodule_divergence.yaml").write_text(synth.entry_yaml(
        "H2_submodule_divergence", "H2", h2_ctx, subs, "submodule_vs_sibling", synth.BASE_COV + ["NNMD", "CasActivity"],
        dict(method="stagewise_bb", parent_directions=["H3a", "H3b"], q=0.05), fin))
    if events is not None:
        h1_ctx = [dict(id=e.context, covariates=e.covariates) for e in {e.context: e for e in events}.values()]
        (reg / "H1_within_complex_GI.yaml").write_text(synth.entry_yaml(
            "H1_within_complex_GI", "H1", h1_ctx, h1m.units_of(events), "member_resolved", synth.BASE_COV,
            dict(method="bh", q=0.05), fin))
    nm = sim.named
    wrn = next(k for k in sim.fmap.members[nm["WRN"]] if symbol(k) == "WRN")
    myb = next(k for k in sim.fmap.members[nm["MYB"]] if symbol(k) == "MYB")
    a1b = next(k for k in sim.fmap.members[nm["EVA"]] if symbol(k) == "ARID1B")
    cdk = {symbol(k): k for k in sim.fmap.members[nm["CDK"]]}
    (reg / "CTRL_PC1.yaml").write_text(synth.ctrl_yaml("CTRL_PC1", "PC1", "MSI_H", "+", [wrn], nm["WRN"],
                                                  "S_C_mean_winsor3", "H3b", synth.BASE_COV, fin))
    (reg / "CTRL_PC20.yaml").write_text(synth.ctrl_yaml("CTRL_PC20", "PC20", "heme", "+", [myb], nm["MYB"],
                                                   "S_C_mean_winsor3", "H3a", synth.H3A_COV, fin, ctx_type="lineage"))
    (reg / "CTRL_PC6.yaml").write_text(synth.ctrl_yaml("CTRL_PC6", "PC6", "H1loss:ARID1A", "+", [a1b], nm["EVA"],
                                                  "member_resolved", "H1", synth.MUT_COV, fin))
    (reg / "CTRL_PC21.yaml").write_text(synth.ctrl_yaml(
        "CTRL_PC21", "PC21", "heme_vs_breast", "+", [f"{cdk['CDK6']} minus {cdk['CDK4']}"], f"{nm['CDK']}.S1",
        "submodule_vs_sibling", "H2", synth.H3A_COV + ["NNMD", "CasActivity"], fin, ctx_type="lineage"))
    (reg / "BENCH_complex_vs_gene.yaml").write_text(synth.entry_yaml(
        "BENCH_complex_vs_gene", "benchmark", [dict(id=c) for c in synth.H3B_CONTEXTS], coarse, "S_C_mean_winsor3",
        synth.BASE_COV, dict(method="none", q=0.05), fin))
    commit_all(root, "synthetic freeze: lockbox, contexts, map, registry")

    def loader(model_ids, genes, source):
        return sim.E.reindex(index=list(model_ids), columns=list(genes) if genes is not None else sim.E.columns)

    return dict(proj=assoc.Project(root=root, outcome_loader=loader), inputs=[rel_ctx, "results/part2/map/membership.tsv"])


# ------------------------------------------------------------------------------------------------ pipeline
def h1_inputs(sim: synth.Synth):
    z = sim.lof * 0
    loss, kind = h1.loss_calls(sim.lof, z, sim.lof, z + 5.0, sim.meta.OncotreeLineage)
    events = h1.enumerate_events(sim.fmap, loss, kind, sim.annot, sim.fold)
    return loss, kind, events


def _scores_for(E, fmap, units, ctxs, annot):
    cache, out = {}, {}
    for c in ctxs:
        rem = cscore.removal_for_context(c)
        if rem not in cache:
            cache[rem] = cscore.unit_scores(E, fmap, units, rem, annot)[0]
        out[c] = cache[rem]
    return out


def _metrics(res, truth, expected=frozenset(), null_units=None, planted=None):
    """False positive: a discovery not in `truth` (nor `expected`); with null_units (H3a one-vs-rest), only
    discoveries on units carrying no planted signal at all count as false (weak true one-vs-rest differences on
    signal-carrying units, |diff| < 0.03, are reported as 'other'). Power over `planted` (default truth)."""
    d = set(map(tuple, res.loc[res.discovery, ["context", "unit"]].values))
    planted = truth if planted is None else planted
    fp = d - truth - set(expected)
    other = set()
    if null_units is not None:
        other = {x for x in fp if x[1] not in null_units}
        fp = fp - other
    tp = d & planted
    return dict(n_discoveries=len(d), true_pos=len(d & truth), false_pos=len(fp), other=len(other),
                n_truth=len(truth), n_planted=len(planted), power=len(tp) / len(planted) if planted else None,
                fdp=len(fp) / max(len(d), 1))


def _h3_metrics(res, sim, key):
    if key == "H3a":
        return _metrics(res, sim.truth["H3a"], null_units=sim.truth["null_units"], planted=sim.truth["H3a_planted"])
    return _metrics(res, sim.truth["H3b"])


def discovery(sim, proj, inputs, n_perm=2000, groups=True, save=True):
    """H3a + H3b discovery runs. Returns dict with results, score frames, covariates."""
    tab = pd.read_parquet(proj.root / REL_CTX)
    E_A = assoc.load_A_effects(proj=proj)
    tabA = tab.loc[E_A.index]
    blocks, strata = tab.PatientID, tab.OncotreeLineage
    out = {}
    for aid, ctxs in (("H3a_lineage_atlas", [f"lineage:{l}" for l in synth.LINEAGES]),
                      ("H3b_driver_family", synth.H3B_CONTEXTS)):
        e, _ = registry.load_entry(proj.registry_dir, aid, root=proj.root)
        Y = _scores_for(E_A, sim.fmap, sim.fmap.coarse(), ctxs, sim.annot)
        covs, _ = covmod.per_context(e, ctxs, tabA)
        res = assoc.run_association(aid, Y, tabA[ctxs], covs, inputs=inputs, proj=proj, blocks=blocks, strata=strata,
                                    n_perm=n_perm, min_context_n=12, min_ref_n=12, release_hits=True, save=save,
                                    perm_seed=4000 if aid.startswith("H3a") else 4100)
        out[aid] = dict(res=res, Y=Y, covs=covs, ctxs=ctxs)
    return E_A, tab, out


def run_pipeline(sim: synth.Synth, root: Path, n_perm: int = 2000, K: int = 1000, B_bench: int = 200) -> dict:
    T, t0 = {}, time.time()
    loss, kind, events = h1_inputs(sim)
    R = build_repo(root, sim, events)
    proj, inputs = R["proj"], R["inputs"]
    T["setup"] = time.time() - t0
    tab = pd.read_parquet(proj.root / REL_CTX)
    E_A = assoc.load_A_effects(proj=proj)
    A = E_A.index
    tabA = tab.loc[A]
    blocks, strata = tab.PatientID, tab.OncotreeLineage
    lossA = loss.loc[A]
    out = dict(truth={k: (sorted(v) if isinstance(v, set) else v) for k, v in sim.truth.items()})

    # ---- controls: one run, kill tally committed before any direction's discovery run
    t = time.time()
    ctab, kill = controls.evaluate_controls(E_A, sim.fmap, tabA, tabA, sim.annot, inputs=inputs, proj=proj)
    synth_commit(proj.root, "controls run (A) + kill tally")
    out["kill_rule"] = {d: v["verdict"] for d, v in kill.items()}
    out["controls_recovered"] = ctab[ctab.level == "C"][["entry", "recovered", "p"]].to_dict("records")
    T["controls"] = time.time() - t

    # ---- H3a / H3b discovery
    t = time.time()
    _, _, disc = discovery(sim, proj, inputs, n_perm=n_perm)
    resA, resB = disc["H3a_lineage_atlas"]["res"], disc["H3b_driver_family"]["res"]
    out["H3a"] = _h3_metrics(resA, sim, "H3a")
    out["H3b"] = _h3_metrics(resB, sim, "H3b")
    T["H3_discovery"] = time.time() - t
    t = time.time()
    rob = {}
    for aid, d in disc.items():
        hits = d["res"][d["res"].discovery]
        grp = {c: (tabA.OncotreePrimaryDisease if c.startswith("lineage:") else tabA.OncotreeLineage) for c in d["ctxs"]}
        rob[aid] = robust.robustness(hits, d["Y"], tabA[d["ctxs"]], d["covs"], grp, E=E_A, members=sim.fmap.members)
    out["robust_loso_frac"] = {k: float(v.loso_robust.mean()) if len(v) else None for k, v in rob.items()}
    T["robustness"] = time.time() - t
    synth_commit(proj.root, "H3a/H3b discovery runs committed")

    # ---- H2: stage-wise inside discovered parents
    t = time.time()
    subs = h2.h2_units(sim.fmap)
    sel = assoc.parent_selection_from({"H3a": resA, "H3b": resB}, {s: sim.fmap.parent[s] for s in subs}, q=0.05)
    ctx2 = [f"lineage:{l}" for l in synth.LINEAGES] + synth.H3B_CONTEXTS
    e2, _ = registry.load_entry(proj.registry_dir, "H2_submodule_divergence", root=proj.root)
    Y2, coefs2 = h2.h2_frames(E_A, sim.fmap, ctx2, cscore.removal_for_context, sim.annot, units=subs)
    cov2, _ = covmod.per_context(e2, ctx2, tabA)
    res2 = assoc.run_association("H2_submodule_divergence", Y2, tabA[ctx2], cov2, inputs=inputs, proj=proj,
                                 blocks=blocks, strata=strata, n_perm=n_perm, min_context_n=12, min_ref_n=12,
                                 parent_selection=sel, release_hits=True, perm_seed=4300)
    out["H2"] = _metrics(res2, sim.truth["H2"], sim.truth.get("H2_expected_sibling", set()),
                         planted=sim.truth.get("H2_planted", set()))
    out["H2"]["n_eligible_tests"] = int(len(res2))
    T["H2"] = time.time() - t

    # ---- H1
    t = time.time()
    e1, _ = registry.load_entry(proj.registry_dir, "H1_within_complex_GI", root=proj.root)
    Y1, ctx1, pairs, comps = h1.assoc_inputs(E_A, events, lossA)
    cov1, _ = covmod.per_context(e1, list(ctx1.columns), tabA)
    res1 = assoc.run_association("H1_within_complex_GI", Y1, ctx1, cov1, inputs=inputs, proj=proj, blocks=blocks,
                                 strata=strata, n_perm=n_perm, min_context_n=12, min_ref_n=12, pairs=pairs,
                                 release_hits=True, perm_seed=4200)
    Ec = E_A.clip(-3, 3)
    ess = Ec.mean()
    genes_cl = sorted(sim.fmap.all_genes(1))
    pair_scores = 10 * Ec[genes_cl].corr()
    cls = h1.classify(events, comps, res1, E_A, lossA, cov1, sim.fmap, sim.annot, ess, pair_scores, K=K)
    correct = n_axes = 0
    calls = {}
    for r in cls.itertuples(index=False):
        exp = sim.truth["H1"].get(r.event, {})
        for ax, col in (("N", "class_complex_wide"), ("P", "class_paralog")):
            if ax in exp:
                n_axes += 1
                correct += getattr(r, col) == exp[ax]
                calls[f"{r.event}|{ax}"] = (getattr(r, col), exp[ax])
    out["H1"] = dict(accuracy=correct / max(n_axes, 1), n_axes=n_axes, calls=calls,
                     collapse_index={r.event: r.collapse_index for r in cls.itertuples(index=False)})
    T["H1"] = time.time() - t
    synth_commit(proj.root, "H2 / H1 discovery runs committed")

    # ---- confirmations: target list committed -> lockbox read -> one-sided L test -> close
    t = time.time()
    Lids = sim.fold.index[sim.fold == "L"]
    tabL = tab.loc[Lids]
    lossL = loss.loc[Lids]
    # BENCH target list (members of H3b discoveries) is committed now; its read must be refused until every
    # confirmatory entry is read and closed
    hb = resB[resB.discovery]
    confirm.write_target_list(proj, "BENCH_complex_vs_gene", sorted({g for u in hb.unit for g in sim.fmap.members[u]}),
                              hits=hb, contexts=set(hb.context))
    gate_before = None
    try:
        assoc.read_lockbox("BENCH_complex_vs_gene", confirm.target_path(proj, "BENCH_complex_vs_gene"), proj=proj)
    except assoc.FirewallError as ex:
        gate_before = str(ex)[:60]
    out["bench_gate_before_confirmations"] = gate_before
    conf = {}

    def run_conf(aid, hits, members_of, score_L, ctx_L, cov_L):
        genes = sorted({g for u in hits.unit for g in members_of(u)})
        confirm.write_target_list(proj, aid, genes, hits=hits, units={u: members_of(u) for u in set(hits.unit)},
                                  contexts=set(hits.context))
        r = confirm.confirm_entry(aid, hits, score_L, ctx_L, cov_L, proj=proj, blocks=tabL.PatientID,
                                  strata=tabL.OncotreeLineage, lineage_contexts={c for c in ctx_L.columns
                                                                                if c.startswith("lineage:")},
                                  n_perm=n_perm, seed=9000)
        synth_commit(proj.root, f"confirmation {aid}")
        confirm.close_entry(proj, aid)
        conf[aid] = dict(n_hits=int(len(hits)), n_confirmed=int(r.confirmed.sum()) if len(r) else 0)
        return r

    # controls
    for aid, e in controls.control_entries(proj).items():
        if e["status"] == "closed":
            continue
        rows = ctab[(ctab.entry == aid) & ctab.beta.notna()]
        hits = pd.DataFrame(dict(context=rows.context, unit=rows.unit, beta=rows.beta, alternative=rows.alternative))
        cu = e["control"]["complex_unit"]
        _, cinfo = controls.complex_readout(e, E_A, sim.fmap, sim.annot)

        def sL(E_L, c, u, e=e, cu=cu, cinfo=cinfo):
            if u == cu:
                return controls.complex_readout(e, E_L, sim.fmap, sim.annot, fit=cinfo)[0]
            return cscore.gene_unit(E_L, u)

        def mem(u, cu=cu):
            if u == cu:
                return sim.fmap.members[u] + (sim.fmap.members[sim.fmap.parent[u]] if sim.fmap.level[u] == 2 else [])
            return [x.strip() for x in u.split(" minus ")]
        ctxs = [c["id"] for c in e["context_family"]["contexts"]]
        covL, _ = covmod.per_context(e, ctxs, tabL)
        run_conf(aid, hits, mem, sL, tabL[ctxs], covL)
    # H3a / H3b
    for aid, d in disc.items():
        hits = d["res"][d["res"].discovery][["context", "unit", "beta", "alternative"]]
        e, _ = registry.load_entry(proj.registry_dir, aid, root=proj.root)
        covL, _ = covmod.per_context(e, d["ctxs"], tabL)
        run_conf(aid, hits, lambda u: sim.fmap.members[u],
                 lambda E_L, c, u: cscore.unit_scores(E_L, sim.fmap, [u], cscore.removal_for_context(c),
                                                      sim.annot)[0][u],
                 tabL[d["ctxs"]], covL)
    # H2
    hits2 = res2[res2.discovery][["context", "unit", "beta", "alternative"]]
    covL2, _ = covmod.per_context(e2, ctx2, tabL)
    run_conf("H2_submodule_divergence", hits2, lambda u: sim.fmap.members[sim.fmap.parent[u]],
             h2.score_L_fn(sim.fmap, coefs2, cscore.removal_for_context, sim.annot), tabL[ctx2], covL2)
    # H1: every non-neutral class call is confirmed on its defining statistic
    called = []
    for r in cls.itertuples(index=False):
        if r.class_paralog == "PARALOG COMPENSATION":
            called.append((f"H1loss:{r.gene}", f"{r.event}|P", r.beta_P))
        if r.class_complex_wide in ("COLLAPSE", "RESIDUAL SENSITISATION"):
            called.append((f"H1loss:{r.gene}", f"{r.event}|N", r.kappa_N))
    hits1 = pd.DataFrame(called, columns=["context", "unit", "beta"]).assign(alternative="two-sided")
    ev = {e.id: e for e in events}
    ctx1L = pd.DataFrame({f"H1loss:{g}": lossL[g] for g in lossL.columns})
    covL1, _ = covmod.per_context(e1, list(ctx1.columns), tabL)
    run_conf("H1_within_complex_GI", hits1, lambda u: ev[u.rsplit("|", 1)[0]].P + ev[u.rsplit("|", 1)[0]].N,
             h1.score_L_fn(events, comps, lossL), ctx1L, covL1)
    out["confirmation"] = conf
    T["confirmations"] = time.time() - t

    # ---- BENCH (gated): complex vs best gene in L
    t = time.time()
    eB, _ = registry.load_entry(proj.registry_dir, "BENCH_complex_vs_gene", root=proj.root)
    hb = resB[resB.discovery]
    if len(hb):
        rb = rob["H3b_driver_family"].set_index(["context", "unit"])
        gst = []
        for r in hb.itertuples(index=False):
            mem = [g for g in sim.fmap.members[r.unit]
                   if g not in cscore.removal_for_context(r.context).removed([g], sim.annot)]
            gst.append(bench.choose_gstar(E_A, mem, tabA[r.context], disc["H3b_driver_family"]["covs"][r.context],
                                          r.beta))
        pairsB = pd.DataFrame(dict(context=hb.context.values, unit=hb.unit.values, beta_A=hb.beta.values, gstar=gst,
                                   breadth=[rb.loc[(c, u), "breadth"] for c, u in zip(hb.context, hb.unit)]))
        tp = confirm.target_path(proj, "BENCH_complex_vs_gene")
        E_L = assoc.read_lockbox("BENCH_complex_vs_gene", tp, proj=proj)   # target list committed above
        covLB, _ = covmod.per_context(eB, synth.H3B_CONTEXTS, tabL)
        bd = bench.bench_delta(pairsB, lambda E, c, u: cscore.unit_scores(E, sim.fmap, [u], cscore.removal_for_context(c),
                                                                           sim.annot)[0][u],
                               E_L.reindex(columns=sim.E.columns), tabL[synth.H3B_CONTEXTS], covLB,
                               blocks_L=tabL.PatientID, B=B_bench)
        out["bench"] = bd["summary"]
    T["bench"] = time.time() - t
    T["total"] = time.time() - t0
    out["runtime_s"] = {k: round(v, 2) for k, v in T.items()}
    out["trials_lines"] = len(proj.trials.read_text().splitlines())
    return out


def synth_commit(root, msg):
    st = registry.git(root, "status", "--porcelain").stdout.strip()
    if st:
        commit_all(root, msg)


# ------------------------------------------------------------------------------------------------ null / power sims
def fdr_sims(n_sims: int = 50, seed0: int = 10_000, null: bool = True, lineage_null: bool = True, effect: float = 0.3,
             n_perm: int = 1000, workdir: Path | None = None) -> pd.DataFrame:
    rows = []
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        for i in range(n_sims):
            sim = synth.simulate(seed0 + i, null=null, lineage_null=lineage_null, effect=effect)
            R = build_repo(Path(td) / f"s{i}", sim)
            _, _, d = discovery(sim, R["proj"], R["inputs"], n_perm=n_perm, save=False)
            row = dict(seed=seed0 + i)
            for aid, key in (("H3a_lineage_atlas", "H3a"), ("H3b_driver_family", "H3b")):
                m = _h3_metrics(d[aid]["res"], sim, key)
                row.update({f"{key}_{k}": v for k, v in m.items()})
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    with tempfile.TemporaryDirectory() as td:
        t = time.time()
        res = run_pipeline(synth.simulate(0), Path(td) / "e2e")
        res["wall_s"] = round(time.time() - t, 1)
    t = time.time()
    nul = fdr_sims(50, null=True, lineage_null=True)
    nul_conf = fdr_sims(50, seed0=20_000, null=True, lineage_null=False)
    mix = fdr_sims(20, seed0=30_000, null=False, lineage_null=False, effect=0.15)
    res["sims"] = dict(
        global_null_50=dict(H3a_mean_fdp=float(nul.H3a_fdp.mean()), H3b_mean_fdp=float(nul.H3b_fdp.mean())),
        H3b_null_with_confounded_lineage_effects_50=dict(H3b_mean_fdp=float(nul_conf.H3b_fdp.mean()),
                                                         H3a_mean_fdp=float(nul_conf.H3a_fdp.mean()),
                                                         H3a_power=float(nul_conf.H3a_power.mean())),
        planted_effect015_20=dict(H3a_mean_fdp=float(mix.H3a_fdp.mean()), H3b_mean_fdp=float(mix.H3b_fdp.mean()),
                                  H3a_power=float(mix.H3a_power.mean()), H3b_power=float(mix.H3b_power.mean())),
        runtime_s=round(time.time() - t, 1))
    txt = json.dumps(res, indent=1, default=str)
    print(txt)
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "synthetic_pipeline.json").write_text(txt + "\n")
