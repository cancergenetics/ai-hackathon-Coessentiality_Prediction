"""Synthetic world for the lockbox driver (scripts/48_lockbox_confirm.py); helper for tests/test_part2_lockbox_driver.py.

part2.synth.simulate plus:
  * planted lockbox pattern: the synthetic effects (MSI_H -> WRN complex, hotspot:KRAS -> KRAS complex, H1 ARID1A ->
    ARID1B compensation, ...) are present in A and L (should confirm); a planted A-only effect (hotspot:BRAF on C0032,
    A models only) is discovered in A and must fail in L;
  * two rare contexts with < 5 L models: LoF:KMT2D (C0030) and fusion:EWSR1::FLI1|ERG (C0031), effects in A and L ->
    discovered in A, not lockbox-testable, rare-context route (LOSO, DEMETER2 on A, leave-one-model-out);
  * a synthetic DEMETER2 (0.8 x CRISPR + noise) served by the outcome loader for source 'demeter2';
  * a GLS comparator map (members partly reshuffled), cluster tiers, coverage tables;
  * the registry entries of the real layout (controls, H3a, H3b, H1, H2, BENCH, BRIDGE-L) with seeds incl. the A3
    lockbox_permutation seeds and the rulings the driver requires, and a frozen lockbox-testable list from the fold.
The A-side discovery runs use the real runner (scripts/44_run_direction.py: run_h3 / run_h1 / run_h2) and
controls.evaluate_controls, committed as on real data. No DepMap data is read.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import part2_synth_pipeline as P  # noqa: E402
from part2 import assoc, controls, registry, synth  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402

RARE = {"LoF:KMT2D": ("C0030", 20, 3, -0.6), "fusion:EWSR1::FLI1|ERG": ("C0031", 18, 4, -0.7)}
FAIL = ("hotspot:BRAF", "C0032", -0.5)          # A-only effect: discovered in A, must not confirm in L
CONFIRM = [("MSI_H", "WRN"), ("hotspot:KRAS", "KRASC")]      # planted in A and L
H3B = synth.H3B_CONTEXTS + list(RARE)
H3B_COV = dict(synth.H3B_COV, **{c: synth.MUT_COV for c in RARE})
SEEDS = {"H3a_lineage_atlas": (4000, 6100), "H3b_driver_family": (4100, 6200), "H1_within_complex_GI": (4200, 6300),
         "H2_submodule_divergence": (4300, 6500)}
E2_LABEL = "DETECTABLE, WITHIN EQUIVALENCE MARGIN (negligible)"
CTRL_INDEX = {"CTRL_PC1": 0, "CTRL_PC6": 5, "CTRL_PC20": 20, "CTRL_PC21": 21}
RULINGS_BENCH = ("estimand_secondary_not_run:\n  - \"A-L correlation vs matched random gene sets\"\n"
                 "  - \"number of distinct complexes among the g*\"\n")       # mirrors A4 (G1)
RULING_H1 = "lockbox:\n  rare_route_fitted_parameters: not evaluable\n"       # mirrors A4 (G3) for H1S


def _script(name):
    mod = "script_" + name.split(".")[0]
    if mod in sys.modules:
        return sys.modules[mod]
    spec = importlib.util.spec_from_file_location(mod, ROOT / "scripts" / name)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod] = m                     # dataclasses need the module registered
    spec.loader.exec_module(m)
    return m


def plant(sim: synth.Synth, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed + 77)
    fold = sim.fold
    A_ids = np.array(sorted(fold.index[fold == "A"]))
    L_ids = np.array(sorted(fold.index[fold == "L"]))
    truth = {}
    for c, (u, nA, nL, d) in RARE.items():
        pick = list(rng.choice(A_ids, nA, replace=False)) + list(rng.choice(L_ids, nL, replace=False))
        col = pd.Series(0.0, index=sim.contexts.index)
        col[pick] = 1.0
        sim.contexts[c] = col
        sim.E.loc[pick, sim.fmap.members[u]] += d
        truth[c] = u
    c, u, d = FAIL
    rows = (sim.contexts[c] == 1) & (fold == "A")
    sim.E.loc[rows.values, sim.fmap.members[u]] += d
    tiers = ["T1", "T2", "T3-neg"]
    sim.fmap.tier = {k: tiers[i % 3] for i, k in enumerate(sorted(sim.fmap.members))}
    return truth


def gls_map(sim: synth.Synth, seed: int = 0) -> FrozenMap:
    rng = np.random.default_rng(seed + 88)
    co = sim.fmap.coarse()
    mem = {f"G{i + 1:04d}": list(sim.fmap.members[u]) for i, u in enumerate(co)}
    keys = list(mem)
    for k in keys:                                    # move ~20% of members to another cluster
        m = mem[k]
        mv = [g for g in m if rng.random() < 0.2]
        if len(m) - len(mv) < 3:
            continue
        for g in mv:
            m.remove(g)
            mem[keys[rng.integers(len(keys))]].append(g)
    mem = {k: v for k, v in mem.items() if 3 <= len(v) <= 60}
    tiers = ["T1", "T2", "T3-neg"]
    return FrozenMap(members=mem, level={k: 1 for k in mem}, parent={k: "" for k in mem},
                     locus={k: False for k in mem}, tier={k: tiers[i % 3] for i, k in enumerate(sorted(mem))})


def _seeds(perm: int | None, lockbox: int | None, extra: dict | None = None) -> str:
    s = ["seeds:"]
    if perm is not None:
        s.append(f"  permutation: {perm}")
    if lockbox is not None:
        s.append(f"  lockbox_permutation: {lockbox}")
    for k, v in (extra or {}).items():
        s.append(f"  {k}: {v}")
    return "\n".join(s) + "\n"


def testable_list(sim, fold, entries: dict, events) -> pd.DataFrame:
    rows = []
    L = fold == "L"
    tab = pd.concat([sim.contexts, sim.lof.add_prefix("H1loss:")], axis=1)
    for k, e in entries.items():
        if e["direction"] == "H1":
            for ev in events:
                col = tab[ev.context]
                rows.append(dict(entry=k, direction="H1", context=ev.context, event=ev.id, role="primary",
                                 n_L=int((col[L] == 1).sum()), ref_L=int((col[L] == 0).sum())))
            continue
        for c in e["context_family"]["contexts"]:
            col = tab[c["id"]]
            rows.append(dict(entry=k, direction=e["direction"], context=c["id"], event="", role="primary",
                             n_L=int((col[L] == 1).sum()), ref_L=int((col[L] == 0).sum())))
    t = pd.DataFrame(rows)
    t["lockbox_testable"] = (t.n_L >= 5) & (t.ref_L >= 5)
    t["route"] = np.where(t.lockbox_testable, "lockbox", "rare-context route")
    return t


def build_world(root: Path, seed: int = 0, n_perm: int = P.N_PERM, K: int = 200, B_bench: int = 30,
                B_bridge: int = 200):
    """Synthetic project with committed A-side runs; returns (World, info)."""
    D = _script("48_lockbox_confirm.py")
    R44 = _script("44_run_direction.py")
    sim = synth.simulate(seed)
    truth = plant(sim, seed)
    loss, kind, events = P.h1_inputs(sim)
    R = P.build_repo(root, sim, events)
    proj = R["proj"]
    fin = ["results/part2/contexts/checksums.json", "results/part2/lockbox/manifest.json",
           "results/part2/map/manifest.json"]
    reg = proj.registry_dir
    # ---- GLS map (frozen with the zavg map)
    gls = gls_map(sim, seed)
    mdir = root / "results" / "part2" / "map"
    gls.to_membership().to_csv(mdir / "gls_membership.tsv", sep="\t", index=False)
    man = json.loads((mdir / "manifest.json").read_text())
    man["files"]["gls_membership.tsv"] = assoc._sha_file(mdir / "gls_membership.tsv")
    (mdir / "manifest.json").write_text(json.dumps(man))
    # ---- registry: H3b with the rare contexts, H2, BENCH, BRIDGE-L, seeds (A3 lockbox seeds) and rulings
    coarse = sim.fmap.coarse()
    h3a_ctx = [dict(id=f"lineage:{l}", type="lineage") for l in synth.LINEAGES]
    h3b_ctx = [dict(id=c, covariates=H3B_COV[c]) for c in H3B]
    bb = dict(method="benjamini_bogomolov", families="contexts", q=0.05)
    (reg / "H3a_lineage_atlas.yaml").write_text(synth.entry_yaml(
        "H3a_lineage_atlas", "H3a", h3a_ctx, coarse, "S_C_mean_winsor3", synth.H3A_COV, bb, fin,
        extra=_seeds(*SEEDS["H3a_lineage_atlas"], {"pseudo_contexts": 4600})))
    (reg / "H3b_driver_family.yaml").write_text(synth.entry_yaml(
        "H3b_driver_family", "H3b", h3b_ctx, coarse, "S_C_mean_winsor3", synth.BASE_COV, bb, fin,
        extra=_seeds(*SEEDS["H3b_driver_family"], {"pseudo_contexts": 4700})))
    from part2 import h2 as h2m, h1 as h1m
    subs = h2m.h2_units(sim.fmap)
    h2_ctx = [dict(id=f"lineage:{l}", type="lineage", covariates=synth.H3A_COV + ["NNMD", "CasActivity"])
              for l in synth.LINEAGES]
    h2_ctx += [dict(id=c, covariates=H3B_COV[c] + ["NNMD", "CasActivity"]) for c in H3B]
    (reg / "H2_submodule_divergence.yaml").write_text(synth.entry_yaml(
        "H2_submodule_divergence", "H2", h2_ctx, subs, "submodule_vs_sibling", synth.BASE_COV + ["NNMD", "CasActivity"],
        dict(method="stagewise_bb", parent_directions=["H3a", "H3b"], q=0.05), fin,
        extra=_seeds(*SEEDS["H2_submodule_divergence"], {"random_splits": 4800})))
    h1_ctx = [dict(id=e.context, covariates=e.covariates) for e in {e.context: e for e in events}.values()]
    (reg / "H1_within_complex_GI.yaml").write_text(synth.entry_yaml(
        "H1_within_complex_GI", "H1", h1_ctx, h1m.units_of(events), "member_resolved", synth.BASE_COV,
        dict(method="bh", q=0.05), fin,
        extra=_seeds(*SEEDS["H1_within_complex_GI"], {"matched_null": 7000}) + RULING_H1))
    for k, i in CTRL_INDEX.items():
        p = reg / f"{k}.yaml"
        p.write_text(p.read_text() + _seeds(4400 + i, 6700 + i))
    L = sim.fold == "L"
    bench_ctx = [c for c in h3a_ctx + h3b_ctx
                 if ((sim.contexts[c["id"]][L] == 1).sum() >= 5) and ((sim.contexts[c["id"]][L] == 0).sum() >= 5)]
    cov_of = {c["id"]: c.get("covariates", synth.H3A_COV) for c in h3a_ctx + h3b_ctx}
    bctx = [dict(c, covariates=cov_of[c["id"]]) for c in bench_ctx]
    (reg / "BENCH_complex_vs_gene.yaml").write_text(synth.entry_yaml(
        "BENCH_complex_vs_gene", "benchmark", bctx, coarse, "S_C_mean_winsor3", synth.BASE_COV,
        dict(method="none", q=0.05), fin,
        extra=_seeds(None, None, {"bootstrap": 5100, "matched_sets": 7100}) + RULINGS_BENCH +
        "multiverse:\n  ci: point estimate only\n"))
    gctx = [dict(c, covariates=cov_of[c["id"]]) for c in h3a_ctx + h3b_ctx]
    (reg / "BRIDGE-L_confirmation_rate.yaml").write_text(synth.entry_yaml(
        "BRIDGE-L_confirmation_rate", "bridge", gctx, gls.coarse(), "S_C_mean_winsor3", synth.BASE_COV, bb, fin,
        extra=_seeds(4500, 6600, {"bootstrap": 5200}) + "multiverse:\n  ci: point estimate only\n"))
    entries = registry.load_all(reg)
    tl = testable_list(sim, sim.fold, entries, events)
    fz = root / "results" / "part2" / "fzx"
    fz.mkdir(parents=True, exist_ok=True)
    tl.to_csv(fz / "lockbox_testable.tsv", sep="\t", index=False)
    P.commit_all(root, "synthetic registry (lockbox driver test): rare contexts, BENCH, BRIDGE-L, A3 seeds")

    base_loader = proj.outcome_loader
    noise = np.random.default_rng(seed + 99).normal(0, 0.1, sim.E.shape)
    D2 = pd.DataFrame(0.8 * sim.E.values + noise, index=sim.E.index, columns=sim.E.columns)

    def loader(model_ids, genes, source):
        if source == "demeter2":
            return D2.reindex(index=list(model_ids), columns=list(genes) if genes is not None else D2.columns)
        return base_loader(model_ids, genes, source)

    proj = assoc.Project(root=root, outcome_loader=loader)
    # ---- A-side runs (the real runner), committed
    tab = pd.read_parquet(root / P.REL_CTX)
    E_A = assoc.load_A_effects(proj=proj)
    A = E_A.index
    tabA = tab.loc[A]
    ctab, kill = controls.evaluate_controls(E_A, sim.fmap, tabA, tabA, sim.annot, inputs=R["inputs"], proj=proj)
    P.synth_commit(root, "controls (A)")
    env = R44.Env(proj=proj, fmap=sim.fmap, ctx=tabA, E_A=E_A, cov=tabA, annot=sim.annot, inputs=R["inputs"],
                  blocks=tabA.PatientID, strata=tabA.OncotreeLineage, homdel={}, strict=False, write=True)
    for aid in ("H3a_lineage_atlas", "H3b_driver_family"):
        R44.run_h3(env, aid, n_perm=n_perm)
    P.synth_commit(root, "H3a / H3b discovery (A)")
    lossA = loss.loc[A]
    genes = sorted(sim.fmap.all_genes(1))
    ps = 10 * E_A.clip(-3, 3)[genes].corr()
    R44.run_h1(env, "H1_within_complex_GI", events, lossA, ps, n_perm=n_perm, K=K)
    cp = root / "results" / "part2" / "directions" / "H1_within_complex_GI" / "classes.tsv"
    cl = pd.read_csv(cp, sep="\t")
    cl["label_final"], cl["label_final_paralog"] = cl.class_complex_wide, cl.class_paralog
    cl["headline_eligible"] = cl.class_complex_wide.isin(["COLLAPSE", "RESIDUAL SENSITISATION"]) | \
        (cl.class_paralog == "PARALOG COMPENSATION")
    sens = cl.class_complex_wide == "RESIDUAL SENSITISATION"       # an E2-style relabel: confirmation cannot promote
    cl.loc[sens, "label_final"] = E2_LABEL
    cl.loc[sens, "headline_eligible"] = False
    cl.to_csv(cp, sep="\t", index=False)
    P.synth_commit(root, "H1 discovery (A) + reporting labels")
    R44.run_h2(env, n_perm=n_perm)
    P.synth_commit(root, "H2 discovery (A)")
    cov_t = {}
    for k, fm in (("zavg", sim.fmap), ("gls", gls)):
        rr = np.random.default_rng(seed + (1 if k == "gls" else 2))
        co = fm.coarse()
        cov_t[k] = pd.DataFrame(dict(cluster_id=co, size=[len(fm.members[u]) for u in co],
                                     stability=rr.random(len(co)), cohesion=rr.random(len(co))))
    world = D.World(proj=proj, fmap=sim.fmap, gls=gls, ctx=tab, loss=tab[[c for c in tab.columns
                                                                                    if c.startswith("H1loss:")]],
                    cov=tab, E_A=E_A, annot=sim.annot, homdel={}, ctrl_homdel={},
                    events={"H1_within_complex_GI": events},
                    inputs_gls=[P.REL_CTX, "results/part2/map/gls_membership.tsv"], pair_scores=ps,
                    coverage_tables=cov_t, n_perm=n_perm, B_bench=B_bench, B_bridge=B_bridge)
    info = dict(sim=sim, truth=truth, events=events, kill=kill, D=D, named=sim.named)
    return world, info
