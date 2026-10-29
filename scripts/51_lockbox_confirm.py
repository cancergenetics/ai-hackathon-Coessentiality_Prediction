#!/usr/bin/env python
"""Lockbox confirmations: the single read (DESIGN v3.1 §0 F4/F5, §2.5-2.6; registry `lockbox` blocks; amendment A3).

Thin driver over the existing machinery; it adds no analysis rule of its own. Order (registry, README "lockbox read
order"): the 23 controls by registry index (seeds.permutation 4400 + index), H3a_lineage_atlas, H3b_driver_family,
H1_within_complex_GI, H1S_silencing_GI, H2_submodule_divergence; then BENCH_complex_vs_gene and
BRIDGE-L_confirmation_rate (refused unless every confirmatory entry has a committed ledger line and is closed).

Per confirmatory entry (confirm.py: write_target_list -> confirm_entry -> close_entry):
  1. A side first (no L value): the hits (controls: every unit with an A estimate; H3a/H3b/H2: the committed
     discoveries; H1/H1S: every evaluable event axis, eligible = the non-neutral class call on its defining statistic,
     beta_P for PARALOG COMPENSATION and kappa_N for COLLAPSE / RESIDUAL SENSITISATION), sorted by (context, unit);
     the route of each hit from the frozen lockbox-testable list (results/part2/fzx/lockbox_testable.tsv; H1/H1S per
     event); for eligible hits that are not lockbox-testable, confirm.rare_route in A (leave-one-lineage-out same
     sign, DEMETER2 on A models same sign, leave-one-model-out same sign; labelled 'not lockbox-testable').
  2. The committed target list (members of every hit unit; H1/H1S: members of every event; H2: the sub-module's
     parent), then confirm.confirm_entry: assoc.read_lockbox (ledger line committed first, one read per entry),
     same score / removals / A-fitted parameters / covariates in L, one-sided test in the A direction, HC3 with
     Satterthwaite df, exact PatientID-blocked Freedman-Lane permutation (10,000; within lineage for non-lineage
     contexts) when < 20 context models are in L; confirmed iff same sign, p < 0.05, >= 5 context and >= 5 reference
     L models; seed = seeds.lockbox_permutation (A3) + hit index. L 95% CIs for every hit, 90% CIs added (H1 neutral
     calls). Only eligible hits on the lockbox route can be confirmed; E2 labels (`label_final`, `headline_eligible`)
     are carried unchanged, so a confirmation never promotes a DETECTABLE-WITHIN-MARGIN event.
  3. Outputs results/part2/lockbox_confirm/<ID>/ (confirmation.tsv, summary.json) and the engine's
     results/part2/assoc/<ID>/confirmation_crispr.json, committed as "Lockbox confirmations (single read): <ID>";
     then confirm.close_entry (status: closed, its own commit).
BENCH: primary pairs = H3a + H3b discoveries in BENCH's (lockbox-testable) contexts; g* = bench.choose_gstar on the
  members after the context's removals, breadth = robust.breadth over C's members (A only); one read of all members of
  the complexes in the primary and secondary pair sets; bench.bench_delta (B 2000, seed 5100); Spearman(Delta,
  breadth); bench.distinct_gstar; multiverse eigengene (A-fitted), g* by largest beta in the A direction, T1/T2 pairs
  (within_line_rank NOT RUN, A3(b)); secondary bench.al_correlation against calib.matched_random_sets (seed 7100).
BRIDGE-L: the registered GLS-map discovery on A (two run_association calls, 23 H3a and 20 H3b contexts, each its own
  BB family at q 0.05, Freedman-Lane re-entry, seed 4500, H3a/H3b removals and covariates), committed; then the
  confirmation of the GLS discoveries exactly as H3a/H3b (one read); bench.bridge_rate (B 2000, seed 5200) on the
  lockbox-testable zavg (H3a/H3b confirmation tables, no new zavg read) and GLS hits; multiverse: unweighted rate,
  matched gene coverage subset (validation.match_coverage, FZ-C decision D5), T1/T2 clusters.
Report: rates per direction (overall, non-locus, beats-all-matched-sets vs not), the H1 per-event table, BENCH and
  BRIDGE-L estimates (results/part2/lockbox_confirm/REPORT.json, rates.tsv, h1_events.tsv).

Registry rulings the entries do not contain are refused up front (preflight) -- nothing is read while a gap is open.
Amendments A4/A5: BENCH secondary A-L random-set correlation NOT RUN, distinct-complex count among g* reported (A
only); BENCH multiverse point estimates, BRIDGE-L multiverse bootstrap CIs; H1S rare-route N axes use the A CRISPR
alignment (`lockbox.rare_route_fitted_parameters`); a gene missing from DEMETER2 makes the rare route NOT EVALUABLE.
Review MUST-1/2 (results/part2/verify/lockbox_code_review.md): refused unless the tracked tree is clean and every input
that decides hits, eligibility or route is committed (lockbox_testable.tsv also against its FZ-X manifest hash); the
target list records each hit's eligibility, route, E2 label and headline flag, and confirm_entry checks the hits
against it; every read matrix is persisted (data_cache/part2/lockbox_reads/, hash committed at once) and `--resume`
finishes an entry from it without calling the loader; every A-side computation and validity check (g*, eigengene
fits, GLS frozen inputs, coverage tables) runs before the corresponding target list is committed. Run the read under
`caffeinate -i` (macOS).

Usage (PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/51_lockbox_confirm.py --dry-check     # A side only: hits, routes, gaps; no L read
    caffeinate -i .venv/bin/python scripts/51_lockbox_confirm.py --run            # the single read
    caffeinate -i .venv/bin/python scripts/51_lockbox_confirm.py --run --resume   # after a crash: persisted reads only
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, bench, calib, confirm, controls, covariates as covmod, cscore, h1, h2, registry, robust  # noqa
from part2 import validation  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402
from part2.registry import RegistryError  # noqa: E402

DIRECTIONS = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI",
              "H2_submodule_divergence")
H3_IDS, H1_IDS, H2_ID = DIRECTIONS[:2], DIRECTIONS[2:4], DIRECTIONS[4]
BENCH_ID, BRIDGE_ID = "BENCH_complex_vs_gene", "BRIDGE-L_confirmation_rate"
N_PERM = 10_000           # registry: exact Freedman-Lane permutation, 10,000 (controls) / the discovery machinery
MIN_A = 12                # BRIDGE-L GLS discovery: H3a/H3b evaluability (>= 12 context and >= 12 reference A models)
SCREEN_Q = 0.1            # BRIDGE-L GLS discovery: screening hits re-tested by permutation (H3a/H3b machinery)
B_BENCH, B_BRIDGE = 2000, 2000   # BENCH / BRIDGE-L estimands: B = 2000
NONNEUTRAL = {"N": ("COLLAPSE", "RESIDUAL SENSITISATION"), "P": ("PARALOG COMPENSATION",)}
OUT_REL = "results/part2/lockbox_confirm"
TESTABLE_REL = "results/part2/fzx/lockbox_testable.tsv"
GLS_MAP = "results/part2/fzx/gls_engine_membership.tsv"
CLUSTERS = {"zavg": "results/part2/construction/zavg_map/clusters.tsv",
            "gls": "results/part2/construction/gls_map/clusters.tsv"}
RARE_LABEL = "not lockbox-testable"
FIT_RULINGS = ("A CRISPR parameters", "DEMETER2 A refit")      # A5 (G3): H1S uses the A CRISPR alignment
D2_NA = "NOT EVALUABLE (DEMETER2 not available)"                # review N5
FZX_MANIFEST = "results/part2/fzx/FZX_manifest.json"


class GapError(RegistryError):
    """A registered step needs a registry ruling that the entries do not contain (nothing is read)."""


def _script(name: str):
    mod = "script_" + name.split(".")[0]
    if mod in sys.modules:
        return sys.modules[mod]
    spec = importlib.util.spec_from_file_location(mod, ROOT / "scripts" / name)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod] = m                     # dataclasses need the module registered
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------------------------------------ world
@dataclass
class World:
    proj: assoc.Project
    fmap: FrozenMap               # zavg_BEF frozen map
    gls: FrozenMap                # GLS comparator map (BRIDGE-L)
    ctx: pd.DataFrame             # all models x context columns (labels only)
    loss: pd.DataFrame            # all models x 'H1loss:<GENE>' (labels only)
    cov: pd.DataFrame             # all models covariate source table (labels, screen metadata)
    E_A: pd.DataFrame             # A models x gene keys (CRISPR)
    annot: cscore.Annotation
    homdel: dict                  # H3/H2/BENCH/BRIDGE context -> hom-del share (removal rule)
    ctrl_homdel: dict             # control context -> hom-del share (controls script rule)
    events: dict                  # H1 entry -> [h1.Event]
    inputs_gls: list              # frozen inputs of the GLS discovery run_association calls
    pair_scores: object = None    # callable -> DataFrame (coherence of matched random sets)
    coverage_tables: dict | None = None   # {'zavg', 'gls'}: level-1 cluster tables (size, stability, cohesion)
    calib_tests: dict = field(default_factory=dict)   # entry -> D8 calibration test table
    n_perm: int = N_PERM
    B_bench: int = B_BENCH
    B_bridge: int = B_BRIDGE
    testable_sha: str | None = None       # FZ-X manifest hash of the frozen lockbox-testable list (real run)
    expect_gls: dict | None = None        # review N13: GLS discovery counts of the rehearsal (determinism check)
    persist: bool = True                  # persist every read matrix (PI crash-recovery decision)
    resume: bool = False                  # finish entries already read, from their persisted matrices only
    log: list = field(default_factory=list)
    _d2: pd.DataFrame | None = None

    @property
    def A(self):
        return self.E_A.index

    def D2_A(self) -> pd.DataFrame:
        """DEMETER2 on A models (rare-context route), through assoc.load_A_effects (A list only)."""
        if self._d2 is None:
            self._d2 = assoc.load_A_effects(source="demeter2", proj=self.proj)
        return self._d2

    def testable(self) -> pd.DataFrame:
        p = self.proj.root / TESTABLE_REL
        registry.committed_blob(self.proj.root, p)                    # MUST-1: committed and clean
        if self.testable_sha is not None and assoc._sha_file(p) != self.testable_sha:
            raise RegistryError(f"{TESTABLE_REL} differs from its FZ-X manifest hash")
        t = pd.read_csv(p, sep="\t", keep_default_na=False)
        t["lockbox_testable"] = t.lockbox_testable.astype(str).str.lower() == "true"
        return t

    def out(self, aid: str | None = None) -> Path:
        p = self.proj.root / OUT_REL
        return p / aid if aid else p

    def say(self, msg: str):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        self.log.append(line)
        print(line, flush=True)


def entry(world: World, aid: str):
    return registry.load_entry(world.proj.registry_dir, aid, root=world.proj.root)


def committed(world: World, path) -> Path:
    """MUST-1: an input that decides hits, eligibility or route must be committed and clean at HEAD."""
    p = Path(path)
    p = p if p.is_absolute() else world.proj.root / p
    registry.committed_blob(world.proj.root, p)
    return p


def require_clean(world: World) -> None:
    r = registry.git(world.proj.root, "status", "--porcelain", "--untracked-files=no", check=False)
    if r.returncode != 0 or r.stdout.strip():
        raise RegistryError(f"tracked files have uncommitted changes (refused): {r.stdout.strip()[:300]}")


def check_A(world: World) -> None:
    """Review S2: the A matrix the driver uses holds A models only."""
    bad = set(world.E_A.index) - set(assoc.a_models(world.proj))
    if bad:
        raise assoc.FirewallError(f"World.E_A holds {len(bad)} non-A models (e.g. {sorted(bad)[:2]})")


def _not_run(items, key: str) -> bool:
    return any(key.lower() in str(s).lower() for s in (items or []))


def primary_contexts(e: dict) -> list:
    return [c["id"] for c in e["context_family"]["contexts"] if c.get("role", "primary") == "primary"]


def lockbox_seed(e: dict) -> int:
    s = (e.get("seeds") or {}).get("lockbox_permutation")
    if s is None:
        raise GapError(f"{e['id']}: seeds.lockbox_permutation is not registered (amendment A3)")
    return int(s)


def confirm_order(proj) -> list:
    allent = registry.load_all(proj.registry_dir)
    ctrls = [k for k, e in allent.items() if e["direction"] == "controls"]
    ctrls = sorted(ctrls, key=lambda k: int(allent[k]["seeds"]["permutation"]))
    return ctrls + [k for k in DIRECTIONS if k in allent]


def committed_table(proj, aid: str, tag: str = "discovery") -> pd.DataFrame:
    """The engine's full table (data_cache, gitignored), checked against its committed summary hash."""
    p = proj.cache_dir / "assoc" / aid / f"{tag}.parquet"
    s = proj.results_dir / "assoc" / aid / f"{tag}_summary.json"
    registry.committed_blob(proj.root, s)
    if assoc._sha_file(p) != json.loads(s.read_text())["full_table_sha256"]:
        raise RuntimeError(f"{aid} {tag} table differs from its committed summary")
    return pd.read_parquet(p)


def rulings(world: World) -> dict:
    """Registry rulings for steps the original entries leave open (read from the committed entries)."""
    out = {}
    eb, _ = entry(world, BENCH_ID)
    out["bench_secondary"] = eb.get("estimand_secondary")
    out["bench_secondary_not_run"] = eb.get("estimand_secondary_not_run")      # A4 (G1)
    out["bench_multiverse_ci"] = (eb.get("multiverse") or {}).get("ci")
    eg, _ = entry(world, BRIDGE_ID)
    out["bridge_multiverse_ci"] = (eg.get("multiverse") or {}).get("ci")
    out["fitted"] = {}
    for k in confirm_order(world.proj) + [BRIDGE_ID]:
        e, _ = entry(world, k)
        out["fitted"][k] = (e.get("lockbox") or {}).get("rare_route_fitted_parameters")
    return out


def _fit_ruling(world: World, aid: str) -> str:
    r = rulings(world)["fitted"].get(aid)
    if r not in FIT_RULINGS:
        raise GapError(f"{aid}: a rare-route hit needs fitted score parameters on DEMETER2 (registry silent: "
                       f"lockbox.rare_route_fitted_parameters, one of {FIT_RULINGS})")
    return r


# ------------------------------------------------------------------------------------------------ hit specifications
@dataclass
class Spec:
    aid: str
    hits: pd.DataFrame            # context, unit, beta, alternative, eligible, lockbox_testable, (+ meta)
    members_of: object            # unit -> gene keys (target list)
    score_A: object               # (context, unit) -> Series over A
    score_L: object               # (E_L, context, unit) -> Series over L
    score_D2: object              # (context, unit) -> Series over DEMETER2 A rows, or None (gene not screened)
    ctx_all: pd.DataFrame         # all models x the hits' contexts
    cov_A: dict
    cov_all: object               # models -> {context: covariates over those L rows} (R6)
    lineage_contexts: set
    extra_genes: list = field(default_factory=list)
    needs_fit: list = field(default_factory=list)   # rare-route hits whose DEMETER2 score needs fitted parameters


def _route(world: World, aid: str, hits: pd.DataFrame, by_event: bool = False) -> pd.Series:
    t = world.testable()
    t = t[t.entry == aid]
    key = (lambda r: (r.context, r.event)) if by_event else (lambda r: (r.context, ""))
    look = {(r.context, r.event if by_event else ""): bool(r.lockbox_testable) for r in t.itertuples(index=False)}
    out = []
    for r in hits.itertuples(index=False):
        k = key(r)
        if k not in look:
            raise RegistryError(f"{aid}: {k} is not in the frozen lockbox-testable list")
        out.append(look[k])
    return pd.Series(out, index=hits.index, dtype=bool)


def l_covariates(world: World, e: dict, ctxs, models) -> dict:
    """R6 fix: the L covariate design is assembled over the L rows actually read, exactly as the A side assembled it
    over A rows (covariates.per_context on the source table of those models; assoc._design then applies the engine's
    conventions, incl. merging FE levels with < 2 models into the reference). The former all-model assembly subset to
    L could leave a single-model lineage as its own group (leverage 1; results/part2/lockbox_confirm/README.md, R6)."""
    return covmod.per_context(e, ctxs, world.cov.reindex(list(models)))[0]


def _cov(world: World, e: dict, ctxs) -> tuple[dict, object]:
    ca, _ = covmod.per_context(e, ctxs, world.cov.reindex(world.A))
    return ca, (lambda models: l_covariates(world, e, ctxs, models))


def _lineage(e: dict, ctxs) -> set:
    items = {c["id"]: c for c in e["context_family"]["contexts"]}
    return {c for c in ctxs if cscore.is_lineage_context(c, items.get(c))}


def _sorted(h: pd.DataFrame) -> pd.DataFrame:
    return h.sort_values(["context", "unit"], kind="mergesort").reset_index(drop=True)


def spec_controls(world: World, aid: str) -> Spec:
    e, _ = entry(world, aid)
    tab = pd.read_csv(committed(world, world.proj.results_dir / "controls" / "controls_summary.tsv"), sep="\t")
    rows = tab[(tab.entry == aid) & tab.beta.notna()]
    hits = _sorted(pd.DataFrame(dict(context=rows.context.values, unit=rows.unit.values,
                                     beta=rows.beta.values.astype(float), alternative=rows.alternative.values,
                                     level=rows.level.values)))
    alt_ok = [(b < 0) == (a == "less") for b, a in zip(hits.beta, hits.alternative)]
    if not all(alt_ok):
        raise RuntimeError(f"{aid}: an A estimate is opposite to the registered one-sided direction")
    hits["eligible"] = True
    ctl = e["control"]
    cu = ctl.get("complex_unit")
    fm, an = world.fmap, world.annot
    ctxs = sorted(set(hits.context))
    hs = world.ctrl_homdel.get(ctxs[0]) if ctxs else None
    has_c = isinstance(cu, str) and cu in fm.members
    cinfo = controls.complex_readout(e, world.E_A, fm, an, hs)[1] if has_c else {}

    def sA(c, u):
        return controls.complex_readout(e, world.E_A, fm, an, hs)[0] if u == cu else cscore.gene_unit(world.E_A, u)

    def sL(E_L, c, u):
        return controls.complex_readout(e, E_L, fm, an, hs, fit=cinfo)[0] if u == cu else cscore.gene_unit(E_L, u)

    def sD2(c, u):
        D = world.D2_A()
        if u != cu:
            try:
                return cscore.gene_unit(D, u)
            except KeyError:
                return None
        fit = None
        if e["score"] == "submodule_vs_sibling":
            fit = cinfo if _fit_ruling(world, aid) == FIT_RULINGS[0] else None
        return controls.complex_readout(e, D, fm, an, hs, fit=fit)[0]

    def mem(u):
        if u == cu:
            return list(fm.members[u]) + (list(fm.members[fm.parent[u]]) if fm.level.get(u) == 2 else [])
        return [x.strip() for x in u.split(" minus ")]

    hits["lockbox_testable"] = _route(world, aid, hits)
    ca, call = _cov(world, e, ctxs)
    needs = [u for u, t in zip(hits.unit, hits.lockbox_testable) if not t and u == cu
             and e["score"] == "submodule_vs_sibling"]
    return Spec(aid, hits, mem, sA, sL, sD2, world.ctx.reindex(columns=ctxs), ca, call, _lineage(e, ctxs),
                needs_fit=needs)


def _h3_scores(world: World, fm: FrozenMap):
    cache = {}

    def removal(c):
        return cscore.removal_for_context(c, world.homdel.get(c))

    def frame(E, c, units, key):
        rem = removal(c)
        k = (key, rem)
        if k not in cache:
            cache[k] = cscore.unit_scores(E, fm, units, rem, world.annot)[0]
        return cache[k]

    return removal, frame


def spec_h3(world: World, aid: str, fm: FrozenMap | None = None, disc: pd.DataFrame | None = None) -> Spec:
    """H3a / H3b (and the GLS hits of BRIDGE-L, fm = GLS map): the committed discoveries."""
    e, _ = entry(world, aid)
    fm = fm or world.fmap
    disc = committed_table(world.proj, aid) if disc is None else disc
    d = disc[disc.discovery.astype(bool)]
    hits = _sorted(d[["context", "unit", "beta", "alternative", "t"]].copy())
    hits["eligible"] = True
    units = fm.coarse()
    removal, frame = _h3_scores(world, fm)

    def sA(c, u):
        return frame(world.E_A, c, units, "A")[u]

    def sL(E_L, c, u):
        return cscore.unit_scores(E_L, fm, [u], removal(c), world.annot)[0][u]

    def sD2(c, u):
        return cscore.unit_scores(world.D2_A(), fm, [u], removal(c), world.annot)[0][u]

    ctxs = sorted(set(hits.context))
    hits["lockbox_testable"] = _route(world, aid, hits)
    ca, call = _cov(world, e, ctxs)
    return Spec(aid, hits, lambda u: list(fm.members[u]), sA, sL, sD2, world.ctx.reindex(columns=ctxs), ca, call,
                _lineage(e, ctxs))


def spec_h2(world: World, aid: str = H2_ID) -> Spec:
    e, _ = entry(world, aid)
    fm = world.fmap
    disc = committed_table(world.proj, aid)
    hits = _sorted(disc[disc.discovery.astype(bool)][["context", "unit", "beta", "alternative", "t"]].copy())
    hits["eligible"] = True
    ctxs = sorted(set(hits.context))
    subs = h2.h2_units(fm)

    def removal(c):
        return cscore.removal_for_context(c, world.homdel.get(c))

    Y2, coefs = h2.h2_frames(world.E_A, fm, ctxs, removal, world.annot, units=subs)
    f_L = h2.score_L_fn(fm, coefs, removal, world.annot)

    def sD2(c, u):
        coef = coefs.get((c, u)) if _fit_ruling(world, aid) == FIT_RULINGS[0] else None
        return h2.residualised(world.D2_A(), fm, u, removal(c), world.annot, coef=coef)[0]

    hits["lockbox_testable"] = _route(world, aid, hits)
    ca, call = _cov(world, e, ctxs)
    needs = list(hits.unit[~hits.lockbox_testable])
    return Spec(aid, hits, lambda u: list(fm.members[fm.parent[u]]), lambda c, u: Y2[c][u], f_L, sD2,
                world.ctx.reindex(columns=ctxs), ca, call, _lineage(e, ctxs), needs_fit=needs)


def spec_h1(world: World, aid: str) -> Spec:
    e, _ = entry(world, aid)
    events = world.events[aid]
    lossA = world.loss.reindex(world.A)
    lossA_g = lossA.rename(columns=lambda c: c.split(":", 1)[1])
    comps = {ev.id: h1.component_scores(world.E_A, ev, lossA_g[ev.gene]) for ev in events}
    disc = committed_table(world.proj, aid)
    rmap = {r.unit: r for r in disc.itertuples(index=False)}
    cls = pd.read_csv(committed(world, world.proj.results_dir / "directions" / aid / "classes.tsv"), sep="\t",
                      keep_default_na=False)
    cmap = {r["event"]: r for r in cls.to_dict("records")}
    rows = []
    for ev in events:
        c = comps[ev.id]
        k = cmap[ev.id]
        for ax in ("N", "P"):
            if ax not in c.y:
                continue
            u = f"{ev.id}|{ax}"
            r = rmap.get(u)
            if r is None or not bool(r.evaluable) or not np.isfinite(r.beta):
                continue
            cl = k["class_complex_wide"] if ax == "N" else k["class_paralog"]
            lab = k.get("label_final" if ax == "N" else "label_final_paralog", cl)
            rows.append(dict(context=ev.context, unit=u, beta=float(r.beta), alternative=r.alternative, t=r.t,
                             event=ev.id, axis=ax, class_engine=cl, label_final=lab,
                             headline_eligible=str(k.get("headline_eligible", "")),
                             eligible=cl in NONNEUTRAL[ax]))
    hits = _sorted(pd.DataFrame(rows))
    ctxs = sorted(set(hits.context))
    ev_of = {ev.id: ev for ev in events}

    def sA(cx, u):
        eid, ax = u.rsplit("|", 1)
        return comps[eid].y[ax]

    def sD2(cx, u):
        eid, ax = u.rsplit("|", 1)
        ev, D = ev_of[eid], world.D2_A()
        if ax == "P":
            return cscore.mean_score(D, ev.P)
        if _fit_ruling(world, aid) == FIT_RULINGS[0]:
            return h1.aligned_score(D, comps[eid].Nprime, comps[eid].b)
        return h1.component_scores(D, ev, lossA_g[ev.gene].reindex(D.index)).y.get("N")

    hits["lockbox_testable"] = _route(world, aid, hits, by_event=True)
    ca, call = _cov(world, e, ctxs)
    needs = list(hits.unit[(~hits.lockbox_testable) & hits.eligible & (hits.axis == "N")])
    genes = sorted({g for ev in events for g in ev.P + ev.N})
    return Spec(aid, hits, lambda u: ev_of[u.rsplit("|", 1)[0]].P + ev_of[u.rsplit("|", 1)[0]].N, sA,
                h1.score_L_fn(events, comps, world.loss), sD2, world.loss.reindex(columns=ctxs), ca, call, set(),
                extra_genes=genes, needs_fit=needs)


def build_spec(world: World, aid: str) -> Spec:
    e, _ = entry(world, aid)
    if e["direction"] == "controls":
        return spec_controls(world, aid)
    if aid in H3_IDS:
        return spec_h3(world, aid)
    if aid in H1_IDS:
        return spec_h1(world, aid)
    if aid == H2_ID:
        return spec_h2(world, aid)
    raise RegistryError(f"no confirmation spec for {aid}")


# ------------------------------------------------------------------------------------------------ confirmation
def _ledger_ids(proj) -> set:
    return {r["analysis_id"] for r in assoc._ledger(proj)}


def rare_routes(world: World, spec: Spec) -> dict:
    """{hit index: confirm.rare_route result} for eligible hits that are not lockbox-testable (A only)."""
    out = {}
    groups = world.cov.OncotreeLineage.reindex(world.A)
    for i, r in spec.hits.iterrows():
        if r.lockbox_testable or not r.eligible:
            continue
        yA = spec.score_A(r.context, r.unit)
        yD = spec.score_D2(r.context, r.unit)
        res = confirm.rare_route(yA, spec.ctx_all[r.context].reindex(world.A), spec.cov_A.get(r.context), groups,
                                 float(r.beta), yD)
        res.pop("loso_full", None)
        res["rare_route_note"] = ""
        if res["demeter2_same_sign"] is None:                  # N5: DEMETER2 missing -> not evaluable, not 'fail'
            res["rare_route_pass"] = None
            res["rare_route_note"] = D2_NA
        out[i] = res
    return out


def _with_ci90(t: pd.DataFrame) -> pd.DataFrame:
    h = sst.t.ppf(0.95, t.df_L.astype(float)) * t.se_L.astype(float)
    t["ci90_L_lo"], t["ci90_L_hi"] = t.beta_L - h, t.beta_L + h
    ci = [c if isinstance(c, (list, tuple)) else [np.nan, np.nan] for c in t.ci95_L]
    t["ci95_L_lo"], t["ci95_L_hi"] = [c[0] for c in ci], [c[1] for c in ci]
    return t.drop(columns=["ci95_L"])


def _promotable(r) -> bool:
    """N4: an E2-relabelled call (label_final differs from the engine class) or a call that is not headline-eligible
    (E2, NOT EVALUABLE, exploratory entry) can never be promoted by a confirmation."""
    lab, cls = getattr(r, "label_final", None), getattr(r, "class_engine", None)
    he = str(getattr(r, "headline_eligible", "True"))
    return not ((isinstance(lab, str) and isinstance(cls, str) and lab != cls) or he == "False")


def _status(r) -> str:
    if not r.eligible:
        return "not eligible (no non-neutral call; L effect reported)"
    if not r.lockbox_testable:
        if r.rare_route_pass is None or (isinstance(r.rare_route_pass, float) and np.isnan(r.rare_route_pass)):
            return f"{RARE_LABEL}: rare-context route {D2_NA}"
        return f"{RARE_LABEL}: rare-context route {'passed' if r.rare_route_pass else 'failed'}"
    if not r.testable_L:
        return "lockbox-testable (frozen list) but < 5 context or reference L models after complete case"
    if r.confirmed_final:
        return "confirmed" if _promotable(r) else "confirmed (not promotable: E2 label or not headline-eligible)"
    return "not confirmed"


def _check_target_hits(world: World, aid: str, ht: pd.DataFrame) -> None:
    """Resume: the committed target list must hold exactly these hits, eligibility and route."""
    tp = committed(world, confirm.target_path(world.proj, aid))
    rec = json.loads(tp.read_text()).get("hits", [])
    want = [(h["context"], h["unit"], float(h["beta_A"]), h["alternative"], h.get("eligible"), h.get("route"))
            for h in rec]
    got = [(str(r.context), str(r.unit), float(r.beta), str(r.alternative), bool(r.eligible), str(r.route))
           for r in ht.itertuples(index=False)]
    if want != got:
        raise RegistryError(f"{aid}: hits / eligibility / route differ from the committed target list")


def run_entry(world: World, spec: Spec, close: bool = True) -> pd.DataFrame:
    proj, aid = world.proj, spec.aid
    require_clean(world)
    check_A(world)
    e, blob = entry(world, aid)
    resuming = aid in _ledger_ids(proj)
    if resuming and not world.resume:
        raise assoc.FirewallError(f"{aid} already has a lockbox ledger line (single read; --resume re-uses the "
                                  f"persisted matrix)")
    if e["status"] != "open":
        raise RegistryError(f"{aid} is not open")
    t0 = time.time()
    seed = lockbox_seed(e)
    rare = rare_routes(world, spec)                                   # A side, before the read
    units = {u: spec.members_of(u) for u in sorted(set(spec.hits.unit))}
    genes = sorted(set(g for v in units.values() for g in v) | set(spec.extra_genes))
    ht = spec.hits.assign(route=np.where(spec.hits.lockbox_testable, "lockbox", RARE_LABEL))
    cols = ["context", "unit", "beta", "alternative", "eligible", "route"] + \
        [c for c in ("label_final", "headline_eligible", "event", "axis") if c in ht]
    h4 = ht[["context", "unit", "beta", "alternative"]]
    if resuming:
        _check_target_hits(world, aid, ht)
        world.say(f"{aid}: RESUME from the persisted read (the loader is not called)")
    else:
        confirm.write_target_list(proj, aid, genes, hits=ht[cols], units=units, contexts=set(spec.hits.context),
                                  msg_suffix=" (single read)")
        world.say(f"{aid}: target list committed ({len(genes)} genes, {len(h4)} hits); reading the lockbox")
    res = confirm.confirm_entry(aid, h4, spec.score_L, spec.ctx_all, spec.cov_all, proj=proj,
                                blocks=world.cov.PatientID, strata=world.cov.OncotreeLineage,
                                lineage_contexts=spec.lineage_contexts, n_perm=world.n_perm, seed=seed, write=True,
                                persist=world.persist and not resuming, resume=resuming)
    out = spec.hits.copy()
    if len(res):
        if list(res.unit) != list(out.unit) or list(res.context) != list(out.context):
            raise RuntimeError(f"{aid}: confirmation rows are not aligned with the hits")
        res = res.drop(columns=["context", "unit", "beta_A"]).rename(columns={"testable": "testable_L"})
        out = pd.concat([out, res.reset_index(drop=True)], axis=1)
        out = _with_ci90(out)
    else:
        for c in ("beta_L", "se_L", "p_L", "confirmed", "testable_L", "same_sign"):
            out[c] = np.nan if c not in ("confirmed", "testable_L", "same_sign") else False
    rr = pd.DataFrame.from_dict(rare, orient="index") if rare else pd.DataFrame()
    for c in ("loso_same_sign", "loso_folds", "lomo_same_sign", "demeter2_same_sign", "rare_route_pass",
              "rare_route_note"):
        out[c] = rr[c].reindex(out.index) if c in rr else np.nan
    out["route"] = np.where(out.lockbox_testable, "lockbox", RARE_LABEL)
    out["confirmed_final"] = out.eligible & out.lockbox_testable & out.confirmed.fillna(False).astype(bool)
    out["rare_route_pass"] = out.rare_route_pass.astype(object)
    out["promotable"] = [bool(r.confirmed_final and _promotable(r)) for r in out.itertuples(index=False)]
    out["status"] = [_status(r) for r in out.itertuples(index=False)]
    d = world.out(aid)
    d.mkdir(parents=True, exist_ok=True)
    out.to_csv(d / "confirmation.tsv", sep="\t", index=False, float_format="%.6g")
    led = [r for r in assoc._ledger(proj) if r["analysis_id"] == aid]
    resumed = dict(resumed=bool(resuming))
    el = out[out.eligible]
    summ = dict(analysis_id=aid, registry_blob=blob, lockbox_seed=seed, n_perm=world.n_perm, n_hits=int(len(out)),
                n_eligible=int(len(el)), n_eligible_lockbox=int(el.lockbox_testable.sum()),
                n_eligible_rare=int((~el.lockbox_testable).sum()), n_confirmed=int(out.confirmed_final.sum()),
                n_rare_pass=int((el.rare_route_pass == True).sum()),  # noqa: E712
                n_target_genes=len(genes), ledger=led, runtime_s=round(time.time() - t0, 1), **resumed,
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (d / "summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=str) + "\n")
    eng = proj.results_dir / "assoc" / aid / "confirmation_crispr.json"
    confirm.commit_paths(proj, [d, eng], f"Lockbox confirmations (single read): {aid}")
    if close:
        confirm.close_entry(proj, aid)
    world.say(f"{aid}: {summ['n_confirmed']}/{summ['n_eligible_lockbox']} confirmed; rare route "
              f"{summ['n_rare_pass']}/{summ['n_eligible_rare']}; closed={close}")
    return out


# ------------------------------------------------------------------------------------------------ late gate
def gate_late(world: World, aid: str) -> bool:
    """BENCH / BRIDGE-L: refused (before any work) unless every confirmatory entry has a committed ledger line and is
    committed with status closed (the same rule assoc.read_lockbox enforces at the read). Returns True when `aid`
    itself already read and --resume finishes it from the persisted matrix."""
    proj = world.proj
    allent = registry.load_all(proj.registry_dir)
    conf = [k for k, e in allent.items() if e["direction"] in registry.CONFIRMATORY]
    if proj.ledger.exists():
        registry.committed_blob(proj.root, proj.ledger)
    done = _ledger_ids(proj)
    missing = [k for k in conf if k not in done]
    open_ = [k for k in conf if allent[k]["status"] != "closed"]
    if missing or open_ or not conf:
        raise assoc.FirewallError(f"{aid} before every confirmation is read and closed: missing {missing[:5]}, "
                                  f"open {open_[:5]}")
    if aid in done:
        if world.resume and allent[aid]["status"] == "open":
            return True
        raise assoc.FirewallError(f"{aid} already has a lockbox ledger line (single read)")
    return False


def _bench_pairs(world: World, e: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    ctxs = primary_contexts(e)
    t = world.testable()
    tb = t[(t.entry == BENCH_ID)]
    ok = {r.context for r in tb.itertuples(index=False) if r.lockbox_testable}
    ctxs = [c for c in ctxs if c in ok]
    tabs = {k: committed_table(world.proj, k) for k in H3_IDS}
    prim = pd.concat([d[d.discovery.astype(bool) & d.context.isin(ctxs)].assign(source=k) for k, d in tabs.items()])
    sec = pd.concat([d[d.evaluable.astype(bool) & d.context.isin(ctxs)].assign(source=k) for k, d in tabs.items()])
    ca, _ = covmod.per_context(e, ctxs, world.cov.reindex(world.A))
    ctxA = world.ctx.reindex(world.A)
    rows = []
    for r in _sorted(prim[["context", "unit", "beta", "source"]]).itertuples(index=False):
        rem = cscore.removal_for_context(r.context, world.homdel.get(r.context))
        mem = world.fmap.members[r.unit]
        removed = set(rem.removed(mem, world.annot))
        use = [g for g in mem if g not in removed and g in world.E_A.columns]
        gstar = bench.choose_gstar(world.E_A, use, ctxA[r.context], ca[r.context], r.beta)
        s = np.sign(r.beta)
        key, gb = None, None
        for g in use:                          # multiverse: g* by largest beta in the A direction (ties: larger t)
            f = robust.fit_one(world.E_A[g].astype(float).clip(-3, 3), ctxA[r.context], ca[r.context])
            if f["ok"] and (key is None or (s * f["beta"], s * f["t"]) > key):
                key, gb = (s * f["beta"], s * f["t"]), g
        br = robust.breadth(world.E_A, [g for g in mem if g in world.E_A.columns], ctxA[r.context], ca[r.context], s)
        rows.append(dict(context=r.context, unit=r.unit, beta_A=float(r.beta), source=r.source, gstar=gstar,
                         gstar_beta=gb, breadth=br["breadth"], breadth_n=br["breadth_n"],
                         tier=world.fmap.tier.get(r.unit, ""), locus=bool(world.fmap.locus.get(r.unit, False))))
    return pd.DataFrame(rows), _sorted(sec[["context", "unit", "beta", "source"]]), ca, {"contexts": ctxs}


def bench_a_side(world: World, e: dict) -> dict:
    """Everything BENCH computes on A, with its validity checks, before the target list is committed (MUST-2b)."""
    pairs, sec_all, covA, info = _bench_pairs(world, e)
    bad = pairs[pairs.gstar.isna() | pairs.gstar_beta.isna() | pairs.breadth.isna()]
    if len(bad):
        raise RuntimeError(f"BENCH: g* / g*(beta) / breadth missing for {len(bad)} pairs, e.g. "
                           f"{bad[['context', 'unit']].head(3).values.tolist()}")
    terc = pd.Series(np.asarray(pd.qcut(pairs.breadth.rank(method="first"), 3, labels=["low", "mid", "high"]))
                     if len(pairs) >= 3 else ["high"] * len(pairs), index=pairs.index)   # N8: fixed once (N11 ties)
    fits = {}
    for r in pairs.itertuples(index=False):
        rem = cscore.removal_for_context(r.context, world.homdel.get(r.context))
        k = (rem, r.unit)
        if k not in fits:
            f = cscore.unit_scores(world.E_A, world.fmap, [r.unit], rem, world.annot, method="eigengene")[2]
            if r.unit not in f:
                raise RuntimeError(f"BENCH: no eigengene A fit for {r.unit} in {r.context}")
            fits[k] = f[r.unit]
    return dict(pairs=pairs, sec_all=sec_all, covA=covA, contexts=info["contexts"], tercile=terc, eig_fits=fits)


def bench_stage(world: World) -> dict:
    resuming = gate_late(world, BENCH_ID)
    require_clean(world)
    check_A(world)
    proj = world.proj
    e, blob = entry(world, BENCH_ID)
    rul = rulings(world)
    sec_spec, sec_not_run = rul["bench_secondary"], rul["bench_secondary_not_run"]
    al_not_run, distinct_not_run = _not_run(sec_not_run, "A-L correlation"), _not_run(sec_not_run, "distinct")
    ci_rule = rul["bench_multiverse_ci"]
    check_bench_rulings(sec_spec, ci_rule, sec_not_run)
    t0 = time.time()
    seeds = e["seeds"]
    a = bench_a_side(world, e)
    pairs, sec_all, ctxs, terc, eig_fits = a["pairs"], a["sec_all"], a["contexts"], a["tercile"], a["eig_fits"]
    sets = {}
    if not al_not_run:
        sec_pairs = pairs[["context", "unit"]] if sec_spec["pair_set"] == "primary" else sec_all[["context", "unit"]]
        units_sec = sorted(set(sec_pairs.unit))
        ess = calib.clip_essentiality(world.E_A)
        ps = world.pair_scores() if callable(world.pair_scores) else world.pair_scores
        sets, set_status = calib.matched_random_sets(world.E_A, world.fmap, units_sec, ess, ps, int(sec_spec["R"]),
                                                     int(seeds["matched_sets"]))
    genes = sorted({g for u in set(pairs.unit) | set(sec_all.unit) for g in world.fmap.members[u]}
                   | {g for v in sets.values() for g in v})
    p4 = pairs.rename(columns={"beta_A": "beta"}).assign(alternative="two-sided")[["context", "unit", "beta",
                                                                                     "alternative"]]
    tp = confirm.target_path(proj, BENCH_ID)
    if resuming:
        want = [(h["context"], h["unit"], float(h["beta_A"]), h["alternative"])
                for h in json.loads(committed(world, tp).read_text()).get("hits", [])]
        if confirm.hit_keys(p4) != want:
            raise RegistryError("BENCH: pairs differ from the committed target list")
        world.say("BENCH: RESUME from the persisted read (the loader is not called)")
        E_L = assoc.resume_lockbox_read(BENCH_ID, proj=proj)
    else:
        confirm.write_target_list(proj, BENCH_ID, genes, hits=p4, contexts=set(ctxs), msg_suffix=" (single read)")
        world.say(f"BENCH: target list committed ({len(genes)} genes, {len(pairs)} primary pairs); reading the lockbox")
        E_L = assoc.read_lockbox(BENCH_ID, tp, proj=proj, persist=world.persist)
    Lidx = E_L.index
    ctxL = world.ctx.reindex(Lidx)
    covL, _ = covmod.per_context(e, ctxs, world.cov.reindex(Lidx))
    blocksL = world.cov.PatientID.reindex(Lidx)

    def rem(c):
        return cscore.removal_for_context(c, world.homdel.get(c))

    def sL(E, c, u):
        return cscore.unit_scores(E, world.fmap, [u], rem(c), world.annot)[0][u]

    B = world.B_bench
    prim = bench.bench_delta(pairs, sL, E_L, ctxL, covL, blocks_L=blocksL, B=B, seed=int(seeds["bootstrap"]),
                             tercile=terc)
    tab = prim["table"]
    rho = sst.spearmanr(tab.delta, tab.breadth, nan_policy="omit")
    out = dict(analysis_id=BENCH_ID, registry_blob=blob, n_L_models=int(len(Lidx)), B=B,
               seed=int(seeds["bootstrap"]), primary=prim["summary"],
               trend_spearman_delta_breadth=dict(rho=float(rho.statistic), p=float(rho.pvalue)))
    out["distinct"] = "NOT RUN" if distinct_not_run else bench.distinct_gstar(pairs, world.fmap)   # A5: reported
    # multiverse (reported, never selected); within_line_rank NOT RUN (A3(b))
    Bm = B if ci_rule == "bootstrap as primary" else 1
    def sL_eig(E, c, u):                         # A-fitted eigengene (fits computed before the read)
        return cscore.unit_scores(E, world.fmap, [u], rem(c), world.annot, method="eigengene",
                                  fits={u: eig_fits[(rem(c), u)]})[0][u]

    mv = {}
    t12 = pairs.tier.isin(["T1", "T2"]).values
    for name, pp, fn, tt in (("eigengene", pairs, sL_eig, terc),
                             ("gstar_largest_beta", pairs.assign(gstar=pairs.gstar_beta), sL, terc),
                             ("T1_T2_pairs", pairs[t12].reset_index(drop=True), sL,
                              terc[t12].reset_index(drop=True))):
        if len(pp):
            r = bench.bench_delta(pp, fn, E_L, ctxL, covL, blocks_L=blocksL, B=Bm, seed=int(seeds["bootstrap"]),
                                  tercile=tt)
            s = r["summary"]
            if Bm == 1:
                s = {k: v for k, v in s.items() if not k.startswith("ci95") and k != "top_tercile_wins"}
            mv[name] = dict(s, n_pairs=int(len(pp)))
    mv["within_line_rank"] = "NOT RUN (amendment A3(b))"
    out["multiverse"] = mv
    out["multiverse_ci"] = ci_rule
    # secondary: A-L correlation vs matched random sets (NOT RUN under A4 G1 / A5)
    out["secondary"] = "NOT RUN (amendments A4 G1, A5)"
    al = None
    if not al_not_run:
        al = _bench_secondary(world, e, ctxs, sets, set_status, sec_spec, sec_pairs, rem, sL, E_L, ctxL, covL)
        out["secondary"] = dict({k: v for k, v in al.items() if k != "table"}, spec=sec_spec,
                                matched_sets_status=set_status.status.value_counts().to_dict(),
                                seed=int(seeds["matched_sets"]))
    out["runtime_s"] = round(time.time() - t0, 1)
    out["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    d = world.out(BENCH_ID)
    d.mkdir(parents=True, exist_ok=True)
    tab.to_csv(d / "bench_pairs.tsv", sep="\t", index=False, float_format="%.6g")
    if al is not None:
        al["table"].to_csv(d / "al_correlation_pairs.tsv", sep="\t", index=False, float_format="%.6g")
    (d / "bench.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=float) + "\n")
    confirm.commit_paths(proj, [d], f"Lockbox confirmations (single read): {BENCH_ID}")
    confirm.close_entry(proj, BENCH_ID, outputs=[d / "bench.json"])
    world.say(f"BENCH: mean Delta {out['primary']['mean']:.3f}; closed")
    return out


def _bench_secondary(world, e, ctxs, sets, set_status, sec_spec, sec_pairs, rem, sL, E_L, ctxL, covL) -> dict:
    """BENCH secondary (only if registered; A4 G1 declared it NOT RUN): bench.al_correlation vs matched sets."""
    ctxA = world.ctx.reindex(world.A)
    covAll, _ = covmod.per_context(e, ctxs, world.cov.reindex(world.A))
    removal_cache = {}

    def sA(c, u):
        k = (rem(c), "A")
        if k not in removal_cache:
            removal_cache[k] = cscore.unit_scores(world.E_A, world.fmap, world.fmap.coarse(), rem(c), world.annot)[0]
        return removal_cache[k][u]

    def setA(c, u, r):
        g = sets.get(f"{u}#r{r}")
        return None if g is None else calib.score_sets(world.E_A, {"s": g}, rem(c), world.annot)["s"]

    def setL(c, u, r):
        return calib.score_sets(E_L, {"s": sets[f"{u}#r{r}"]}, rem(c), world.annot)["s"]

    return bench.al_correlation(sec_pairs, sA, lambda c, u: sL(E_L, c, u), ctxA, ctxL, covAll, covL, setA, setL,
                                R=int(sec_spec["R"]), method=str(sec_spec["correlation"]))


def check_bench_rulings(sec_spec, ci_rule, not_run=None) -> None:
    if ci_rule not in ("bootstrap as primary", "point estimate only"):
        raise GapError("BENCH multiverse: `multiverse.ci` ruling missing ('bootstrap as primary' | "
                       "'point estimate only')")
    if _not_run(not_run, "A-L correlation"):      # A4 (G1) / A5: the A-L random-set correlation is NOT RUN
        return
    need = ("correlation", "R", "pair_set", "matching", "summary")
    if not isinstance(sec_spec, dict) or any(k not in sec_spec for k in need):
        raise GapError(f"BENCH secondary: registry ruling `estimand_secondary` with {need} is missing")
    if sec_spec["correlation"] not in ("pearson", "spearman") or sec_spec["pair_set"] not in ("primary", "secondary"):
        raise GapError("BENCH secondary: correlation must be pearson|spearman and pair_set primary|secondary")
    if "matched_random_sets" not in str(sec_spec["matching"]):
        raise GapError("BENCH secondary: only calib.matched_random_sets is implemented as the matching rule")


def gls_discovery(world: World) -> pd.DataFrame:
    """The registered GLS-map discovery pass on A (two run_association calls), committed before the BRIDGE-L read."""
    proj = world.proj
    e, _ = entry(world, BRIDGE_ID)
    if proj.trials.exists() and any(json.loads(l).get("analysis_id") == BRIDGE_ID
                                    for l in proj.trials.read_text().splitlines() if l.strip()):
        raise RegistryError("BRIDGE-L: trials.jsonl already holds a GLS discovery line (single pass)")
    units = world.gls.coarse()
    if registry.set_hash(units) != e["complex_set_hash"]:
        raise RegistryError("GLS coarse clusters do not match the BRIDGE-L complex_set_hash")
    reg = set(primary_contexts(e))
    removal, frame = _h3_scores(world, world.gls)
    A = world.A
    res_all, paths = [], []
    for tag, src in (("discovery_H3a", H3_IDS[0]), ("discovery_H3b", H3_IDS[1])):
        ctxs = primary_contexts(entry(world, src)[0])
        if not set(ctxs) <= reg:
            raise RegistryError(f"BRIDGE-L context family lacks {sorted(set(ctxs) - reg)[:3]}")
        Y = {c: frame(world.E_A, c, units, "A") for c in ctxs}
        covs, _ = covmod.per_context(e, ctxs, world.cov.reindex(A))
        res = assoc.run_association(BRIDGE_ID, Y, world.ctx.reindex(A)[ctxs], covs, inputs=world.inputs_gls,
                                    proj=proj, blocks=world.cov.PatientID.reindex(A),
                                    strata=world.cov.OncotreeLineage.reindex(A), n_perm=world.n_perm,
                                    perm_screen_q=SCREEN_Q, perm_seed=int(e["seeds"]["permutation"]),
                                    min_context_n=MIN_A, min_ref_n=MIN_A, release_hits=True, tag=tag)
        res_all.append(res.assign(family=tag))
        paths.append(proj.results_dir / "assoc" / BRIDGE_ID / f"{tag}_summary.json")
        world.say(f"BRIDGE-L GLS {tag}: {int(res.discovery.sum())} discoveries / {int(res.evaluable.sum())} tests")
    res = pd.concat(res_all, ignore_index=True)
    if world.expect_gls is not None:                   # N13: same deterministic pass as the rehearsal (REHEARSAL.md)
        got = {tag: int(r.discovery.sum()) for tag, r in zip(("discovery_H3a", "discovery_H3b"), res_all)}
        if got != world.expect_gls:
            raise RuntimeError(f"BRIDGE-L GLS discovery counts {got} differ from the rehearsal {world.expect_gls}; "
                               f"stopped before the BRIDGE-L target list")
    d = world.out(BRIDGE_ID)
    d.mkdir(parents=True, exist_ok=True)
    res[res.discovery].to_csv(d / "gls_discoveries.tsv", sep="\t", index=False, float_format="%.6g")
    confirm.commit_paths(proj, [proj.trials, d] + paths,
                         "BRIDGE-L: registered GLS-map discovery on A (before its lockbox read)")
    return res


def bridge_a_side(world: World, e: dict) -> dict:
    """GLS frozen inputs and the D5 coverage matching, checked before any BRIDGE-L work (MUST-2b)."""
    man = list(dict.fromkeys(["results/part2/lockbox/manifest.json"] + list(e["frozen_inputs"])))
    assoc.check_frozen(world.proj, world.inputs_gls, man)
    if world.coverage_tables is None:
        raise GapError("BRIDGE-L multiverse: no coverage tables for the matched gene coverage subset")
    mz, mg, cov = validation.match_coverage(world.coverage_tables["zavg"], world.coverage_tables["gls"])
    return dict(keep={("zavg", u) for u in mz.cluster_id} | {("gls", u) for u in mg.cluster_id}, coverage=cov)


def _gls_trials(proj) -> list:
    if not proj.trials.exists():
        return []
    return [json.loads(l) for l in proj.trials.read_text().splitlines()
            if l.strip() and json.loads(l).get("analysis_id") == BRIDGE_ID]


def bridge_stage(world: World) -> dict:
    resuming = gate_late(world, BRIDGE_ID)
    require_clean(world)
    check_A(world)
    proj = world.proj
    e, blob = entry(world, BRIDGE_ID)
    ci_rule = rulings(world)["bridge_multiverse_ci"]
    if ci_rule not in ("bootstrap as primary", "point estimate only"):
        raise GapError("BRIDGE-L multiverse: `multiverse.ci` ruling missing")
    t0 = time.time()
    a = bridge_a_side(world, e)
    if world.resume and _gls_trials(proj):          # the GLS pass already committed: re-use it, never re-run
        res = pd.concat([committed_table(proj, BRIDGE_ID, tag).assign(family=tag)
                         for tag in ("discovery_H3a", "discovery_H3b")], ignore_index=True)
    else:
        res = gls_discovery(world)
    spec = spec_h3(world, BRIDGE_ID, fm=world.gls, disc=res)
    conf = run_entry(world, spec, close=False)
    rows = []
    for k in H3_IDS:
        z = pd.read_csv(world.out(k) / "confirmation.tsv", sep="\t")
        z = z[z.route == "lockbox"]
        rows.append(pd.DataFrame(dict(map="zavg", source=k, context=z.context, unit=z.unit,
                                      confirmed=z.confirmed_final.astype(bool))))
    g = conf[conf.route == "lockbox"]
    rows.append(pd.DataFrame(dict(map="gls", source="BRIDGE-L", context=g.context, unit=g.unit,
                                  confirmed=g.confirmed_final.astype(bool))))
    hits = pd.concat(rows, ignore_index=True)
    fm = {"zavg": world.fmap, "gls": world.gls}
    hits["size"] = [len(fm[m].members[u]) for m, u in zip(hits["map"], hits.unit)]
    hits["tier"] = [fm[m].tier.get(u, "") for m, u in zip(hits["map"], hits.unit)]
    B, seed = world.B_bridge, int(e["seeds"]["bootstrap"])
    out = dict(analysis_id=BRIDGE_ID, registry_blob=blob, B=B, seed=seed,
               n_hits={m: int((hits["map"] == m).sum()) for m in ("zavg", "gls")},
               primary=bench.bridge_rate(hits, B=B, seed=seed))
    Bm = B if ci_rule == "bootstrap as primary" else 0
    mv = dict(unweighted=bench.bridge_unweighted(hits, B=Bm, seed=seed))
    keep = a["keep"]
    sub = hits[[(m, u) in keep for m, u in zip(hits["map"], hits.unit)]]
    mv["matched_gene_coverage"] = dict(_bridge_mv(sub, Bm, seed), coverage_matching=a["coverage"],
                                       n_hits={m: int((sub["map"] == m).sum()) for m in ("zavg", "gls")})
    sub = hits[hits.tier.isin(["T1", "T2"])]
    mv["T1_T2_clusters"] = dict(_bridge_mv(sub, Bm, seed), n_hits={m: int((sub["map"] == m).sum()) for m in ("zavg", "gls")})
    out["multiverse"], out["multiverse_ci"] = mv, ci_rule
    out["runtime_s"] = round(time.time() - t0, 1)
    out["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    d = world.out(BRIDGE_ID)
    hits.to_csv(d / "bridge_hits.tsv", sep="\t", index=False)
    (d / "bridge.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=float) + "\n")
    confirm.commit_paths(proj, [d], f"Lockbox confirmations (single read): {BRIDGE_ID} estimate")
    confirm.close_entry(proj, BRIDGE_ID, outputs=[d / "bridge.json",
                                                   proj.results_dir / "assoc" / BRIDGE_ID / "confirmation_crispr.json"])
    world.say(f"BRIDGE-L: diff {out['primary']['diff']:.3f} {out['primary']['ci95']}; closed")
    return out


def _bridge_mv(hits: pd.DataFrame, B: int, seed: int) -> dict:
    if hits["map"].nunique() < 2:
        return dict(note="one map has no hits in this subset")
    if B > 0:
        return bench.bridge_rate(hits, B=B, seed=seed)
    r = bench.bridge_rate(hits, B=1, seed=seed)
    return {k: v for k, v in r.items() if k not in ("ci95", "supports_bridge")}


# ------------------------------------------------------------------------------------------------ report
def _flags(world: World, aid: str, tab: pd.DataFrame) -> pd.Series:
    """D8 'beats all matched sets' (H3a/H3b) / 'beats all its random splits' (H2): |t| of the discovery above the
    largest |t| of its unit's sets in that context; NaN when no set exists or no calibration table is given."""
    ct = world.calib_tests.get(aid)
    if ct is None or "t" not in tab:
        return pd.Series(np.nan, index=tab.index)
    m = ct[ct.calibration.isin(["matched_random_set", "random_split"])]
    mx = m.assign(at=m.t.abs()).groupby(["context", "unit"]).at.max()
    v = [mx.get((c, u), np.nan) for c, u in zip(tab.context, tab.unit)]
    return pd.Series([np.nan if not np.isfinite(x) else bool(abs(t) > x) for t, x in zip(tab.t, v)],
                     index=tab.index, dtype=object)


def _rate(df: pd.DataFrame) -> dict:
    lb = df[df.eligible & df.lockbox_testable]
    rr = df[df.eligible & ~df.lockbox_testable]
    return dict(n_eligible=int(df.eligible.sum()), n_lockbox=int(len(lb)), n_confirmed=int(lb.confirmed_final.sum()),
                rate=float(lb.confirmed_final.mean()) if len(lb) else None, n_rare=int(len(rr)),
                n_rare_pass=int((rr.rare_route_pass == True).sum()),  # noqa: E712
                n_rare_not_evaluable=int(rr.rare_route_pass.isna().sum()))


def report(world: World) -> dict:
    out, rates = {}, []
    order = confirm_order(world.proj) + ([BRIDGE_ID] if (world.out(BRIDGE_ID) / "confirmation.tsv").exists() else [])
    for aid in order:
        p = world.out(aid) / "confirmation.tsv"
        if not p.exists():
            continue
        t = pd.read_csv(p, sep="\t")
        for c in ("eligible", "lockbox_testable", "confirmed_final"):
            t[c] = t[c].astype(str).str.lower() == "true"
        t["rare_route_pass"] = t.rare_route_pass.map({True: True, False: False, "True": True, "False": False})
        fm = world.gls if aid == BRIDGE_ID else world.fmap
        t["locus"] = [bool(fm.locus.get(str(u).split("|")[0], False)) for u in t.unit]
        t["beats_all"] = _flags(world, aid, t)
        r = dict(overall=_rate(t), non_locus=_rate(t[~t.locus]))
        if t.beats_all.notna().any():
            r["beats_all_matched_sets"] = _rate(t[t.beats_all == True])  # noqa: E712
            r["not_beats_all"] = _rate(t[t.beats_all == False])  # noqa: E712
            r["no_sets"] = _rate(t[t.beats_all.isna()])
        out[aid] = r
        for k, v in r.items():
            rates.append(dict(entry=aid, subset=k, **v))
        if aid in H1_IDS:
            cols = ["event", "axis", "class_engine", "label_final", "headline_eligible", "eligible", "route",
                    "beta", "beta_L", "ci95_L_lo", "ci95_L_hi", "ci90_L_lo", "ci90_L_hi", "p_L", "p_L_perm",
                    "n_context_L", "n_ref_L", "confirmed_final", "promotable", "rare_route_pass", "rare_route_note",
                    "status"]
            t[[c for c in cols if c in t]].rename(columns={"beta": "beta_A"}).to_csv(
                world.out() / f"{aid}_events.tsv", sep="\t", index=False, float_format="%.6g")
    for aid, f in ((BENCH_ID, "bench.json"), (BRIDGE_ID, "bridge.json")):
        p = world.out(aid) / f
        if p.exists():
            out[aid] = json.loads(p.read_text())
    d = world.out()
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rates).to_csv(d / "rates.tsv", sep="\t", index=False, float_format="%.6g")
    (d / "REPORT.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=float) + "\n")
    return out


# ------------------------------------------------------------------------------------------------ preflight / run
def preflight(world: World) -> dict:
    """A side only (no L value is touched): entries, order, seeds, frozen list coverage, hits, routes, rulings."""
    proj = world.proj
    gaps, entries = [], {}
    order = confirm_order(proj)
    led = _ledger_ids(proj)
    for fn in (require_clean, check_A):
        try:
            fn(world)
        except (RegistryError, assoc.FirewallError) as ex:
            gaps.append(str(ex))
    for aid in order + [BENCH_ID, BRIDGE_ID]:
        e, blob = entry(world, aid)
        info = dict(direction=e["direction"], status=e["status"], blob=blob)
        try:
            info["lockbox_seed"] = lockbox_seed(e) if aid != BENCH_ID else None
        except GapError as ex:
            gaps.append(str(ex))
        if aid in led and not (world.resume and e["status"] == "open"):
            gaps.append(f"{aid}: a lockbox ledger line already exists")
        if confirm.target_path(proj, aid).exists() and not (world.resume and aid in led):
            gaps.append(f"{aid}: a target list already exists")
        if aid in order:
            try:
                sp = build_spec(world, aid)
            except (RegistryError, RuntimeError, assoc.FreezeError, KeyError) as ex:
                gaps.append(f"{aid}: {ex}")
                entries[aid] = info
                continue
            h = sp.hits
            info.update(n_hits=int(len(h)), n_eligible=int(h.eligible.sum()),
                        n_rare=int((h.eligible & ~h.lockbox_testable).sum()),
                        rare_contexts=sorted(set(h.context[h.eligible & ~h.lockbox_testable])))
            fit_needed = [u for u in sp.needs_fit
                          if not bool(h.lockbox_testable[h.unit == u].all()) and bool(h.eligible[h.unit == u].any())]
            if fit_needed:
                try:
                    info["rare_route_fitted_parameters"] = _fit_ruling(world, aid)
                except GapError as ex:
                    gaps.append(f"{ex} [hits: {fit_needed[:8]}]")
                info["needs_fitted_parameters"] = sorted(fit_needed)
            reg_fit = (e.get("lockbox") or {}).get("rare_route_fitted_hits")      # review N14
            if reg_fit is not None and sorted(reg_fit) != sorted(fit_needed):
                gaps.append(f"{aid}: registered rare_route_fitted_hits {sorted(reg_fit)} differ from the hits that "
                            f"need fitted parameters {sorted(fit_needed)}")
        entries[aid] = info
    rul = rulings(world)
    try:
        check_bench_rulings(rul["bench_secondary"], rul["bench_multiverse_ci"], rul["bench_secondary_not_run"])
    except GapError as ex:
        gaps.append(str(ex))
    if rul["bridge_multiverse_ci"] not in ("bootstrap as primary", "point estimate only"):
        gaps.append("BRIDGE-L multiverse: `multiverse.ci` ruling missing ('bootstrap as primary' | "
                    "'point estimate only')")
    eb, _ = entry(world, BRIDGE_ID)
    if registry.set_hash(world.gls.coarse()) != eb["complex_set_hash"]:
        gaps.append("GLS map does not match the BRIDGE-L complex_set_hash")
    try:                                              # MUST-2b / S3: late A-side failures surface here, before any read
        b = bridge_a_side(world, eb)
        entries[BRIDGE_ID]["coverage_matching"] = b["coverage"]
    except Exception as ex:                           # noqa: BLE001 -- every failure is a gap
        gaps.append(f"BRIDGE-L A side: {type(ex).__name__}: {ex}")
    try:
        ba = bench_a_side(world, entry(world, BENCH_ID)[0])
        entries[BENCH_ID].update(n_pairs=int(len(ba["pairs"])), n_secondary_pairs=int(len(ba["sec_all"])),
                                 n_eigengene_fits=len(ba["eig_fits"]),
                                 tercile_counts=ba["tercile"].value_counts().to_dict())
    except Exception as ex:                           # noqa: BLE001
        gaps.append(f"BENCH A side: {type(ex).__name__}: {ex}")
    return dict(order=order + [BENCH_ID, BRIDGE_ID], entries=entries, gaps=gaps, rulings=rul)


def run_all(world: World) -> dict:
    pf = preflight(world)
    if pf["gaps"]:
        raise GapError("preflight: " + " | ".join(pf["gaps"]))
    t0 = time.time()
    for aid in confirm_order(world.proj):
        if world.resume and confirm.is_closed(world.proj, aid):
            continue
        run_entry(world, build_spec(world, aid))
    b = bench_stage(world) if not (world.resume and confirm.is_closed(world.proj, BENCH_ID)) else \
        json.loads((world.out(BENCH_ID) / "bench.json").read_text())
    g = bridge_stage(world) if not (world.resume and confirm.is_closed(world.proj, BRIDGE_ID)) else \
        json.loads((world.out(BRIDGE_ID) / "bridge.json").read_text())
    rep = report(world)
    confirm.commit_paths(world.proj, [world.out() / "REPORT.json", world.out() / "rates.tsv"] +
                         sorted(world.out().glob("*_events.tsv")),
                         "Lockbox confirmations (single read): reporting layer")
    return dict(report=rep, bench=b, bridge=g, runtime_s=round(time.time() - t0, 1),
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


# ------------------------------------------------------------------------------------------------ R6 relook
def _finite(v) -> bool:
    return v is not None and isinstance(v, (int, float)) and np.isfinite(v)


def relook(world: World, ruling: str, mode: str) -> pd.DataFrame:
    """PI ruling R6: a disclosed, NON-CONFIRMATORY relook through the resume path (assoc.resume_lockbox_read with
    relook=<ruling>; never the loader; no ledger line -- it is not a read). The L covariate design is the fixed one
    (assembled over the L rows read). mode 'check': every committed successful L fit (beta_L finite) is recomputed and
    returned next to the committed values (reproduction check); 'refit': only the lockbox-route fits that failed.
    Seeds, hit order, blocks and strata are those of the single read (lockbox_permutation + hit index)."""
    if mode not in ("check", "refit"):
        raise ValueError(mode)
    proj, rows = world.proj, []
    for aid in confirm_order(proj) + [BRIDGE_ID]:
        e, _ = entry(world, aid)
        if aid == BRIDGE_ID:
            res = pd.concat([committed_table(proj, BRIDGE_ID, tag).assign(family=tag)
                             for tag in ("discovery_H3a", "discovery_H3b")], ignore_index=True)
            spec = spec_h3(world, BRIDGE_ID, fm=world.gls, disc=res)
        else:
            spec = build_spec(world, aid)
        eng = json.loads(committed(world, proj.results_dir / "assoc" / aid / "confirmation_crispr.json").read_text())
        eh = eng["hits"]
        if [(h["context"], h["unit"]) for h in eh] != list(zip(spec.hits.context, spec.hits.unit)):
            raise RuntimeError(f"{aid}: committed confirmation rows differ from the hits")
        if mode == "check":
            sel = [i for i, h in enumerate(eh) if _finite(h["beta_L"])]
        else:
            sel = [i for i, h in enumerate(eh) if not _finite(h["beta_L"]) and bool(spec.hits.lockbox_testable.iloc[i])
                   and bool(spec.hits.eligible.iloc[i])]
        if not sel:
            continue
        E_L = assoc.resume_lockbox_read(aid, proj=proj, relook=ruling)
        covL = spec.cov_all(E_L.index)
        seed = lockbox_seed(e)
        for i in sel:
            r = spec.hits.iloc[i]
            st = None if r.context in spec.lineage_contexts else world.cov.OncotreeLineage
            y = spec.score_L(E_L, r.context, r.unit)
            tt = confirm.lockbox_test(y, spec.ctx_all[r.context], covL.get(r.context), float(r.beta),
                                      blocks=world.cov.PatientID, strata=st, n_perm=world.n_perm, seed=seed + i)
            h = eh[i]
            row = dict(entry=aid, hit_index=i, context=r.context, unit=r.unit, beta_A=float(r.beta),
                       eligible=bool(r.eligible), lockbox_testable=bool(r.lockbox_testable))
            for k in ("beta_L", "se_L", "p_L", "p_L_perm", "n_L", "n_context_L", "confirmed"):
                row[f"{k}_committed"] = h.get(k)
                row[f"{k}_relook"] = tt[k]
            lo, hi = tt["ci95_L"]
            row.update(ci95_L_lo_relook=lo, ci95_L_hi_relook=hi, testable_L_relook=tt["testable"],
                       same_sign_relook=tt["same_sign"])
            rows.append(row)
        world.say(f"R6 relook ({mode}) {aid}: {len(sel)} fits")
    return pd.DataFrame(rows)


def check_identical(chk: pd.DataFrame, tol: float = 1e-10) -> dict:
    """Reproduction check: relook estimates and p-values equal the committed ones (relative tolerance for float
    round-off only)."""
    out = dict(n=int(len(chk)))
    bad = np.zeros(len(chk), bool)
    for k in ("beta_L", "se_L", "p_L"):
        a = chk[f"{k}_committed"].astype(float).values
        b = chk[f"{k}_relook"].astype(float).values
        d = np.abs(a - b) > tol * np.maximum(1.0, np.abs(a))
        d |= np.isnan(a) != np.isnan(b)
        bad |= d
        out[f"max_abs_diff_{k}"] = float(np.nanmax(np.abs(a - b))) if len(a) else 0.0
    flips = chk.confirmed_committed.astype(bool) != (chk.testable_L_relook & chk.same_sign_relook &
                                                     (chk.p_L_relook.astype(float) < 0.05))
    out.update(n_differ=int(bad.sum()), identical=bool(not bad.any()),
               n_differ_by_entry=chk[bad].entry.value_counts().to_dict(),
               n_confirmed_status_differs=int(flips.sum()))
    return out


def real_world(outcome_loader=None) -> World:
    """The real inputs. outcome_loader: only for the dress rehearsal in a throw-away clone (scripts/verify/
    lockbox_rehearsal.py); the real run uses the project loader."""
    R44 = _script("44_run_direction.py")
    R34 = _script("34_run_controls.py")
    env = R44.real_env()
    root = env.proj.root
    if outcome_loader is not None:
        env.proj = assoc.Project(root=root, outcome_loader=outcome_loader)
    for rel in (R44.H1_EVENTS, R44.H1_SPLIT, GLS_MAP, FZX_MANIFEST, *CLUSTERS.values()):
        registry.committed_blob(root, root / rel)                     # MUST-1
    ctx_full = pd.read_parquet(root / R44.CTX)
    loss = pd.read_parquet(root / R44.H1_LOSS)
    cov = covmod.source_table(ctx_full, ctx_full.index)
    gls = FrozenMap.from_membership(pd.read_csv(root / GLS_MAP, sep="\t"))
    cl = {k: pd.read_csv(root / p, sep="\t") for k, p in CLUSTERS.items()}
    cl = {k: v[v.level == 1] for k, v in cl.items()}
    ct = {}
    for k in ("H3a_lineage_atlas", "H3b_driver_family", "H2_submodule_divergence"):
        p = root / "results" / "part2" / "directions" / k / "calibration_D8_tests.tsv.gz"
        if p.exists():
            registry.committed_blob(root, p)
            ct[k] = pd.read_csv(p, sep="\t")
    man = json.loads((root / FZX_MANIFEST).read_text())["files"]
    return World(proj=env.proj, fmap=env.fmap, gls=gls, ctx=ctx_full, loss=loss, cov=cov, E_A=env.E_A,
                 annot=env.annot, homdel=env.homdel, ctrl_homdel=R34.homdel_shares(ctx_full, env.E_A.index),
                 events={k: R44.load_events(k, root) for k in H1_IDS}, inputs_gls=[R44.CTX, GLS_MAP, R44.HOMDEL],
                 pair_scores=lambda: R44.load_pair_scores(root), coverage_tables=cl, calib_tests=ct,
                 testable_sha=man[TESTABLE_REL], expect_gls=dict(discovery_H3a=1206, discovery_H3b=139))


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--dry-check", action="store_true")
    g.add_argument("--run", action="store_true")
    ap.add_argument("--resume", action="store_true", help="with --run: finish entries already read, from their "
                                                          "persisted matrices; the loader is never called for them")
    g.add_argument("--relook", metavar="RULING", help="PI-ruled non-confirmatory relook through the resume path "
                                                      "(R6); needs --mode and --out")
    ap.add_argument("--mode", choices=["check", "refit"])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    t0 = time.time()
    world = real_world()
    world.resume = bool(a.resume)
    if a.relook:
        if not (a.mode and a.out):
            raise SystemExit("--relook needs --mode and --out")
        tab = relook(world, a.relook, a.mode)
        tab.to_csv(a.out, sep="\t", index=False, float_format="%.10g")
        res = dict(mode=a.mode, n=int(len(tab)), runtime_s=round(time.time() - t0, 1),
                   peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        if a.mode == "check":
            res.update(check_identical(tab))
        print(json.dumps(res, indent=1, default=str))
        return
    if a.dry_check:
        pf = preflight(world)
        pf["dry_runtime_s"] = round(time.time() - t0, 1)
        pf["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        print(json.dumps(pf, indent=1, default=str))
        return
    out = run_all(world)
    print(json.dumps(dict(runtime_s=out["runtime_s"], peak_rss_bytes=out["peak_rss_bytes"]), indent=1))


if __name__ == "__main__":
    main()
