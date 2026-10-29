"""A-side discovery run of ONE registered direction (DESIGN v3.1 §2.5-2.6; registry order H3a_lineage_atlas ->
H3b_driver_family -> H1_within_complex_GI -> H1S_silencing_GI -> H2_submodule_divergence). Single look: refuses to
run when results/part2/directions/<ID>/summary.json exists or trials.jsonl already holds a discovery line for <ID>.

Thin runner: it loads the committed registry entry and the frozen inputs, calls the engine and writes outputs. It
adds no analysis logic. Engine calls:
  H3a / H3b   cscore.unit_scores (context removals) -> robust.pseudo_context_calibration (registry `calibration`;
              pre-declared response: raw p < 0.01 also re-tested) -> assoc.run_association (BB, Freedman-Lane
              permutation re-entry)
  H1 / H1S    h1.assoc_inputs (frozen FZ-X event list) -> assoc.run_association (pooled BH, permutation re-entry)
              -> h1.classify (TOST, P-vs-N contrast, breadth, matched random-set null K = 1000, seed 7000 + index)
  H2          assoc.parent_selection_from (committed H3a / H3b discovery tables) -> h2.h2_frames ->
              assoc.run_association (stagewise_bb)
  --sensitivity SMARCB1_CN (after H1; D7, not a trial): h1.classify for Z0109|SMARCB1 with the registry's CN-type
              covariate set (h1.COVARIATES["CN"]: + Ploidy, WGD, TP53_altered); q / discovery flags are those of the
              primary H1 run. No trial line.
  --dry-check  real registry, manifests, contexts, covariates, map and A index; synthetic N(-0.5,1) outcome in place
              of CRISPRGeneEffect (column names and model IDs from the cache sidecar), random pair scores; the
              trial ledger and every writer are disabled. No dependency value is read.

Usage (from the project root, PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/44_run_direction.py --id H3a_lineage_atlas [--dry-check]
    .venv/bin/python scripts/44_run_direction.py --id H1_within_complex_GI --sensitivity SMARCB1_CN
Outputs (results/part2/directions/<ID>/): summary.json, tests.tsv.gz (every test, 95% CI), discoveries.tsv,
calibration.json (H3a/H3b), classes.tsv (H1/H1S), plus the engine's results/part2/assoc/<ID>/discovery_summary.json,
data_cache/part2/assoc/<ID>/discovery.parquet (gitignored) and one trials.jsonl line.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import resource
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from part2 import assoc, covariates as covmod, cscore, h1, h2, registry, robust, stats  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402

MAP = "results/part2/fzx/zavg_engine_membership.tsv"
CTX = "data_cache/part2/contexts/context_table.parquet"
HOMDEL = "results/part2/fzx/lof_homdel_shares.tsv"
H1_LOSS = "data_cache/part2/fzx/h1_loss_table.parquet"
H1_EVENTS = "results/part2/fzx/h1_events.tsv"
H1_SPLIT = "results/part2/fzx/h1_event_split_A2.tsv"
ZAVG_A_BLOCK = "data_cache/part2/blocks/zavg_BEF_bd28f78d612bef9f_S8c264ca513f8.npz"   # FZ-C zavg_BEF S x S on A
A_HASH = "bd28f78d612bef9f"
ORDER = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI",
         "H2_submodule_divergence")

N_PERM = 10_000          # registry `statistic`: 10,000 permutations, GPD tail only when p_perm < 1e-4
SCREEN_Q = 0.1           # registry: screening hits (adjusted q < 0.1) re-tested by permutation
MIN_N = 12               # registry: >= 12 context and >= 12 reference A models
K_NULL = 1000            # H1 matched null: K = 1000 random member sets
N_PSEUDO = 100           # registry `calibration`: 100 pseudo-contexts per direction
CAL_THRESHOLD = 0.075    # pre-declared response threshold
CAL_RAW_P = 0.01         # pre-declared response: every real test with p < 0.01 re-tested
DEFAULT_SEEDS = {"H3a": dict(permutation=4000, pseudo_contexts=4600), "H3b": dict(permutation=4100, pseudo_contexts=4700),
                 "H1": dict(permutation=4200, matched_null=7000), "H1S": dict(permutation=4200, matched_null=7000),
                 "H2": dict(permutation=4300, random_splits=4800)}   # synthetic entries carry no `seeds`


@dataclass
class Env:
    proj: assoc.Project
    fmap: FrozenMap
    ctx: pd.DataFrame            # A x context columns
    E_A: pd.DataFrame            # A x gene keys
    cov: pd.DataFrame            # A x covariate source columns
    annot: cscore.Annotation
    inputs: list
    blocks: pd.Series            # PatientID over A
    strata: pd.Series            # OncotreeLineage over A
    homdel: dict = field(default_factory=dict)
    strict: bool = True          # real run: registry must carry seeds; outputs refused twice
    write: bool = True


# ------------------------------------------------------------------------------------------------ helpers
def entry_of(env: Env, aid: str):
    return registry.load_entry(env.proj.registry_dir, aid, root=env.proj.root)


def primary_contexts(entry: dict) -> list:
    return [c["id"] for c in entry["context_family"]["contexts"] if c.get("role", "primary") == "primary"]


def seeds_of(env: Env, entry: dict) -> dict:
    s = entry.get("seeds")
    if s is None:
        if env.strict:
            raise registry.RegistryError(f"{entry['id']} has no `seeds`")
        s = DEFAULT_SEEDS[entry["direction"]]
    return dict(s)


def outdir_of(env: Env, aid: str) -> Path:
    return env.proj.results_dir / "directions" / aid


def guard_single_look(env: Env, aid: str):
    if (outdir_of(env, aid) / "summary.json").exists():
        raise registry.RegistryError(f"{aid}: summary.json exists; the discovery run is a single look")
    if env.proj.trials.exists():
        for l in env.proj.trials.read_text().splitlines():
            r = json.loads(l) if l.strip() else {}
            if r.get("analysis_id") == aid and r.get("tag") == "discovery":
                raise registry.RegistryError(f"{aid}: trials.jsonl already holds a discovery line (single look)")


def with_ci(res: pd.DataFrame, fmap: FrozenMap, level: float = 0.95) -> pd.DataFrame:
    out = res.copy()
    lo, hi = h1.ci(out["beta"].values.astype(float), out["se"].values.astype(float), out["df"].values.astype(float),
                   level)
    out["ci95_lo"], out["ci95_hi"] = lo, hi
    qcol = "q_adj_perm" if "q_adj_perm" in out else "q_adj"
    out["q_final"] = out[qcol]
    base = out.unit.astype(str).str.split("|").str[0]
    out["tier"] = [fmap.tier.get(u, fmap.tier.get(fmap.parent.get(u, ""), "")) for u in base]
    out["locus"] = [bool(fmap.locus.get(u, False)) for u in base]
    return out


def tsv_gz(df: pd.DataFrame, path: Path):
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, float_format="%.6g")
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
        gz.write(buf.getvalue().encode())


def kill_verdict(env: Env, direction: str):
    p = env.proj.results_dir / "controls" / "kill_rule.json"
    if not p.exists():
        return None
    return json.loads(p.read_text()).get(direction, {}).get("verdict")


def write_common(env: Env, aid: str, entry: dict, blob: str, res: pd.DataFrame, extra: dict, t0: float) -> dict:
    """tests.tsv.gz, discoveries.tsv, summary.json (counts, runtime, peak RSS, the trial line)."""
    r = with_ci(res, env.fmap)
    ev = r[r.evaluable]
    disc = r[r.discovery].assign(_abs_t=lambda d: -d.t.abs()).sort_values(["p_final", "_abs_t"]).drop(columns="_abs_t")
    cols = ["context", "unit", "n", "n_context", "beta", "se", "df", "ci95_lo", "ci95_hi", "t", "p", "p_perm",
            "p_final", "q_final", "tier", "locus", "covariates_dropped"]
    cols = [c for c in cols if c in disc.columns]
    trial = json.loads(env.proj.trials.read_text().splitlines()[-1]) if env.proj.trials.exists() else {}
    trial = trial if trial.get("analysis_id") == aid else {}
    summ = dict(analysis_id=aid, direction=entry["direction"], registry_blob=blob, head=trial.get("head"),
                standing=entry.get("standing", "confirmatory"), headline=entry.get("headline", True),
                kill_rule_verdict=kill_verdict(env, entry["direction"]), fdr=entry["fdr_family"],
                m_evaluable=int(len(ev)), n_tests=int(len(r)), n_not_evaluable=int((~r.evaluable).sum()),
                n_families=int(ev.context.nunique()), registered_m=entry.get("family_m", {}).get("value"),
                n_p_lt_005=int((ev.p < 0.05).sum()), n_reject_analytic=int(r.reject.sum()),
                n_screened_perm=int(r.p_perm.notna().sum()), n_discovery=int(r.discovery.sum()),
                n_discovery_nonlocus=int((r.discovery & ~r.locus).sum()),
                discoveries_per_context={k: int(v) for k, v in r[r.discovery].context.value_counts().items()},
                trial_line=trial, runtime_s=round(time.time() - t0, 1),
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    summ.update(extra)
    if env.write:
        od = outdir_of(env, aid)
        od.mkdir(parents=True, exist_ok=True)
        tsv_gz(r, od / "tests.tsv.gz")
        disc[cols].to_csv(od / "discoveries.tsv", sep="\t", index=False, float_format="%.6g")
        (od / "summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=str) + "\n")
    return summ


# ------------------------------------------------------------------------------------------------ calibration
def pseudo_calibration(env: Env, entry: dict, ctxs: list, Y: dict, covs: dict, seed: int) -> dict:
    """Registry `calibration`: N_PSEUDO pseudo-contexts per direction, pseudo-context j matched (size per lineage
    stratum, PatientID blocks; lineage contexts: random sets of the same size) to context j mod M in registry order,
    seed + context index; each tested against every unit with the same HC3 machinery (robust.pseudo_context_
    calibration). Pooled fraction p < 0.05 vs expected 0.05; response fires iff > CAL_THRESHOLD. BB rejections over
    the pseudo tests (families = pseudo-contexts) are reported (expected ~0)."""
    items = {c["id"]: c for c in entry["context_family"]["contexts"]}
    M = len(ctxs)
    rows, allp, allfam = [], [], []
    for k, c in enumerate(ctxs):
        n = len(range(k, N_PSEUDO, M))
        if n == 0:
            continue
        lin = cscore.is_lineage_context(c, items.get(c))
        alt = {"+": "less", "-": "greater", "two-sided": "two-sided"}[registry.contexts_of(entry)[c]["sign"]]
        st = None if lin else env.strata
        r = robust.pseudo_context_calibration(Y[c], env.ctx[c], covs[c], n=n, seed=seed + k, strata=st,
                                              blocks=env.blocks, alternative=alt, threshold=CAL_THRESHOLD)
        # same draws, p-values kept for the BB report (identical to the engine's loop; checked below)
        pcs = robust.pseudo_contexts(env.ctx[c], n, seed + k, strata=st, blocks=env.blocks)
        ps_c = []
        for j in pcs:
            X, _, ok = assoc._design(pcs[j].reindex(Y[c].index), covs[c].reindex(Y[c].index))
            p = h1.fit_columns(X, Y[c].values[ok], alt).p.values
            ps_c.append(p)
            allp.append(p)
            allfam.append(np.full(len(p), f"{c}|{j}", dtype=object))
        pc = np.concatenate(ps_c)
        pc = pc[np.isfinite(pc)]
        if len(pc) != r["n_tests"] or not np.isclose(np.mean(pc < 0.05), r["frac_p_lt_005"]):
            raise RuntimeError(f"pseudo-context replicate differs from the engine for {c}")
        rows.append(dict(context=c, lineage_context=lin, n_pseudo=n, seed=seed + k, n_tests=r["n_tests"],
                         frac_p_lt_005=r["frac_p_lt_005"], frac_p_lt_001=r["frac_p_lt_001"]))
    p = np.concatenate(allp)
    fam = np.concatenate(allfam)
    fin = np.isfinite(p)
    bb = stats.benjamini_bogomolov(p[fin], fam[fin], float(entry["fdr_family"]["q"]))
    frac = float(np.mean(p[fin] < 0.05))
    return dict(n_pseudo=int(sum(r["n_pseudo"] for r in rows)), n_tests=int(fin.sum()), frac_p_lt_005=frac,
                expected=0.05, frac_p_lt_001=float(np.mean(p[fin] < 0.01)), threshold=CAL_THRESHOLD,
                response_fires=bool(frac > CAL_THRESHOLD), bb_rejections=int(bb["reject"].sum()),
                bb_families_selected=int(bb["R"]), per_context=rows,
                not_run=["random gene sets matched for size, essentiality and coherence",
                         "dispersed never-expressed sets",
                         "segment-matched null for CN contexts (H3b)"],
                not_run_reason="no engine implementation; reported-only negative controls that do not enter "
                               "discovery (the pre-declared response depends on the pseudo-context fraction only)")


# ------------------------------------------------------------------------------------------------ H3a / H3b
def run_h3(env: Env, aid: str, n_perm: int = N_PERM) -> dict:
    t0 = time.time()
    entry, blob = entry_of(env, aid)
    if env.strict:
        guard_single_look(env, aid)
    ctxs = primary_contexts(entry)
    units = env.fmap.coarse()
    if registry.set_hash(units) != entry["complex_set_hash"]:
        raise registry.RegistryError("coarse units do not match the registered complex_set_hash")
    cache, Y, rem_info = {}, {}, {}
    for c in ctxs:
        rem = cscore.removal_for_context(c, env.homdel.get(c))
        if rem not in cache:
            cache[rem] = cscore.unit_scores(env.E_A, env.fmap, units, rem, env.annot)
        Y[c] = cache[rem][0]
        info = cache[rem][1]
        rem_info[c] = dict(n_units_with_removal=int((info.n_removed > 0).sum()),
                           n_units_not_evaluable=int((~info.evaluable).sum()))
    covs, dropped = covmod.per_context(entry, ctxs, env.cov)
    sd = seeds_of(env, entry)
    cal = pseudo_calibration(env, entry, ctxs, Y, covs, int(sd["pseudo_contexts"]))
    res = assoc.run_association(aid, Y, env.ctx[ctxs], covs, inputs=env.inputs, proj=env.proj, blocks=env.blocks,
                                strata=env.strata, n_perm=n_perm, perm_screen_q=SCREEN_Q,
                                perm_screen_p=CAL_RAW_P if cal["response_fires"] else None,
                                perm_seed=int(sd["permutation"]), min_context_n=MIN_N, min_ref_n=MIN_N,
                                release_hits=True, tag="discovery", save=env.write)
    extra = dict(n_perm=n_perm, seeds=sd, calibration=dict((k, v) for k, v in cal.items() if k != "per_context"),
                 covariates_dropped_by_rule={c: v for c, v in dropped.items() if v}, removals=rem_info,
                 covariates_dropped_by_design={c: sorted(set(";".join(res[res.context == c].covariates_dropped)
                                                             .split(";")) - {""}) for c in ctxs})
    if env.write:
        od = outdir_of(env, aid)
        od.mkdir(parents=True, exist_ok=True)
        (od / "calibration.json").write_text(json.dumps(cal, indent=1, sort_keys=True, default=float) + "\n")
    summ = write_common(env, aid, entry, blob, res, extra, t0)
    return dict(res=res, Y=Y, covs=covs, summary=summ, calibration=cal)


# ------------------------------------------------------------------------------------------------ H1 / H1S
def load_events(aid: str, root: Path) -> list:
    ev = pd.read_csv(root / H1_EVENTS, sep="\t")
    sp = pd.read_csv(root / H1_SPLIT, sep="\t")
    keep = set(sp.event[sp.entry == aid])
    split = lambda s: [x for x in str(s).split(";") if x and x != "nan"]  # noqa: E731
    out = []
    for r in ev.to_dict("records"):
        if r["event"] not in keep:
            continue
        out.append(h1.Event(index=int(r["index"]), id=r["event"], unit=r["unit"], gene=r["gene"], link=r["link"],
                            type=r["type"], P=split(r["P"]), N=split(r["N"]), n_A=int(r["n_A"]),
                            ref_A=int(r["ref_A"]), n_L=int(r["n_L"]), ref_L=int(r["ref_L"]),
                            testable=bool(r["testable"]), locus=bool(r["locus"]), covariates=split(r["covariates"])))
    if len(out) != len(keep):
        raise RuntimeError(f"{aid}: {len(keep)} events in the split, {len(out)} in h1_events.tsv")
    return out


def load_loss_A(root: Path, A) -> pd.DataFrame:
    lt = pd.read_parquet(root / H1_LOSS).reindex(A)
    lt.columns = [c.split(":", 1)[1] for c in lt.columns]
    return lt


def load_pair_scores(root: Path) -> pd.DataFrame:
    """Coherence for the H1 matched null: zavg_BEF on A (FZ-C S x S block), keyed by DepMap column keys."""
    from part2 import data, scores
    blk = scores.load_block(root / ZAVG_A_BLOCK)
    if blk.meta["model_set_hash"] != A_HASH:
        raise RuntimeError("zavg block is not the A block")
    keys = data.protein_genes().depmap_column.values[blk.s_index]
    return pd.DataFrame(blk.vals, index=keys, columns=keys)


def run_h1(env: Env, aid: str, events: list, loss_A: pd.DataFrame, pair_scores: pd.DataFrame,
           n_perm: int = N_PERM, K: int = K_NULL) -> dict:
    t0 = time.time()
    entry, blob = entry_of(env, aid)
    if env.strict:
        guard_single_look(env, aid)
    Y1, ctx1, pairs, comps = h1.assoc_inputs(env.E_A, events, loss_A)
    if set(ctx1.columns) != set(primary_contexts(entry)):
        raise registry.RegistryError("H1 loss columns differ from the registered contexts")
    items = {c["id"]: c for c in entry["context_family"]["contexts"]}
    for e in events:                                       # frozen event covariates == registry context covariates
        if list(items[e.context].get("covariates", entry["covariates"])) != list(e.covariates):
            raise registry.RegistryError(f"{e.id}: event covariates differ from the registry")
    cov1, dropped = covmod.per_context(entry, list(ctx1.columns), env.cov)
    sd = seeds_of(env, entry)
    res = assoc.run_association(aid, Y1, ctx1, cov1, inputs=env.inputs, proj=env.proj, blocks=env.blocks,
                                strata=env.strata, n_perm=n_perm, perm_screen_q=SCREEN_Q,
                                perm_seed=int(sd["permutation"]), min_context_n=MIN_N, min_ref_n=MIN_N, pairs=pairs,
                                release_hits=True, tag="discovery", save=env.write)
    t1 = time.time()
    ess = env.E_A.clip(-cscore.CLIP, cscore.CLIP).mean()
    cls = h1.classify(events, comps, res, env.E_A, loss_A, cov1, env.fmap, env.annot, ess, pair_scores,
                      q=float(entry["fdr_family"]["q"]), K=K, seed0=int(sd["matched_null"]))
    cls.insert(1, "event_index", [e.index for e in events])
    calls = dict(complex_wide=cls.class_complex_wide.value_counts().to_dict(),
                 paralog=cls.class_paralog.value_counts().to_dict())
    extra = dict(n_perm=n_perm, K_null=K, seeds=sd, n_events=len(events), class_counts=calls,
                 classify_runtime_s=round(time.time() - t1, 1),
                 covariates_dropped_by_rule={c: v for c, v in dropped.items() if v},
                 calibration=dict(note="H1/H1S calibration is the per-event matched random-set null (classes.tsv: "
                                       "N_/P_null_q025, q975, status, sets); the entry has no pseudo-context field"))
    if env.write:
        od = outdir_of(env, aid)
        od.mkdir(parents=True, exist_ok=True)
        cls.to_csv(od / "classes.tsv", sep="\t", index=False, float_format="%.6g")
    summ = write_common(env, aid, entry, blob, res, extra, t0)
    return dict(res=res, classes=cls, comps=comps, covs=cov1, summary=summ)


def run_sensitivity_cn(env: Env, aid: str, gene: str, events: list, loss_A: pd.DataFrame,
                       pair_scores: pd.DataFrame, K: int = K_NULL, label: str = "SMARCB1_CN") -> pd.DataFrame:
    """D7 pre-declared sensitivity (not a trial): the event(s) of `gene` re-classified with the CN-type covariate
    set of the registry (h1.COVARIATES['CN']); q / discovery flags from the committed primary run."""
    entry, _ = entry_of(env, aid)
    sd = seeds_of(env, entry)
    od = outdir_of(env, aid)
    sp = od / f"sensitivity_{label}.tsv"
    if env.strict and sp.exists():
        raise registry.RegistryError(f"{sp} exists (single look)")
    res = pd.read_parquet(env.proj.cache_dir / "assoc" / aid / "discovery.parquet")
    summ = json.loads((env.proj.results_dir / "assoc" / aid / "discovery_summary.json").read_text())
    if assoc._sha_file(env.proj.cache_dir / "assoc" / aid / "discovery.parquet") != summ["full_table_sha256"]:
        raise RuntimeError("primary H1 table differs from its committed summary")
    evs = [e for e in events if e.gene == gene]
    rows = []
    for setname, names in (("primary (registered)", None), ("CN-covariate sensitivity", h1.COVARIATES["CN"])):
        covs = {}
        for e in evs:
            nm = names if names is not None else covmod.covariates_for(entry, e.context)[0]
            bad = covmod.context_functions(e.context)
            covs[e.context] = covmod.assemble([n for n in nm if n not in bad], env.cov)
        comps = {e.id: h1.component_scores(env.E_A, e, loss_A[e.gene]) for e in evs}
        cls = h1.classify(evs, comps, res, env.E_A, loss_A, covs, env.fmap, env.annot,
                          env.E_A.clip(-cscore.CLIP, cscore.CLIP).mean(), pair_scores,
                          q=float(entry["fdr_family"]["q"]), K=K, seed0=int(sd["matched_null"]))
        for e, r in zip(evs, cls.to_dict("records")):
            Z = covs[e.context]
            ok = Z.notna().all(axis=1) & loss_A[e.gene].reindex(Z.index).notna()
            nm = names if names is not None else covmod.covariates_for(entry, e.context)[0]
            r.update(covariate_set=setname, covariates=";".join(nm), n_complete=int(ok.sum()),
                     n_loss_complete=int((ok & (loss_A[e.gene].reindex(Z.index) == 1)).sum()),
                     note="sensitivity analysis (D7), not a trial; q and discovery flags from the primary run")
            rows.append(r)
    out = pd.DataFrame(rows)
    if env.write:
        od.mkdir(parents=True, exist_ok=True)
        out.to_csv(sp, sep="\t", index=False, float_format="%.6g")
    return out


# ------------------------------------------------------------------------------------------------ H2
def committed_discovery(env: Env, aid: str) -> pd.DataFrame:
    p = env.proj.cache_dir / "assoc" / aid / "discovery.parquet"
    summ_p = env.proj.results_dir / "assoc" / aid / "discovery_summary.json"
    registry.committed_blob(env.proj.root, summ_p)
    if assoc._sha_file(p) != json.loads(summ_p.read_text())["full_table_sha256"]:
        raise RuntimeError(f"{aid} discovery table differs from its committed summary")
    return pd.read_parquet(p)


def run_h2(env: Env, aid: str = "H2_submodule_divergence", n_perm: int = N_PERM, parents: dict | None = None) -> dict:
    t0 = time.time()
    entry, blob = entry_of(env, aid)
    if env.strict:
        guard_single_look(env, aid)
    ctxs = primary_contexts(entry)
    subs = h2.h2_units(env.fmap)
    if registry.set_hash(subs) != entry["complex_set_hash"]:
        raise registry.RegistryError("H2 sub-modules do not match the registered complex_set_hash")
    if parents is None:
        parents = {"H3a": committed_discovery(env, "H3a_lineage_atlas"),
                   "H3b": committed_discovery(env, "H3b_driver_family")}
    sel = assoc.parent_selection_from(parents, {s: env.fmap.parent[s] for s in subs}, q=float(entry["fdr_family"]["q"]))
    removal_of = lambda c: cscore.removal_for_context(c, env.homdel.get(c))  # noqa: E731
    Y2, coefs = h2.h2_frames(env.E_A, env.fmap, ctxs, removal_of, env.annot, units=subs)
    cov2, dropped = covmod.per_context(entry, ctxs, env.cov)
    sd = seeds_of(env, entry)
    res = assoc.run_association(aid, Y2, env.ctx[ctxs], cov2, inputs=env.inputs, proj=env.proj, blocks=env.blocks,
                                strata=env.strata, n_perm=n_perm, perm_screen_q=SCREEN_Q,
                                perm_seed=int(sd["permutation"]), min_context_n=MIN_N, min_ref_n=MIN_N,
                                parent_selection=sel, release_hits=True, tag="discovery", save=env.write)
    fam = sorted({(c, p) for c, p in sel["selected"]})
    extra = dict(n_perm=n_perm, seeds=sd, n_eligible_families=len(fam),
                 eligible_families=[f"{c}|{p}" for c, p in fam],
                 levels={c: v for c, v in sel["level"].items() if any(c == f[0] for f in fam)},
                 covariates_dropped_by_rule={c: v for c, v in dropped.items() if v},
                 calibration=dict(not_run=["random member splits of the same parents (seed 4800)",
                                           "H3 pseudo-contexts through the stage-wise machinery"],
                                  not_run_reason="no engine implementation; reported-only negative controls that "
                                                 "do not enter discovery"))
    summ = write_common(env, aid, entry, blob, res, extra, t0)
    return dict(res=res, Y=Y2, covs=cov2, summary=summ, selection=sel)


# ------------------------------------------------------------------------------------------------ real env
def real_env(dry: bool = False) -> Env:
    proj = assoc.Project()
    root = proj.root
    fmap = FrozenMap.from_membership(pd.read_csv(root / MAP, sep="\t"))
    ctx_full = pd.read_parquet(root / CTX)
    A = assoc.a_models(proj)
    if dry:
        side = json.loads(next((root / "data_cache" / "part2").glob("CRISPRGeneEffect.*.f32.json")).read_text())
        cols = side["columns"]
        rng = np.random.default_rng(0)
        Aset = set(A)
        ids = [m for m in side["models"] if m in Aset]
        E_A = pd.DataFrame((rng.standard_normal((len(ids), len(cols))) - 0.5).astype(np.float32), index=ids,
                           columns=cols)   # N(-0.5, 1): H1 baselines |b_j| >= 0.1, so N axes are exercised
    else:
        E_A = assoc.load_A_effects(proj=proj)
    if registry.set_hash(E_A.index) != A_HASH or not set(E_A.index) <= set(A):
        raise assoc.FirewallError("effect matrix is not the A model set")
    Aidx = E_A.index
    cov = covmod.source_table(ctx_full, Aidx)
    hd = pd.read_csv(root / HOMDEL, sep="\t")
    return Env(proj=proj, fmap=fmap, ctx=ctx_full.reindex(Aidx), E_A=E_A, cov=cov,
               annot=cscore.load_annotation(E_A.columns), inputs=[CTX, MAP, HOMDEL], blocks=cov.PatientID,
               strata=cov.OncotreeLineage, homdel=dict(zip(hd.context, hd.homdel_share)), strict=not dry,
               write=not dry)


def dry_check(aid: str) -> dict:
    """Manifests and hashes, placeholders, A hash, synthetic end-to-end pass and a runtime estimate."""
    t0 = time.time()
    env = real_env(dry=True)
    entry, blob = entry_of(env, aid)
    txt = (env.proj.registry_dir / f"{aid}.yaml").read_text()
    out = dict(analysis_id=aid, registry_blob=blob, status=entry["status"], pending=entry.get("pending"),
               placeholder_strings=[p for p in registry.PLACEHOLDERS if f": {p}" in txt or f"- {p}" in txt],
               A_models=len(env.E_A), A_hash=registry.set_hash(env.E_A.index),
               seeds=entry.get("seeds"), q=entry["fdr_family"])
    inputs = list(env.inputs) + ([H1_LOSS, H1_EVENTS, H1_SPLIT] if aid.startswith("H1") else [])
    out["frozen_inputs"] = assoc.check_frozen(env.proj, inputs, ["results/part2/lockbox/manifest.json"]
                                              + list(entry["frozen_inputs"]))
    blk = env.proj.root / ZAVG_A_BLOCK
    out["zavg_A_block_sha256"] = assoc._sha_file(blk) if aid.startswith("H1") else None
    missing = [c for c in primary_contexts(entry) if not c.startswith("H1loss:") and c not in env.ctx.columns]
    out["missing_context_columns"] = missing
    env.inputs = inputs
    calls = []
    orig_append, orig_save = assoc._append_jsonl, assoc.save_outputs
    assoc._append_jsonl = lambda path, rec: calls.append(rec)
    assoc.save_outputs = lambda *a, **k: {}
    try:
        if aid in ("H3a_lineage_atlas", "H3b_driver_family"):
            r = run_h3(env, aid, n_perm=N_PERM)
            out["calibration_synthetic"] = {k: r["calibration"][k] for k in ("n_pseudo", "n_tests", "frac_p_lt_005",
                                                                              "response_fires", "bb_rejections")}
        elif aid.startswith("H1"):
            events = load_events(aid, env.proj.root)
            loss_A = load_loss_A(env.proj.root, env.E_A.index)
            genes = sorted(env.fmap.all_genes(1))
            rng = np.random.default_rng(1)
            ps = pd.DataFrame(rng.standard_normal((len(genes), len(genes))), index=genes, columns=genes)
            r = run_h1(env, aid, events, loss_A, ps, n_perm=N_PERM, K=K_NULL)
            out["n_events"] = len(events)
        else:
            r = run_h2(env, aid, n_perm=N_PERM)
            out["n_eligible_families"] = r["summary"]["n_eligible_families"]
        s = r["summary"]
        out.update(m_evaluable=s["m_evaluable"], n_tests=s["n_tests"], registered_m=s["registered_m"],
                   n_families=s["n_families"], n_screened_synthetic=s["n_screened_perm"],
                   n_discovery_synthetic=s["n_discovery"], trial_line_preview={k: calls[-1][k] for k in (
                       "models_hash", "n_models", "n_units", "n_tests", "n_evaluable", "inputs")} if calls else None)
        # runtime estimate: one Freedman-Lane permutation (10,000) on the real design of the first context
        res = r["res"]
        ev = res[res.evaluable]
        c0 = ev.context.iloc[0]
        Z = r["covs"][c0]
        X, _, ok = assoc._design(env.ctx[c0].reindex(env.E_A.index) if not c0.startswith("H1loss:") else
                                 load_loss_A(env.proj.root, env.E_A.index)[c0.split(":", 1)[1]],
                                 Z.reindex(env.E_A.index))
        ids = env.E_A.index[ok]
        y = np.random.default_rng(2).standard_normal(len(ids))
        tp = time.time()
        assoc.permutation_p(y, X, ids, "two-sided", N_PERM, np.random.default_rng(3), env.blocks,
                            env.strata.reindex(ids).fillna("NA").astype(str).values)
        out["perm_s_per_hit_within_lineage"] = round(time.time() - tp, 2)
        tp = time.time()
        assoc.permutation_p(y, X, ids, "two-sided", N_PERM, np.random.default_rng(3), env.blocks, None)
        out["perm_s_per_hit_unrestricted"] = round(time.time() - tp, 2)
        out["dry_runtime_s"] = round(time.time() - t0, 1)
        out["peak_rss_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    finally:
        assoc._append_jsonl, assoc.save_outputs = orig_append, orig_save
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", required=True, choices=ORDER)
    ap.add_argument("--dry-check", action="store_true")
    ap.add_argument("--sensitivity", choices=["SMARCB1_CN"])
    a = ap.parse_args()
    if a.dry_check:
        print(json.dumps(dry_check(a.id), indent=1, default=str))
        return
    env = real_env()
    t0 = time.time()
    if a.sensitivity:
        if a.id != "H1_within_complex_GI":
            raise SystemExit("SMARCB1_CN is an H1 sensitivity")
        env.inputs = env.inputs + [H1_LOSS, H1_EVENTS, H1_SPLIT]
        assoc.check_frozen(env.proj, env.inputs, ["results/part2/lockbox/manifest.json"] +
                           list(entry_of(env, a.id)[0]["frozen_inputs"]))
        out = run_sensitivity_cn(env, a.id, "SMARCB1", load_events(a.id, env.proj.root),
                                 load_loss_A(env.proj.root, env.E_A.index), load_pair_scores(env.proj.root))
        print(out.to_string())
    elif a.id in ("H3a_lineage_atlas", "H3b_driver_family"):
        print(json.dumps(run_h3(env, a.id)["summary"], indent=1, default=str))
    elif a.id.startswith("H1"):
        env.inputs = env.inputs + [H1_LOSS, H1_EVENTS, H1_SPLIT]
        r = run_h1(env, a.id, load_events(a.id, env.proj.root), load_loss_A(env.proj.root, env.E_A.index),
                   load_pair_scores(env.proj.root))
        print(json.dumps(r["summary"], indent=1, default=str))
    else:
        print(json.dumps(run_h2(env, a.id)["summary"], indent=1, default=str))
    print(f"wall {time.time() - t0:.1f} s; peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
