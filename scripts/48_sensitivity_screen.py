"""Pre-declared screen-quality sensitivity of H3a / H3b (RESEARCH_TREE D8; results/part2/directions/SENSITIVITY_SCREEN.md).

Sensitivity analysis, NOT a trial: no line is written to results/part2/trials.jsonl (checked by hash before / after),
the registered outputs are not touched, and every output is labelled `sensitivity_screen_quality`.

What changes against the registered run (scripts/44_run_direction.py run_h3): each context's covariate list is the
registered list (covariates.covariates_for, §2.2 rule) plus NNMD (-> CRISPRConfounders ScreenNNMD) and CasActivity,
appended last, assembled by the engine (covariates.assemble) from the same per-ModelID source table
(covariates.source_table). Complete case, as the engine's design builder does. Everything else is the registered
procedure, through the same engine calls and the runner's own helpers: scores (cscore.unit_scores with each context's
removals), pseudo-context calibration (run_direction.pseudo_calibration, registry seed, 100 pseudo-contexts, the
pre-declared response raw p < 0.01 iff the fraction exceeds 0.075), assoc.run_association (HC3 + Satterthwaite,
>= 12 / >= 12, BB families = contexts at q = 0.05, screen q < 0.1, Freedman-Lane 10,000 permutations, PatientID blocks,
registry permutation seed + row, GPD tail only below 1e-4, BB re-run).

The engine's run record (which run_association appends to `Project.trials`) goes to the sensitivity folder
(`engine_record.jsonl`) through a Project whose `trials` path points there; the engine itself is unchanged.

Reports, against the committed primary run (data_cache/part2/assoc/<ID>/discovery.parquet, SHA-256 checked against
the committed discovery_summary.json): discoveries surviving / lost / new per context, the pseudo-context rate before
and after, complete-case dropout per lineage and per context, and the Spearman correlation between unit essentiality
(mean S_C over A) and each H3a lineage beta (the review's lymphoid "scale effect"), before and after.

Usage (project root, PART1_DIR set, PYTHONHASHSEED=0):
    .venv/bin/python scripts/48_sensitivity_screen.py --id H3a_lineage_atlas|H3b_driver_family [--dry-check]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as ss

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, covariates as covmod, cscore, registry  # noqa: E402

spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
R = importlib.util.module_from_spec(spec)
sys.modules["run_direction"] = R
spec.loader.exec_module(R)

LABEL = "sensitivity_screen_quality"
EXTRA = ("NNMD", "CasActivity")             # registry names; NNMD -> ScreenNNMD (covariates.ALIASES)
IDS = ("H3a_lineage_atlas", "H3b_driver_family")


@dataclass
class SensitivityProject(assoc.Project):
    """The real project, except that the engine's run record goes to `record` instead of the trial ledger."""
    record: Path | None = None

    @property
    def trials(self):
        if self.record is None:
            raise RuntimeError("SensitivityProject needs a record path")
        return Path(self.record)


def _sha(p: Path) -> str | None:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def sensitivity_covariates(entry: dict, ctxs, table: pd.DataFrame, extra=EXTRA) -> tuple[dict, dict]:
    """{context: registered covariates + extra (appended last)} and {context: names}."""
    covs, names = {}, {}
    for c in ctxs:
        nm, _ = covmod.covariates_for(entry, c)
        nm = list(nm) + [e for e in extra if e not in nm]
        covs[c] = covmod.assemble(nm, table)
        names[c] = nm
    return covs, names


def dropout(ctx: pd.DataFrame, cov_prim: dict, cov_sens: dict, lineage: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Complete-case counts per context (context / reference models, registered vs sensitivity design) and the models
    lost by adding the covariates, per OncotreeLineage over all A models."""
    rows, lost_any = [], pd.Series(False, index=ctx.index)
    for c in cov_sens:
        x = ctx[c]
        okp = x.notna() & cov_prim[c].reindex(ctx.index).notna().all(axis=1)
        oks = x.notna() & cov_sens[c].reindex(ctx.index).notna().all(axis=1)
        lost = okp & ~oks
        lost_any |= lost
        rows.append(dict(context=c, n_context_primary=int((okp & (x == 1)).sum()),
                         n_ref_primary=int((okp & (x == 0)).sum()), n_context_sens=int((oks & (x == 1)).sum()),
                         n_ref_sens=int((oks & (x == 0)).sum()), dropped_context=int((lost & (x == 1)).sum()),
                         dropped_ref=int((lost & (x == 0)).sum())))
    per_ctx = pd.DataFrame(rows)
    lin = lineage.reindex(ctx.index).fillna("NA")
    per_lin = pd.DataFrame(dict(n_A=lin.value_counts(), dropped=lin[lost_any].value_counts())).fillna(0).astype(int)
    per_lin["frac_dropped"] = (per_lin.dropped / per_lin.n_A).round(4)
    per_lin = per_lin.sort_values("n_A", ascending=False).rename_axis("OncotreeLineage").reset_index()
    return per_ctx, per_lin


def compare(prim: pd.DataFrame, sens: pd.DataFrame) -> pd.DataFrame:
    """Per context: evaluable m, discoveries primary / sensitivity, surviving (both), lost, new."""
    key = ["context", "unit"]
    a = prim[key + ["evaluable", "discovery"]].rename(columns={"evaluable": "ev_p", "discovery": "d_p"})
    b = sens[key + ["evaluable", "discovery"]].rename(columns={"evaluable": "ev_s", "discovery": "d_s"})
    m = a.merge(b, on=key, how="outer")
    for c in ("ev_p", "d_p", "ev_s", "d_s"):
        m[c] = m[c].fillna(False).astype(bool)
    g = m.groupby("context")
    out = pd.DataFrame(dict(m_primary=g.ev_p.sum(), m_sens=g.ev_s.sum(), disc_primary=g.d_p.sum(),
                            disc_sens=g.d_s.sum(), surviving=g.apply(lambda d: int((d.d_p & d.d_s).sum())),
                            lost=g.apply(lambda d: int((d.d_p & ~d.d_s).sum())),
                            new=g.apply(lambda d: int((~d.d_p & d.d_s).sum())))).astype(int)
    out["surviving_frac"] = (out.surviving / out.disc_primary.where(out.disc_primary > 0)).round(4)
    return out.rename_axis("context").reset_index()


def scale_rho(res: pd.DataFrame, ess: pd.Series) -> pd.DataFrame:
    """Spearman rho between unit essentiality (mean S_C over A) and beta, per context, over its evaluable tests."""
    rows = []
    for c, g in res[res.evaluable].groupby("context"):
        r = ss.spearmanr(g.beta.values, ess.reindex(g.unit).values, nan_policy="omit")
        rows.append(dict(context=c, n=int(len(g)), rho=float(r.statistic), p=float(r.pvalue)))
    return pd.DataFrame(rows, columns=["context", "n", "rho", "p"])


def run(env, aid: str, primary: pd.DataFrame, primary_cal: dict, outdir: Path, n_perm: int = R.N_PERM) -> dict:
    t0 = time.time()
    if aid not in IDS:
        raise ValueError(f"{aid}: the screen-quality sensitivity is H3a / H3b only")
    if env.write and (outdir / "summary.json").exists():
        raise registry.RegistryError(f"{outdir}/summary.json exists: the sensitivity runs once")
    entry, blob = R.entry_of(env, aid)
    ctxs = R.primary_contexts(entry)
    units = env.fmap.coarse()
    if registry.set_hash(units) != entry["complex_set_hash"]:
        raise registry.RegistryError("coarse units do not match the registered complex_set_hash")
    cache, Y = {}, {}
    for c in ctxs:
        rem = cscore.removal_for_context(c, env.homdel.get(c))
        if rem not in cache:
            cache[rem] = cscore.unit_scores(env.E_A, env.fmap, units, rem, env.annot)
        Y[c] = cache[rem][0]
    cov_prim, _ = covmod.per_context(entry, ctxs, env.cov)
    covs, names = sensitivity_covariates(entry, ctxs, env.cov)
    per_ctx_drop, per_lin_drop = dropout(env.ctx[ctxs], cov_prim, covs, env.strata)
    sd = R.seeds_of(env, entry)
    cal = R.pseudo_calibration(env, entry, ctxs, Y, covs, int(sd["pseudo_contexts"]))
    rec = outdir / "engine_record.jsonl" if env.write else Path(tempfile.mkdtemp()) / "engine_record.jsonl"
    proj_s = SensitivityProject(root=env.proj.root, outcome_loader=env.proj.outcome_loader, record=rec)
    trials_before = _sha(env.proj.trials)
    res = assoc.run_association(aid, Y, env.ctx[ctxs], covs, inputs=env.inputs, proj=proj_s, blocks=env.blocks,
                                strata=env.strata, n_perm=n_perm, perm_screen_q=R.SCREEN_Q,
                                perm_screen_p=R.CAL_RAW_P if cal["response_fires"] else None,
                                perm_seed=int(sd["permutation"]), min_context_n=R.MIN_N, min_ref_n=R.MIN_N,
                                release_hits=True, tag=LABEL, save=env.write)
    if _sha(env.proj.trials) != trials_before:
        raise RuntimeError("the trial ledger changed during a sensitivity run")
    r = R.with_ci(res, env.fmap)
    cmp_ = compare(primary, r)
    lin_ctx = [c for c in ctxs if cscore.is_lineage_context(c)]
    ess = next(iter(cache.values()))[0].mean() if aid.startswith("H3a") else None
    rho = None
    if ess is not None:
        rb, ra = scale_rho(primary[primary.context.isin(lin_ctx)], ess), scale_rho(r[r.context.isin(lin_ctx)], ess)
        rho = rb.merge(ra, on="context", suffixes=("_primary", "_sens"))
    ev = r[r.evaluable]
    disc_p = primary[primary.discovery]
    both = disc_p.merge(r[r.discovery][["context", "unit"]], on=["context", "unit"])
    summ = dict(analysis_id=aid, label=LABEL, standing="sensitivity analysis (pre-declared, D8); not a trial",
                registry_blob=blob, head=R.assoc._head(env.proj), covariates=names, added=list(EXTRA),
                confounder_columns={"NNMD": "ScreenNNMD", "CasActivity": "CasActivity"},
                n_perm=n_perm, seeds=sd, m_evaluable=int(len(ev)), m_evaluable_primary=int(primary.evaluable.sum()),
                n_families=int(ev.context.nunique()), n_screened_perm=int(r.p_perm.notna().sum()),
                n_discovery=int(r.discovery.sum()), n_discovery_nonlocus=int((r.discovery & ~r.locus).sum()),
                n_discovery_primary=int(len(disc_p)), n_surviving=int(len(both)),
                n_lost=int(len(disc_p) - len(both)), n_new=int(r.discovery.sum() - len(both)),
                surviving_frac=round(len(both) / len(disc_p), 4) if len(disc_p) else None,
                pseudo_context=dict(primary=primary_cal.get("frac_p_lt_005"), sensitivity=cal["frac_p_lt_005"],
                                    change=(cal["frac_p_lt_005"] - primary_cal["frac_p_lt_005"])
                                    if primary_cal.get("frac_p_lt_005") is not None else None,
                                    primary_bb_rejections=primary_cal.get("bb_rejections"),
                                    sensitivity_bb_rejections=cal["bb_rejections"],
                                    response_fires=cal["response_fires"], threshold=R.CAL_THRESHOLD),
                dropout_complete_case=dict(models_dropped_any_context=int(per_lin_drop.dropped.sum()),
                                           per_lineage=per_lin_drop.to_dict("records")),
                trial_ledger_sha256_unchanged=trials_before,
                runtime_s=round(time.time() - t0, 1),
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    if rho is not None:
        lym = rho[rho.context == "lineage:Lymphoid"]
        summ["scale_effect"] = dict(definition="Spearman rho between mean S_C over A (unit essentiality) and the "
                                               "lineage beta, over the context's evaluable units",
                                    lymphoid=lym.to_dict("records")[0] if len(lym) else None)
    if env.write:
        outdir.mkdir(parents=True, exist_ok=True)
        R.tsv_gz(r.assign(label=LABEL), outdir / "tests.tsv.gz")
        cols = ["context", "unit", "n", "n_context", "beta", "se", "df", "ci95_lo", "ci95_hi", "t", "p", "p_perm",
                "p_final", "q_final", "tier", "locus", "covariates_dropped"]
        d = r[r.discovery].assign(_a=lambda x: -x.t.abs()).sort_values(["p_final", "_a"])
        d[cols].to_csv(outdir / "discoveries.tsv", sep="\t", index=False, float_format="%.6g")
        cmp_.to_csv(outdir / "comparison_per_context.tsv", sep="\t", index=False)
        per_ctx_drop.to_csv(outdir / "dropout_per_context.tsv", sep="\t", index=False)
        per_lin_drop.to_csv(outdir / "dropout_per_lineage.tsv", sep="\t", index=False)
        if rho is not None:
            rho.to_csv(outdir / "scale_effect_rho.tsv", sep="\t", index=False, float_format="%.4g")
        (outdir / "calibration.json").write_text(json.dumps(dict(cal, label=LABEL), indent=1, sort_keys=True,
                                                            default=float) + "\n")
        (outdir / "summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=str) + "\n")
    return dict(res=r, summary=summ, comparison=cmp_, dropout_context=per_ctx_drop, dropout_lineage=per_lin_drop,
                rho=rho, calibration=cal)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", required=True, choices=IDS)
    ap.add_argument("--dry-check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    env = R.real_env(dry=a.dry_check)
    entry, _ = R.entry_of(env, a.id)
    assoc.check_frozen(env.proj, env.inputs, ["results/part2/lockbox/manifest.json"] + list(entry["frozen_inputs"]))
    od = env.proj.results_dir / "directions" / a.id / LABEL
    if a.dry_check:   # synthetic outcome: covariates, dropout and m only; no output, no primary table read
        prim = pd.DataFrame(columns=["context", "unit", "evaluable", "discovery", "beta"])
        out = run(env, a.id, prim.astype({"evaluable": bool, "discovery": bool}), {}, od)
        s = out["summary"]
        print(json.dumps({k: s[k] for k in ("m_evaluable", "n_families", "covariates", "dropout_complete_case",
                                            "runtime_s", "peak_rss_bytes")}, indent=1, default=str))
        print(out["dropout_context"].to_string())
        return
    prim = R.with_ci(R.committed_discovery(env, a.id), env.fmap)
    pcal = json.loads((env.proj.results_dir / "directions" / a.id / "calibration.json").read_text())
    out = run(env, a.id, prim, pcal, od)
    s = out["summary"]
    print(json.dumps({k: v for k, v in s.items() if k not in ("dropout_complete_case", "covariates")}, indent=1,
                     default=str))
    print(out["comparison"].to_string())
    if out["rho"] is not None:
        print(out["rho"].to_string())
    print(f"wall {time.time() - t0:.1f} s; peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
