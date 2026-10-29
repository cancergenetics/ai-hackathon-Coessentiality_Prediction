#!/usr/bin/env python
"""R-3 (referee M2(e)): DEMETER2 orthogonal check of the R-2 corrected leads (and, as a supplement, every H3a
strongest finding). POST HOC, DESCRIPTIVE, no FDR (results/part2/revision/INTENT.md).

DEMETER2 on A models only (assoc.load_A_effects(source="demeter2")); the same complex score as the discovery and the
rare-context route (cscore.unit_scores with the context's removals); the same HC3 / Satterthwaite design rules with the
context's registered covariates as available on those models. Reported per hit: n, n_context, beta_D2, t, same sign as
beta_A, one-sided p in beta_A's direction and two-sided p.
Needs results/part2/revision/r2_strongest_h3a.tsv (scripts/59). Outputs: results/part2/revision/r3_*.
Usage (PART1_DIR set, PYTHONHASHSEED=0): .venv/bin/python scripts/60_rev_demeter2.py
"""
from __future__ import annotations

import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sst

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import assoc, covariates as covmod, cscore, revision as rv, stats  # noqa: E402

OUT = ROOT / "results" / "part2" / "revision"
AID = "H3a_lineage_atlas"
MIN_N, LOW_N = 3, 12


def main():
    t0 = time.time()
    spec = importlib.util.spec_from_file_location("run_direction", ROOT / "scripts" / "44_run_direction.py")
    R = importlib.util.module_from_spec(spec)
    sys.modules["run_direction"] = R
    spec.loader.exec_module(R)
    env = R.real_env()
    entry, _ = R.entry_of(env, AID)
    hits = pd.read_csv(OUT / "r2_strongest_h3a.tsv", sep="\t")
    D2 = assoc.load_A_effects(source="demeter2", proj=env.proj)
    if not set(D2.index) <= set(env.E_A.index):
        raise assoc.FirewallError("DEMETER2 rows outside A")
    ctxs = sorted(set(hits.context))
    covs, _ = covmod.per_context(entry, ctxs, env.cov)
    rows = []
    for h in hits.itertuples(index=False):
        rem = cscore.removal_for_context(h.context, env.homdel.get(h.context))
        Y, info, _ = cscore.unit_scores(D2, env.fmap, [h.unit], rem, env.annot)
        y = Y[h.unit]
        n_mem_d2 = int(info.n_used.iloc[0]) if "n_used" in info else np.nan
        f = rv.fit_one_hc3(y, env.ctx[h.context], covs[h.context], min_n=MIN_N)
        alt = "greater" if h.beta_A > 0 else "less"
        p_dir = float(stats.pvalue_t(f["t"], f["df"], alt)) if np.isfinite(f["t"]) else np.nan
        p_two = float(stats.pvalue_t(f["t"], f["df"], "two-sided")) if np.isfinite(f["t"]) else np.nan
        rows.append(dict(rank=h.rank, context=h.context, unit=h.unit, lead_corrected=h.lead_corrected,
                         tier=h.tier, beta_A=h.beta_A, t_A=h.t_A, residual_z_nl=h.residual_z_nl,
                         n_members_used_D2=n_mem_d2, n_D2=f["n"], n_context_D2=f["n_context"],
                         n_covariates_D2=f["n_covariates"], beta_D2=f["beta"], se_D2=f["se"], t_D2=f["t"],
                         df_D2=f["df"], same_sign=bool(np.sign(f["beta"]) == np.sign(h.beta_A))
                         if np.isfinite(f["beta"]) else None, p_dir_D2=p_dir, p_two_sided_D2=p_two,
                         low_n=bool(f["n_context"] < LOW_N), members=h.members))
    t = pd.DataFrame(rows)

    def rates(d):
        ev = d[d.beta_D2.notna()]
        return dict(n_hits=int(len(d)), n_evaluable=int(len(ev)),
                    same_sign=int(ev.same_sign.astype(bool).sum()),
                    same_sign_rate=float(ev.same_sign.astype(bool).mean()) if len(ev) else None,
                    p_dir_lt_005=int((ev.p_dir_D2 < 0.05).sum()),
                    p_dir_lt_005_rate=float((ev.p_dir_D2 < 0.05).mean()) if len(ev) else None,
                    p_two_sided_lt_005_same_sign=int(((ev.p_two_sided_D2 < 0.05) & ev.same_sign.astype(bool)).sum()),
                    spearman_beta_A_vs_D2=float(sst.spearmanr(ev.beta_A, ev.beta_D2)[0]) if len(ev) > 2 else None,
                    median_n_context_D2=float(ev.n_context_D2.median()) if len(ev) else None,
                    n_low_n=int(ev.low_n.sum()))
    leads = t[t.lead_corrected]
    summ = dict(label="POST HOC, DESCRIPTIVE (revision R-3); A models only; no FDR",
                n_A_models_with_D2=int(len(D2)), n_D2_genes=int(D2.shape[1]),
                leads=rates(leads), leads_not_low_n=rates(leads[~leads.low_n]), all_207=rates(t),
                non_leads=rates(t[~t.lead_corrected]),
                by_context_leads={c: rates(g) for c, g in leads.groupby("context")},
                runtime_s=round(time.time() - t0, 1),
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    t.sort_values(["lead_corrected", "rank"], ascending=[False, True]).to_csv(
        OUT / "r3_demeter2_per_hit.tsv", sep="\t", index=False, float_format="%.5g")
    (OUT / "r3_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if k != "by_context_leads"}, indent=1, default=float))
    print(leads.assign(az=leads.residual_z_nl.abs()).sort_values("az", ascending=False)[
        ["context", "unit", "beta_A", "beta_D2", "n_D2", "n_context_D2", "same_sign", "p_dir_D2"]].to_string(index=False))


if __name__ == "__main__":
    main()
