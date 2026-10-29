"""Registered per-hit robustness and the paper's strongest-findings rule (DESIGN v3.1 §2.4-2.5; registry
`success_criterion` and `multiverse.robustness` of H3a_lineage_atlas, H3b_driver_family, H1_within_complex_GI).

Rule (fixed in advance; DESIGN §2.5 line "Paper's strongest findings"):
    discovery q < 0.05 AND lockbox-confirmed (or the rare-context route, labelled) AND LOSO-robust (same sign,
    p < 0.05) AND not a locus cluster; ranked by the A statistic x tier.
Applied with the R7 reporting rule (PI, 2026-10-29): the PRIMARY lockbox result decides; a defect-sensitive call
(primary and corrected relook differ) cannot support a claim. H1: only non-neutral class calls are eligible, and the
E2 amendment bars promotion of DETECTABLE-WITHIN-MARGIN events (`headline_eligible`).

LOSO unit, as registered:
    H3a   leave-one-OncotreePrimaryDisease-out  (H3a_lineage_atlas.yaml `multiverse.robustness`)
    H3b   leave-one-lineage-out (OncotreeLineage) (H3b_driver_family.yaml: "leave-one-lineage-out (LOSO)")
    H1    leave-one-lineage-out (OncotreeLineage) (registry abbreviation as expanded in every entry that spells it
          out, incl. the H1S A5 amendment; same unit as the H1 rare route already run in A)
The fold and leave-top-2-out mechanics are the pre-registered engine functions robust.loso / robust.leave_top2_out
(every evaluable fold same sign as the full A fit and p < 0.05 under the discovery test's alternative; folds leaving
< 2 context or < 2 reference models, or a singular design, are skipped and counted). Leave-top-2-out and breadth are
reported per hit; only LOSO enters the rule (the registry's `use`: "reported per hit, never selected"; the
success_criterion names LOSO only).

Ranking: the registry gives no numeric tier weights, so "A statistic x tier" is applied as tier class first (T1, then
T2, then T3 = T3-neg / T3-nc, i.e. "complex" before "co-dependency module", DESIGN §2.3), then |t_A| descending.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2 import robust

LOSO_UNIT = {"H3a": "OncotreePrimaryDisease", "H3b": "OncotreeLineage", "H1": "OncotreeLineage"}
PASS_CALLS = ("confirmed", "not lockbox-testable: rare-context route passed")
TIER_CLASS = {"T1": 1, "T2": 2, "T3-neg": 3, "T3-nc": 3}


def _bool(s: pd.Series) -> pd.Series:
    return s.map(lambda v: v if isinstance(v, (bool, np.bool_)) else str(v).strip() == "True").astype(bool)


def h3_candidates(disc: pd.DataFrame, calls: pd.DataFrame, q: float = 0.05) -> pd.DataFrame:
    """Every H3 discovery with the rule's non-robustness criteria. disc: discoveries.tsv (context, unit, beta, t,
    q_final, tier, locus); calls: primary_relabelled_R7 rows of the entry (context, unit, primary_call,
    defect_sensitive)."""
    x = disc.merge(calls[["context", "unit", "primary_call", "defect_sensitive"]], on=["context", "unit"],
                   how="left", validate="one_to_one")
    x["locus"] = _bool(x.locus)
    x["defect_sensitive"] = _bool(x.defect_sensitive.fillna(False))
    x["discovery_q_lt_q"] = x.q_final.astype(float) < q
    x["lockbox_pass"] = x.primary_call.isin(PASS_CALLS)
    x["route"] = np.where(x.primary_call.eq(PASS_CALLS[1]), "rare-context route (not lockbox-testable)", "lockbox")
    x["candidate"] = x.discovery_q_lt_q & x.lockbox_pass & ~x.defect_sensitive & ~x.locus
    return x


def h1_candidates(disc: pd.DataFrame, calls: pd.DataFrame, events: pd.DataFrame, headline: pd.DataFrame,
                  q: float = 0.05) -> pd.DataFrame:
    """The H1 eligible (non-neutral) calls with the rule's non-robustness criteria. disc: H1 discoveries.tsv
    (unit = event|axis); calls: h1_events_primary_vs_corrected_R7.tsv; events: h1_events.tsv (locus);
    headline: confirmation table rows (unit, headline_eligible) carrying the E2 label unchanged."""
    c = calls[_bool(calls.eligible)].copy()
    c["unit"] = c.event + "|" + c.axis
    c["context"] = "H1loss:" + c.event.str.split("|").str[1]
    x = c.merge(disc[["unit", "beta", "t", "q_final", "tier"]], on="unit", how="left", validate="one_to_one")
    x = x.merge(events[["event", "locus"]], on="event", how="left", validate="many_to_one")
    x = x.merge(headline[["unit", "headline_eligible"]], on="unit", how="left", validate="one_to_one")
    x["locus"] = _bool(x.locus)
    x["defect_sensitive"] = _bool(x.defect_sensitive)
    x["headline_eligible"] = _bool(x.headline_eligible.fillna(False))
    x["discovery_q_lt_q"] = x.q_final.astype(float) < q
    x["lockbox_pass"] = x.primary_call.isin(PASS_CALLS)
    x["route"] = np.where(x.primary_call.eq(PASS_CALLS[1]), "rare-context route (not lockbox-testable)", "lockbox")
    x["candidate"] = (x.discovery_q_lt_q & x.lockbox_pass & ~x.defect_sensitive & ~x.locus & x.headline_eligible)
    return x


def hit_robustness(y: pd.Series, ctx: pd.Series, cov: pd.DataFrame | None, groups: pd.Series,
                   alternative: str = "two-sided") -> dict:
    """Full A fit, registered LOSO and leave-top-2-out for one hit (engine functions, A models only)."""
    full = robust.fit_one(y, ctx, cov, alternative)
    out = dict(beta_A_refit=full["beta"], t_A_refit=full["t"], p_A_refit=full["p"], n_A=full["n"],
               n_context_A=full["n_context"])
    out.update(robust.loso(y, ctx, cov, groups, alternative))
    out.update(robust.leave_top2_out(y, ctx, cov, alternative))
    out["top2_same_sign"] = bool(np.isfinite(out["top2_beta"]) and np.sign(out["top2_beta"]) == np.sign(full["beta"]))
    return out


def apply_rule(x: pd.DataFrame) -> pd.DataFrame:
    """strongest = candidate AND loso_robust; ranked by tier class, then |t_A| (descending). Adds tier_class, rank."""
    x = x.copy()
    x["loso_robust"] = _bool(x.loso_robust.fillna(False))
    x["strongest"] = _bool(x.candidate) & x.loso_robust
    x["tier_class"] = x.tier.map(TIER_CLASS).fillna(9).astype(int)
    x["abs_t_A"] = x.t.astype(float).abs()
    x = x.sort_values(["strongest", "tier_class", "abs_t_A", "context", "unit"],
                      ascending=[False, True, False, True, True], kind="mergesort").reset_index(drop=True)
    x["rank"] = np.where(x.strongest, np.cumsum(x.strongest), 0)
    return x


# ------------------------------------------------------------------------------------------------ post hoc annotations
# Descriptive only (PI ruling after the robustness run): computed from committed A-side files, never used as filters.

# Control context -> hit contexts it overlaps by definition (same defining alteration, or a superset / subset of the
# same lines). Exact equality is handled separately. Written from the registry context definitions.
CONTROL_CONTEXT_OVERLAP = {
    "heme": {"lineage:Lymphoid", "lineage:Myeloid"},
    "FOCAD_del_or_MSI_H": {"MSI_H", "H1loss:FOCAD"},
    "RAS_hotspot": {"hotspot:KRAS", "hotspot:NRAS"},
    "BRAF_vs_RAS_hotspot": {"hotspot:BRAF", "hotspot:KRAS", "hotspot:NRAS"},
    "melanoma": {"lineage:Skin"},
    "PTEN_LoF_noPIK3CAhs": {"LoF:PTEN", "H1loss:PTEN"},
    "PIK3CA_hotspot_noPTENLoF": {"hotspot:PIK3CA"},
    "female_or_LoY": {"LoY", "sex_female"},
    "LoF:ARID1A": {"H1loss:ARID1A"},
    "LoF:SMARCA4": {"H1loss:SMARCA4"},
    "LoF:STAG2": {"H1loss:STAG2"},
    "LoF:RB1": {"H1loss:RB1"},
    "VPS4B_loss": {"H1loss:VPS4B"},
    "SMARCB1_LoF_or_SS18SSX": {"H1loss:SMARCB1", "fusion:SS18::SSX"},
    "9p21_loss": {"H1loss:CDKN2A", "H1loss:MTAP"},
    "TP53_altered": {"H1loss:TP53"},
    "NRF2_active": {"H1loss:KEAP1"},
}


def control_recovery(context: str, unit: str, controls: pd.DataFrame) -> str:
    """'readout (PCx)' if (context, complex unit) is a pre-registered complex-level control readout, counting the
    control context's by-definition equivalents (CONTROL_CONTEXT_OVERLAP); 'same unit, other control context (PCx:
    ctx)' if the unit is a control's complex readout in a non-overlapping context; '' otherwise. controls: rows of
    controls_summary.tsv with level == 'C' (design_id, context, unit)."""
    base = unit.split("|")[0]
    exact, other = [], []
    for r in controls.itertuples(index=False):
        if r.unit != base:
            continue
        if r.context == context or context in CONTROL_CONTEXT_OVERLAP.get(r.context, set()):
            exact.append(r.design_id)
        else:
            other.append(f"{r.design_id}: {r.context}")
    if exact:
        return f"readout ({', '.join(sorted(set(exact)))})"
    if other:
        return f"same unit, other control context ({'; '.join(sorted(set(other)))})"
    return ""


def unit_essentiality(gene_ess: pd.Series, members: dict) -> pd.Series:
    """Unit baseline essentiality as in the calibration matching (calib.matched_random_sets): mean over the unit's
    members of the gene's mean clipped A effect (calib.clip_essentiality)."""
    return pd.Series({u: float(np.nanmean(gene_ess.reindex(list(m)).values.astype(float))) if len(m) else np.nan
                      for u, m in members.items()})


def scale_effect(tests: pd.DataFrame, ess: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per context: OLS of beta on unit essentiality over ALL evaluable tested units (descriptive, not a test).
    Returns (per-context table: n_units, slope, intercept, spearman_rho, resid_sd; per-test table: context, unit,
    essentiality, fitted, residual_beta, residual_z, scale_effect_concordant = sign(beta) == sign(fitted))."""
    from scipy import stats as sst
    t = tests[_bool(tests.evaluable)][["context", "unit", "beta"]].copy()
    t["beta"] = t.beta.astype(float)
    t["essentiality"] = t.unit.map(ess).astype(float)
    t = t[np.isfinite(t.beta) & np.isfinite(t.essentiality)]
    ctx_rows, out = [], []
    for c, g in t.groupby("context", sort=True):
        x, y = g.essentiality.values, g.beta.values
        b, a = np.polyfit(x, y, 1)
        fit = a + b * x
        res = y - fit
        sd = float(np.std(res, ddof=2))
        ctx_rows.append(dict(context=c, n_units=len(g), slope=float(b), intercept=float(a),
                             spearman_rho=float(sst.spearmanr(x, y)[0]), resid_sd=sd))
        out.append(g.assign(fitted=fit, residual_beta=res, residual_z=res / sd,
                            scale_effect_concordant=np.sign(y) == np.sign(fit)))
    return pd.DataFrame(ctx_rows), pd.concat(out, ignore_index=True)
