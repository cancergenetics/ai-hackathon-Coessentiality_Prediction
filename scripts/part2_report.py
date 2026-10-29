"""Shared, deterministic derivations for the Part 2 discovery figures and tables (Figs 3-4, Tables S3-S5).

Reads committed result files only (results/part2/...). No dataset, data_cache, lockbox or hu.MAP file is read.
Nothing here is a new test: every quantity is a re-tabulation of committed run outputs.

  discoveries(ID)      H3a / H3b primary discoveries with
                         beats_all      |t| > max |t| of the unit's matched random sets in that context (D8 flag;
                                        reproduces 502 / 1,145 and 103 / 161 of calibration_D8.json)
                         outside_beta   post hoc 95% filter (SPECIFICITY.md; specificity_discoveries.tsv.gz)
                         cons           kept under the GPD sensitivity (discoveries_report.tsv: discovery_cons)
                         screen         kept under the screen-quality sensitivity (sensitivity_screen_quality/)
                         control        registered control readout of the same unit (controls_summary.tsv)
  never_expressed(ID)  per-context fraction p < 0.05 of the dispersed never-expressed sets
  h2_ranks()           rank of each H2 discovery's |t| among its own 20 random member splits (1 = beats all)
  h1_events()          the 24 H1 events with engine class, reported label and member composition
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import fig_style as S

DIR = "results/part2/directions"
H3A, H3B, H1, H1S, H2 = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI",
                         "H1S_silencing_GI", "H2_submodule_divergence")
MEMBERSHIP = "results/part2/fzx/zavg_engine_membership.tsv"
H1_EVENTS = "results/part2/fzx/h1_events.tsv"


def _tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(S.ROOT / rel, sep="\t")


def members() -> pd.Series:
    """unit id -> list of member gene symbols (frozen FZ-X engine membership; attached genes excluded)."""
    m = _tsv(MEMBERSHIP)
    m = m[m.role == "member"].copy()
    m["sym"] = m.gene.str.split(" ").str[0]
    return m.groupby("unit_id").sym.apply(lambda s: sorted(s))


def gene_label(unit: str, mem: pd.Series, k: int = 4) -> str:
    g = mem.get(unit, [])
    return ", ".join(g[:k]) + (f" +{len(g) - k}" if len(g) > k else "")


def short_ctx(c: str) -> str:
    return (c.replace("lineage:", "").replace("Peripheral Nervous System", "PNS")
            .replace("Ovary/Fallopian Tube", "Ovary").replace("Bladder/Urinary Tract", "Bladder")
            .replace("fusion:EWSR1::FLI1|ERG", "EWSR1 fusion").replace("Esophagus/Stomach", "Esoph./Stomach"))


def discoveries(ID: str) -> pd.DataFrame:
    d = _tsv(f"{DIR}/{ID}/discoveries.tsv")
    cal = _tsv(f"{DIR}/{ID}/calibration_D8_tests.tsv.gz")
    mrs = cal[cal.calibration == "matched_random_set"]
    mx = mrs.assign(a=mrs.t.abs()).groupby(["context", "unit"]).a.max().rename("max_abs_t_matched")
    d = d.join(mx, on=["context", "unit"])
    d["has_sets"] = d.max_abs_t_matched.notna()
    d["beats_all"] = d.t.abs() > d.max_abs_t_matched
    spec = _tsv(f"{DIR}/specificity_discoveries.tsv.gz")
    spec = spec[spec.analysis_id == ID][["context", "unit", "outside_beta"]]
    d = d.merge(spec, on=["context", "unit"], how="left")
    d["outside_beta"] = d.outside_beta.astype(object).eq(True)
    rep = _tsv(f"{DIR}/{ID}/discoveries_report.tsv")[["context", "unit", "discovery_cons", "p_perm_report"]]
    d = d.merge(rep, on=["context", "unit"], how="left").rename(columns={"discovery_cons": "cons"})
    sens = _tsv(f"{DIR}/{ID}/sensitivity_screen_quality/discoveries.tsv")[["context", "unit"]]
    sens = sens.assign(screen=True)
    d = d.merge(sens, on=["context", "unit"], how="left")
    d["screen"] = d.screen.astype(object).eq(True)
    ctl = _tsv(S.CONTROLS)
    ctl = ctl[ctl.level == "C"]
    exact = set(zip(ctl.context, ctl.unit))
    cu = set(ctl.unit)
    d["control"] = [("control readout" if (c, u) in exact else "control unit, other context" if u in cu else "")
                    for c, u in zip(d.context, d.unit)]
    d["focus"] = d.beats_all & ~d.locus
    return d


def never_expressed(ID: str) -> pd.Series:
    cal = _tsv(f"{DIR}/{ID}/calibration_D8_tests.tsv.gz")
    ne = cal[cal.calibration == "never_expressed"]
    return (ne.p < 0.05).groupby(ne.context).mean()


def matched_rates(ID: str) -> pd.DataFrame:
    j = S.load_json(f"{DIR}/{ID}/calibration_D8.json")["matched_random_sets"]
    rows = [(c, v["real"]["frac_p_lt_005"], v["random"]["frac_p_lt_005"], v["real"]["n_tests"])
            for c, v in j["per_context"].items()]
    return pd.DataFrame(rows, columns=["context", "real", "random", "n_real"]).set_index("context")


def h2_ranks() -> pd.DataFrame:
    d = _tsv(f"{DIR}/{H2}/discoveries.tsv")
    cal = _tsv(f"{DIR}/{H2}/calibration_D8_tests.tsv.gz")
    rs = cal[cal.calibration == "random_split"]
    g = {k: v.t.to_numpy() for k, v in rs.groupby(["context", "unit"])}
    ranks, nsplit, nsame = [], [], []
    for _, r in d.iterrows():
        s = g[(r.context, r.unit)]
        ranks.append(int((np.abs(s) >= abs(r.t)).sum()) + 1)
        nsplit.append(len(s))
        nsame.append(int((np.sign(s) == np.sign(r.t)).sum()))
    d["rank_vs_splits"] = ranks
    d["n_splits"] = nsplit
    d["splits_same_sign"] = nsame
    spec = _tsv(f"{DIR}/specificity_discoveries.tsv.gz")
    spec = spec[spec.analysis_id == H2][["context", "unit", "outside_beta", "outside_t"]]
    d = d.merge(spec, on=["context", "unit"], how="left")
    d["parent"] = d.unit.str.split(".").str[0]
    return d


def h1_events() -> pd.DataFrame:
    c = _tsv(f"{DIR}/{H1}/classes.tsv")
    ev = _tsv(H1_EVENTS)[["event", "P", "N"]]
    c = c.merge(ev, on="event", how="left")
    disc = _tsv(f"{DIR}/{H1}/discoveries.tsv")
    disc["event"] = disc.unit.str.rsplit("|", n=1).str[0]
    disc["axis"] = disc.unit.str.rsplit("|", n=1).str[1]
    for ax in ("N", "P"):
        sub = disc[disc.axis == ax].set_index("event")
        c[f"disc_{ax}"] = c.event.map(sub.q_final).notna()
        c[f"q_{ax}"] = c.event.map(sub.q_final)
        c[f"ci95_lo_{ax}"] = c.event.map(sub.ci95_lo)
        c[f"ci95_hi_{ax}"] = c.event.map(sub.ci95_hi)
        c[f"n_context_{ax}"] = c.event.map(sub.n_context)
    return c


def jdump(obj, rel: str) -> None:
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    (S.ROOT / rel).write_text(json.dumps(obj, indent=1, sort_keys=True, default=conv) + "\n")
