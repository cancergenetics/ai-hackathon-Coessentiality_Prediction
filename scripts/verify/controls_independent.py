"""Independent recomputation of selected pre-registered controls (controls review).

Written by the independent reviewer. Does NOT import src/part2 (no assoc.py, no cscore.py).
Inputs read: registry definitions (by hand, encoded below), the frozen map
results/part2/fzx/zavg_engine_membership.tsv, the context table, A_models.txt, CRISPRGeneEffect.csv
(restricted to A models and the needed columns), plus the raw covariate sources required by DESIGN 2.5
(CRISPRConfounders.csv for library multi-hot and ScreenType3DO; Model.csv for ModelType) and the Ensembl-116
paralogue table (only to check which cluster members are paralogs of the context gene).

Estimand (registry): S_C = mean over the frozen cluster's members (role == member, after the registered
removals) of CRISPRGeneEffect clipped to [-3, 3] (nanmean); OLS of S_C on context (1 vs 0; NaN excluded) +
registered covariates; HC3 SE; one-sided p (alternative 'less' for design sign '+', 'greater' for '-').
Lineage FE: one dummy per OncotreeLineage level (full); a variant with levels of < 2 models merged into the
reference is also reported (columns *_merge1).
df: Bell-McCaffrey/Satterthwaite approximation for the HC3 contrast under a homoskedastic working model;
the residual-df p is printed alongside.

Usage: python scripts/verify/controls_independent.py   (run from the repo/worktree root)
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.getcwd()
GE = os.path.join(ROOT, "datasets/DepMAP/CRISPRGeneEffect.csv")
CONF = os.path.join(ROOT, "datasets/DepMAP/CRISPRConfounders.csv")
MODEL = os.path.join(ROOT, "datasets/DepMAP/Model.csv")
PARA = os.path.join(ROOT, "datasets/Paralogs/ensembl116_human_paralogues.tsv")
CTX = os.path.join(ROOT, "data_cache/part2/contexts/context_table.parquet")
MAP = os.path.join(ROOT, "results/part2/fzx/zavg_engine_membership.tsv")
AMOD = os.path.join(ROOT, "results/part2/lockbox/A_models.txt")
LBX = os.path.join(ROOT, "results/part2/lockbox/lockbox_models.txt")
SUMMARY = os.path.join(ROOT, "results/part2/controls/controls_summary.tsv")

BASE = ["library", "ScreenType3DO", "ModelType"]
MUT = ["log_burden", "MSIScore"]

# Encoded from docs/registry/CTRL_*.yaml (context id, design sign, gene units, complex unit, readout, removals,
# covariates).  'paralog_of' => member_resolved paralog component (H1 machinery).
CONTROLS = {
    "PC1": dict(context="MSI_H", sign="+", genes=["WRN (7486)"], unit="Z0157", remove=[],
                cov=["lineage_FE"] + BASE),
    "PC19a": dict(context="PTEN_LoF_noPIK3CAhs", sign="+", genes=["PIK3CB (5291)"], unit="Z0077",
                  remove=["PTEN"], remove_paralogs_of=["PTEN"], cov=["lineage_FE"] + BASE + MUT),
    "PC6": dict(context="LoF:ARID1A", sign="+", genes=["ARID1B (57492)"], unit="Z0109", remove=["ARID1A"],
                paralog_of="ARID1A", cov=["lineage_FE"] + BASE + MUT),
    "PC2": dict(context="TP53_altered", sign="+", genes=["MDM2 (4193)", "PPM1D (8493)", "USP7 (7874)"],
                unit="Z0166", remove=["TP53"], remove_paralogs_of=["TP53"], cov=["lineage_FE"] + BASE + MUT,
                # registry sign_mapping: DESIGN '+' in TP53 WT; column codes altered = 1, so registered
                # sign on TP53_altered is '-' (alternative = greater)
                alternative="greater"),
    "PC17": dict(context="melanoma", sign="+", genes=["SOX10 (6663)", "MITF (4286)"], unit="Z0018", remove=[],
                 cov=BASE),
    "PC15": dict(context="SMARCB1_LoF_or_SS18SSX", sign="+", genes=["BRD9 (65980)"], unit="Z0186",
                 remove=["SMARCB1", "SS18", "SS18L1"], remove_paralogs_of=["SMARCB1", "SS18"],
                 cov=["lineage_FE"] + BASE + MUT),
}


def sym(g):
    return g.split(" (")[0]


def ols_hc3(y, d, Z):
    """OLS y ~ 1 + d + Z.  Returns beta_d, HC3 SE, t, Satterthwaite (Bell-McCaffrey) df, residual df."""
    X = np.column_stack([np.ones(len(y)), d, Z]) if Z.shape[1] else np.column_stack([np.ones(len(y)), d])
    # drop exactly collinear nuisance columns (estimand for d is unaffected; keeps leverages exact)
    q, r = np.linalg.qr(X)
    keep = np.abs(np.diag(r)) > 1e-8 * np.abs(np.diag(r)).max()
    assert keep[1], "context column collinear with covariates"
    X = X[:, keep]
    n, p = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    H_diag = np.einsum("ij,jk,ik->i", X, XtXi, X)
    a = (XtXi @ X.T)[1]  # c'(X'X)^-1 x_i for the context coefficient
    one_m_h = 1.0 - H_diag
    sing = one_m_h < 1e-10  # leverage-1 rows (singleton FE levels): a_i == 0 and e_i == 0
    w = np.where(sing, 0.0, a / np.where(sing, 1.0, one_m_h))
    var = np.sum(w ** 2 * e ** 2)
    se = np.sqrt(var)
    t = b[1] / se
    # Bell-McCaffrey / Satterthwaite df for V = e' D e, e = M y, D = diag(w^2)
    M = np.eye(n) - X @ XtXi @ X.T
    MDM = (M * (w ** 2)[None, :]) @ M
    df_s = np.trace(MDM) ** 2 / np.sum(MDM * MDM.T)
    return b[1], se, t, df_s, n - p, n


def main():
    A = [l.strip() for l in open(AMOD) if l.strip()]
    L = set(l.strip() for l in open(LBX) if l.strip())
    assert not (set(A) & L), "A and L overlap"
    mp = pd.read_csv(MAP, sep="\t")
    para = pd.read_csv(PARA, sep="\t", usecols=["Gene name", "Human paralogue associated gene name"])
    ctx = pd.read_parquet(CTX).reindex(A)
    conf = pd.read_csv(CONF, index_col=0).reindex(A)
    model = pd.read_csv(MODEL, index_col=0, usecols=["ModelID", "ModelType"]).reindex(A)

    needed = set()
    members = {}
    for k, c in CONTROLS.items():
        m = mp[(mp.unit_id == c["unit"]) & (mp.role == "member")].gene.tolist()
        rm = set(c.get("remove", []))
        for g in c.get("remove_paralogs_of", []):
            rm |= set(para.loc[para["Gene name"] == g, "Human paralogue associated gene name"].dropna())
        m_used = [g for g in m if sym(g) not in rm]
        removed = [g for g in m if sym(g) in rm]
        if "paralog_of" in c:
            pg = set(para.loc[para["Gene name"] == c["paralog_of"], "Human paralogue associated gene name"])
            m_used = [g for g in m_used if sym(g) in pg]
        members[k] = (m_used, removed, len(m))
        needed |= set(m_used) | set(c["genes"])
    hdr = pd.read_csv(GE, nrows=0).columns
    idcol = hdr[0]
    ge = pd.read_csv(GE, usecols=[idcol] + sorted(needed), index_col=0)
    ge = ge.reindex(A)  # A models only; lockbox rows never enter
    assert set(ge.index) <= set(A)

    lib_cols = [c for c in conf.columns if c.startswith("Library")]
    summ = pd.read_csv(SUMMARY, sep="\t")
    rows = []
    for k, c in CONTROLS.items():
        alt = c.get("alternative", "less" if c["sign"] == "+" else "greater")
        cov = pd.DataFrame(index=A)
        if "lineage_FE" in c["cov"]:
            cov = cov.join(pd.get_dummies(ctx["OncotreeLineage"], prefix="lin", drop_first=True, dtype=float))
        if "library" in c["cov"]:
            cov = cov.join(conf[lib_cols])
        if "ScreenType3DO" in c["cov"]:
            cov["ScreenType3DO"] = conf["ScreenType3DO"]
        if "ModelType" in c["cov"]:
            cov = cov.join(pd.get_dummies(model["ModelType"], prefix="mt", drop_first=True, dtype=float))
        for v in ("log_burden", "MSIScore"):
            if v in c["cov"]:
                cov[v] = ctx[v]
        m_used, removed, n_mem = members[k]
        units = [(g, "G", ge[g]) for g in c["genes"]]
        sc = ge[m_used].clip(-3, 3).mean(axis=1, skipna=True)
        units.append((c["unit"], "C", sc))
        for name, lvl, y in units:
            df = pd.concat([y.rename("y"), ctx[c["context"]].rename("d"), cov], axis=1).dropna()
            Z = df.drop(columns=["y", "d"])
            beta, se, t, dfs, dfr, n = ols_hc3(df.y.values, df.d.values, Z.values)
            f = stats.t.cdf if alt == "less" else stats.t.sf
            p_s, p_r = f(t, dfs), f(t, dfr)
            # variant: indicator levels with < 2 models in the analysed sample merged into the reference
            # (singleton lineages); explains the ~1e-3 relative differences vs the producer
            small = [col for col in Z.columns if col.startswith(("lin_", "mt_")) and min(Z[col].sum(), len(Z) - Z[col].sum()) < 2]
            b2, se2, t2, dfs2, _, _ = ols_hc3(df.y.values, df.d.values, Z.drop(columns=small).values)
            p_s2 = f(t2, dfs2)
            prod = summ[(summ.design_id == k) & (summ.unit == name) & (summ.level == lvl)]
            pr = prod.iloc[0] if len(prod) else None
            rows.append(dict(control=k, unit=name, level=lvl, n=n, n_ctx=int(df.d.sum()), beta=beta, se=se, t=t,
                             df_satt=dfs, p_satt=p_s, p_resid_df=p_r, alt=alt,
                             beta_merge1=b2, se_merge1=se2, p_merge1=p_s2, n_merged_levels=len(small),
                             prod_n=None if pr is None else pr.n, prod_beta=None if pr is None else pr.beta,
                             prod_se=None if pr is None else pr.beta / pr.t,
                             prod_t=None if pr is None else pr.t, prod_p=None if pr is None else pr.p,
                             members_used=";".join(sym(g) for g in m_used) if lvl == "C" else "",
                             removed=";".join(sym(g) for g in removed) if lvl == "C" else ""))
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(out.drop(columns=["members_used", "removed"]).to_string(float_format=lambda x: f"{x:.4g}"))
    print()
    for r in out[out.level == "C"].itertuples():
        print(f"{r.control} {r.unit}: members used = {r.members_used}; removed = {r.removed or 'none'}")
    os.makedirs(os.path.join(ROOT, "results/part2/verify"), exist_ok=True)
    out.to_csv(os.path.join(ROOT, "results/part2/verify/controls_independent.tsv"), sep="\t", index=False)
    kill_rule_check(summ)
    return 0


def kill_rule_check(summ):
    """Recompute the DESIGN 2.6 kill-rule tally from the registry (direction, counted, complex_unit, sign) and the
    complex-level rows of controls_summary.tsv; compare with kill_rule.json."""
    import glob
    import json
    import re

    def field(txt, block, key):  # 'key: value' inside a top-level block (no YAML lib in the venv)
        b = re.search(rf"^{block}:\n((?:[ #].*\n|\n)*)", txt, re.M).group(1)
        return re.search(rf"^  {key}: (.*)$", b, re.M).group(1).strip().strip('"')

    tally = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "docs/registry/CTRL_*.yaml"))):
        txt = open(f).read()
        e = dict(id=re.search(r"^id: (\S+)", txt, re.M).group(1))
        kr = dict(direction=field(txt, "kill_rule", "direction"), counted=field(txt, "kill_rule", "counted") == "true")
        ctl = {k: field(txt, "control", k) for k in ("design_id", "design_sign", "complex_unit")}
        if not kr.get("counted"):
            continue
        d = kr["direction"]
        t = tally.setdefault(d, dict(controls=[], evaluable=0, recovered=0))
        t["controls"].append(ctl["design_id"])
        cu = str(ctl.get("complex_unit", ""))
        if not cu.startswith(("Z", "G")) or " " in cu:
            continue  # pre-declared 'not evaluable at complex level'
        r = summ[(summ.entry == e["id"]) & (summ.level == "C") & (summ.unit == cu)]
        assert len(r) == 1, e["id"]
        r = r.iloc[0]
        if not (r.n_context >= 12):
            continue
        t["evaluable"] += 1
        alt = "less" if ctl["design_sign"] == "+" else "greater"
        if e["id"] == "CTRL_PC2":
            alt = "greater"  # registry sign_mapping: TP53_altered coded 1, registered sign '-'
        assert r.alternative == alt, (e["id"], r.alternative, alt)
        ok_sign = r.beta < 0 if alt == "less" else r.beta > 0
        t["recovered"] += int(ok_sign and r.p < 0.05)
    kr_json = json.load(open(os.path.join(ROOT, "results/part2/controls/kill_rule.json")))
    print("\nkill rule (independent tally vs kill_rule.json):")
    for d, t in sorted(tally.items()):
        frac = t["recovered"] / t["evaluable"] if t["evaluable"] else None
        verdict = "n/a" if frac is None else ("passes" if frac >= 0.5 else "KILLED")
        j = kr_json[d]
        same = (sorted(j["controls"]) == sorted(t["controls"]) and j["evaluable"] == t["evaluable"]
                and j["recovered"] == t["recovered"])
        print(f"  {d}: counted={len(t['controls'])} evaluable={t['evaluable']} recovered={t['recovered']} "
              f"-> {verdict}; matches kill_rule.json: {same}")


if __name__ == "__main__":
    sys.exit(main())
