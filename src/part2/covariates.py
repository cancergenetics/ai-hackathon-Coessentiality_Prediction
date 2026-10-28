"""Covariate assembly per context (DESIGN v3.1 §2.2 covariate rule, §2.5 model).

Registry covariate names -> numeric design columns, from a covariate source table (models x raw columns):

    lineage_FE      OncotreeLineage one-hot (missing lineage -> NaN row); levels with < 2 analysed models are merged
                    into the reference level by the design builder (assoc._design)
    library         Library* multi-hot columns of CRISPRConfounders
    ScreenType3DO   CRISPRConfounders ScreenType3DO
    ModelType       Model.csv ModelType one-hot (collinear columns, e.g. identical to ScreenType3DO, are dropped by
                    the design builder and recorded)
    log_burden, MSIScore, Ploidy, WGD, TP53_altered       context-table columns (numeric)
    NNMD -> ScreenNNMD, CasActivity, doubling_time -> ScreenDoublingTime, GrowthPattern -> GrowthPattern* one-hot

Covariate-of-context rule: no covariate may be a function of the context. ``covariates_for`` resolves the list of a
registry context (its own `covariates` override the entry list) and drops (and reports) any covariate that is a
function of the context: lineage_FE for lineage contexts; MSIScore / log_burden for MSI-H-based contexts; Ploidy for
WGD contexts; TP53_altered for TP53 contexts (and for H1 loss of TP53).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2.cscore import is_lineage_context

ALIASES = {"NNMD": "ScreenNNMD", "doubling_time": "ScreenDoublingTime"}
FE_NAMES = {"lineage_FE": "OncotreeLineage", "ModelType": "ModelType"}
PREFIX_NAMES = {"library": "Library", "GrowthPattern": "GrowthPattern"}


def context_functions(cid: str, item: dict | None = None) -> set:
    """Covariates that are functions of the context `cid` (DESIGN §2.2)."""
    out = set()
    if is_lineage_context(cid, item):
        out.add("lineage_FE")
    if cid in ("MSI_H", "FOCAD_del_or_MSI_H") or "MSI" in cid.replace("inclMSI", ""):
        out |= {"MSIScore", "log_burden"}
    if cid.startswith("WGD"):
        out |= {"Ploidy", "Ploidy_anyprofile", "WGD", "WGD_anyprofile"}
    if "TP53" in cid:
        out.add("TP53_altered")
    if cid.startswith("H1loss:") and cid.split(":", 1)[1] == "TP53":
        out.add("TP53_altered")
    return out


def covariates_for(entry: dict, cid: str) -> tuple[list, list]:
    """(covariate names after the rule, names dropped by the rule) for registry context `cid`."""
    items = entry["context_family"]["contexts"]
    item = next((c for c in items if c["id"] == cid), None) if isinstance(items, list) else None
    names = list((item or {}).get("covariates", entry["covariates"]))
    bad = context_functions(cid, item)
    return [n for n in names if n not in bad], [n for n in names if n in bad]


def assemble(names, table: pd.DataFrame) -> pd.DataFrame:
    """Numeric design columns (models x columns) for covariate `names` from the source `table`."""
    parts = []
    for n in names:
        if n in FE_NAMES:
            src = table[FE_NAMES[n]]
            d = pd.get_dummies(src, prefix=n, prefix_sep=":", dtype=float)
            d = d.loc[:, sorted(d.columns)]
            if d.shape[1]:
                d = d.iloc[:, 1:]                       # first level = reference
            d.loc[src.isna().values, :] = np.nan
            parts.append(d)
        elif n in PREFIX_NAMES:
            cols = sorted(c for c in table.columns if c.startswith(PREFIX_NAMES[n]))
            if not cols:
                raise KeyError(f"covariate {n}: no {PREFIX_NAMES[n]}* columns in the covariate table")
            d = table[cols].astype(float)
            if n == "GrowthPattern" and d.shape[1] > 1:
                d = d.iloc[:, 1:]
            parts.append(d.add_prefix(f"{n}:"))
        else:
            col = ALIASES.get(n, n)
            if col not in table.columns:
                raise KeyError(f"covariate {n} ({col}) not in the covariate table")
            parts.append(table[[col]].astype(float).rename(columns={col: n}))
    if not parts:
        return pd.DataFrame(index=table.index)
    return pd.concat(parts, axis=1)


def per_context(entry: dict, contexts, table: pd.DataFrame) -> tuple[dict, dict]:
    """{context: covariate DataFrame} and {context: dropped names} for run_association."""
    out, dropped = {}, {}
    for c in contexts:
        names, bad = covariates_for(entry, c)
        out[c] = assemble(names, table)
        dropped[c] = bad
    return out, dropped
