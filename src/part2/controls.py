"""Pre-registered controls evaluator (DESIGN v3.1 §2.7; docs/registry/CTRL_*.yaml). One run, after FZ-C + FZ-X.

For every CTRL entry: gene-level units (raw effect; 'A minus B' contrasts) and, when the filled
``control.complex_unit`` is a frozen unit ID, the complex-level readout of the entry's `score`:
  S_C_mean_winsor3      the unit's primary score after the control's removals
  member_resolved       paralog component (H1 machinery): mean of the unit's members that are paralogs of the lost
                        gene (PC13: X-linked gametologs of chrY genes); must contain the target, else not evaluable
  submodule_vs_sibling  the target's sub-module residualised on its sibling complement (A-fitted OLS; H2 machinery)
Each entry goes through assoc.run_association (one-sided HC3 test in the registered direction, no multiplicity
correction). Recovered (kill rule) = complex-level readout one-sided p < 0.05 in A, evaluable iff the complex unit
exists, the readout is non-empty and the context has >= 12 A models after complete-case. Kill rule per direction:
dropped iff evaluable > 0 and recovered / evaluable < 0.5; 'n/a (exploratory)' with 0 evaluable.

CONTROL_SPECS is the machine-readable form of each entry's `removed_members` / `readout_C` text (the registry text is
the specification; tests/test_part2_controls.py checks every listed gene against that text).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from part2 import assoc, covariates as covmod, registry
from part2.cmap import FrozenMap, symbol
from part2.cscore import Annotation, Removal, gene_unit, mean_score, unit_scores
from part2.h2 import residualised

KILL_THRESHOLD = 0.5
MIN_A = 12

CONTROL_SPECS = {
    "PC1": dict(removal=Removal()),
    "PC2": dict(removal=Removal(genes=("TP53",))),
    "PC4": dict(removal=Removal(genes=("KRAS", "NRAS"))),
    "PC5": dict(removal=Removal(genes=("BRAF", "KRAS", "NRAS"))),
    "PC6": dict(lost=("ARID1A",), target="ARID1B", cn_if_homdel=True),
    "PC7": dict(lost=("SMARCA4",), target="SMARCA2", cn_if_homdel=True),
    "PC8": dict(lost=("STAG2",), target="STAG1", cn_if_homdel=True),
    "PC9": dict(removal=Removal(genes=("MTAP", "CDKN2A", "CDKN2B"), cn_genes=("MTAP",), paralogs=False)),
    "PC10": dict(removal=Removal(genes=("FOCAD",), cn_genes=("FOCAD",), paralogs=False)),
    "PC11": dict(removal=Removal(genes=("KEAP1", "CUL3"), retain=("NFE2L2",))),
    "PC12": dict(lost=("VPS4B",), target="VPS4A", removal=Removal(genes=("VPS4B",), cn_genes=("VPS4B",),
                                                                   paralogs=False)),
    "PC13": dict(lost="chrY", target="DDX3X", removal=Removal(ygenes=True)),
    "PC15": dict(removal=Removal(genes=("SMARCB1", "SS18"), retain=("BRD9",))),
    "PC17": dict(removal=Removal()),
    "PC18": dict(removal=Removal()),
    "PC19a": dict(removal=Removal(genes=("PTEN",))),
    "PC19b": dict(removal=Removal(retain=("PIK3CA",))),
    "PC20": dict(removal=Removal()),
    "PC21": dict(removal=Removal()),
    "PC22": dict(removal=Removal(genes=("SMARCB1", "SS18"))),
}


def control_entries(proj) -> dict:
    out = {}
    for p in sorted(proj.registry_dir.glob("CTRL_*.yaml")):
        e, _ = registry.load_entry(proj.registry_dir, p.stem, root=proj.root)
        out[p.stem] = e
    return out


def _removal(spec: dict, homdel_share: float | None) -> Removal:
    if "removal" in spec:
        return spec["removal"]
    lost = spec["lost"]
    cn = lost if (spec.get("cn_if_homdel") and (homdel_share or 0) >= 0.5) else ()
    return Removal(genes=lost, cn_genes=cn, paralogs=False)        # paralogs are the readout, not removed


def complex_readout(entry: dict, E: pd.DataFrame, fmap: FrozenMap, annot: Annotation, homdel_share=None,
                    fit: dict | None = None) -> tuple[pd.Series, dict]:
    """Complex-level readout over E's rows. fit: A-fitted parameters (submodule_vs_sibling) for lockbox rows."""
    ctl = entry["control"]
    unit = ctl["complex_unit"]
    spec = CONTROL_SPECS[ctl["design_id"]]
    nan = pd.Series(np.nan, index=E.index)
    if unit not in fmap.members:
        return nan, dict(evaluable=False, reason=str(unit))
    rem = _removal(spec, homdel_share)
    score = entry["score"]
    if score == "S_C_mean_winsor3":
        Y, info, _ = unit_scores(E, fmap, [unit], rem, annot)
        ok = bool(info.evaluable.iloc[0])
        return Y[unit], dict(evaluable=ok, reason="" if ok else "< 2 members after removals",
                             removed=info.removed.iloc[0], n_used=int(info.n_used.iloc[0]))
    if score == "member_resolved":
        mem = [g for g in fmap.members[unit] if g not in set(rem.removed(fmap.members[unit], annot))]
        if spec["lost"] == "chrY":
            P = [g for g in mem if annot.x_gametolog(symbol(g))]
        else:
            P = [g for g in mem if any(annot.is_paralog(symbol(g), l) for l in spec["lost"])]
        if spec["target"] not in {symbol(g) for g in P}:
            return nan, dict(evaluable=False, reason=f"paralog component lacks {spec['target']}", P=P)
        return mean_score(E, P), dict(evaluable=True, reason="", P=";".join(P))
    if score == "submodule_vs_sibling":
        if fmap.level.get(unit) != 2:
            return nan, dict(evaluable=False, reason="complex unit is not a sub-module")
        y, cf, _ = residualised(E, fmap, unit, rem, annot, coef=None if fit is None else fit.get("coef"))
        if cf is None:
            return nan, dict(evaluable=False, reason="< 2 members in sub-module or sibling complement")
        return y, dict(evaluable=True, reason="", coef=cf)
    raise ValueError(f"control score {score!r} not supported")


def evaluate_controls(E_A: pd.DataFrame, fmap: FrozenMap, contexts_A: pd.DataFrame, cov_table: pd.DataFrame,
                      annot: Annotation, *, inputs=(), proj=None, homdel_share: dict | None = None,
                      write: bool = True) -> tuple[pd.DataFrame, dict]:
    """Runs every CTRL entry once. Returns (per-unit table, kill tally); writes results/part2/controls/."""
    proj = proj or assoc.Project()
    outdir = proj.results_dir / "controls"
    if write and (outdir / "kill_rule.json").exists():
        raise registry.RegistryError("controls already evaluated (results/part2/controls/kill_rule.json exists); "
                                     "the controls script runs once")
    entries = control_entries(proj)
    rows = []
    for aid, e in entries.items():
        ctl = e["control"]
        ctxs = [c["id"] for c in e["context_family"]["contexts"] if c.get("role", "primary") == "primary"]
        cid = ctxs[0]
        hs = (homdel_share or {}).get(cid)
        cols = {u: (gene_unit(E_A, u) if all(x.strip() in E_A.columns for x in u.split(" minus "))
                    else pd.Series(np.nan, index=E_A.index)) for u in ctl["gene_units"]}
        cinfo = dict(evaluable=False, reason="gene level only")
        cu = ctl.get("complex_unit")
        if isinstance(cu, str) and cu in fmap.members:
            y, cinfo = complex_readout(e, E_A, fmap, annot, hs)
            cols[cu] = y
        Y = pd.DataFrame(cols, index=E_A.index)
        covs, dropped = covmod.per_context(e, ctxs, cov_table.reindex(E_A.index))
        res = assoc.run_association(aid, Y, contexts_A[ctxs].reindex(E_A.index), covs, inputs=inputs, proj=proj,
                                    tag="controls")
        kr = e.get("kill_rule", {})
        for r in res.itertuples(index=False):
            is_c = isinstance(cu, str) and r.unit == cu
            ev = bool(is_c and cinfo["evaluable"] and r.evaluable and r.n_context >= MIN_A)
            rows.append(dict(entry=aid, design_id=ctl["design_id"], direction=kr.get("direction", "none"),
                             counted=bool(kr.get("counted", False)), context=r.context, unit=r.unit,
                             level="C" if is_c else "G", n=r.n, n_context=r.n_context, beta=r.beta, t=r.t, p=r.p,
                             alternative=r.alternative, evaluable_kill=ev, recovered=bool(ev and r.p < 0.05),
                             gene_level_p_lt_005=bool((not is_c) and r.evaluable and r.p < 0.05),
                             note=cinfo.get("reason", "") if is_c else "",
                             covariates_dropped_by_rule=";".join(dropped[r.context])))
        if not (isinstance(cu, str) and cu in fmap.members):
            rows.append(dict(entry=aid, design_id=ctl["design_id"], direction=kr.get("direction", "none"),
                             counted=bool(kr.get("counted", False)), context=cid, unit=str(cu), level="C",
                             evaluable_kill=False, recovered=False, note="not evaluable at complex level"))
    tab = pd.DataFrame(rows)
    kill = kill_tally(tab)
    if write:
        outdir.mkdir(parents=True, exist_ok=True)
        tab.to_csv(outdir / "controls_summary.tsv", sep="\t", index=False)
        (outdir / "kill_rule.json").write_text(json.dumps(kill, indent=1, sort_keys=True) + "\n")
    return tab, kill


def kill_tally(tab: pd.DataFrame) -> dict:
    """Per direction: evaluable complex-level controls, recovered, fraction, verdict."""
    out = {}
    c = tab[(tab.level == "C") & tab.counted]
    for d in ("H3a", "H3b", "H1", "H2"):
        s = c[c.direction == d]
        n_ev, n_rec = int(s.evaluable_kill.sum()), int(s.recovered.sum())
        if n_ev == 0:
            verdict = "n/a (0 evaluable controls; direction exploratory)"
        elif n_rec / n_ev < KILL_THRESHOLD:
            verdict = "killed by controls"
        else:
            verdict = "passes"
        out[d] = dict(n_controls=int(len(s)), evaluable=n_ev, recovered=n_rec,
                      fraction=(n_rec / n_ev) if n_ev else None, verdict=verdict,
                      thin_coverage=bool(0 < n_ev <= 2),
                      controls=sorted(s.design_id.unique().tolist()))
    return out
