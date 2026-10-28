"""H2 sub-module divergence (DESIGN v3.1 §2.6; docs/registry/H2_submodule_divergence.yaml).

Score: for sub-module s of parent P and context c (after c's removals, applied to s and to its sibling complement),
    S_s, S_sib = primary mean scores (clip +/-3, 50% rule) of s and of P minus s
    Y_s = S_s - (a_s + b_s S_sib),  (a_s, b_s) by OLS over all A models (complete-case)
The A-fitted (a_s, b_s) are applied unchanged to lockbox rows (``apply``).

Testing is stage-wise inside parents discovered by the parent directions (H3a for lineage contexts, H3b for driver
contexts): assoc.parent_selection_from(...) gives the eligible (c, P) families and the level q R_x / m_x;
assoc.run_association with fdr method stagewise_bb tests only those sub-modules.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from part2.cmap import FrozenMap
from part2.cscore import Annotation, Removal, mean_score


def eligible_parents(fmap: FrozenMap, min_children: int = 2) -> list:
    return [p for p in fmap.coarse() if len(fmap.submodules(p)) >= min_children]


def h2_units(fmap: FrozenMap, min_children: int = 2) -> list:
    """Sub-modules of parents with >= 2 sub-modules (the registered unit set)."""
    return sorted(s for p in eligible_parents(fmap, min_children) for s in fmap.submodules(p))


def residualised(E: pd.DataFrame, fmap: FrozenMap, sub: str, removal: Removal | None = None,
                 annot: Annotation | None = None, coef: tuple | None = None, min_members: int = 2):
    """(Y_s over E's rows, (a, b), raw (S_s, S_sib)). coef=None fits (a, b) on E's rows (A)."""
    rm = set(removal.removed(fmap.members[fmap.parent[sub]], annot or Annotation())) if removal else set()
    s_mem = [g for g in fmap.members[sub] if g not in rm and g in E.columns]
    sib = [g for g in fmap.sibling_complement(sub) if g not in rm and g in E.columns]
    if len(s_mem) < min_members or len(sib) < min_members:
        return pd.Series(np.nan, index=E.index), None, None
    Ss, Sb = mean_score(E, s_mem), mean_score(E, sib)
    if coef is None:
        ok = Ss.notna() & Sb.notna()
        b, a = np.polyfit(Sb[ok].values, Ss[ok].values, 1)
        coef = (float(a), float(b))
    a, b = coef
    return Ss - (a + b * Sb), coef, (Ss, Sb)


def h2_frames(E_A: pd.DataFrame, fmap: FrozenMap, contexts, removal_of, annot: Annotation, units=None):
    """{context: A x units residualised scores} and {(context, unit): (a, b)}; removal_of(context) -> Removal."""
    units = h2_units(fmap) if units is None else units
    Y, coefs = {}, {}
    for c in contexts:
        rem = removal_of(c)
        cols = {}
        for s in units:
            y, cf, _ = residualised(E_A, fmap, s, rem, annot)
            cols[s] = y
            coefs[(c, s)] = cf
        Y[c] = pd.DataFrame(cols, index=E_A.index)
    return Y, coefs


def score_L_fn(fmap: FrozenMap, coefs: dict, removal_of, annot: Annotation):
    def f(E_L, context, unit):
        cf = coefs.get((context, unit))
        if cf is None:                                  # not evaluable in A: never refit on lockbox rows
            return pd.Series(np.nan, index=E_L.index)
        y, _, _ = residualised(E_L, fmap, unit, removal_of(context), annot, coef=cf)
        return y
    return f
