"""FZ-LB lockbox draw (DESIGN v3.1 §2.1): L = 25% of U, seed 3000, PatientID-blocked, multilabel-stratified.

Labels (model labels only; no outcome is read): OncotreeLineage one-hot (missing lineage -> 'lineage:NA') plus the
genotype list fixed in DESIGN §2.1 (``GENOTYPE_LABELS``; primary context definitions of contexts.py, NaN counted as
label-absent).

Algorithm: iterative stratification (Sechidis, Tsoumakas & Vlahavas 2011, ECML-PKDD) on patient *blocks*, two folds
with ratios (L, A) = (0.25, 0.75):
  1. Blocks = PatientID groups (a model without PatientID is its own block). A block's label vector is the per-label
     count of its models; its size is its model count.
  2. Desired counts: per label l with n_l positive models, the L fold wants t_l and A wants n_l - t_l, where
     t_l = round_half_up(0.25 n_l), except for genotype labels with 17 <= n_l < 20, for which t_l = 5 (">= 5 L models
     where n_U permits": n_U >= 17 is the smallest size that also leaves >= 12 models in A, the confirmatory size).
     Total desired: L = round_half_up(0.25 |U|), A = |U| - that.
  3. Repeat until every block is assigned: take the label with the fewest positive models among unassigned blocks
     (>0; ties -> label order after a seeded shuffle); visit the unassigned blocks carrying it in seeded random order;
     put each into the fold with the largest remaining desired count for that label (ties -> larger remaining total
     desired size -> seeded coin); subtract the block's label counts and size from that fold's desired counts.
     Every model carries a lineage label, so all blocks are assigned.
The draw is a deterministic function of (U, labels, PatientID, seed): blocks and labels are sorted before the
seeded permutations, and only numpy's default_rng(seed) is used.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# DESIGN §2.1 genotype list -> context_table columns (primary definitions)
GENOTYPE_LABELS = {
    "MSI-H": "MSI_H",
    "LoF ARID1A": "LoF:ARID1A", "LoF SMARCA4": "LoF:SMARCA4", "LoF STAG2": "LoF:STAG2", "LoF SMARCB1": "LoF:SMARCB1",
    "LoF PTEN": "LoF:PTEN", "LoF RB1": "LoF:RB1", "LoF KMT2D": "LoF:KMT2D", "LoF CREBBP": "LoF:CREBBP",
    "LoF EP300": "LoF:EP300",
    "NRF2-active": "NRF2_active",
    "KMT2A-r": "fusion:KMT2A-r",
    "9p21/MTAP hom-del": "9p21_loss",
    "VPS4B loss": "VPS4B_loss",
    "KRAS hotspot": "hotspot:KRAS", "NRAS hotspot": "hotspot:NRAS", "BRAF hotspot": "hotspot:BRAF",
    "PIK3CA hotspot": "hotspot:PIK3CA",
}
FRAC_L = 0.25
SEED = 3000
MIN_L = 5
MIN_A_CONFIRM = 12


def round_half_up(x: float) -> int:
    return int(np.floor(x + 0.5))


def target_L(n: int, genotype: bool, frac: float = FRAC_L) -> int:
    t = round_half_up(frac * n)
    if genotype and t < MIN_L and n >= MIN_L + MIN_A_CONFIRM:
        t = MIN_L
    return t


def label_matrix(ctx: pd.DataFrame) -> pd.DataFrame:
    """Models x labels (0/1 int) from a context table (needs OncotreeLineage + GENOTYPE_LABELS columns)."""
    lin = ctx["OncotreeLineage"].fillna("NA").astype(str)
    lab = pd.get_dummies(lin, prefix="lineage", prefix_sep=":").astype(int)
    for name, col in GENOTYPE_LABELS.items():
        lab[name] = (ctx[col].values == 1).astype(int)
    return lab


def iterative_stratification(labels: pd.DataFrame, blocks: pd.Series, frac: float = FRAC_L, seed: int = SEED,
                             genotype_cols=()) -> pd.Series:
    """Returns a Series ModelID -> 'L' / 'A'. See module docstring."""
    labels = labels.sort_index()
    blocks = blocks.reindex(labels.index)
    blocks = blocks.where(blocks.notna(), pd.Series("__model__" + labels.index, index=labels.index))
    rng = np.random.default_rng(seed)
    bl_names = np.array(sorted(blocks.unique()))
    B = labels.groupby(blocks).sum().reindex(bl_names)           # blocks x labels counts
    size = blocks.value_counts().reindex(bl_names).values
    C = B.values.astype(np.int64)
    lab_names = list(labels.columns)
    gset = set(genotype_cols)
    n_l = C.sum(0)
    tL = np.array([target_L(int(n), name in gset, frac) for n, name in zip(n_l, lab_names)])
    desired = np.vstack([tL, n_l - tL]).astype(float)              # fold 0 = L, fold 1 = A
    nL = round_half_up(frac * len(labels))
    desired_tot = np.array([nL, len(labels) - nL], float)
    lab_order = rng.permutation(len(lab_names))                    # tie-break order of labels
    rank = np.empty_like(lab_order)
    rank[lab_order] = np.arange(len(lab_order))
    assigned = np.full(len(bl_names), -1)
    while (assigned < 0).any():
        un = assigned < 0
        rem = C[un].sum(0)
        cand = np.flatnonzero(rem > 0)
        if len(cand) == 0:                                         # blocks without labels (not expected)
            l = None
            todo = np.flatnonzero(un)
        else:
            m = rem[cand].min()
            ties = cand[rem[cand] == m]
            l = ties[np.argmin(rank[ties])]
            todo = np.flatnonzero(un & (C[:, l] > 0))
        todo = todo[rng.permutation(len(todo))]
        for b in todo:
            if l is not None:
                d = desired[:, l]
                best = np.flatnonzero(d == d.max())
            else:
                best = np.array([0, 1])
            if len(best) > 1:
                dt = desired_tot[best]
                best = best[dt == dt.max()]
            f = int(best[0] if len(best) == 1 else best[rng.integers(len(best))])
            assigned[b] = f
            desired[f] -= C[b]
            desired_tot[f] -= size[b]
    fold = pd.Series(np.where(assigned == 0, "L", "A"), index=bl_names)
    return pd.Series(fold.reindex(blocks.values).values, index=labels.index, name="fold")


def stratum_table(labels: pd.DataFrame, fold: pd.Series, genotype_cols=(), frac: float = FRAC_L) -> pd.DataFrame:
    gset = set(genotype_cols)
    rows = []
    for c in labels.columns:
        x = labels[c]
        n = int(x.sum())
        rows.append(dict(stratum=c, kind="genotype" if c in gset else "lineage", n_U=n,
                         n_A=int(x[fold == "A"].sum()), n_L=int(x[fold == "L"].sum()),
                         target_L=target_L(n, c in gset, frac)))
    t = pd.DataFrame(rows)
    t["dev_L"] = t.n_L - t.target_L
    return t
