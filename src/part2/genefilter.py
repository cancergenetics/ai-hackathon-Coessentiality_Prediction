"""Dependency-bearing gene filter (DESIGN v2 §2.3; computed once on the construction set A, then fixed).

S = protein genes (Part 1 genes.tsv) with
  * effect < thr (-0.5) in >= min_frac (1%) of the given models (denominator: all given models; NaN never counts),
  * non-NaN in >= min_obs (95%) of the given models,
  * not Part 1 low-coverage (genes.tsv low_coverage == 0).
Optional (round-2 algorithms SHOULD 8, off by default): also admit genes with effect > +pos_thr in >= min_frac of
models (growth suppressors such as TP53 / CDKN1A).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FilterConfig:
    thr: float = -0.5
    min_frac: float = 0.01
    min_obs: float = 0.95
    exclude_low_coverage: bool = True
    include_positive: bool = False
    pos_thr: float = 0.5


def dependency_bearing(Xg: np.ndarray, low_coverage: np.ndarray | None = None, cfg: FilterConfig = FilterConfig()):
    """Xg: models x genes (already restricted to the gene universe). Returns (mask, per-gene stats dict)."""
    Xg = np.asarray(Xg, dtype=np.float64)
    n = Xg.shape[0]
    obs = (~np.isnan(Xg)).sum(0) / n
    with np.errstate(invalid="ignore"):
        neg = (Xg < cfg.thr).sum(0) / n
        pos = (Xg > cfg.pos_thr).sum(0) / n
    keep = (neg >= cfg.min_frac)
    if cfg.include_positive:
        keep |= pos >= cfg.min_frac
    keep &= obs >= cfg.min_obs
    if cfg.exclude_low_coverage and low_coverage is not None:
        keep &= ~np.asarray(low_coverage, bool)
    return keep, dict(frac_below=neg, frac_above=pos, frac_observed=obs)


def selective_mask(frac_below: np.ndarray, lo: float = 0.05, hi: float = 0.90) -> np.ndarray:
    """'Selective' genes for the split-half criterion: lo <= fraction of models below -0.5 < hi."""
    return (frac_below >= lo) & (frac_below < hi)


def dependency_bearing_for_models(model_ids, cfg: FilterConfig = FilterConfig(), return_stats: bool = False):
    """Filter over Part 1 genes.tsv for a DepMap model set (use the construction set A)."""
    from part2 import data
    ge = data.load_gene_effect()
    genes = data.protein_genes()
    Xg = ge.subset(model_ids)[:, genes.depmap_index.values]
    mask, st = dependency_bearing(Xg, genes.low_coverage.values, cfg)
    return (mask, st) if return_stats else mask
