"""Label-free selection of (z0, k) and stability reporting (DESIGN v2 §2.3; round-2 algorithms MUST 4, SHOULD 9).

Split-half recovery: R PatientID-blocked, lineage-stratified splits of A; each half = one full run (S x S block) and
the same consensus/persistence pipeline (5 seeds, single run). A half-1 cluster of size 5-30 (any level, not a
locus cluster when arm labels are given) is recovered if some half-2 cluster (any level) has Jaccard >= 0.5.
Criterion = fraction of selective genes (5-90% of A models below -0.5) lying in recovered clusters, averaged over
both directions and the R splits. Reported next to the same statistic on degree-rewired half graphs. Ties go to the
smaller graph (higher z0, then lower k).

Stability by size bin: per final cluster, mean best-match Jaccard against the stage-1 runs at its gamma; the chance
baseline is the same quantity for clusters built from degree-preserving rewires of the run graphs (each edge keeps
its per-run presence and weight).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd

from part2 import consensus as C
from part2.graph import build_graph, rewire_runs

SIZE_BINS = (3, 5, 10, 20, 61)


@dataclass(frozen=True)
class SplitHalfConfig:
    size_lo: int = 5
    size_hi: int = 30
    jaccard: float = 0.5
    exclude_locus: bool = True


def recovered_coverage(map1: C.ClusterMap, map2: C.ClusterMap, selective: np.ndarray, cfg=SplitHalfConfig(),
                       arm=None, band=None, chrom=None, pos=None):
    """Fraction of selective genes covered by map1 clusters (size_lo..size_hi) recovered in map2."""
    from part2.locus import locus_flag
    cov = np.zeros(len(selective), bool)
    others = [c["members"] for c in map2.clusters]
    n_eval = n_rec = 0
    for c in map1.clusters:
        m = c["members"]
        if not (cfg.size_lo <= len(m) <= cfg.size_hi):
            continue
        if cfg.exclude_locus and arm is not None and locus_flag(m, arm, band, chrom=chrom, pos=pos)["locus"]:
            continue
        n_eval += 1
        if any(C.best_jaccard_sets(m, o) >= cfg.jaccard for o in others):
            cov[m] = True
            n_rec += 1
    return float((cov & selective).sum() / max(selective.sum(), 1)), n_eval, n_rec


def split_half_statistic(half_blocks: list[tuple[np.ndarray, np.ndarray]], z0: float, k: int, selective,
                         ccfg: C.ClusterConfig = C.ClusterConfig(), scfg=SplitHalfConfig(), cap: float = 3.0,
                         arm=None, band=None, rewire_seed: int | None = 4000, chrom=None, pos=None) -> dict:
    """half_blocks: [(block_half1, block_half2), ...] S x S score blocks per split (same S)."""
    rows = []
    for si, (b1, b2) in enumerate(half_blocks):
        g1, g2 = build_graph(b1, k, z0, cap), build_graph(b2, k, z0, cap)
        variants = [("observed", g1, g2)]
        if rewire_seed is not None:
            (r1,), (r2,) = rewire_runs([g1], rewire_seed + 2 * si), rewire_runs([g2], rewire_seed + 2 * si + 1)
            variants.append(("rewired", r1, r2))
        for name, h1, h2 in variants:
            m1 = C.build_map([h1], ccfg, b1 if name == "observed" else None, z0)
            m2 = C.build_map([h2], ccfg, b2 if name == "observed" else None, z0)
            a, ne1, nr1 = recovered_coverage(m1, m2, selective, scfg, arm, band, chrom, pos)
            b, ne2, nr2 = recovered_coverage(m2, m1, selective, scfg, arm, band, chrom, pos)
            nl = dataclasses.replace(scfg, exclude_locus=False)
            a2 = recovered_coverage(m1, m2, selective, nl)[0]
            b2 = recovered_coverage(m2, m1, selective, nl)[0]
            rows.append(dict(split=si, variant=name, cov_1to2=a, cov_2to1=b, stat=(a + b) / 2,
                             stat_locus_included=(a2 + b2) / 2,
                             n_eval=ne1 + ne2, n_recovered=nr1 + nr2, edges_h1=h1.m, edges_h2=h2.m,
                             n_clusters_h1=len(m1.clusters), n_clusters_h2=len(m2.clusters)))
    df = pd.DataFrame(rows)
    obs = df[df.variant == "observed"]
    rew = df[df.variant == "rewired"]
    return dict(z0=z0, k=k, stat=float(obs.stat.mean()), stat_sd=float(obs.stat.std(ddof=0)),
                stat_locus_included=float(obs.stat_locus_included.mean()),
                rewired_stat=float(rew.stat.mean()) if len(rew) else np.nan, per_split=df)


def select_z0k(results: list[dict]) -> dict:
    """Highest split-half statistic; ties (equal to 1e-12) to the smaller graph: higher z0, then lower k."""
    return sorted(results, key=lambda d: (-round(d["stat"], 12), -d["z0"], d["k"]))[0]


def stability_by_size(cmap: C.ClusterMap, levels: list[C.Level], bins=SIZE_BINS) -> pd.DataFrame:
    rows = []
    jidx = {lv.gamma: j for j, lv in enumerate(levels)}
    for c in cmap.clusters:
        lv = levels[jidx[c["gamma"]]]
        _, stab = C.support(c["members"], lv, None, C.ClusterConfig())
        rows.append(dict(id=c["id"], level=c["level"], size=len(c["members"]), stability=stab))
    df = pd.DataFrame(rows, columns=["id", "level", "size", "stability"])
    df["size_bin"] = pd.cut(df["size"], bins=list(bins), right=False)
    return df


def stability_report(graphs, cfg: C.ClusterConfig = C.ClusterConfig(), ref_scores=None, z0=0.0, rewire_seed=4000):
    """Observed vs degree-preserving-rewire chance baseline, mean stability per size bin."""
    lv = C.build_levels(graphs, cfg)
    cm = C.build_map(graphs, cfg, ref_scores, z0, levels=lv)
    obs = stability_by_size(cm, lv)
    rg = rewire_runs(graphs, rewire_seed)
    lvr = C.build_levels(rg, cfg)
    cmr = C.build_map(rg, cfg, None, z0, levels=lvr)
    rew = stability_by_size(cmr, lvr)
    agg = lambda d: d.groupby("size_bin", observed=False).stability.agg(["count", "mean"])
    out = agg(obs).join(agg(rew), lsuffix="_observed", rsuffix="_rewired")
    return out, cm, cmr
