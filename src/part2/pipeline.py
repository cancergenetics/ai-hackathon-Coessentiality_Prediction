"""Glue: score blocks of a set of runs -> graphs -> consensus/persistence map (same code for toy and DepMap)."""
from __future__ import annotations

import numpy as np

from part2 import consensus as C
from part2.graph import build_graph
from part2.scores import zavg_block


def matrix_runs(X, S_mask, n_boot: int, frac: float, seed0: int, config=None, blocks=None):
    """Full-data block + n_boot row subsamples (frac of rows, by blocks if given, seeds seed0..). In-memory matrix
    version used on toy data (DepMap runs use scores.score_block with PatientID-blocked model lists)."""
    gidx = np.arange(X.shape[1])
    full = zavg_block(X, gidx, S_mask, K_all=0, config=config)[0]
    boots = []
    for b in range(n_boot):
        rng = np.random.default_rng(seed0 + b)
        rows = np.sort(rng.choice(len(X), int(round(frac * len(X))), replace=False))
        boots.append(zavg_block(np.ascontiguousarray(X[rows]), gidx, S_mask, K_all=0, config=config)[0])
    return full, boots


def build(full_block, boot_blocks, z0: float, k: int, ccfg: C.ClusterConfig = C.ClusterConfig(), cap: float = 3.0):
    graphs = [build_graph(b, k, z0, cap) for b in boot_blocks]
    levels = C.build_levels(graphs, ccfg)
    cmap = C.build_map(graphs, ccfg, full_block, z0, levels=levels)
    return cmap, graphs, levels
