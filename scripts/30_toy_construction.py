#!/usr/bin/env python
"""Toy validation of the Part 2 construction pipeline (known answers; no DepMap labels or contexts).

Scenarios (src/part2/toy.py): planted (complexes 3-40, nested sub-modules, shared subunit, pan-essential core,
never-essential genes with arm-block noise, lineage low-rank, library batch), null (no complexes), armblock (no
complexes; strong CN-like arm blocks among dependency-bearing genes).
Same pipeline as DepMap, small settings: locked CoessentialityModel with Config(k=10, d=64) (toy has 300 lines,
3,000 genes), gene filter, B bootstraps of 80% of lines, mutual top-k graph, consensus + persistence + top-down
nesting, cohesion, locus flag. Reports ARI / best-match F1, persistent-cluster counts, (z0, k) split-half selection
(with rewired baseline), stability by size bin vs degree-preserving rewire.
Writes results/part2/infra/toy_results.json and toy_results.md.

  .venv/bin/python scripts/30_toy_construction.py [--quick]
"""
import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from part2 import consensus as C  # noqa: E402
from part2 import genefilter, locus, paths, pipeline, selection, toy  # noqa: E402
from part2.graph import build_graph  # noqa: E402
from part2.part1 import Config  # noqa: E402
from part2.scores import zavg_block  # noqa: E402

TOY_CFG = Config(k=10, d=64, n_null=200_000)   # small settings for 300 lines x 3,000 genes (DepMap uses the locked Config())


def evaluate(T, S, cmap):
    sidx = np.flatnonzero(S)
    loc = np.full(len(S), -1)
    loc[sidx] = np.arange(len(sidx))
    arm, band = T.arm[sidx], T.band[sidx]
    flags = [locus.locus_flag(c["members"], arm, band)["locus"] for c in cmap.clusters]
    out = dict(n_level1=len(cmap.level(1)), n_level2=len(cmap.level(2)), n_locus=int(sum(flags)),
               n_nonlocus=int(len(flags) - sum(flags)), sizes=sorted(int(c["size"]) for c in cmap.clusters),
               counts={k: v for k, v in cmap.counts.items()})
    if T.complexes:
        truth = [loc[c][loc[c] >= 0] for c in T.complexes]
        subs = [loc[c][loc[c] >= 0] for c in T.submodules]
        pg = np.unique(np.concatenate(truth))
        out.update(ARI_level1=toy.adjusted_rand(truth, cmap.labels(len(sidx), 1), pg),
                   F1_level1=toy.best_match_f1(truth, [c["members"] for c in cmap.level(1)]).round(3).tolist(),
                   F1_any_level=toy.best_match_f1(truth, [c["members"] for c in cmap.clusters]).round(3).tolist(),
                   F1_submodules_any_level=toy.best_match_f1(subs, [c["members"] for c in cmap.clusters]).round(3).tolist(),
                   planted_sizes=[len(t) for t in truth],
                   frac_clustered_genes_planted=float(np.isin(np.concatenate([c["members"] for c in cmap.clusters]), pg).mean()) if cmap.clusters else None)
    return out


def main():
    C.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer bootstraps / splits (tests)")
    ap.add_argument("--out", default=str(paths.RESULTS / "infra"))
    a = ap.parse_args()
    B, R = (4, 1) if a.quick else (8, 2)
    seeds = (0, 1, 2) if a.quick else C.SEEDS
    z0, k = 5.0, 10
    res = dict(settings=dict(config=dataclasses.asdict(TOY_CFG), n_boot=B, frac=0.8, boot_seeds=f"1000..{1000 + B - 1}",
                             z0=z0, k=k, leiden_seeds=list(seeds), split_half_R=R), scenarios={})
    t0 = time.perf_counter()
    for scen in ("planted", "null", "armblock"):
        T = toy.simulate(scen, seed=0)
        S, st = genefilter.dependency_bearing(T.X)
        full, boots = pipeline.matrix_runs(T.X, S, B, 0.8, 1000, TOY_CFG)
        graphs = [build_graph(b, k, z0) for b in boots]
        sc = dict(n_S=int(S.sum()), never_in_S=int(((T.gene_class == "never") & S).sum()), edges_per_run=[g.m for g in graphs])
        for mode in ("jaccard", "containment"):
            for ng in (11, 13):
                cfg = C.ClusterConfig(seeds=seeds, support_mode=mode, gammas=C.gamma_grid(n=ng))
                cm = C.build_map(graphs, cfg, full, z0, levels=C.build_levels(graphs, cfg))
                sc[f"{mode}_grid{ng}"] = evaluate(T, S, cm)
        # selection: split-half over (z0, k) grid on PatientID-like halves (here: random halves of lines)
        sidx = np.flatnonzero(S)
        sel = genefilter.selective_mask(st["frac_below"][sidx])
        halves = []
        for r in range(R):
            perm = np.random.default_rng(3000 + r).permutation(len(T.X))
            h1, h2 = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])
            halves.append(tuple(zavg_block(np.ascontiguousarray(T.X[h]), np.arange(T.X.shape[1]), S, 0, TOY_CFG)[0]
                                for h in (h1, h2)))
        grid = []
        for zz in (4.0, 5.0, 6.0):
            for kk in (10, 20):
                cfg = C.ClusterConfig(seeds=seeds)
                rr = selection.split_half_statistic(halves, zz, kk, sel, cfg, arm=T.arm[sidx], band=T.band[sidx])
                row = dict(z0=zz, k=kk, stat=rr["stat"], rewired_stat=rr["rewired_stat"],
                           stat_locus_included=rr["stat_locus_included"])
                if T.complexes:
                    gr = [build_graph(b, kk, zz) for b in boots]
                    cm = C.build_map(gr, cfg, full, zz, levels=C.build_levels(gr, cfg))
                    row["ARI_level1_full_pipeline"] = evaluate(T, S, cm)["ARI_level1"]
                grid.append(row)
        chosen = selection.select_z0k([dict(r) for r in grid])
        sc["selection_grid"] = grid
        sc["selected"] = dict(z0=chosen["z0"], k=chosen["k"], stat=chosen["stat"])
        # dyad control: the planted complexes cut into pairs score 0 on the criterion by construction
        if T.complexes:
            loc = np.full(len(S), -1)
            loc[sidx] = np.arange(len(sidx))
            dy = C.ClusterMap(clusters=[dict(members=loc[c][loc[c] >= 0][i:i + 2], id=j, level=1)
                                        for j, c in enumerate(T.complexes) for i in range(0, len(c) - 1, 2)])
            sc["dyad_map_criterion"] = selection.recovered_coverage(dy, dy, sel)[0]
            tr = C.ClusterMap(clusters=[dict(members=loc[c][loc[c] >= 0], id=j, level=1) for j, c in enumerate(T.complexes)])
            sc["truth_map_criterion"] = selection.recovered_coverage(tr, tr, sel)[0]
        # stability by size vs degree-preserving rewire
        stab, _, _ = selection.stability_report(graphs, C.ClusterConfig(seeds=seeds), full, z0)
        sc["stability_by_size"] = {str(i): {kk: (None if not np.isfinite(v) else float(v)) for kk, v in r.items()}
                                   for i, r in stab.iterrows()}
        res["scenarios"][scen] = sc
        print(scen, json.dumps({k: v for k, v in sc.items() if k in ("n_S", "selected", "jaccard_grid11", "dyad_map_criterion")}, default=str)[:1500], flush=True)
    res["runtime_s"] = round(time.perf_counter() - t0, 1)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    name = "toy_results_quick" if a.quick else "toy_results"
    (out / f"{name}.json").write_text(json.dumps(res, indent=1, default=str))
    print("written", out / f"{name}.json", res["runtime_s"], "s")


if __name__ == "__main__":
    main()
