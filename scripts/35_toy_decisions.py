#!/usr/bin/env python
"""Toy check of the pre-freeze FZ-C algorithm decisions (results/part2/construction/DECISIONS.md).

Planted truth (src/part2/toy.py: planted / null / armblock scenarios, toy seeds 0-2), the same pipeline as DepMap at
toy scale (locked CoessentialityModel with Config(k=10, d=64); 20 row-subsamples of 80%, seeds 1000-1019; Leiden
SEEDS 0-4; z0 = 5, k = 10). No DepMap data, label or context is read.

D1 support  : 'jaccard' (round-2 wording) vs 'containment' (>= 75% of the members in one stage-1 cluster).
D2 placement: 'greedy' (coarse-first, trim, drop) vs 'nested' (persistence-first, containment = nesting, partial
              overlap > 50% resolved to the more persistent cluster).
D3 attachment: held-out planted members (one per planted complex of size >= 5, removed from S before clustering)
              must attach to the cluster of their own complex; never/pan/selective genes must not attach; random
              size-matched member sets on the null scenario must collect ~0 attachments.
Writes results/part2/construction/toy_decisions.json.

  .venv/bin/python scripts/35_toy_decisions.py [--seeds 0 1 2] [--B 20]
"""
import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from part2 import attach, genefilter, locus, paths, toy  # noqa: E402
from part2 import consensus as C  # noqa: E402
from part2.graph import build_graph  # noqa: E402
from part2.part1 import Config  # noqa: E402
from part2.scores import zavg_block  # noqa: E402

TOY_CFG = Config(k=10, d=64, n_null=200_000)
VARIANTS = {"V0_jaccard_greedy": ("jaccard", "greedy"), "V1_containment_greedy": ("containment", "greedy"),
            "V2_containment_nested": ("containment", "nested"), "V3_jaccard_nested": ("jaccard", "nested")}
Z0, K = 5.0, 10


def blocks(X, S, B, seed0=1000, frac=0.8):
    G = X.shape[1]
    allg = np.ones(G, bool)
    full = zavg_block(X, np.arange(G), allg, K_all=0, config=TOY_CFG)[0]         # all genes x all genes
    boots = []
    for b in range(B):
        rng = np.random.default_rng(seed0 + b)
        rows = np.sort(rng.choice(len(X), int(round(frac * len(X))), replace=False))
        boots.append(zavg_block(np.ascontiguousarray(X[rows]), np.arange(G), S, K_all=0, config=TOY_CFG)[0])
    return full, boots


def evaluate(T, sidx, cm):
    loc = np.full(T.X.shape[1], -1)
    loc[sidx] = np.arange(len(sidx))
    arm, band = T.arm[sidx], T.band[sidx]
    cl = cm.clusters
    flags = [locus.locus_flag(c["members"], arm, band)["locus"] for c in cl]
    out = dict(n_level1=len(cm.level(1)), n_level2=len(cm.level(2)), n_locus=int(sum(flags)),
               n_nonlocus=int(len(flags) - sum(flags)), sizes_level1=sorted(int(len(c["members"])) for c in cm.level(1)),
               counts={k: v for k, v in cm.counts.items()})
    if not T.complexes:
        out["n_clusters_total"] = len(cl)
        return out
    truth = [loc[c][loc[c] >= 0] for c in T.complexes]
    subs = [loc[c][loc[c] >= 0] for c in T.submodules]
    f1_any = toy.best_match_f1(truth, [c["members"] for c in cl])
    f1_l1 = toy.best_match_f1(truth, [c["members"] for c in cm.level(1)])
    sizes = np.array([len(t) for t in truth])
    pg = np.unique(np.concatenate(truth))
    # false clusters: < 50% of members from any single planted complex
    false = 0
    for c in cl:
        best = max(len(np.intersect1d(c["members"], t)) for t in truth) / len(c["members"])
        false += best < 0.5
    out.update(planted_sizes=sizes.tolist(), F1_any=f1_any.round(3).tolist(), F1_level1=f1_l1.round(3).tolist(),
               recovered_any_F1ge0p8=int((f1_any >= 0.8).sum()), recovered_level1_F1ge0p8=int((f1_l1 >= 0.8).sum()),
               recovered_small_3to4=int((f1_any[sizes <= 4] >= 0.8).sum()), n_small_3to4=int((sizes <= 4).sum()),
               submodules_F1_any=toy.best_match_f1(subs, [c["members"] for c in cl]).round(3).tolist(),
               nested40_parent_F1_level1=float(f1_l1[int(np.argmax(sizes))]),
               ARI_level1=toy.adjusted_rand(truth, cm.labels(len(sidx), 1), pg), false_clusters=int(false))
    if len(truth) > len(toy.SIZES):                         # 'hier' extras: 3 parents (+9 subs), 2 overlapping pairs
        n0 = len(toy.SIZES)
        out.update(hier_parents_F1_level1=f1_l1[n0:n0 + 3].round(3).tolist(),
                   hier_parents_recovered_level1=int((f1_l1[n0:n0 + 3] >= 0.8).sum()),
                   hier_subs_recovered_any=int((np.array(out["submodules_F1_any"][2:]) >= 0.8).sum()),
                   hier_overlap_F1_any=f1_any[n0 + 3:].round(3).tolist(),
                   hier_overlap_recovered_any=int((f1_any[n0 + 3:] >= 0.8).sum()))
    return out


def attachment_eval(T, S_used, held, full, cm, z0):
    """Attach every non-S gene (incl. held-out planted members) to the clusters (all levels)."""
    sidx = np.flatnonzero(S_used)
    ridx = np.flatnonzero(~S_used)
    rect = full[np.ix_(ridx, sidx)]
    cl = [c["members"] for c in cm.clusters]
    if not cl:
        return dict(n_clusters=0)
    df, summ = attach.attach(rect, cl, z0, row_ids=ridx, cluster_ids=np.arange(len(cl)))
    loc = np.full(T.X.shape[1], -1)
    loc[sidx] = np.arange(len(sidx))
    truth = [loc[c][loc[c] >= 0] for c in T.complexes]
    correct = wrong = 0
    missed = []
    for h in held:
        ts = [i for i, c in enumerate(T.complexes) if h in c]
        got = df[df.gene == h]
        ok_cl = [j for j, c in enumerate(cl) if any(len(np.intersect1d(c, truth[t])) / len(c) >= 0.5 for t in ts)]
        if got.cluster.isin(ok_cl).any():
            correct += 1
        else:
            missed.append(int(h))
        wrong += int((~got.cluster.isin(ok_cl)).sum())
    cls = T.gene_class[df.gene.values]
    heldset = set(map(int, held))
    other = [g for g in df.gene.values if int(g) not in heldset]
    oc = T.gene_class[np.array(other, int)] if other else np.array([])
    return dict(summary=summ, n_heldout=len(held), heldout_attached_correct=correct,
                heldout_wrong_cluster_pairs=wrong, heldout_missed=missed,
                other_attached_genes_by_class={k: int((oc == k).sum()) for k in np.unique(oc)} if len(oc) else {},
                n_attached_pairs=int(len(df)), attached_classes={k: int((cls == k).sum()) for k in np.unique(cls)})


def summarise(res) -> str:
    """Markdown table: per variant, sums over toy seeds."""
    rows = []
    for name in VARIANTS:
        acc = dict(rec=0, n=0, small=0, nsmall=0, false=0, sub=0, nsub=0, null=0, arm_nonlocus=0, arm=0,
                   hp=0, hs=0, ho=0)
        for r in res["runs"]:
            e = r["variants"][name]
            if r["scenario"] in ("planted", "hier"):
                acc["rec"] += e["recovered_any_F1ge0p8"]
                acc["n"] += len(e["planted_sizes"])
                acc["small"] += e["recovered_small_3to4"]
                acc["nsmall"] += e["n_small_3to4"]
                acc["false"] += e["false_clusters"]
                acc["sub"] += int((np.array(e["submodules_F1_any"]) >= 0.8).sum())
                acc["nsub"] += len(e["submodules_F1_any"])
            if r["scenario"] == "hier":
                acc["hp"] += e["hier_parents_recovered_level1"]
                acc["ho"] += e["hier_overlap_recovered_any"]
            if r["scenario"] == "null":
                acc["null"] += e["n_level1"] + e["n_level2"]
            if r["scenario"] == "armblock":
                acc["arm_nonlocus"] += e["n_nonlocus"]
                acc["arm"] += e["n_level1"] + e["n_level2"]
        rows.append(f"| {name} | {acc['rec']}/{acc['n']} | {acc['small']}/{acc['nsmall']} | {acc['sub']}/{acc['nsub']} | "
                    f"{acc['hp']}/{3 * len(res['settings']['toy_seeds'])} | {acc['ho']}/{4 * len(res['settings']['toy_seeds'])} | "
                    f"{acc['false']} | {acc['null']} | {acc['arm_nonlocus']} of {acc['arm']} |")
    head = ("| Variant | Planted complexes recovered (F1 >= 0.8, any level; planted + hier) | of which size 3-4 | "
            "Sub-modules recovered | hier parents at level 1 | hier overlapping complexes | False clusters "
            "(< 50% from one complex) | Null-scenario clusters | Arm-block non-locus clusters |\n"
            "|---|---|---|---|---|---|---|---|---|\n")
    att = []
    for r in res["runs"]:
        e = r["variants"]["V2_containment_nested"].get("attachment")
        if e and "summary" in e:
            att.append(f"| {r['seed']} | {r['scenario']} | {e['n_heldout']} | {e['heldout_attached_correct']} | "
                       f"{e['heldout_wrong_cluster_pairs']} | {e['other_attached_genes_by_class']} | "
                       f"{e['summary']['n_tests']} |")
        if "null_random_sets_attachment" in r:
            x = r["null_random_sets_attachment"]
            att.append(f"| {r['seed']} | null, 14 random size-matched sets | 0 | - | - | {x['n_attached']} attached | {x['n_tests']} |")
    ahead = ("| Toy seed | Scenario | Held-out members | Attached to own complex | Wrong-cluster pairs | Other genes "
             "attached (by class) | Tests |\n|---|---|---|---|---|---|---|\n")
    return head + "\n".join(rows) + "\n\nAttachment (D3, V2 map):\n\n" + ahead + "\n".join(att) + "\n"


def main():
    C.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--B", type=int, default=20)
    ap.add_argument("--scenarios", nargs="+", default=["planted", "hier", "null", "armblock"])
    ap.add_argument("--out", default=str(paths.RESULTS / "construction" / "toy_decisions.json"))
    ap.add_argument("--summarise-only", action="store_true")
    a = ap.parse_args()
    if a.summarise_only:
        md = summarise(json.loads(Path(a.out).read_text()))
        Path(a.out).with_suffix(".md").write_text(md)
        print(md)
        return
    t0 = time.perf_counter()
    res = dict(settings=dict(config=str(TOY_CFG), B=a.B, frac=0.8, boot_seeds=f"1000..{1000 + a.B - 1}",
                             leiden_seeds=list(C.SEEDS), z0=Z0, k=K, toy_seeds=a.seeds, variants=VARIANTS), runs=[])
    for seed in a.seeds:
        for scen in a.scenarios:
            T = toy.simulate(scen, seed=seed)
            S0, _ = genefilter.dependency_bearing(T.X)
            held = []
            if T.complexes:                       # D3: hold one member of each planted complex of size >= 5 out of S
                shared = set(np.intersect1d(T.complexes[5], T.complexes[6]).tolist())
                for c in T.complexes:
                    cs = [g for g in c if S0[g] and g not in shared]
                    if len(c) >= 5 and cs:
                        held.append(int(cs[-1]))
            S = S0.copy()
            S[held] = False
            sidx = np.flatnonzero(S)
            full, boots = blocks(T.X, S, a.B)
            ref = full[np.ix_(sidx, sidx)]
            graphs = [build_graph(b, K, Z0) for b in boots]
            row = dict(seed=seed, scenario=scen, n_S=int(S.sum()), n_heldout=len(held), variants={})
            lv = {}
            for name, (sup, plc) in VARIANTS.items():
                cfg = C.ClusterConfig(support_mode=sup, placement=plc)
                if sup not in lv:
                    lv[sup] = C.build_levels(graphs, cfg)
                cm = C.build_map(graphs, cfg, ref, Z0, levels=lv[sup])
                ev = evaluate(T, sidx, cm)
                if name == "V2_containment_nested":
                    ev["attachment"] = attachment_eval(T, S, held, full, cm, Z0)
                row["variants"][name] = ev
            if scen == "null":                     # D3 FDR check: random size-matched member sets in S
                rng = np.random.default_rng(4000 + seed)
                rnd = [np.sort(rng.choice(len(sidx), s, replace=False)) for s in toy.SIZES]
                ridx = np.flatnonzero(~S)
                df, summ = attach.attach(full[np.ix_(ridx, sidx)], rnd, Z0)
                row["null_random_sets_attachment"] = dict(summ, n_attached=int(len(df)))
            res["runs"].append(row)
            v = row["variants"]
            print(seed, scen, {k: (e.get("recovered_any_F1ge0p8"), e.get("false_clusters"), e["n_level1"],
                                   e["n_level2"], e["n_locus"], e.get("hier_parents_recovered_level1"),
                                   e.get("hier_subs_recovered_any"), e.get("hier_overlap_recovered_any"))
                              for k, e in v.items()}, flush=True)
    res["runtime_s"] = round(time.perf_counter() - t0, 1)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    out.with_suffix(".md").write_text(summarise(json.loads(out.read_text())))
    print("written", out, res["runtime_s"], "s")


if __name__ == "__main__":
    main()
