#!/usr/bin/env python
"""FZ-C: (z0, k) selection by split-half recovery and the zavg_BEF / GLS consensus-persistence maps on A.

DESIGN v3.1 §2.3 + docs/reviews/round2_algorithms.md + results/part2/construction/DECISIONS.md (D1-D5).
Inputs: cached S x S blocks of the construction runs (scripts/34_run_A_blocks.sh), the gene filter on A, HGNC
location + Ensembl coordinates (locus flag). No physical-interaction or protein-abundance data, context label, hu.MAP file or lockbox model.

  select : split-half statistic for (z0, k) in {4,5,6} x {10,20}, R = 4 splits (seeds 2000-2003), observed and
           degree-rewired; zavg graphs, and GLS graphs edge-count-matched to the zavg half graph at the same (z0, k).
           zavg's (z0, k) = highest statistic, ties -> higher z0, then lower k. GLS's own preference is reported only.
  build  : zavg map on the 20 subsample graphs (seeds 1000-1019) x Leiden SEEDS at the selected (z0, k), cohesion on
           the full-A block; GLS map with per-run edge-count-matched graphs (same S, k), cohesion on the GLS A block
           rank-mapped to the zavg A block (D4); locus flags; attachment of non-S genes to the zavg map (D3).
Writes results/part2/construction/{selection,zavg_map,gls_map}/ and data_cache/part2/construction/ (levels).

  PYTHONHASHSEED=0 .venv/bin/python scripts/36_build_maps.py --stage all --n-proc 4
"""
import argparse
import dataclasses
import hashlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from part2 import attach, construct, data, genefilter, locus, paths, runsets, scores, selection  # noqa: E402
from part2 import consensus as C  # noqa: E402
from part2.graph import build_graph  # noqa: E402

OUT = paths.RESULTS / "construction"
CACHE = paths.CACHE / "construction"
GRID = [(z0, k) for z0 in (4.0, 5.0, 6.0) for k in (10, 20)]
# FZ-C decisions D1 (containment support) and D2 (nested placement); everything else = DESIGN / round-2 defaults
CCFG = C.ClusterConfig(support_mode="containment", placement="nested")
SCFG = selection.SplitHalfConfig()


def gene_context():
    """S mask, S-local gene table (symbol, entrez, essentiality class), locus arrays."""
    S = runsets.S_mask()
    _, st = genefilter.dependency_bearing_for_models(list(runsets.A()), return_stats=True)
    genes = data.protein_genes()
    sidx = np.flatnonzero(S)
    fb = st["frac_below"][sidx]
    ess = np.where(fb >= 0.9, "pan", np.where(fb >= 0.05, "selective", "other"))
    loc = locus.gene_locations(genes).iloc[sidx]
    g = pd.DataFrame(dict(gene_index=sidx, symbol=genes.symbol.values[sidx], entrez=genes.entrez.values[sidx],
                          frac_below_A=fb, ess_class=ess, location=loc.location.values, arm=loc.arm.values,
                          band=loc.band.values, chrom=loc.chrom.values, pos=loc.pos.values))
    return S, g


# ------------------------------------------------------------------------------------------ selection
def _sel_job(args):
    method, z0, k, si, gctx = args
    S, g = gctx
    sel = genefilter.selective_mask(g.frac_below_A.values)
    arm, band, chrom, pos = g.arm.values, g.band.values, g.chrom.values, g.pos.values
    t = time.perf_counter()
    za = construct.load(f"half{runsets.HALF_SEEDS[si]}a", "zavg_BEF", S).vals
    zb = construct.load(f"half{runsets.HALF_SEEDS[si]}b", "zavg_BEF", S).vals
    g1, g2 = build_graph(za, k, z0), build_graph(zb, k, z0)
    extra = {}
    if method == "gls":
        ga = construct.load(f"half{runsets.HALF_SEEDS[si]}a", "gls", S).vals
        gb = construct.load(f"half{runsets.HALF_SEEDS[si]}b", "gls", S).vals
        (h1, i1), (h2, i2) = construct.gls_matched(ga, g1, k), construct.gls_matched(gb, g2, k)
        za, zb = construct.rank_map(ga, za), construct.rank_map(gb, zb)     # cohesion scale (D4)
        del ga, gb
        g1, g2 = h1, h2
        extra = dict(matched_h1=i1, matched_h2=i2)
    rows = selection.split_half_rows(si, g1, g2, za, zb, z0, sel, CCFG, SCFG, arm, band, 4000, chrom, pos)
    for r in rows:
        r.update(method=method, z0=z0, k=k, runtime_s=round(time.perf_counter() - t, 1), **extra)
    print(f"sel {method} z0={z0} k={k} split={si} {time.perf_counter() - t:.0f}s "
          f"stat={rows[0]['stat']:.4f} peak_rss={scores.peak_rss_gb():.2f}GB", flush=True)
    return rows


def stage_select(gctx, n_proc):
    t0 = time.perf_counter()
    jobs = [(m, z0, k, si, gctx) for m in ("zavg_BEF", "gls") for z0, k in GRID for si in range(len(runsets.HALF_SEEDS))]
    with mp.get_context("spawn").Pool(n_proc) as pool:
        res = pool.map(_sel_job, jobs, chunksize=1)
    rows = [r for rr in res for r in rr]
    df = pd.DataFrame(rows)
    d = OUT / "selection"
    d.mkdir(parents=True, exist_ok=True)
    df.drop(columns=[c for c in df.columns if c.startswith("matched_")]).to_csv(d / "split_half_per_split.tsv", sep="\t", index=False)
    summ = {}
    for m in ("zavg_BEF", "gls"):
        out = []
        for z0, k in GRID:
            sub = df[(df.method == m) & (df.z0 == z0) & (df.k == k)]
            s = selection.summarise_split_half(z0, k, sub.to_dict("records"))
            s.pop("per_split")
            obs = sub[sub.variant == "observed"]
            s.update(mean_edges_half=float((obs.edges_h1 + obs.edges_h2).mean() / 2),
                     mean_clusters_half=float((obs.n_clusters_h1 + obs.n_clusters_h2).mean() / 2),
                     mean_eval_clusters=float(obs.n_eval.mean() / 2), mean_recovered_clusters=float(obs.n_recovered.mean() / 2))
            out.append(s)
        ch = selection.select_z0k([dict(x) for x in out])
        summ[m] = dict(grid=out, selected=dict(z0=ch["z0"], k=ch["k"], stat=ch["stat"]))
    summ["rule"] = ("highest mean split-half statistic (fraction of selective genes in recovered size 5-30 non-locus "
                    "clusters, Jaccard >= 0.5, both directions, R = 4); ties (1e-12) -> higher z0, then lower k")
    summ["primary"] = dict(zavg_BEF=summ["zavg_BEF"]["selected"],
                           gls="matched to the zavg selection (same S, k; per-run edge count); GLS's own split-half "
                               "preference is reported as a secondary variant only")
    summ["runtime_s"] = round(time.perf_counter() - t0, 1)
    (d / "selection.json").write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps({m: summ[m]["selected"] for m in ("zavg_BEF", "gls")}), summ["runtime_s"], "s", flush=True)
    return summ


# ------------------------------------------------------------------------------------------ maps
def cluster_table(cm, g, prefix, ref, levels):
    stab = selection.stability_by_size(cm, levels).set_index("id").stability
    rows, memb = [], []
    l1 = sorted(cm.level(1), key=lambda c: c["id"])
    name = {c["id"]: f"{prefix}{i + 1:04d}" for i, c in enumerate(l1)}
    for p in l1:
        kids = sorted([c for c in cm.level(2) if c["parent"] == p["id"]], key=lambda c: (-len(c["members"]), int(c["members"][0])))
        for j, c in enumerate(kids):
            name[c["id"]] = f"{name[p['id']]}.{j + 1}"
    for c in sorted(cm.clusters, key=lambda c: name[c["id"]]):
        m = np.asarray(c["members"])
        lf = locus.locus_flag(m, g.arm.values, g.band.values, chrom=g.chrom.values, pos=g.pos.values)
        ess = g.ess_class.values[m]
        rows.append(dict(cluster_id=name[c["id"]], level=c["level"],
                         parent=name[c["parent"]] if c["level"] == 2 else "", size=len(m),
                         gamma=round(c["gamma"], 6), j_start=c["j_start"], j_end=c["j_end"], n_adjacent=c["n_adjacent"],
                         support=round(c["support"], 4), stability=round(float(stab[c["id"]]), 4),
                         cohesion=round(C.cohesion(m, ref), 4), origin=c.get("origin", ""),
                         n_pan=int((ess == "pan").sum()), n_selective=int((ess == "selective").sum()),
                         locus=lf["locus"], top_arm=lf["top_arm"], arm_share=round(lf["arm_share"], 3),
                         window_share=round(lf["window_share"], 3), top_band=lf["top_band"],
                         members=";".join(g.symbol.values[m]), entrez=";".join(map(str, g.entrez.values[m]))))
        for i in m:
            memb.append(dict(cluster_id=name[c["id"]], level=c["level"], symbol=g.symbol.values[i],
                             entrez=int(g.entrez.values[i]), gene_index=int(g.gene_index.values[i]),
                             ess_class=g.ess_class.values[i]))
    return pd.DataFrame(rows), pd.DataFrame(memb), name


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_manifest(d: Path):
    files = sorted(p for p in d.iterdir() if p.is_file() and p.name != "MANIFEST.sha256")
    (d / "MANIFEST.sha256").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in files))


def stage_build(gctx, sel, n_proc):
    S, g = gctx
    z0, k = sel["zavg_BEF"]["selected"]["z0"], sel["zavg_BEF"]["selected"]["k"]
    names = [f"boot{s}" for s in runsets.BOOT_SEEDS]
    CACHE.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    zg, zmeta = construct.zavg_graphs(names, k, z0, S)
    t_g = time.perf_counter() - t0
    blkA = construct.load("A", "zavg_BEF", S)
    refA = blkA.vals
    results = {}
    for method in ("zavg_BEF", "gls"):
        t = time.perf_counter()
        if method == "zavg_BEF":
            graphs, gmeta, ref, prefix, tg = zg, zmeta, refA, "Z", t_g
        else:
            graphs, gmeta = construct.gls_graphs(names, k, zg, S)
            glsA = construct.load("A", "gls", S)
            ref = construct.rank_map(glsA.vals, refA)
            gls_meta_A = dict(model_set_hash=glsA.meta["model_set_hash"])
            del glsA
            prefix, tg = "G", time.perf_counter() - t
        cm, levels, tim = construct.build_map(graphs, ref, z0, CCFG, n_proc)
        tab, memb, name = cluster_table(cm, g, prefix, ref, levels)
        d = OUT / ("zavg_map" if method == "zavg_BEF" else "gls_map")
        d.mkdir(parents=True, exist_ok=True)
        tab.to_csv(d / "clusters.tsv", sep="\t", index=False)
        memb.to_csv(d / "membership.tsv", sep="\t", index=False)
        np.savez_compressed(CACHE / f"{method}_levels.npz", gammas=np.array([lv.gamma for lv in levels]),
                            **{f"runs_{j}": lv.runs.astype(np.int32) for j, lv in enumerate(levels)})
        att_summ = None
        t_att = 0.0
        if method == "zavg_BEF":
            ta = time.perf_counter()
            rv, ridx = blkA.rect
            cl = [np.asarray(c["members"]) for c in cm.clusters]
            ids = np.array([name[c["id"]] for c in cm.clusters])
            genes = data.protein_genes()
            df, att_summ = attach.attach(rv, cl, z0, row_ids=ridx, cluster_ids=ids)
            df.insert(1, "symbol", genes.symbol.values[df.gene.values])
            df.insert(2, "entrez", genes.entrez.values[df.gene.values])
            # D3 locus guard: attachments to locus clusters, or of a gene on the cluster's dominant arm / within
            # 10 Mb of a member, are flagged and excluded from H1 event definitions
            gl = locus.gene_locations(genes)
            ct = tab.set_index("cluster_id")
            mem_of = {name[c["id"]]: np.asarray(c["members"]) for c in cm.clusters}
            gi = df.gene.values.astype(int)
            garm, gch, gpos = gl.arm.values[gi], gl.chrom.values[gi], gl.pos.values[gi]
            near = []
            for i, cid in enumerate(df.cluster.values):
                m = mem_of[cid]
                same = (g.chrom.values[m] == gch[i]) & (g.pos.values[m] >= 0) & (gpos[i] >= 0)
                near.append(bool((np.abs(g.pos.values[m][same] - gpos[i]) <= locus.WINDOW_BP).any()))
            df["gene_arm"] = garm
            df["cluster_locus"] = ct.locus.reindex(df.cluster.values).values.astype(bool)
            df["same_arm_as_cluster"] = (garm != "") & (garm == ct.top_arm.reindex(df.cluster.values).values)
            df["within_10Mb_of_member"] = near
            df["locus_attach"] = df.cluster_locus | df.same_arm_as_cluster | df.within_10Mb_of_member
            df = df.rename(columns={"gene": "gene_index"})
            df.to_csv(d / "attachments.tsv", sep="\t", index=False, float_format="%.4g")
            t_att = time.perf_counter() - ta
        params = dict(method=method, z0=z0, k=k, weight_cap=3.0, cluster_config=dataclasses.asdict(CCFG),
                      split_half_config=dataclasses.asdict(SCFG), selection=sel["zavg_BEF"]["selected"],
                      gls_matching=("same S and k as zavg; per run: top mutual-top-k GLS pairs up to the zavg edge "
                                    "count of the same model set; weights quantile-matched; cohesion on the GLS A "
                                    "block rank-mapped onto the zavg A block (D4)") if method == "gls" else None,
                      S=dict(n=int(S.sum()), hash=scores.mask_hash(S), rule="effect < -0.5 in >= 1% of A, non-NaN "
                             ">= 95% of A, not Part 1 low-coverage"),
                      A=dict(n_models=len(runsets.A()), model_set_hash=data.model_set_hash(list(runsets.A()))),
                      runs=gmeta, reference_block=dict(run="A", model_set_hash=blkA.meta["model_set_hash"],
                                                       null_moments=blkA.meta.get("null_moments")),
                      attachment=dict(rule="x >= ceil(|C|/2) members with score >= z0 and BH q < 0.01 over all "
                                      "gene x cluster tests; hypergeometric size-matched null; clusters = all "
                                      "levels; genes = protein genes outside S", summary=att_summ) if att_summ else None,
                      counts=cm.counts, seeds=dict(leiden=list(CCFG.seeds), bootstraps=list(runsets.BOOT_SEEDS),
                                                   split_half=list(runsets.HALF_SEEDS), rewire=4000),
                      code_commit=scores.code_commit())
        if method == "gls":
            params["reference_block"]["gls_A_model_set_hash"] = gls_meta_A["model_set_hash"]
        (d / "parameters.json").write_text(json.dumps(params, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
        rt = dict(graphs_s=round(tg, 1), **tim, attachment_s=round(t_att, 1), total_s=round(time.perf_counter() - t + (t_g if method == "zavg_BEF" else 0), 1),
                  peak_rss_gb_process=round(scores.peak_rss_gb(), 2), n_proc=n_proc)
        (d / "runtime.json").write_text(json.dumps(rt, indent=1))
        write_manifest(d)
        l1 = tab[tab.level == 1]
        results[method] = dict(n_level1=int(len(l1)), n_level2=int((tab.level == 2).sum()),
                               genes_level1=int(l1["size"].sum()), locus_level1=int(l1.locus.sum()), runtime=rt)
        print(method, json.dumps(results[method]), flush=True)
    return results


def main():
    C.ensure_hashseed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["select", "build", "all"], default="all")
    ap.add_argument("--n-proc", type=int, default=4)
    a = ap.parse_args()
    gctx = gene_context()
    if a.stage in ("select", "all"):
        sel = stage_select(gctx, a.n_proc)
    else:
        sel = json.loads((OUT / "selection" / "selection.json").read_text())
    if a.stage in ("build", "all"):
        stage_build(gctx, sel, a.n_proc)


if __name__ == "__main__":
    main()
