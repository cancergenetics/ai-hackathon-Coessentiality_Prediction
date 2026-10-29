"""GPD-tail sensitivity and reporting tables (PI ruling R3; pre-declared in results/part2/directions/SENSITIVITY_GPD.md).
Sensitivity analysis only: no trial line, no registered output changed.

For each direction with a committed discovery run: p_cons = max(p_final, analytic p); the entry's FDR procedure is
re-run on p_cons with assoc._fdr; overlap with the primary discoveries is reported. Also writes
discoveries_report.tsv (ranked by |t|; p_perm shown as '<1e-4' whenever the GPD branch was used).
Usage: .venv/bin/python scripts/45_gpd_sensitivity.py [ID ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from part2 import assoc, h1, h2, registry  # noqa: E402
from part2.cmap import FrozenMap  # noqa: E402

IDS = ("H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI", "H2_submodule_divergence")
GPD_BELOW = 1e-4
MAP = "results/part2/fzx/zavg_engine_membership.tsv"


def committed_table(proj, aid) -> pd.DataFrame:
    summ = proj.results_dir / "assoc" / aid / "discovery_summary.json"
    registry.committed_blob(proj.root, summ)
    p = proj.cache_dir / "assoc" / aid / "discovery.parquet"
    if assoc._sha_file(p) != json.loads(summ.read_text())["full_table_sha256"]:
        raise RuntimeError(f"{aid}: discovery table differs from the committed summary")
    return pd.read_parquet(p)


def p_report(r) -> str:
    if not np.isfinite(r.p_perm):
        return "n/a (not permutation re-tested)"
    return "<1e-4" if r.p_perm < GPD_BELOW else f"{r.p_perm:.4g}"


def run(proj, aid, fmap, parents=None) -> dict:
    entry, _ = registry.load_entry(proj.registry_dir, aid, root=proj.root)
    ff = entry["fdr_family"]
    res = committed_table(proj, aid)
    res["evaluable"] = res.evaluable.astype(bool)
    sel = None
    if ff["method"] == "stagewise_bb":
        res["eligible"] = res.eligible.astype(bool)
        subs = h2.h2_units(fmap)
        sel = assoc.parent_selection_from(parents, {s: fmap.parent[s] for s in subs}, q=float(ff["q"]))
    res["p_cons"] = np.fmax(res.p_final.astype(float), res.p.astype(float))
    adj, rej, nsel = assoc._fdr(res, "p_cons", ff, sel)
    res["q_cons"], res["discovery_cons"] = adj, rej
    prim, cons = res.discovery.astype(bool), res.discovery_cons.astype(bool)
    per_ctx = {c: dict(primary=int(prim[g.index].sum()), conservative=int(cons[g.index].sum()),
                       both=int((prim & cons)[g.index].sum()))
               for c, g in res.groupby("context") if prim[g.index].any() or cons[g.index].any()}
    out = dict(analysis_id=aid, label="sensitivity analysis (GPD tail; PI ruling R3); not a trial",
               fdr=ff, n_primary=int(prim.sum()), n_conservative=int(cons.sum()), both=int((prim & cons).sum()),
               primary_only=int((prim & ~cons).sum()), conservative_only=int((~prim & cons).sum()),
               n_families_selected_conservative=nsel,
               n_gpd_used_among_primary=int((prim & (res.p_perm < GPD_BELOW)).sum()),
               n_pperm_below_analytic_among_primary=int((prim & (res.p_perm < res.p)).sum()), per_context=per_ctx)
    if aid.startswith("H1"):
        out["axes_changing"] = res.loc[prim != cons, ["unit"]].assign(
            primary=prim[prim != cons], conservative=cons[prim != cons]).to_dict("records")
    d = res[prim].copy()
    lo, hi = h1.ci(d["beta"].values.astype(float), d["se"].values.astype(float), d["df"].values.astype(float), 0.95)
    d["ci95_lo"], d["ci95_hi"] = lo, hi
    d["p_perm_report"] = [p_report(r) for r in d.itertuples(index=False)]
    d["q_registered"] = d.q_adj_perm if "q_adj_perm" in d else d.q_adj
    base = d.unit.astype(str).str.split("|").str[0]
    d["tier"] = [fmap.tier.get(u, fmap.tier.get(fmap.parent.get(u, ""), "")) for u in base]
    d["locus"] = [bool(fmap.locus.get(u, False)) for u in base]
    d = d.assign(_a=-d.t.abs()).sort_values(["_a", "context", "unit"]).drop(columns="_a")
    cols = ["context", "unit", "n", "n_context", "beta", "ci95_lo", "ci95_hi", "t", "p", "p_perm_report",
            "q_registered", "discovery_cons", "tier", "locus"]
    od = proj.results_dir / "directions" / aid
    d[cols].to_csv(od / "discoveries_report.tsv", sep="\t", index=False, float_format="%.6g")
    (od / "gpd_sensitivity.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=str) + "\n")
    return out, res


def main():
    proj = assoc.Project()
    fmap = FrozenMap.from_membership(pd.read_csv(proj.root / MAP, sep="\t"))
    ids = sys.argv[1:] or [i for i in IDS if (proj.results_dir / "assoc" / i / "discovery_summary.json").exists()]
    summ, tabs = {}, {}
    for aid in ids:
        parents = None
        if aid.startswith("H2"):
            parents = {"H3a": committed_table(proj, "H3a_lineage_atlas"), "H3b": committed_table(proj, "H3b_driver_family")}
        o, tabs[aid] = run(proj, aid, fmap, parents)
        summ[aid] = {k: v for k, v in o.items() if k not in ("per_context", "axes_changing", "fdr")}
        print(json.dumps(summ[aid]))
    p = proj.results_dir / "directions" / "gpd_sensitivity_summary.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(summ)
    p.write_text(json.dumps(old, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
