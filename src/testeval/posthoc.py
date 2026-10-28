"""Post hoc reference (secondary, not a candidate): PC-removed Pearson ``spec_k20_tinf_std`` exactly as Part 1 defined
it for its review-requested reference (final_eval.py --posthoc-reference, look '2-posthoc'): gene-standardised,
top-20 PCs removed, no whitening (tau = inf), cosine = results/layers/manifest.json (screen_v1) column 64.

Part 1's predict_pairs.py does not expose this layer, so it is computed here with Part 1's own LayerEngine on the
raw CSV (Part 1 predict.load_matrix: hash-checked against the layer manifests), Part 1's ID mapping
(predict.map_ids) and Part 1's unmapped rule (score = min(mapped) - 1, mapped = 0). Output and sidecar follow
predict_pairs.py's format. The frozen config and the genome-wide null moments are asserted equal to the manifest.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from testeval import part1  # noqa: F401  (puts the locked clone on sys.path)
from src import layers, predict  # noqa: E402
from src.layers import INF, LayerEngine, Spectral  # noqa: E402

NAME = "spec_k20_tinf_std"
CFG = Spectral(20, INF, standardize=True)
MANIFEST = predict.RESULTS / "layers" / "manifest.json"


def manifest_column() -> dict:
    cols = {c["name"]: c for c in json.loads(MANIFEST.read_text())["columns"]}
    c = cols[NAME]
    if c["config"] != CFG.to_dict():
        raise RuntimeError(f"manifest config {c['config']} != frozen {CFG.to_dict()}")
    return c


def run(pairs_path, out_path, record_path=None) -> dict:
    t0 = time.perf_counter()
    ref = manifest_column()
    pairs = predict.read_pairs(pairs_path)
    X, genes, msrc = predict.load_matrix(cache=False)
    eng = LayerEngine(X)
    del X
    null = layers.null_pairs(eng)
    mo = layers.moments(eng.compute(CFG, *null))
    d_mean = abs(mo["null_mean"] - ref["null_mean"]) / ref["null_sd"]
    d_sd = abs(mo["null_sd"] - ref["null_sd"]) / ref["null_sd"]
    if d_mean > predict.NULL_TOL or d_sd > predict.NULL_TOL:
        raise RuntimeError(f"null moments differ from manifest: {mo} vs {ref}")
    ia, ib, umap = predict.map_ids(pairs, genes)
    ok = (ia >= 0) & (ib >= 0) & (ia != ib)
    s = np.full(len(pairs), np.nan)
    if ok.any():
        lo, hi = np.minimum(ia[ok], ib[ok]), np.maximum(ia[ok], ib[ok])
        s[ok] = eng.compute(CFG, lo, hi)
    mapped = ok & np.isfinite(s)
    floor = (float(s[mapped].min()) if mapped.any() else 0.0) - 1.0
    s[~mapped] = floor
    out = pairs[["uniprot_a", "uniprot_b"]].copy()
    out["score"] = s
    out["mapped"] = mapped.astype(np.int8)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, na_rep="NA")
    rec = dict(tool="src/testeval/posthoc.py (Part 1 LayerEngine, locked clone)", method=NAME, config=CFG.to_dict(),
               role="post hoc reference (Part 1 review B1), secondary, not a candidate",
               part1=part1.lock_state(),
               inputs=dict(pairs=dict(path=str(pairs_path), sha256=predict.sha256_file(pairs_path)), matrix=msrc),
               output=dict(path=str(out_path), sha256=predict.sha256_file(out_path)),
               null=dict(**mo, manifest_null_mean=ref["null_mean"], manifest_null_sd=ref["null_sd"],
                         rel_diff_mean=d_mean, rel_diff_sd=d_sd),
               counts=dict(n_pairs=len(pairs), n_mapped=int(mapped.sum()), n_unmapped=int((~mapped).sum()),
                           unmapped_floor=floor),
               unmapped_rule="score = min(mapped scores in this file) - 1, mapped = 0 (Part 1 DESIGN §3.10)",
               runtime=dict(total_s=round(time.perf_counter() - t0, 2)))
    rp = Path(record_path) if record_path else out_path.with_name(out_path.name + ".run.json")
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text(json.dumps(rec, indent=1, default=float))
    return rec
