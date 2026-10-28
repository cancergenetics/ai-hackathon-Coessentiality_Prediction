#!/usr/bin/env python
"""Build the Part 2 context table (DESIGN v3.1 §2.2) over U from DepMap Omics + Model.csv (+ Ensembl/HGNC loci,
CCLE proteomics for the LoF functional check). No gene-effect / dependency value is read.

Outputs
  data_cache/part2/contexts/   (gitignored)  context_table.parquet, lof_matrix.parquet, cn_<event>.parquet,
                                             segment_matrix.parquet
  results/part2/contexts/      (committed)   context_summary.tsv, lof_functional_check.tsv, segments_ge8.tsv,
                                             build_log.json (CN-scale check, LoY threshold, fusion counts, ...),
                                             checksums.json (SHA-256 of every output + content hashes + inputs)

  PYTHONHASHSEED=0 .venv/bin/python scripts/32_build_contexts.py
"""
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import consensus, contexts, data, paths  # noqa: E402

INPUTS = [paths.OMICS_PROFILES, paths.OMICS_DAMAGING, paths.OMICS_HOTSPOT, paths.OMICS_CN, paths.OMICS_EXPR,
          paths.OMICS_SIGNATURES, paths.OMICS_FUSIONS, paths.MODEL_CSV, paths.HGNC_FILE, paths.GENE_COORDS,
          paths.CCLE_PROTEOMICS_XLSX]


def git_head():
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def main():
    consensus.ensure_hashseed()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    t0 = time.time()
    U = data.universe_model_ids()
    assert len(U) == 1345, len(U)
    b = contexts.build(U)
    t_build = time.time() - t0

    cache, res = paths.CONTEXT_CACHE, paths.CONTEXT_RESULTS
    cache.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)
    out = {"context_table.parquet": b.table, "lof_matrix.parquet": b.lof, "segment_matrix.parquet": b.segment_matrix}
    for k, v in b.cn_events.items():
        out[f"cn_{k}.parquet"] = v
    files, content = {}, {}
    for name, df in out.items():
        df = df.copy()
        df.columns = df.columns.astype(str)
        p = cache / name
        df.to_parquet(p)
        files[str(p.relative_to(paths.PROJECT))] = data.sha256_file(p)
        content[name] = contexts.content_hash(df)

    summ = contexts.summary(b)
    summ.to_csv(res / "context_summary.tsv", sep="\t", index=False)
    b.functional.to_csv(res / "lof_functional_check.tsv", sep="\t", index=False, float_format="%.4g")
    # committed: hom-del and high-amp segments with >= 8 U models (loss / amp segments are many and short; full table
    # in the cache: segments.tsv)
    b.segments.to_csv(cache / "segments.tsv", sep="\t", index=False)
    seg8 = b.segments[(b.segments.n_events >= 8) & b.segments.event.isin(["homdel", "highamp"])].drop(columns="genes")
    seg8.to_csv(res / "segments_ge8.tsv", sep="\t", index=False)
    blog = dict(b.log, runtime_s=round(time.time() - t0, 1), build_s=round(t_build, 1), code_commit=git_head())
    (res / "build_log.json").write_text(json.dumps(blog, indent=1, default=str) + "\n")
    files[str((cache / "segments.tsv").relative_to(paths.PROJECT))] = data.sha256_file(cache / "segments.tsv")
    for name in ("context_summary.tsv", "lof_functional_check.tsv", "segments_ge8.tsv", "build_log.json"):
        files[str((res / name).relative_to(paths.PROJECT))] = data.sha256_file(res / name)
    ck = dict(files=files, content_sha256=content, universe_hash=data.model_set_hash(U), n_U=len(U),
              inputs={str(p.relative_to(paths.PROJECT)): data.sha256_file(p) for p in INPUTS},
              code_commit=git_head())
    ck["files"].pop(str((res / "build_log.json").relative_to(paths.PROJECT)))   # runtime varies; not frozen
    (res / "checksums.json").write_text(json.dumps(ck, indent=1, sort_keys=True) + "\n")
    logging.info("done in %.1f s; %d context columns", time.time() - t0, len(summ))


if __name__ == "__main__":
    main()
