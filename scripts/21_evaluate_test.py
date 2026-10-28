"""Evaluate scored pair sets under docs/TEST_PROTOCOL.md.

  --mode fit     fit decision rules + calibrators on ALL hu.MAP train pairs (scores from
                 scripts/20_score_test_pairs.sh train) -> results/test_eval/fitted_rules.json (no test file is read)
  --mode test    THE SINGLE TEST LOOK. Refuses unless docs/TEST_PROTOCOL.md and fitted_rules.json are committed and
                 src/, scripts/ are clean, and no earlier look is logged in results/test_eval/looks.jsonl. Logs the
                 look BEFORE reading labels, writes results/test_eval/test_results.json (+ figures) with the
                 protocol commit hash, then logs completion.
  --mode dryrun  development data only (Part 1 dev / IFT split of the train pairs):
                 A  reproduction of Part 1's IFT look 1 with Part 1's definitions (HT-AP etc.),
                 B  the full protocol pipeline on the IFT (rules fitted on Part 1 dev pairs),
                 C  the test-mode code path on all train pairs (in-sample; mechanics, runtime, expected CI width)
                 -> results/test_eval_dryrun/
Re-executes itself with PYTHONHASHSEED=0 (Part 1 Louvain tie-breaking).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

if os.environ.get("PYTHONHASHSEED") != "0":
    os.execve(sys.executable, [sys.executable, *sys.argv], dict(os.environ, PYTHONHASHSEED="0"))

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from testeval import evaluate as E  # noqa: E402
from testeval import paths  # noqa: E402
from testeval import part1  # noqa: E402
from testeval.protocol import append_jsonl, git_state, read_jsonl, require_frozen, sha256_file  # noqa: E402

T = {}


def tick(name, t0):
    T[name] = round(time.perf_counter() - t0, 1)
    print(f"  [{name}] {T[name]} s", flush=True)


def jdefault(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def dump(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=jdefault))


def score_inputs(score_dir: Path) -> dict:
    out = {}
    for m in E.METHODS:
        f = score_dir / f"{m}.tsv"
        out[m] = dict(path=str(f.relative_to(paths.PROJECT)), sha256=sha256_file(f), predictor=E.PREDICTOR_ARG[m])
    return out


def load_matrix():
    from src import predict
    X, genes, src = predict.load_matrix(cache=False)
    return X, genes, src


def p1_assignment() -> pd.DataFrame:
    import hashlib
    import io
    b = paths.P1_ASSIGNMENT.read_bytes()
    if hashlib.sha256(b).hexdigest() != paths.P1_ASSIGNMENT_SHA.read_text().split()[0]:
        raise RuntimeError("Part 1 assignment hash mismatch")
    return pd.read_csv(io.BytesIO(b))


def ht_weights_fn(assign: pd.DataFrame, side: str):
    """Part 1 design weights for a cluster subsample: within-cluster 1, between-cluster 1/f (f = share of mapped
    proteins on this side); pairs with an unmapped protein 1/f (Part 1 final_eval §3.10 rule)."""
    m = assign[assign.mapped]
    f = float((m.side == side).mean())
    c = assign.set_index("protein").cluster

    def fn(t, mapped):
        ca, cb = c.reindex(t.a).values, c.reindex(t.b).values
        both = np.logical_and.reduce(list(mapped.values()))
        w = np.where((ca == cb) & both, 1.0, 1.0 / f)
        return w
    fn.f = f
    return fn


def figures(res: dict, out_dir: Path, tag: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for m, d in res["metrics"].items():
        rel = d.get("calibration", {}).get("isotonic", {}).get("reliability_log16")
        if not rel:
            continue
        x = [max(r["mean_pred"], 1e-5) for r in rel]
        y = [max(r["obs_rate"], 1e-5) for r in rel]
        axes[0].plot(x, y, marker="o", ms=3, lw=1, label=m)
        relp = d["calibration"]["platt"]["reliability_log16"]
        axes[1].plot([max(r["mean_pred"], 1e-5) for r in relp], [max(r["obs_rate"], 1e-5) for r in relp],
                     marker="o", ms=3, lw=1, label=m)
    for ax, t in zip(axes, ("isotonic (fit on train pairs)", "Platt (fit on train pairs)")):
        ax.plot([1e-5, 1], [1e-5, 1], color="grey", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("mean predicted probability (16 log-spaced bins)")
        ax.set_ylabel("observed co-complex rate")
        ax.set_title(t, fontsize=10)
    axes[0].legend(fontsize=7)
    fig.suptitle(f"Reliability: {tag} (values < 1e-5 drawn at 1e-5)", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"reliability_{tag}.png", dpi=130)
    plt.close(fig)


def subset_ess(es, ess, keep):
    return None if ess is None else (ess[0][keep], ess[1][keep])


def run_protocol(es, rules, ess, mc_seeds, label, n_boot=E.N_BOOT):
    """Primary (ITT) + secondary analysis sets, identical for dry run and test."""
    out = {}
    t0 = time.perf_counter()
    cl = E.annotate_clusters(es)
    tick(f"{label}: clusters", t0)
    t0 = time.perf_counter()
    out["itt"] = E.evaluate(es, rules, cl, None, n_boot=n_boot, mc_seeds=mc_seeds)
    tick(f"{label}: ITT (primary) metrics + bootstraps", t0)
    keep = es.mapped_all
    sub = es.subset(keep, es.name + " mapped-only")
    clm = dict(cl, ca=cl["ca"][keep], cb=cl["cb"][keep], comp_a=cl["comp_a"][keep], comp_b=cl["comp_b"][keep])
    t0 = time.perf_counter()
    out["mapped_only"] = E.evaluate(sub, rules, clm, subset_ess(es, ess, keep), n_boot=n_boot, full=False)
    tick(f"{label}: mapped-only", t0)
    t = es.table
    if "in_train_pos" in t:
        for name, k in (("excl_train_overlap", ~(t.in_train_pos.values | t.in_train_neg.values)),
                        ("both_proteins_unseen_in_train", ~(t.a_in_train.values | t.b_in_train.values))):
            if k.all() or not (es.y & k).any():
                out[name] = dict(note="identical to ITT" if k.all() else "no positives", n=int(k.sum()))
                continue
            s2 = es.subset(k, f"{es.name} {name}")
            cl2 = dict(cl, ca=cl["ca"][k], cb=cl["cb"][k], comp_a=cl["comp_a"][k], comp_b=cl["comp_b"][k])
            t0 = time.perf_counter()
            out[name] = E.evaluate(s2, rules, cl2, None, n_boot=n_boot, full=False)
            tick(f"{label}: {name}", t0)
    return out, cl


# --------------------------------------------------------------------------------------------------- modes
def mode_fit():
    t0 = time.perf_counter()
    sd = paths.SCORES / "train"
    fit = E.load_set("hu.MAP train (all pairs)", paths.TRAIN_POS, paths.TRAIN_NEG, sd)
    rules = E.fit_rules(fit)
    rules.update(inputs=dict(train_pos=sha256_file(paths.TRAIN_POS), train_neg=sha256_file(paths.TRAIN_NEG),
                             scores=score_inputs(sd)), audit=fit.audit, part1=part1.assert_locked(),
                 note="fitted on train pairs only; no test file read. Frozen with the protocol before the look.")
    dump(rules, paths.FITTED_RULES)
    tick("fit rules on train pairs", t0)
    print(json.dumps({m: dict(f1max_threshold=r["f1max"]["threshold"], f1_fit=r["f1max"]["f1"])
                      for m, r in rules["methods"].items()}, indent=1))


def mode_test(relook_reason=None):
    st = require_frozen()
    g = git_state(paths=(str(paths.FITTED_RULES.relative_to(paths.PROJECT)),))
    if g["dirty"]:
        raise SystemExit("fitted_rules.json is not committed / modified")
    looks = [r for r in read_jsonl(paths.LOOKS) if r.get("look") == "test"]
    done = [r for r in looks if r.get("status") == "done"]
    if looks and not relook_reason:
        raise SystemExit(f"a test look is already logged ({len(looks)} records): single look only; a re-run needs "
                         "--relook-reason (logged; non-confirmatory if a completed look exists, TEST_PROTOCOL §8.3)")
    if any(r.get("protocol_commit") != st["protocol_last_commit"] for r in looks):
        raise SystemExit("the protocol changed since the first logged look: not allowed (TEST_PROTOCOL §8.3)")
    p1 = part1.assert_locked()
    rules = json.loads(paths.FITTED_RULES.read_text())
    rules_sha = sha256_file(paths.FITTED_RULES)
    sd = paths.SCORES / "test"
    head = dict(look="test", status="started", protocol_commit=st["protocol_last_commit"], head=st["head"],
                protocol_sha256=st["protocol_sha256"], fitted_rules_sha256=rules_sha, relook_reason=relook_reason,
                confirmatory=not done, n_prior_records=len(looks))
    append_jsonl(paths.LOOKS, head)                          # logged BEFORE any label is read
    t0 = time.perf_counter()
    pos, neg = paths.test_files()
    es = E.load_set("hu.MAP test", pos, neg, sd, train_paths=(paths.TRAIN_POS, paths.TRAIN_NEG))
    X, genes, msrc = load_matrix()
    ess = E.essential_flags(es, X, genes)
    del X
    tick("load test pairs + scores", t0)
    res, cl = run_protocol(es, rules, ess, E.SEEDS[1:], "test")
    out = dict(protocol=dict(path="docs/TEST_PROTOCOL.md", commit=st["protocol_last_commit"], head=st["head"],
                             sha256=st["protocol_sha256"]),
               fitted_rules=dict(path=str(paths.FITTED_RULES.relative_to(paths.PROJECT)), sha256=rules_sha),
               part1=p1, inputs=dict(test_pos=sha256_file(pos), test_neg=sha256_file(neg), scores=score_inputs(sd),
                                     matrix=msrc),
               audit=es.audit, confirmatory=not done, relook_reason=relook_reason, prior_look_records=looks,
               results=res, runtime_s=T)
    dump(out, paths.RESULTS / "test_results.json")
    figures(res["itt"], paths.RESULTS / "figures", "test_itt")
    append_jsonl(paths.LOOKS, dict(head, status="done", result="results/test_eval/test_results.json",
                                   verdict=res["itt"].get("verdict", {}).get("candidates")))


def compare(ours: dict, ref: dict, keys=None) -> dict:
    out = {}
    for k, v in ref.items():
        if keys and k not in keys:
            continue
        if isinstance(v, (int, float)) and not isinstance(v, bool) and k in ours and ours[k] is not None:
            out[k] = dict(ours=ours[k], part1=v, abs_diff=abs(float(ours[k]) - float(v)))
    return out


def mode_dryrun(n_boot, parts="ABCI"):
    out_dir = paths.DRYRUN
    d = paths.CACHE / "dryrun"
    assign = p1_assignment()
    X, genes, msrc = load_matrix()
    summary = dict(part1=part1.assert_locked(), matrix=msrc)

    # ---------------------------------------------------------------- A: reproduction of Part 1 IFT look 1
    t0 = time.perf_counter()
    wfn_ift = ht_weights_fn(assign, "ift")
    ift = E.load_set("Part 1 IFT", d / "ift_positive_pairs.txt", d / "ift_negative_pairs.txt",
                     paths.SCORES / "dryrun_ift", weights_fn=wfn_ift)
    ess_ift = E.essential_flags(ift, X, genes)
    if "A" in parts:
        _part_a(ift, ess_ift, assign, wfn_ift, out_dir, summary, t0)
    if "B" in parts:
        _part_b(assign, d, wfn_ift, ess_ift, n_boot, out_dir, summary)
    if "C" in parts:
        _part_c(X, genes, n_boot, out_dir, summary)
    if "I" in parts:
        summary["I_batch_invariance"] = _part_i()
    summary["runtime_s"] = T
    dump(summary, out_dir / f"summary{'' if parts == 'ABCI' else '_' + parts}.json")


def _part_i() -> dict:
    """I: the score of a pair must not depend on which file it was scored in (IFT file vs all-train file): max |diff|
    over the IFT pairs mapped in both, and equality of the mapped flags, per method."""
    from testeval import pairs as P
    out = {}
    for m in E.METHODS:
        a = P.read_scores(paths.SCORES / "dryrun_ift" / f"{m}.tsv")
        b = P.read_scores(paths.SCORES / "train" / f"{m}.tsv")
        j = a.merge(b, on=["a", "b"], suffixes=("_ift", "_train"), validate="one_to_one")
        mm = j.mapped_ift.values
        out[m] = dict(n_pairs=int(len(a)), n_matched=int(len(j)), mapped_flags_equal=bool((j.mapped_ift == j.mapped_train).all()),
                      max_abs_diff_mapped=float(np.abs(j.score_ift.values[mm] - j.score_train.values[mm]).max()))
    return out


def _part_a(ift, ess_ift, assign, wfn_ift, out_dir, summary, t0):
    from testeval.part1 import metrics as p1m
    keep = ift.mapped_all
    sc = ift.subset(keep, "Part 1 IFT scorable")
    a = assign.set_index("protein")
    ca, cb = a.cluster.reindex(sc.table.a).values, a.cluster.reindex(sc.table.b).values
    oa, ob = a.component.reindex(sc.table.a).values, a.component.reindex(sc.table.b).values
    m_ = assign[assign.mapped]
    largest = int(m_.loc[m_.component_size.idxmax(), "component"])
    f = wfn_ift.f
    w = p1m.ht_weights(ca, cb, f)
    look = json.loads(paths.P1_LOOK01.read_text())
    ph = json.loads(paths.P1_POSTHOC.read_text())
    rep = dict(f_ift=f, f_ift_part1=look["ht_f"], counts=dict(pos=int(sc.y.sum()), neg=int((~sc.y).sum()),
                                                            unmapped_pos=int((ift.y & ~keep).sum()),
                                                            unmapped_neg=int((~ift.y & ~keep).sum())),
               counts_part1=look["counts"], metrics={}, bootstrap={})
    ys = sc.y
    for m in E.METHODS:
        s = sc.scores[m]
        suite = p1m.metric_suite(ys, s, w, ca, cb, oa, ob, ess_ift[0][keep], ess_ift[1][keep], largest_component=largest)
        yu = ift.y[~keep]
        suite["ap_ht_with_unmapped_min_score"] = p1m.weighted_ap(
            np.r_[ys, yu], np.r_[s, np.full(len(yu), s.min() - 1.0)], np.r_[w, np.full(len(yu), 1.0 / f)])
        r = E.Ranked(ys, s)
        suite["harness_ranked_ap_ht"] = r.ap(w)
        suite["harness_ranked_auroc_ht"] = r.auroc(w)
        ref = look["metrics"].get(m) or (ph["metrics"].get(m) if m == E.REFERENCE else None)
        rep["metrics"][m] = dict(ours=suite, vs_part1=compare(suite, ref) if ref else None,
                                 part1_source="look01" if m in look["metrics"] else "look 2-posthoc")
        if ref:
            rep["metrics"][m]["vs_part1"]["harness_ranked_ap_ht"] = dict(
                ours=suite["harness_ranked_ap_ht"], part1=ref["ap_ht"], abs_diff=abs(suite["harness_ranked_ap_ht"] - ref["ap_ht"]))
    look_scorers = list(look["metrics"])
    diffs = [tuple(k.split("-", 1)) for k in look["bootstrap"]["diff"]]
    boot = p1m.cluster_bootstrap(ys, {m: sc.scores[m] for m in look_scorers}, ca, cb, f=f, n_boot=1000,
                                 seed=look["bootstrap"]["seed"], diffs=diffs)
    rep["bootstrap"] = dict(ap={m: compare(boot["ap"][m], look["bootstrap"]["ap"][m]) for m in look_scorers},
                            diff={k: compare(boot["diff"][k], look["bootstrap"]["diff"][k]) for k in boot["diff"]})
    boot_ph = p1m.cluster_bootstrap(ys, {E.REFERENCE: sc.scores[E.REFERENCE], E.PRIMARY: sc.scores[E.PRIMARY]}, ca, cb,
                                    f=f, n_boot=1000, seed=0, diffs=[(E.REFERENCE, E.PRIMARY)])
    rep["bootstrap"]["posthoc_diff"] = compare(boot_ph["diff"][f"{E.REFERENCE}-{E.PRIMARY}"],
                                               ph["bootstrap"]["diff"][f"{E.REFERENCE}-{E.PRIMARY}"])
    # harness cluster-weight generator == Part 1 cluster_bootstrap (same stream, same replicates)
    fns = {m: E._ap_fn(ys, sc.scores[m]) for m in look_scorers}
    hb = E.B.run(E.B.cluster_weights(ca, cb, 1000, 0), w, fns,
                 [(f"{x}|ap", f"{y}|ap") for x, y in diffs])
    rep["bootstrap"]["harness_generator_vs_part1"] = {
        m: dict(ci_low_abs_diff=abs(hb["stat"][f"{m}|ap"]["ci_low"] - boot["ap"][m]["ci_low"]),
                ci_high_abs_diff=abs(hb["stat"][f"{m}|ap"]["ci_high"] - boot["ap"][m]["ci_high"])) for m in look_scorers}
    maxdiff = max(v["abs_diff"] for m in rep["metrics"].values() if m["vs_part1"] for v in m["vs_part1"].values())
    maxboot = max(v["abs_diff"] for grp in (rep["bootstrap"]["ap"], rep["bootstrap"]["diff"]) for x in grp.values()
                  for v in x.values())
    rep["max_abs_diff_point_metrics"] = maxdiff
    rep["max_abs_diff_bootstrap"] = maxboot
    tick("A: IFT look-1 reproduction", t0)
    rep["runtime_s"] = T["A: IFT look-1 reproduction"]
    dump(rep, out_dir / "A_ift_reproduction.json")
    summary["A"] = dict(ap_ht={m: dict(ours=rep["metrics"][m]["ours"]["ap_ht"],
                                       part1=(rep["metrics"][m]["vs_part1"] or {}).get("ap_ht", {}).get("part1"))
                               for m in E.METHODS},
                        max_abs_diff_point_metrics=maxdiff, max_abs_diff_bootstrap=maxboot)



def _part_b(assign, d, wfn_ift, ess_ift, n_boot, out_dir, summary):
    # ---------------------------------------------------------------- B: full protocol on IFT (fit on dev)
    t0 = time.perf_counter()
    wfn_dev = ht_weights_fn(assign, "dev")
    dev = E.load_set("Part 1 dev", d / "dev_positive_pairs.txt", d / "dev_negative_pairs.txt",
                     paths.SCORES / "dryrun_dev", weights_fn=wfn_dev)
    rules_dev = E.fit_rules(dev)
    tick("B: fit rules on dev pairs", t0)
    ift2 = E.load_set("Part 1 IFT", d / "ift_positive_pairs.txt", d / "ift_negative_pairs.txt",
                      paths.SCORES / "dryrun_ift", weights_fn=wfn_ift,
                      train_paths=(d / "dev_positive_pairs.txt", d / "dev_negative_pairs.txt"))
    resB, _ = run_protocol(ift2, rules_dev, ess_ift, E.SEEDS[1:], "B", n_boot)
    dump(dict(design="IFT pairs as the evaluation set (HT design weights 1/f for between-cluster and unmapped "
                     "pairs, f = %.4f), rules fitted on Part 1 dev pairs (HT weights, f_dev = %.4f); bootstrap units "
                     "= clusters rebuilt from IFT positives (protocol rule)" % (wfn_ift.f, wfn_dev.f),
              fitted_rules=rules_dev, audit=ift2.audit, results=resB), out_dir / "B_ift_protocol.json")
    figures(resB["itt"], out_dir / "figures", "dryrun_B_ift_itt")
    summary["B"] = {k: {m: resB["itt"]["metrics"][m][k] for m in E.METHODS} for k in ("ap", "auroc", "macro_cluster_ap")}
    summary["B"]["verdict"] = resB["itt"].get("verdict", {}).get("candidates")



def _part_c(X, genes, n_boot, out_dir, summary):
    # ---------------------------------------------------------------- C: test-mode code path on all train pairs
    t0 = time.perf_counter()
    tr = E.load_set("hu.MAP train (all pairs, in-sample)", paths.TRAIN_POS, paths.TRAIN_NEG, paths.SCORES / "train",
                    train_paths=(paths.TRAIN_POS, paths.TRAIN_NEG))
    rules_tr = E.fit_rules(tr)
    ess_tr = E.essential_flags(tr, X, genes)
    tick("C: load train + fit rules", t0)
    resC, clC = run_protocol(tr, rules_tr, ess_tr, E.SEEDS[1:], "C", n_boot)
    dump(dict(design="all train pairs at natural prevalence, weights 1 (test-mode code path); rules fitted on the "
                     "same pairs -> calibration / thresholds are IN-SAMPLE here (mechanics only); AP / bootstrap "
                     "width are a proxy for the expected test precision", audit=tr.audit, results=resC),
         out_dir / "C_train_as_test.json")
    figures(resC["itt"], out_dir / "figures", "dryrun_C_train_itt")
    summary["C"] = {k: {m: resC["itt"]["metrics"][m][k] for m in E.METHODS} for k in ("ap", "auroc", "macro_cluster_ap")}
    summary["C"]["delta_ap_ci"] = {k: v for k, v in resC["itt"]["bootstrap"]["cluster"]["diff"].items()
                                   if k.startswith(E.PRIMARY + "|ap")}
    summary["C"]["cluster_stats"] = resC["itt"]["cluster_stats"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["fit", "test", "dryrun"])
    ap.add_argument("--relook-reason", default=None, help="only for a logged, non-confirmatory re-run")
    ap.add_argument("--n-boot", type=int, default=E.N_BOOT)
    ap.add_argument("--parts", default="ABCI", help="dry run parts to run (A, B, C, I)")
    a = ap.parse_args()
    t0 = time.perf_counter()
    if a.mode == "fit":
        mode_fit()
    elif a.mode == "test":
        if a.n_boot != E.N_BOOT:
            raise SystemExit("the test look uses the pre-registered B = 2000")
        mode_test(a.relook_reason)
    else:
        mode_dryrun(a.n_boot, a.parts)
    print(f"total {round(time.perf_counter() - t0, 1)} s")


if __name__ == "__main__":
    main()
