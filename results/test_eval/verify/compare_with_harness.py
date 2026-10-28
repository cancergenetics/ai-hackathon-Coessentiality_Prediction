"""Compare verifier_results.json (committed before test_results.json was opened) with the harness's test_results.json.
Writes results/test_eval/verify/comparison.json. Tolerance 1e-9 (point metrics and CI endpoints), verdicts identical."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
V = json.loads((ROOT / "results/test_eval/verify/verifier_results.json").read_text())
H = json.loads((ROOT / "results/test_eval/test_results.json").read_text())
TOL = 1e-9
METHODS = ["zavg_BEF", "B", "gls_signed", "spearman_signed", "pearson_signed", "spec_k20_tinf_std"]
rows = []


def add(group, name, v, h):
    if isinstance(v, (bool, str)) or isinstance(h, (bool, str)) or v is None or h is None:
        ok = v == h
        d = None
    else:
        v, h = float(v), float(h)
        d = 0.0 if (np.isnan(v) and np.isnan(h)) else abs(v - h)
        ok = d <= TOL
    rows.append(dict(group=group, name=name, verifier=v, harness=h, abs_diff=d, ok=bool(ok)))


# audit
a, ha = V["audit"], H["audit"]
for k, hk in (("itt_pos", "analysed_pos"), ("itt_neg", "analysed_neg"), ("conflicts", "conflicts"),
              ("pos_duplicates", "pos_duplicate_lines"), ("neg_duplicates", "neg_duplicate_lines"),
              ("pos_self", "pos_self_lines"), ("neg_self", "neg_self_lines")):
    add("audit", k, a[k], ha[hk])
add("audit", "test_pairs_in_train", a["test_pairs_in_train"], sum(ha["train_overlap"].values()))
for m in METHODS:
    add("audit", f"{m} mapped pos", V["floors"][m]["n_mapped_pos"], ha["mapped_counts"][m]["pos"])
    add("audit", f"{m} mapped neg", V["floors"][m]["n_mapped_neg"], ha["mapped_counts"][m]["neg"])

KMAP = {"ap": "ap", "auroc": "auroc", "ap_no_largest": "ap_no_largest_component", "macro": "macro_cluster_ap",
        "f1_b": "f1_rule_b", "prec_b": "precision_rule_b", "rec_b": "recall_rule_b", "brier_iso": "brier_isotonic"}
CAL = {"brier": "brier", "bss_train_prev": "bss_vs_fit_prevalence", "bss_test_prev": "bss_vs_observed_prevalence",
       "citl": "calibration_in_the_large", "cal_intercept": "calibration_intercept", "cal_slope": "calibration_slope",
       "ece_equal_mass20": "ece_mass20", "ece_equal_width15": "ece_width15"}

for vs, hs in (("ITT", "itt"), ("S1", "mapped_only")):
    vv, hh = V["sets"][vs], H["results"][hs]
    g = f"{vs} clusters"
    vi, hc = vv["info"], hh["cluster_stats"]
    for k, hk in (("n_proteins", "n_proteins"), ("n_components", "n_components"), ("n_clusters", "n_clusters"),
                  ("n_singleton_clusters", "n_singleton_clusters"), ("components_split", "components_split"),
                  ("largest_component", "largest_component"), ("largest_component_size", "largest_component_size")):
        add(g, k, vi[k], hc[k])
    add(g, "n_pos", vi["n_pos"], hh["n_pos"])
    add(g, "n_neg", vi["n_neg"], hh["n_neg"])
    add(g, "bootstrap units", vv["bootstrap"]["n_clusters"], hh["bootstrap"]["cluster"]["n_units"])
    for m in METHODS:
        g = f"{vs} point"
        vm, hm = vv["metrics"][m], hh["metrics"][m]
        for k in ("ap", "auroc", "macro_cluster_ap", "n_clusters_macro", "ap_no_largest_component", "ap_no_top10_clusters",
                  "ap_hard_neg", "P@R0.1", "P@R0.25", "P@R0.5"):
            add(g, f"{m} {k}", vm[k], hm[k])
        add(g, f"{m} n_hard_neg", vm["n_hard_neg"], hm["part1_metric_suite"]["n_hard_neg"])
        dr = hm["decision_rules"]
        for k in ("threshold", "tp", "fp", "precision", "recall", "f1"):
            add(g, f"{m} rule_b {k}", vm["rule_b"][k], dr["rule_b_f1max_dev"][k])
        for k in ("threshold", "K", "tp", "fp", "precision", "recall", "f1"):
            add(g, f"{m} rule_c {k}", vm["rule_c"][k], dr["rule_c_topk"][k])
        if m == "zavg_BEF":
            for k in ("threshold", "tp", "fp", "precision", "recall", "f1"):
                add(g, f"{m} rule_a(D1) {k}", vm["rule_a_nonblind_D1"][k], dr["rule_a_z3"][k])
        for cal in ("isotonic", "platt"):
            for k, hk in CAL.items():
                add(f"{vs} calibration", f"{m} {cal} {k}", vm[cal][k], hm["calibration"][cal][hk])
        st = hh["bootstrap"]["cluster"]["stat"]
        for k, hk in KMAP.items():
            b = vv["bootstrap"]["per_method"][m][k]
            h = st[f"{m}|{hk}"]
            add(f"{vs} bootstrap", f"{m} {hk} point", b["point"], h["point"])
            add(f"{vs} bootstrap", f"{m} {hk} ci_low", b["ci"][0], h["ci_low"])
            add(f"{vs} bootstrap", f"{m} {hk} ci_high", b["ci"][1], h["ci_high"])
            add(f"{vs} bootstrap", f"{m} {hk} boot_sd", b["boot_sd"], h["boot_sd"])
    for hk, h in hh["bootstrap"]["cluster"]["diff"].items():
        left, right = hk.split(" - ")
        a_, ka = left.split("|")
        b_, _ = right.split("|")
        vk = {v_: k_ for k_, v_ in KMAP.items()}[ka]
        mine = vv["bootstrap"]["diff"].get(f"{vk}:{a_}-{b_}")
        if mine is None:
            rows.append(dict(group=f"{vs} bootstrap diff", name=hk, verifier=None, harness=h["diff"], abs_diff=None,
                             ok=False, note="not computed by verifier"))
            continue
        for k, vval in (("diff", mine["diff"]), ("ci_low", mine["ci"][0]), ("ci_high", mine["ci"][1]),
                        ("rel_gain", mine["rel_gain"]), ("rel_ci_low", mine["rel_ci"][0]), ("rel_ci_high", mine["rel_ci"][1]),
                        ("p_two_sided", mine["p_two_sided"])):
            add(f"{vs} bootstrap diff", f"{hk} {k}", vval, h[k])
    # §7.1 verdict
    lv = hh["locked_verdict"]
    for pk, p in vv["verdict"]["pairs"].items():
        add(f"{vs} verdict", f"{pk} passes", p["passes"], lv["pairs"][pk]["passes"])
        add(f"{vs} verdict", f"{pk} gate_ii", p["gate_ii"], lv["pairs"][pk]["gate_ii_ift"])
        add(f"{vs} verdict", f"{pk} rel_gain", p["rel_gain"], lv["pairs"][pk]["rel_gain"])
        for k, hk in (("macro_cluster_ap", "macro_cluster_ap"), ("ap_no_largest_component", "ap_ht_no_largest_component")):
            add(f"{vs} verdict", f"{pk} NI {k}", p["ni"][k], lv["pairs"][pk]["non_inferiority"][hk]["passed"])
    for c, cv in vv["verdict"]["candidates"].items():
        hc_ = lv["candidates"][c]
        add(f"{vs} verdict", f"{c} clearly_beats", cv["clearly_beats"], hc_["clearly_beats_ift"])
        add(f"{vs} verdict", f"{c} status", cv["status"], hc_["status"])
        add(f"{vs} verdict", f"{c} confirmed", cv["confirmed"], hc_["confirmed"])
    # §7.2 statistic
    cs = hh["ci_statistic"]
    for pk, p in vv["verdict"]["pairs"].items():
        c, b = pk.split("-", 1)
        hp = cs["pairs"][f"{c} vs {b}"]
        add(f"{vs} CI statistic", f"{pk} ci_excludes_0", p["ci_excludes_0"], hp["ci_excludes_0"])
        add(f"{vs} CI statistic", f"{pk} borderline", p["borderline"], hp["borderline"])
        for k, hk in (("macro", "macro_cluster_ap"), ("ap_no_largest", "ap_no_largest_component")):
            r = vv["bootstrap"]["ratio"][f"{k}:{c}/{b}"]
            hr = hp["non_degradation_ratio_ci"][hk]
            add(f"{vs} CI statistic", f"{pk} ratio {hk}", r["ratio"], hr["ratio"])
            add(f"{vs} CI statistic", f"{pk} ratio {hk} ci_low", r["ci"][0], hr["ci_low"])
            add(f"{vs} CI statistic", f"{pk} ratio {hk} ci_high", r["ci"][1], hr["ci_high"])
    for c in ("zavg_BEF", "B"):
        add(f"{vs} CI statistic", f"{c} ci_excludes_0 vs all", vv["verdict"]["candidates"][c]["ci_excludes_0_all"],
            cs["ci_excludes_0_vs_all_baselines"][c])
    for b, p in vv["verdict"]["holm_zavg_dAP"].items():
        add(f"{vs} CI statistic", f"Holm zavg_BEF vs {b}", p, cs["holm_adjusted_p_primary"][f"zavg_BEF vs {b}"])

# C strata point AP (ITT clusters/universe)
for k in ("C1", "C2", "C3"):
    hs = H["results"][f"stratum_{k}"]
    add("strata", f"{k} n_pos", V["sets"]["ITT"]["metrics"]["zavg_BEF"][f"n_pos_{k}"], hs["n_pos"])
    add("strata", f"{k} n_neg", V["sets"]["ITT"]["metrics"]["zavg_BEF"][f"n_neg_{k}"], hs["n_neg"])
    for m in METHODS:
        add("strata", f"{k} {m} ap", V["sets"]["ITT"]["metrics"][m][f"ap_{k}"], hs["metrics"][m]["ap"])
# summary
for s, vs in (("S1", "S1"), ("ITT", "ITT")):
    for c in ("zavg_BEF", "B"):
        add("summary", f"locked_verdict_{s} {c}", V["sets"][vs]["verdict"]["candidates"][c]["confirmed"] if c == "zavg_BEF"
            else V["sets"][vs]["verdict"]["candidates"][c]["clearly_beats"], H["summary"][f"locked_verdict_{s}"][c])

num = [r for r in rows if r["abs_diff"] is not None]
bad = [r for r in rows if not r["ok"]]
groups = {}
for r in rows:
    g = groups.setdefault(r["group"], dict(n=0, n_fail=0, max_abs_diff=0.0))
    g["n"] += 1
    g["n_fail"] += int(not r["ok"])
    if r["abs_diff"] is not None:
        g["max_abs_diff"] = max(g["max_abs_diff"], r["abs_diff"])
out = dict(tolerance=TOL, n_compared=len(rows), n_numeric=len(num), n_fail=len(bad),
           max_abs_diff=max(r["abs_diff"] for r in num), by_group=groups, failures=bad, rows=rows)
(ROOT / "results/test_eval/verify/comparison.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
