"""The FZ-R registry draft (docs/registry/*.yaml) against DESIGN v3.1 §2.2, §2.5-2.7 and the frozen FZ-LB counts.

Checks: every entry parses and validates; the expected set of entries exists; context IDs are columns of the frozen
context table summary and their n_A / n_L / lockbox_testable equal results/part2/lockbox/context_counts.tsv (H1 loss
contexts, added at FZ-X: results/part2/fzx/h1_contexts.tsv); the H3a and H3b families follow their inclusion rules;
the covariate-of-context rule; the control-to-direction kill map; signs; gene-only control hashes; placeholder /
pending bookkeeping, refusal of pending entries and fill-only edits. After FZ-X every placeholder is filled: the
filled values are checked against the frozen FZ-X lists (results/part2/fzx/; tests/test_part2_fzx.py re-derives them).
No dependency data is read."""
import copy
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part2 import registry  # noqa: E402
from part2.registry import RegistryError  # noqa: E402

REG = ROOT / "docs" / "registry"
COUNTS = pd.read_csv(ROOT / "results/part2/lockbox/context_counts.tsv", sep="\t", index_col=0)
SUMMARY = pd.read_csv(ROOT / "results/part2/contexts/context_summary.tsv", sep="\t", index_col=0)
FZX = ROOT / "results" / "part2" / "fzx"
FZX_MANIFEST = "results/part2/fzx/FZX_manifest.json"
FZC_MANIFEST = "results/part2/construction/FZC_manifest.json"
H1CTX = pd.read_csv(FZX / "h1_contexts.tsv", sep="\t", index_col=0)
CTRL_UNITS = pd.read_csv(FZX / "control_units.tsv", sep="\t", index_col="design_id", keep_default_na=False)
ZAVG_UNITS = (FZX / "units_zavg_coarse.txt").read_text().split()
ZAVG_ALL_UNITS = set(pd.read_csv(FZX / "zavg_engine_membership.tsv", sep="\t").unit_id)

PCS = [f"PC{i}" for i in range(1, 19)] + ["PC19a", "PC19b", "PC20", "PC21", "PC22"]
DIRECTION_IDS = ["H3a_lineage_atlas", "H3b_driver_family", "H1_within_complex_GI", "H1S_silencing_GI",
                 "H2_submodule_divergence",
                 "BENCH_complex_vs_gene", "BRIDGE-L_confirmation_rate"]
# DESIGN §2.6 control -> direction map (kill rule); 'none' = gene-level only
KILL_MAP = {"H3a": ["PC17", "PC20"],
            "H3b": ["PC1", "PC2", "PC4", "PC5", "PC9", "PC10", "PC11", "PC13", "PC18", "PC19a", "PC19b"],
            "H1": ["PC6", "PC7", "PC8", "PC12", "PC15"], "H2": ["PC21", "PC22"], "none": ["PC3", "PC14", "PC16"]}
# DESIGN §2.7 expected signs ('+' = more dependent in context) and the registered sign on the context column
DESIGN_SIGN = {**{p: "+" for p in PCS}, "PC3": "-", "PC5": "-", "PC14": "-"}
FLIPPED = {"PC2", "PC3"}       # DESIGN context is TP53 WT; the column TP53_altered codes altered = 1
H3B_DESIGN = ["hotspot:KRAS", "hotspot:NRAS", "hotspot:BRAF", "hotspot:PIK3CA", "TP53_altered", "MSI_H", "WGD",
              "9p21_loss", "MYC_amp", "MYCN_amp", "LoY", "sex_female", "NRF2_active", "LoF:ARID1A", "LoF:SMARCA4",
              "LoF:STAG2", "LoF:PTEN", "LoF:RB1", "LoF:KMT2D", "LoF:CREBBP", "fusion:EWSR1::FLI1|ERG",
              "fusion:KMT2A-r", "fusion:BCR::ABL1", "fusion:PAX3|7::FOXO1", "fusion:SS18::SSX"]
MUT_CONTEXT_FUNCTIONS = {"MSI_H": {"MSIScore", "log_burden"}, "FOCAD_del_or_MSI_H": {"MSIScore", "log_burden"},
                         "WGD": {"Ploidy", "Ploidy_anyprofile"}, "WGD_anyprofile": {"Ploidy", "Ploidy_anyprofile"},
                         "TP53_altered": {"TP53_altered"}}
LINEAGE_CONTEXTS = {"melanoma", "heme", "heme_vs_breast"}


@pytest.fixture(scope="module")
def entries():
    return registry.load_all(REG)


def primary(e):
    cs = e["context_family"]["contexts"]
    return [] if cs == "FZ-X" else [c for c in cs if c.get("role", "primary") == "primary"]


def all_ctx(e):
    cs = e["context_family"]["contexts"]
    return [] if cs == "FZ-X" else cs


def covs_of(e, c):
    return c.get("covariates", e["covariates"])


# ------------------------------------------------------------------ whole registry
def test_expected_entries(entries):
    assert set(entries) == {f"CTRL_{p}" for p in PCS} | set(DIRECTION_IDS)
    assert all(e["status"] == "open" for e in entries.values())


def test_files_are_ascii_and_parse_strictly():
    for p in REG.glob("*.yaml"):
        txt = p.read_text()
        assert txt.isascii(), p.name
        assert registry.validate(registry.load_yaml(txt), p.stem) == [], p.name


def test_common_fields(entries):
    for k, e in entries.items():
        assert e["success_criterion"].strip() and e["stopping_rule"].strip(), k
        assert isinstance(e.get("family_m"), dict) and e["family_m"]["formula"].strip(), k
        assert "value" in e["family_m"], k
        assert isinstance(e.get("multiverse"), dict) and e["multiverse"], k
        assert isinstance(e.get("lockbox"), dict) and str(e["lockbox"].get("target_list", "")).endswith(f"{k}.json"), k
        fi = e["frozen_inputs"]
        assert fi[:2] == ["results/part2/contexts/checksums.json", "results/part2/lockbox/manifest.json"], k
        assert FZX_MANIFEST in fi, k                               # every entry depends on the FZ-X freeze
        assert "FZ-X" not in fi and "FZ-C" not in fi, k            # filled at FZ-X
        gene_only = e["direction"] == "controls" and e["control"]["level"] == "G"
        assert fi[2:] == ([FZX_MANIFEST] if gene_only else [FZC_MANIFEST, FZX_MANIFEST]), k


def test_contexts_exist_and_counts_match_fzlb(entries):
    for k, e in entries.items():
        for c in all_ctx(e):
            cid = c["id"]
            if cid.startswith("H1loss:"):                          # FZ-X loss columns (not in the FZ-LB table)
                assert k in ("H1_within_complex_GI", "H1S_silencing_GI") and cid in H1CTX.index, (k, cid)
                r = H1CTX.loc[cid]
            else:
                assert cid in SUMMARY.index and cid in COUNTS.index, (k, cid)
                r = COUNTS.loc[cid]
            assert c["n_A"] == r.n_A and c["n_L"] == r.n_L, (k, cid)
            if "ref_A" in c:
                assert c["ref_A"] == r.ref_A and c["ref_L"] == r.ref_L, (k, cid)
            if "lockbox_testable" in c:
                assert c["lockbox_testable"] == bool(r.n_L >= 5 and r.ref_L >= 5), (k, cid)
            assert c.get("role", "primary") in ("primary", "sensitivity"), (k, cid)


def test_covariate_of_context_rule(entries):
    for k, e in entries.items():
        for c in all_ctx(e):
            cv = set(covs_of(e, c))
            assert not cv & MUT_CONTEXT_FUNCTIONS.get(c["id"], set()), (k, c["id"], cv)
            if c["id"].startswith("lineage:") or c["id"] in LINEAGE_CONTEXTS:
                assert "lineage_FE" not in cv, (k, c["id"])


def test_placeholders_are_pending_and_block_running(entries, tmp_path):
    for k, e in entries.items():
        ph = registry.placeholders(e)
        assert sorted(e.get("pending", [])) == sorted(ph), k
        if ph:
            with pytest.raises(RegistryError, match="placeholders"):
                registry.load_entry(REG, k, require_committed=False)
            assert registry.load_entry(REG, k, require_committed=False, allow_pending=True)[0]["id"] == k


# ------------------------------------------------------------------ controls
def test_controls_kill_map_signs_and_levels(entries):
    seen = {}
    for pc in PCS:
        e = entries[f"CTRL_{pc}"]
        assert e["direction"] == "controls" and e["control"]["design_id"] == pc
        kr = e["kill_rule"]
        seen.setdefault(kr["direction"], []).append(pc)
        assert kr["counted"] is (kr["direction"] != "none"), pc
        assert e["control"]["design_sign"] == DESIGN_SIGN[pc], pc
        flip = {"+": "-", "-": "+"}
        assert e["sign"] == (flip[DESIGN_SIGN[pc]] if pc in FLIPPED else DESIGN_SIGN[pc]), pc
        assert e["fdr_family"]["method"] == "none" and e["fdr_family"]["q"] == 0.05, pc
        gene_only = kr["direction"] == "none"
        assert (e["control"]["level"] == "G") is gene_only, pc
        if gene_only:
            assert e["score"] == "gene"
            assert e["complex_set_hash"] == registry.set_hash(e["control"]["gene_units"]), pc
            assert e["family_m"]["value"] == len(e["control"]["gene_units"]), pc
        else:
            assert "not evaluable at complex level" in e["control"]["complex_unit_rule"], pc
            assert kr["level"] == "C", pc
            # FZ-X fill: the frozen unit (or the pre-declared 'not evaluable' label), hash and m from control_units.tsv
            cu = e["control"]["complex_unit"]
            r = CTRL_UNITS.loc[pc]
            assert cu == r.unit, pc
            assert cu == "not evaluable at complex level" or cu in ZAVG_ALL_UNITS, pc
            units = e["control"]["gene_units"] + ([cu] if cu != "not evaluable at complex level" else [])
            assert e["complex_set_hash"] == registry.set_hash(units) == r.complex_set_hash, pc
            ev = str(r.complex_evaluable) == "True"
            assert e["family_m"]["value"] == len(e["control"]["gene_units"]) + int(ev), pc
        assert len(primary(e)) == 1, pc
    assert {d: sorted(v) for d, v in seen.items()} == {d: sorted(v) for d, v in KILL_MAP.items()}


def test_control_scores(entries):
    for pc in ("PC6", "PC7", "PC8", "PC12", "PC13"):
        assert entries[f"CTRL_{pc}"]["score"] == "member_resolved", pc
    for pc in ("PC21", "PC22"):
        assert entries[f"CTRL_{pc}"]["score"] == "submodule_vs_sibling", pc
        assert {"NNMD", "CasActivity"} <= set(entries[f"CTRL_{pc}"]["covariates"]), pc


# ------------------------------------------------------------------ directions
def test_h3a_family(entries):
    e = entries["H3a_lineage_atlas"]
    lin = [c for c in COUNTS.index if c.startswith("lineage:")]
    want = sorted([c for c in lin if COUNTS.loc[c, "n_A"] >= 12] + ["melanoma"])
    assert sorted(c["id"] for c in primary(e)) == want and len(want) == 23
    excl = {x["id"]: x["n_A"] for x in e["context_family"]["excluded"]}
    assert excl == {c: int(COUNTS.loc[c, "n_A"]) for c in lin if COUNTS.loc[c, "n_A"] < 12}
    assert e["family_m"]["n_families"] == 23
    assert e["kill_rule"]["controls"] == KILL_MAP["H3a"]
    assert e["fdr_family"] == {"method": "benjamini_bogomolov", "families": "contexts", "q": 0.05}
    assert e["complex_set_hash"] == registry.set_hash(ZAVG_UNITS) and e["score"] == "S_C_mean_winsor3"
    assert e["family_m"]["n_units"] == len(ZAVG_UNITS) == 367


def test_h3b_family(entries):
    e = entries["H3b_driver_family"]
    inc = [c["id"] for c in primary(e)]
    exc = {x["id"]: x["n_A"] for x in e["context_family"]["excluded"]}
    assert sorted(inc + list(exc)) == sorted(H3B_DESIGN)
    assert all(COUNTS.loc[c, "n_A"] >= 12 for c in inc) and all(COUNTS.loc[c, "n_A"] < 12 for c in exc)
    assert exc == {c: int(COUNTS.loc[c, "n_A"]) for c in exc}
    assert e["family_m"]["n_families"] == len(inc) == 20
    assert e["lockbox"]["testable_contexts"] == sum(COUNTS.loc[c, "n_L"] >= 5 for c in inc)
    assert e["kill_rule"]["controls"] == KILL_MAP["H3b"]
    for c in primary(e):                                      # DESIGN §2.5 covariate sets by context type
        cv = set(covs_of(e, c))
        assert {"lineage_FE", "library", "ScreenType3DO", "ModelType"} <= cv, c["id"]
        if c["id"].startswith(("hotspot:", "LoF:")) or c["id"] in ("TP53_altered", "NRF2_active"):
            assert {"log_burden", "MSIScore"} <= cv, c["id"]
        if c["id"] in ("9p21_loss", "MYC_amp", "MYCN_amp"):
            assert {"Ploidy", "WGD", "TP53_altered"} <= cv, c["id"]
        assert c.get("removed_members"), c["id"]
    for c in all_ctx(e):
        if c.get("role") == "sensitivity":
            assert c["id"] not in inc


@pytest.mark.parametrize("eid,loss_class,units_file", [
    ("H1_within_complex_GI", "genetic", "units_h1_components_A2.txt"),
    ("H1S_silencing_GI", "silencing", "units_h1s_components_A2.txt")])
def test_h1_and_h1s_event_split(entries, eid, loss_class, units_file):
    """FZ-X event list split by amendments A1 + A2: H1 = genetic events (>= 50% of A loss models with a LoF or hom-del
    call), H1S = silencing events; contexts, events, hash and m of each entry come from the frozen lists."""
    e = entries[eid]
    cf = e["context_family"]
    share = (H1CTX.n_A_LoF + H1CTX.n_A_homdel) / H1CTX.n_A
    cls = pd.Series(["genetic" if x >= 0.5 else "silencing" for x in share], index=H1CTX.index)
    want = [c for c in H1CTX.index if cls[c] == loss_class]
    assert [c["id"] for c in cf["contexts"]] == want
    ev = pd.read_csv(FZX / "h1_events.tsv", sep="\t")
    assert sorted(x for c in cf["contexts"] for x in c["events"]) == \
        sorted(ev.event[cls.reindex("H1loss:" + ev.gene).values == loss_class])
    assert all(c["loss_class"] == loss_class for c in cf["contexts"])
    for c in cf["contexts"]:
        r = H1CTX.loc[c["id"]]
        assert c["type"] == r.type and c["covariates"] == r.covariates.split(";"), c["id"]
        assert c["lockbox_testable"] == bool(c["n_L"] >= 5 and c["ref_L"] >= 5), c["id"]
    units = (FZX / units_file).read_text().split()
    assert e["complex_set_hash"] == registry.set_hash(units) and e["family_m"]["value"] == len(units)
    assert [a["id"] for a in e["amendments"]] == ["A1", "A2", "E1"]           # E1: score_rule erratum (R1)


def test_h1s_is_exploratory_copy_of_h1(entries):
    h, s = entries["H1_within_complex_GI"], entries["H1S_silencing_GI"]
    assert s["standing"] == "exploratory" and s["headline"] is False and "standing" not in h
    assert s["direction"] == "H1S" and h["direction"] == "H1" and s["status"] == "open"   # own label (A2)
    assert "THIS entry only" in s["fdr_rule"]
    for k in ("score", "score_rule", "covariates", "covariate_rule", "statistic", "classes", "matched_null",
              "fdr_family", "paralog_definition", "units", "seeds", "multiverse", "frozen_inputs",
              "implementation_note"):
        assert s[k] == h[k], k
    assert s["family_m"]["formula"] == h["family_m"]["formula"]
    assert "never enter any kill tally or main claim" in s["kill_rule"]["rule"]
    assert "not enter the kill rule" in s["success_criterion"] and "exploratory" in s["hypothesis"].lower()


def test_h1_entry(entries):
    e = entries["H1_within_complex_GI"]
    cf = e["context_family"]
    assert "enumeration_rule" in cf
    rule = cf["enumeration_rule"]
    for needle in (">= 12 loss", "functional check", "hom-del", "silencing", "attached", "Lockbox-testable"):
        assert needle in rule, needle
    assert "min(%id target, %id query) >= 20" in e["paralog_definition"]["rule"]
    for cls in ("PARALOG COMPENSATION", "RESIDUAL SENSITISATION", "COLLAPSE", "NEUTRAL", "INCONCLUSIVE"):
        assert cls in e["classes"], cls
    assert "[-0.1, 0.1]" in e["classes"] and "2.5%" in e["classes"]
    assert e["score"] == "member_resolved" and e["fdr_family"]["method"] == "bh"
    assert e["kill_rule"]["controls"] == KILL_MAP["H1"]
    assert "pending" not in e and registry.placeholders(e) == []


def test_h2_entry(entries):
    e = entries["H2_submodule_divergence"]
    h3 = [c["id"] for c in primary(entries["H3a_lineage_atlas"])] + \
         [c["id"] for c in primary(entries["H3b_driver_family"])]
    assert [c["id"] for c in primary(e)] == h3
    assert e["fdr_family"]["method"] == "stagewise_bb" and e["fdr_family"]["parent_directions"] == ["H3a", "H3b"]
    assert e["score"] == "submodule_vs_sibling"
    for c in primary(e):
        assert {"NNMD", "CasActivity"} <= set(covs_of(e, c)), c["id"]
    assert e["kill_rule"]["controls"] == KILL_MAP["H2"]


def test_bench_and_bridge_run_last(entries):
    b, g = entries["BENCH_complex_vs_gene"], entries["BRIDGE-L_confirmation_rate"]
    assert b["direction"] == "benchmark" and g["direction"] == "bridge"
    for e in (b, g):
        assert any("H1/H2" in r for r in e["run_after"])
    testable = [c["id"] for c in primary(entries["H3a_lineage_atlas"]) + primary(entries["H3b_driver_family"])
                if c["lockbox_testable"]]
    assert [c["id"] for c in primary(b)] == testable
    assert len(primary(g)) == 43


def test_every_kill_map_control_is_registered(entries):
    for d, pcs in KILL_MAP.items():
        for pc in pcs:
            assert entries[f"CTRL_{pc}"]["kill_rule"]["direction"] == d


# ------------------------------------------------------------------ schema extension
def _pending_version(e: dict, fields: dict) -> dict:
    """The FZ-R form of a filled entry: the given fields set back to their placeholders, listed in `pending`."""
    old = copy.deepcopy(e)
    for path, ph in fields.items():
        node, keys = old, re.findall(r"[^.\[\]]+|\[\d+\]", path)
        for kk in keys[:-1]:
            node = node[int(kk[1:-1])] if kk.startswith("[") else node[kk]
        last = keys[-1]
        if last.startswith("["):
            node[int(last[1:-1])] = ph
        else:
            node[last] = ph
    old["pending"] = sorted(fields)
    return old


H3A_PENDING = {"complex_set_hash": "FZ-C", "family_m.n_units": "FZ-C", "family_m.value": "FZ-X",
               "frozen_inputs[2]": "FZ-C", "frozen_inputs[3]": "FZ-X"}


def _base():
    """The H3a entry in its FZ-R (pending) form, rebuilt from the filled file."""
    return _pending_version(registry.load_yaml((REG / "H3a_lineage_atlas.yaml").read_text()), H3A_PENDING)


def test_pending_must_match_placeholders():
    e = _base()
    e2 = copy.deepcopy(e)
    e2["pending"] = e2["pending"][:-1]
    assert any("pending" in x for x in registry.validate(e2, e2["id"]))
    e3 = copy.deepcopy(e)
    e3["complex_set_hash"] = "FZ-Q"
    assert any("complex_set_hash" in x for x in registry.validate(e3, e3["id"]))
    e4 = copy.deepcopy(e)
    e4["status"] = "closed"
    assert any("closed" in x for x in registry.validate(e4, e4["id"]))
    entries_h1 = _pending_version(registry.load_yaml((REG / "H1_within_complex_GI.yaml").read_text()),
                                  {"context_family.contexts": "FZ-X", "complex_set_hash": "FZ-X",
                                   "family_m.value": "FZ-X", "frozen_inputs[2]": "FZ-C", "frozen_inputs[3]": "FZ-X"})
    assert registry.validate(entries_h1, entries_h1["id"]) == []
    e5 = copy.deepcopy(entries_h1)
    del e5["context_family"]["enumeration_rule"]
    assert any("enumeration_rule" in x for x in registry.validate(e5, entries_h1["id"]))


def test_fill_only_changes():
    old = _base()
    new = copy.deepcopy(old)
    new["complex_set_hash"] = "0123456789abcdef"
    new["family_m"]["n_units"] = 812
    new["family_m"]["value"] = 18000
    new["frozen_inputs"][2] = "results/part2/complexes/manifest.json"
    new["frozen_inputs"][3] = "results/part2/fzx/manifest.json"
    del new["pending"]
    assert registry.validate(new, new["id"]) == []
    assert registry.fill_only_changes(old, new) == []
    bad = copy.deepcopy(new)
    bad["fdr_family"]["q"] = 0.1
    assert registry.fill_only_changes(old, bad) == ["fdr_family.q"]
    bad2 = copy.deepcopy(new)
    bad2["context_family"]["contexts"] = bad2["context_family"]["contexts"][:-1]
    assert registry.fill_only_changes(old, bad2)
    # H1: the context list placeholder may be replaced by a list
    h1 = _pending_version(registry.load_yaml((REG / "H1_within_complex_GI.yaml").read_text()),
                          {"context_family.contexts": "FZ-X", "complex_set_hash": "FZ-X", "family_m.value": "FZ-X",
                           "frozen_inputs[2]": "FZ-C", "frozen_inputs[3]": "FZ-X"})
    assert registry.validate(h1, h1["id"]) == []
    h1n = copy.deepcopy(h1)
    h1n["context_family"]["contexts"] = [{"id": "H1loss:ARID1A", "definition": "toy", "sign": "two-sided"}]
    assert registry.fill_only_changes(h1, h1n) == []
