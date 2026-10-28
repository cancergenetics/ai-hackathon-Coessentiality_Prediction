"""Toy simulation with known answers for the construction pipeline (DESIGN v2 §3; round-1 algorithms 'Toy').

Lines x genes gene-effect matrix with
  * planted complexes of sizes 3-40 driven by a selective dependency factor (few lines strongly dependent),
    the 40-gene complex made of two sub-modules (25 + 15) with their own factors (nesting),
    one shared subunit (member of two complexes), one pan-essential complex (core machine, mean -1.5),
  * pan-essential singletons, selective singletons (dependency-bearing, no partner), a few hub genes,
  * never-essential genes (60%) carrying chromosome-arm block noise (CN-like), weaker arm noise on all genes,
  * lineage structure of rank n_lineages (gene loadings on lineage indicators) and a library batch (rank 1).
Scenarios: 'planted' (all of the above), 'null' (no complexes, hubs or shared factors), 'armblock' (no complexes;
strong arm-block noise on a set of dependency-bearing genes -> co-dependency that is pure locus artefact),
'hier' (FZ-C decision D2: 'planted' plus three cohesive parents of 3 x 8-gene sub-modules with a dominant parent
factor, and two pairs of 10-gene complexes sharing 3 subunits; the extra truth is appended to complexes /
submodules; the 'planted' truth sets are identical, later random draws differ).
Everything is a deterministic function of (scenario, seed, sizes).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SIZES = (3, 3, 4, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40)


@dataclass
class Toy:
    X: np.ndarray                     # lines x genes float32
    complexes: list                   # list of gene-index arrays (level-1 truth)
    submodules: list                  # list of gene-index arrays (level-2 truth, inside the 40-complex)
    arm: np.ndarray                   # str per gene, e.g. '7p'
    band: np.ndarray                  # str per gene, e.g. '7p2'
    gene_class: np.ndarray            # 'never', 'pan', 'complex', 'selective', 'hub', 'armdep'
    lineage: np.ndarray               # per line
    info: dict = field(default_factory=dict)


def simulate(scenario: str = "planted", n_lines: int = 300, n_genes: int = 3000, seed: int = 0,
             sizes=SIZES, n_lineages: int = 6, n_arms: int = 40, bands_per_arm: int = 5,
             frac_never: float = 0.6, frac_pan: float = 0.08, noise: float = 0.2, arm_sd: float = 0.12,
             factor_strength: float = 1.0, p_dependent: float = 0.15) -> Toy:
    rng = np.random.default_rng(seed)
    L, G = n_lines, n_genes
    lineage = rng.integers(0, n_lineages, L)
    library = (rng.random(L) < 0.3).astype(float)
    arm_id = rng.integers(0, n_arms, G)
    band_id = rng.integers(0, bands_per_arm, G)
    arm = np.array([f"{a // 2 + 1}{'pq'[a % 2]}" for a in arm_id])
    band = np.array([f"{arm[g]}{band_id[g] + 1}" for g in range(G)])
    arm_factor = rng.standard_normal((L, n_arms))            # CN-like per line x arm

    X = noise * rng.standard_normal((L, G))
    cls = np.full(G, "selective", dtype=object)
    perm = rng.permutation(G)
    n_never, n_pan = int(frac_never * G), int(frac_pan * G)
    never, pan = perm[:n_never], perm[n_never:n_never + n_pan]
    rest = perm[n_never + n_pan:]
    cls[never], cls[pan] = "never", "pan"
    X[:, pan] += -1.5
    X[:, never] *= 0.5
    # arm-block noise: strong on never-essential genes, weak on everything
    load = np.where(cls == "never", arm_sd, 0.3 * arm_sd) * rng.uniform(0.5, 1.0, G)
    X += arm_factor[:, arm_id] * load
    # lineage (rank n_lineages) + library batch (rank 1)
    lin = np.eye(n_lineages)[lineage]
    X += lin @ (0.15 * rng.standard_normal((n_lineages, G)) * np.where(cls == "never", 0.3, 1.0))
    X += np.outer(library, 0.05 * rng.standard_normal(G))

    def selective_factor():
        dep = rng.random(L) < p_dependent
        return factor_strength * (-1.0 * dep + 0.3 * rng.standard_normal(L))

    complexes, subs = [], []
    pos = 0
    if scenario in ("planted", "hier"):
        need = sum(sizes)
        pool = list(rest[: need + 10])
        for si, sz in enumerate(sizes):
            members = np.array(pool[pos:pos + sz])
            pos += sz
            if sz == 30:                                    # pan-essential core machine
                f = 0.4 * rng.standard_normal(L)
                X[:, members] += -1.5 - X[:, members].mean(0) * 0 + np.outer(f, rng.uniform(0.7, 1.0, sz))
                cls[members] = "complex"
            elif sz == 40:                                  # two sub-modules + shared parent factor
                fp = selective_factor()
                for part in (members[:25], members[25:]):
                    fs = 0.45 * selective_factor()
                    X[:, part] += -0.3 + np.outer(fp, rng.uniform(0.7, 1.0, len(part))) + \
                        np.outer(fs, rng.uniform(0.6, 0.9, len(part)))
                    subs.append(np.sort(part))
                cls[members] = "complex"
            else:
                f = selective_factor()
                X[:, members] += -0.3 + np.outer(f, rng.uniform(0.6, 1.0, sz))
                cls[members] = "complex"
            complexes.append(np.sort(members))
        # shared subunit: last member of the size-6 complex also loads on the size-8 complex's factor
        c6, c8 = complexes[sizes.index(6)], complexes[sizes.index(8)]
        shared = int(c6[-1])
        f8 = X[:, c8].mean(1) - X[:, c8].mean()
        X[:, shared] += 0.6 * f8
        complexes[sizes.index(8)] = np.sort(np.r_[c8, shared])
        # hubs: weak loading on three complex factors each
        hubs = np.array(pool[pos:pos + 3])
        for h in hubs:
            for ci in rng.choice(len(complexes), 3, replace=False):
                cm = complexes[ci]
                X[:, h] += 0.3 * (X[:, cm].mean(1) - X[:, cm].mean())
        cls[hubs] = "hub"
        if scenario == "hier":
            extra = rest[need + 10:]
            q = 0
            for _ in range(3):                              # cohesive parents with 3 sub-modules each
                members = np.array(extra[q:q + 24])
                q += 24
                fp = selective_factor()
                for i in range(3):
                    part = members[8 * i: 8 * i + 8]
                    fs = 0.4 * selective_factor()
                    X[:, part] += -0.3 + np.outer(fp, rng.uniform(0.75, 1.0, 8)) + np.outer(fs, rng.uniform(0.6, 0.9, 8))
                    subs.append(np.sort(part))
                cls[members] = "complex"
                complexes.append(np.sort(members))
            for _ in range(2):                              # overlapping complexes: 10 + 10 sharing 3 subunits
                members = np.array(extra[q:q + 17])
                q += 17
                a, b = members[:10], members[7:]
                fa, fb = selective_factor(), selective_factor()
                X[:, members[:7]] += -0.3 + np.outer(fa, rng.uniform(0.6, 1.0, 7))
                X[:, members[10:]] += -0.3 + np.outer(fb, rng.uniform(0.6, 1.0, 7))
                X[:, members[7:10]] += -0.3 + np.outer(0.7 * fa + 0.7 * fb, rng.uniform(0.8, 1.0, 3))
                cls[members] = "complex"
                complexes.extend([np.sort(a), np.sort(b)])
        sel_rest = np.array([g for g in rest if cls[g] == "selective"])
    else:
        sel_rest = rest
    # selective singletons: private dependency factor
    for g in sel_rest:
        X[:, g] += -0.3 + selective_factor() * rng.uniform(0.6, 1.0)
    if scenario == "armblock":
        # strong CN-like arm blocks among dependency-bearing genes: groups of genes on the same arm
        arms_used = rng.choice(n_arms, 8, replace=False)
        for a in arms_used:
            g = sel_rest[arm_id[sel_rest] == a]
            X[:, g] += np.outer(arm_factor[:, a], rng.uniform(0.6, 0.9, len(g)))
            cls[g] = "armdep"
    elif scenario not in ("planted", "null", "hier"):
        raise ValueError(scenario)
    return Toy(X.astype(np.float32), complexes, subs, arm, band, cls.astype(str), lineage,
               info=dict(scenario=scenario, seed=seed, n_lines=L, n_genes=G))


# ------------------------------------------------------------------------------------------ evaluation
def adjusted_rand(truth_sets, found_labels: np.ndarray, genes: np.ndarray) -> float:
    """ARI over `genes` (planted-complex genes): truth = complex id, found = cluster id or singleton."""
    from sklearn.metrics import adjusted_rand_score
    t = np.full(len(found_labels), -1)
    for i, s in enumerate(truth_sets):
        t[s] = i
    f = found_labels.copy()
    nxt = f.max() + 1 if len(f) else 0
    for g in np.flatnonzero(f < 0):
        f[g] = nxt
        nxt += 1
    return float(adjusted_rand_score(t[genes], f[genes]))


def best_match_f1(truth_sets, found_sets) -> np.ndarray:
    out = []
    for s in truth_sets:
        best = 0.0
        for c in found_sets:
            i = len(np.intersect1d(s, c))
            if i:
                best = max(best, 2 * i / (len(s) + len(c)))
        out.append(best)
    return np.array(out)
