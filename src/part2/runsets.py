"""Model sets of the FZ-C construction runs (DESIGN v3.1 §2.1, §2.3). Every set is a subset of the construction set
A (results/part2/lockbox/A_models.txt, verified against the committed FZ-LB manifest by ``assoc.a_models``); no
other model list (all models, U, lockbox L) can be produced here.

  A            all 1,007 A models (reference block: cohesion, attachment, tiers)
  boot<seed>   PatientID-blocked 80% subsample of A, seeds 1000-1019 (consensus runs)
  half<seed>a/b  PatientID-blocked, OncotreeLineage-stratified 50/50 split of A, seeds 2000-2003 (split-half selection)
"""
from __future__ import annotations

from functools import lru_cache

BOOT_SEEDS = tuple(range(1000, 1020))
HALF_SEEDS = tuple(range(2000, 2004))
BOOT_FRAC = 0.8


@lru_cache(maxsize=1)
def A() -> tuple[str, ...]:
    from part2.assoc import a_models
    return tuple(sorted(a_models()))


def boot(seed: int) -> list[str]:
    from part2 import data
    return data.blocked_subsample(list(A()), BOOT_FRAC, seed)


def half(seed: int) -> tuple[list[str], list[str]]:
    from part2 import data
    ids = list(A())
    lin = data.load_models().reindex(ids)["OncotreeLineage"]
    return data.blocked_split(ids, 0.5, seed, strata=lin)


def run_names() -> list[str]:
    return ["A"] + [f"boot{s}" for s in BOOT_SEEDS] + [f"half{s}{h}" for s in HALF_SEEDS for h in "ab"]


def models(name: str) -> list[str]:
    if name == "A":
        return list(A())
    if name.startswith("boot"):
        s = int(name[4:])
        if s not in BOOT_SEEDS:
            raise ValueError(f"bootstrap seed {s} not in {BOOT_SEEDS}")
        return boot(s)
    if name.startswith("half"):
        s, h = int(name[4:-1]), name[-1]
        if s not in HALF_SEEDS or h not in "ab":
            raise ValueError(name)
        return half(s)[0 if h == "a" else 1]
    raise ValueError(f"unknown run set {name!r}; allowed: A, boot<1000-1019>, half<2000-2003>a|b")


def S_mask():
    """Dependency-bearing gene set S computed once on A (DESIGN v3.1 §2.3)."""
    from part2 import genefilter
    return genefilter.dependency_bearing_for_models(list(A()))
