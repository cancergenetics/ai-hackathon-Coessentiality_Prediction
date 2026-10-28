"""Frozen complex map as used by the association engine (DESIGN v3.1 §2.3-2.4; FZ-C).

The construction code (consensus.ClusterMap) holds clusters as dicts with local member indices over the gene set S.
At FZ-C the map is written as a membership table (one row per (unit, gene)) with columns

    unit_id   'C0007' (coarse, level 1) or 'C0007.S2' (sub-module, level 2)
    level     1 | 2
    parent    parent unit_id for level 2, '' for level 1
    gene      DepMap column key 'SYMBOL (ENTREZ)'
    role      'member' | 'attached' (attachment rule; attached genes define H1 events, never score members)
    locus     bool (locus flag of the unit)
    tier      'T1' | 'T2' | 'T3-neg' | 'T3-nc' | ''

``FrozenMap.from_cluster_map`` converts a ClusterMap + the S gene keys; ``to_membership`` / ``from_membership``
round-trip the table (TSV/parquet is the caller's choice; its SHA-256 goes into the FZ-C manifest).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


def symbol(key: str) -> str:
    """'SYMBOL (ENTREZ)' -> 'SYMBOL' (keys without ' (' are returned unchanged)."""
    return str(key).split(" (")[0]


@dataclass
class FrozenMap:
    members: dict                                  # unit_id -> list of gene keys (fixed members)
    level: dict                                    # unit_id -> 1 | 2
    parent: dict = field(default_factory=dict)     # unit_id -> parent unit_id ('' for coarse)
    attached: dict = field(default_factory=dict)   # coarse unit_id -> list of attached gene keys
    locus: dict = field(default_factory=dict)
    tier: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- views
    def coarse(self) -> list[str]:
        return sorted(u for u, l in self.level.items() if l == 1)

    def submodules(self, parent: str | None = None) -> list[str]:
        return sorted(u for u, l in self.level.items() if l == 2 and (parent is None or self.parent.get(u) == parent))

    def units_of_gene(self, gene: str, level: int = 1) -> list[str]:
        return sorted(u for u, m in self.members.items() if self.level[u] == level and gene in m)

    def units_of_symbol(self, sym: str, level: int = 1) -> list[str]:
        return sorted(u for u, m in self.members.items() if self.level[u] == level and sym in {symbol(g) for g in m})

    def sibling_complement(self, sub: str) -> list[str]:
        """Members of the parent not in the sub-module (incl. parent members in no sub-module)."""
        p = self.parent[sub]
        s = set(self.members[sub])
        return [g for g in self.members[p] if g not in s]

    def all_genes(self, level: int = 1) -> set:
        return {g for u, m in self.members.items() if self.level[u] == level for g in m}

    # ---------------------------------------------------------------- io
    def to_membership(self) -> pd.DataFrame:
        rows = []
        for u in sorted(self.members):
            for g in self.members[u]:
                rows.append(dict(unit_id=u, level=self.level[u], parent=self.parent.get(u, ""), gene=g, role="member",
                                 locus=bool(self.locus.get(u, False)), tier=self.tier.get(u, "")))
            for g in self.attached.get(u, []):
                rows.append(dict(unit_id=u, level=self.level[u], parent=self.parent.get(u, ""), gene=g,
                                 role="attached", locus=bool(self.locus.get(u, False)), tier=self.tier.get(u, "")))
        return pd.DataFrame(rows, columns=["unit_id", "level", "parent", "gene", "role", "locus", "tier"])

    @classmethod
    def from_membership(cls, df: pd.DataFrame) -> "FrozenMap":
        need = {"unit_id", "level", "gene"}
        if not need <= set(df.columns):
            raise ValueError(f"membership table needs columns {sorted(need)}")
        fm = cls(members={}, level={})
        for u, g in df.groupby("unit_id", sort=True):
            lev = int(g.level.iloc[0])
            role = g.role if "role" in g else pd.Series("member", index=g.index)
            fm.members[u] = list(g.gene[role == "member"])
            fm.level[u] = lev
            par = g.parent.iloc[0] if "parent" in g else ""
            fm.parent[u] = "" if pd.isna(par) else str(par)
            att = list(g.gene[role == "attached"])
            if att:
                fm.attached[u] = att
            if "locus" in g:
                fm.locus[u] = bool(g.locus.iloc[0])
            if "tier" in g:
                fm.tier[u] = "" if pd.isna(g.tier.iloc[0]) else str(g.tier.iloc[0])
        bad = [u for u in fm.members if fm.level[u] == 2 and fm.parent[u] not in fm.members]
        if bad:
            raise ValueError(f"sub-modules without a parent unit: {bad[:3]}")
        return fm

    @classmethod
    def from_cluster_map(cls, cm, gene_keys, locus: dict | None = None, tier: dict | None = None) -> "FrozenMap":
        """consensus.ClusterMap (local member indices) + gene keys of S -> FrozenMap. IDs: C%04d / C%04d.S%d."""
        keys = list(gene_keys)
        fm = cls(members={}, level={})
        cid = {}
        for c in cm.level(1):
            u = f"C{int(c['id']):04d}"
            cid[c["id"]] = u
            fm.members[u] = [keys[i] for i in sorted(c["members"])]
            fm.level[u], fm.parent[u] = 1, ""
        k = {}
        for c in cm.level(2):
            p = cid[c["parent"]]
            k[p] = k.get(p, 0) + 1
            u = f"{p}.S{k[p]}"
            fm.members[u] = [keys[i] for i in sorted(c["members"])]
            fm.level[u], fm.parent[u] = 2, p
        fm.locus = dict(locus or {})
        fm.tier = dict(tier or {})
        return fm
