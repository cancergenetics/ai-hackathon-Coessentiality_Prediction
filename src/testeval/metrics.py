"""Point metrics for the test evaluation (docs/TEST_PROTOCOL.md §4-5).

AP, HT/plain weighted AP, precision at recall, macro per-cluster AP and the Part 1 metric suite are Part 1's own
functions (``testeval.part1.metrics``). Added here, each tested against sklearn / Part 1 on toys:
  * ``Ranked``: scores sorted once (descending, mergesort) with sklearn tie blocks (Part 1 ``PresortedScores``
    ordering) giving O(n) weighted AP (identical formula to Part 1), AUROC (ties count 1/2, = sklearn
    roc_auc_score with sample_weight) and confusion counts at any threshold, for fast bootstraps.
  * decision rules: F1-maximising threshold (ties -> the highest such threshold), fixed thresholds, top-k.
Thresholds always mean "predicted positive iff score >= t".
"""
from __future__ import annotations

import numpy as np

from testeval.part1 import metrics as p1m


class Ranked:
    def __init__(self, y, s):
        y = np.asarray(y).astype(bool)
        s = np.asarray(s, dtype=np.float64)
        if np.isnan(s).any():
            raise ValueError("NaN scores")
        self.order = np.argsort(-s, kind="mergesort")
        self.ys = y[self.order]
        self.ss = s[self.order]
        self.block_end = np.flatnonzero(np.r_[self.ss[1:] != self.ss[:-1], True])
        self.block_score = self.ss[self.block_end]
        self._yf = self.ys.astype(np.float64)
        self._nf = 1.0 - self._yf

    def _cum(self, w):
        if w is None:
            return np.cumsum(self._yf)[self.block_end], np.cumsum(self._nf)[self.block_end]
        ws = np.asarray(w, dtype=np.float64)[self.order]
        return np.cumsum(ws * self._yf)[self.block_end], np.cumsum(ws * self._nf)[self.block_end]

    @staticmethod
    def _ap(tp, fp) -> float:
        P = tp[-1]
        if P <= 0:
            return np.nan
        prec = tp / np.maximum(tp + fp, 1e-300)
        drec = np.diff(np.r_[0.0, tp / P])
        return float(np.sum(drec * prec))

    @staticmethod
    def _auroc(tp, fp) -> float:
        P, N = tp[-1], fp[-1]
        if P <= 0 or N <= 0:
            return np.nan
        dtp = np.diff(np.r_[0.0, tp])
        dfp = np.diff(np.r_[0.0, fp])
        fp_above = np.r_[0.0, fp[:-1]]
        # a positive in a block is ranked above all negatives below the block, ties count 1/2
        wins = np.sum(dtp * (N - fp_above - dfp) + 0.5 * dtp * dfp)
        return float(wins / (P * N))

    def _confusion(self, tp, fp, t) -> dict:
        P, N = tp[-1], fp[-1]
        k = np.searchsorted(-self.block_score, -t, side="right") - 1      # last block with score >= t
        TP, FP = (float(tp[k]), float(fp[k])) if k >= 0 else (0.0, 0.0)
        return _prf(TP, FP, P, N)

    def ap(self, w=None) -> float:
        return self._ap(*self._cum(w))

    def auroc(self, w=None) -> float:
        return self._auroc(*self._cum(w))

    def confusion(self, t: float, w=None) -> dict:
        return self._confusion(*self._cum(w), t)

    def summary(self, w=None, t: float | None = None) -> dict:
        """AP, AUROC and (if t is given) the confusion at t, from one pass over the sorted weights."""
        tp, fp = self._cum(w)
        out = dict(ap=self._ap(tp, fp), auroc=self._auroc(tp, fp))
        if t is not None:
            out["confusion"] = self._confusion(tp, fp, t)
        return out

    def f1_max(self, w=None) -> dict:
        tp, fp = self._cum(w)
        P = tp[-1]
        f1 = 2 * tp / np.maximum(tp + fp + P, 1e-300)
        k = int(np.argmax(f1))                      # first maximum = highest threshold attaining it
        return dict(threshold=float(self.block_score[k]), **_prf(float(tp[k]), float(fp[k]), P, fp[-1]))

    def kth_score(self, k: int) -> float:
        """Score of the k-th ranked pair (1-based); thresholding at it predicts >= k pairs (whole tie block)."""
        k = int(min(max(k, 1), len(self.ss)))
        return float(self.ss[k - 1])


def _prf(TP, FP, P, N) -> dict:
    prec = TP / (TP + FP) if TP + FP > 0 else float("nan")
    rec = TP / P if P > 0 else float("nan")
    f1 = 2 * TP / (TP + FP + P) if TP + FP + P > 0 else float("nan")
    return dict(tp=TP, fp=FP, fn=float(P - TP), tn=float(N - FP), n_pred_pos=TP + FP, precision=prec, recall=rec,
                f1=f1)


def ap(y, s, w=None) -> float:
    """Part 1 weighted_ap (sklearn average_precision_score with sample_weight)."""
    return p1m.weighted_ap(y, s, w)


def auroc(y, s, w=None) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s, sample_weight=w))


def top_cluster_ids(y, ca, cb, n=10) -> list:
    """Part 1 scripts/21_phaseB_v2.top_clusters: the n clusters with most within-cluster positives (stable order)."""
    y = np.asarray(y).astype(bool)
    ca, cb = np.asarray(ca), np.asarray(cb)
    c, n_pos = np.unique(ca[y & (ca == cb)], return_counts=True)
    return [int(x) for x in c[np.argsort(-n_pos, kind="stable")][:n]]
