"""Evaluation metrics. All functions take numpy arrays and ignore masked-out entries."""

from __future__ import annotations

import numpy as np


def brier(p: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean((p - y) ** 2)) if len(y) else float("nan")


def log_loss(p: np.ndarray, y: np.ndarray, eps: float = 1e-6) -> float:
    p = np.clip(p, eps, 1 - eps)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p))) if len(y) else float("nan")


def auroc(p: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUROC; NaN if only one class is present."""
    pos, neg = p[y == 1], p[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(len(order))
    ranks[order] = np.arange(1, len(order) + 1)
    # average ranks for ties
    allv = np.concatenate([pos, neg])
    for v in np.unique(allv):
        idx = allv == v
        if idx.sum() > 1:
            ranks[idx] = ranks[idx].mean()
    return float((ranks[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def ece(p: np.ndarray, y: np.ndarray, bins: int = 10) -> float:
    """Expected calibration error with equal-width bins."""
    if not len(y):
        return float("nan")
    edges = np.linspace(0, 1, bins + 1)
    idx = np.clip(np.digitize(p, edges) - 1, 0, bins - 1)
    return float(sum(abs(p[idx == b].mean() - y[idx == b].mean()) * (idx == b).mean() for b in range(bins) if (idx == b).any()))


def sens_exceeds_fpr(p: np.ndarray, y: np.ndarray, threshold: float = 0.5) -> dict[str, float]:
    """Proposition 8.6: a positive forecast raises the risk exactly when sensitivity > FPR."""
    yhat = p >= threshold
    D, N = (y == 1).sum(), (y == 0).sum()
    tp, fp = (yhat & (y == 1)).sum(), (yhat & (y == 0)).sum()
    sens = tp / D if D else float("nan")
    fpr = fp / N if N else float("nan")
    return {"sensitivity": float(sens), "false_positive_rate": float(fpr),
            "positive_raises_risk": bool(D * fp < tp * N)}


def mae(p: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(p - y))) if len(y) else float("nan")
