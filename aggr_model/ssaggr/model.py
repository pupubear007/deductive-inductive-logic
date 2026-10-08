"""Sparse hierarchical regression of aggressiveness on pathogen expression (trained by FISTA).

For library n of isolate i on host h (time point t),

    y_ih  ~  a_{h,t} + g_{h,t} * c_n  +  x_n . w  [ + x_n . v_h ]

* y_ih is the isolate's aggressiveness on host h (log sAUDPC, centred within host): one value
  per isolate and host, shared by that isolate's libraries.
* c_n is the colonization covariate (logit pathogen transcript share), unpenalised.
* w are shared gene weights: the **fixed program** hypothesis (same genes, same direction on
  every host; a more aggressive isolate runs the program harder).
* v_h are host-specific deviations: the **host-responsive** hypothesis (the gene-to-aggressiveness
  map is retuned on each host). They carry a larger L1 penalty, so they are used only where the
  data demand it.

L1 penalties give exact zeros (FISTA with soft-thresholding), so the non-zero genes are the
candidate aggressiveness determinants. Library weights make every (isolate, host) cell count
equally, so isolates with more libraries do not dominate: the isolate, not the library, is the
unit of evidence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Design:
    X: np.ndarray          # n x G standardized expression
    group: np.ndarray      # n, intercept group (host x time point) index
    host: np.ndarray       # n, host index
    cov: np.ndarray        # n, colonization covariate (0 if unused)
    y: np.ndarray          # n, target
    weight: np.ndarray     # n, library weights (sum to 1)
    n_groups: int
    n_hosts: int


def build_matrix(d: Design, kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Augmented design [group intercepts | group covariate slopes | X | X * 1{host=h} ...]
    and a per-column penalty multiplier (0 = unpenalised, 1 = shared, ratio = host-specific)."""
    n, G = d.X.shape
    oh = np.zeros((n, d.n_groups))
    oh[np.arange(n), d.group] = 1.0
    blocks = [oh, oh * d.cov[:, None], d.X]
    pen = [np.zeros(d.n_groups), np.zeros(d.n_groups), np.ones(G)]
    if kind == "responsive" and d.n_hosts > 1:
        for h in range(d.n_hosts):
            blocks.append(d.X * (d.host == h)[:, None])
            pen.append(np.full(G, -1.0))  # placeholder, replaced by ratio in fit()
    return np.hstack(blocks), np.concatenate(pen)


def lambda_max(d: Design) -> float:
    """Smallest L1 weight that sets every gene weight to zero (after fitting intercepts)."""
    r = d.y.copy()
    for g in range(d.n_groups):
        m = d.group == g
        if m.any():
            r[m] -= np.average(d.y[m], weights=d.weight[m])
    return float(np.max(np.abs(d.X.T @ (d.weight * r)))) + 1e-12


def fit(d: Design, lam: float, kind: str = "fixed", host_ratio: float = 2.0, ridge: float = 1e-3,
        max_iter: int = 400, tol: float = 1e-6, beta0: np.ndarray | None = None, l1_ratio: float = 0.5) -> np.ndarray:
    """Elastic net: minimise 0.5 * sum_n w_n (y_n - D_n b)^2
         + lam * sum_j pen_j [ l1_ratio |b_j| + (1 - l1_ratio)/2 b_j^2 ] + ridge/2 |b_pen|^2
    by FISTA with soft-thresholding (exact zeros for unselected genes). The L2 part keeps groups of
    co-expressed genes together instead of picking one of them at random, which matters because
    aggressiveness-associated genes are strongly co-regulated."""
    D, pen = build_matrix(d, kind)
    pen = np.where(pen < 0, host_ratio, pen)
    w, y = d.weight, d.y
    sw = np.sqrt(w)[:, None]
    l2 = ridge + lam * (1.0 - l1_ratio) * pen
    L = float(np.linalg.norm(D * sw, 2) ** 2) + float(l2.max()) + 1e-12
    step = 1.0 / L
    beta = np.zeros(D.shape[1]) if beta0 is None or beta0.shape[0] != D.shape[1] else beta0.copy()
    z, t = beta.copy(), 1.0
    thr = step * lam * l1_ratio * pen
    rmask = (pen > 0).astype(float) * l2
    for _ in range(max_iter):
        grad = D.T @ (w * (D @ z - y)) + rmask * z
        b_new = z - step * grad
        b_new = np.sign(b_new) * np.maximum(np.abs(b_new) - thr, 0.0)
        t_new = (1 + (1 + 4 * t * t) ** 0.5) / 2
        z = b_new + ((t - 1) / t_new) * (b_new - beta)
        if np.max(np.abs(b_new - beta)) < tol:
            beta = b_new
            break
        beta, t = b_new, t_new
    return beta


def predict(d: Design, beta: np.ndarray, kind: str) -> np.ndarray:
    D, _ = build_matrix(d, kind)
    if D.shape[1] != beta.shape[0]:
        raise ValueError("coefficient vector does not match the design")
    return D @ beta


def gene_weights(beta: np.ndarray, d: Design, kind: str) -> dict[str, np.ndarray]:
    """Split the coefficient vector: shared w and, for the responsive model, v_h per host."""
    G = d.X.shape[1]
    off = 2 * d.n_groups
    out = {"w": beta[off:off + G]}
    if kind == "responsive" and d.n_hosts > 1:
        for h in range(d.n_hosts):
            out[f"v{h}"] = beta[off + G * (h + 1): off + G * (h + 2)]
    return out
