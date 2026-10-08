"""Sparse hierarchical regression of aggressiveness on pathogen expression (weighted elastic net).

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

L1 penalties give exact zeros, so the non-zero genes are the
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


def _split(d: Design, kind: str, host_ratio: float):
    D, pen = build_matrix(d, kind)
    pen = np.where(pen < 0, host_ratio, pen)
    un = pen == 0
    return D[:, un], D[:, ~un], pen[~un], un


def _project_out(U: np.ndarray, w: np.ndarray, *arrs):
    """Weighted residuals after regressing on the unpenalised columns U (Frisch-Waugh-Lovell): the
    penalised problem on the residuals has exactly the same solution for the gene weights."""
    UtW = U.T * w
    G = np.linalg.pinv(UtW @ U)
    return [A - U @ (G @ (UtW @ A)) for A in arrs], G, UtW


def lambda_max(d: Design, l1_ratio: float = 0.5) -> float:
    """Smallest penalty at which every gene weight is zero (after the unpenalised terms)."""
    U, X, pen, _ = _split(d, "fixed", 1.0)
    (yr,), _, _ = _project_out(U, d.weight, d.y)
    return float(np.max(np.abs(X.T @ (d.weight * yr)))) / max(l1_ratio, 1e-3) + 1e-12


def fit_path(d: Design, lams: list[float], kind: str = "fixed", host_ratio: float = 2.0,
             l1_ratio: float = 0.5, tol: float = 1e-4, max_iter: int = 5000) -> list[np.ndarray]:
    """Weighted elastic net along a decreasing path of penalties (warm-started coordinate descent):

        minimise 0.5 * sum_n w_n (y_n - D_n b)^2 + lam * sum_j pen_j [l1_ratio |b_j| + (1 - l1_ratio)/2 b_j^2]

    over gene weights (pen_j = 1 shared, = host_ratio host-specific); group intercepts and the
    colonization slope are unpenalised and profiled out exactly. The per-column penalty factor is
    applied by rescaling columns (exact for the L1 part; the L2 part then scales with pen_j^2). The
    L2 part keeps co-expressed genes together instead of picking one of them at random.
    Returns one full coefficient vector per penalty, in the order given."""
    import warnings
    from sklearn.linear_model import enet_path
    U, X, pen, un = _split(d, kind, host_ratio)
    w = d.weight
    (Xr, yr), G, UtW = _project_out(U, w, X, d.y)
    sw = np.sqrt(w * len(w))  # sample weights folded into the rows (sklearn scales by 1/n)
    Xs = (Xr / pen[None, :]) * sw[:, None]
    ys = yr * sw
    order = np.argsort(lams)[::-1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, coefs, _ = enet_path(Xs, ys, l1_ratio=l1_ratio, alphas=np.asarray(lams, float)[order],
                                tol=tol, max_iter=max_iter)
    out: list[np.ndarray] = [None] * len(lams)  # type: ignore[list-item]
    for k, j in enumerate(order):
        b_pen = coefs[:, k] / pen
        beta = np.zeros(len(un))
        beta[un] = G @ (UtW @ (d.y - X @ b_pen))
        beta[~un] = b_pen
        out[j] = beta
    return out


def fit(d: Design, lam: float, kind: str = "fixed", host_ratio: float = 2.0, ridge: float = 1e-3,
        max_iter: int = 5000, tol: float = 1e-4, beta0: np.ndarray | None = None, l1_ratio: float = 0.5) -> np.ndarray:
    """Single penalty; see fit_path. ``ridge`` and ``beta0`` are kept for API compatibility."""
    return fit_path(d, [lam], kind, host_ratio, l1_ratio, tol, max(max_iter, 1000))[0]


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
