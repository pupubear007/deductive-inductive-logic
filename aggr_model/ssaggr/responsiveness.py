"""Objective 2 index: how much an isolate's pathogen transcriptome changes from host to host.

For isolate i at one time point, with log-normalised expression of shared genes,

    R2_host(i) = SS_between hosts / SS_total     (within isolate i, Euclidean, summed over genes)

which is the PERMANOVA R^2 of host within that isolate. Colonization (logit pathogen transcript
share) is regressed out per gene and host first, so a difference in fungal load cannot
masquerade as host-responsive expression.

Caveat for the current data (CH4 limitations): the soybean and sunflower experiments were run and
quantified separately, so host and batch are confounded. The absolute R2_host therefore includes
batch; only differences *between isolates* are interpretable.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np
import pandas as pd
from scipy import stats


def residualize(expr: pd.DataFrame, meta: pd.DataFrame, cov: str = "transcript_share") -> pd.DataFrame:
    """Per gene and host, remove the linear effect of logit(cov) across all inoculated libraries."""
    if cov not in meta.columns:
        return expr
    from .normalize import logit
    out = expr.copy()
    for h in meta["host"].unique():
        cols = meta.index[(meta["host"] == h)]
        c = logit(meta.loc[cols, cov].to_numpy(float))
        X = np.column_stack([np.ones_like(c), c - c.mean()])
        Y = expr[cols].to_numpy(float).T
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        out[cols] = (Y - X[:, 1:] @ beta[1:]).T
    return out


def host_responsiveness(expr: pd.DataFrame, meta: pd.DataFrame, fdr: float = 0.05) -> pd.DataFrame:
    rows = []
    for iso, g in meta.groupby("isolate"):
        hosts = g["host"].unique()
        if len(hosts) < 2:
            continue
        X = expr[g.index].to_numpy(float)  # genes x libraries
        grand = X.mean(1, keepdims=True)
        sst = float(((X - grand) ** 2).sum())
        ssb = 0.0
        for h in hosts:
            cols = (g["host"] == h).to_numpy()
            ssb += cols.sum() * float(((X[:, cols].mean(1, keepdims=True) - grand) ** 2).sum())
        n_de = 0
        if len(hosts) == 2:
            a = X[:, (g["host"] == hosts[0]).to_numpy()]
            b = X[:, (g["host"] == hosts[1]).to_numpy()]
            if a.shape[1] > 1 and b.shape[1] > 1:
                _, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
                p = np.nan_to_num(p, nan=1.0)
                n_de = int((stats.false_discovery_control(p) < fdr).sum())
        rows.append({"isolate": iso, "R2_host": ssb / sst if sst > 0 else np.nan, "n_host_DE_genes": n_de,
                     "n_libraries": int(len(g))})
    return pd.DataFrame(rows).set_index("isolate")


def exact_spearman(x: pd.Series, y: pd.Series, max_perm: int = 50000, seed: int = 0) -> dict[str, float]:
    """Spearman rho with an exact (or Monte Carlo) permutation p-value over isolates.
    With 6 isolates there are only 720 orderings, so p cannot go below 1/720."""
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    if n < 3:
        return {"rho": float("nan"), "p_perm": float("nan"), "n": n}
    rx, ry = d.iloc[:, 0].rank().to_numpy(), d.iloc[:, 1].rank().to_numpy()
    rho = float(np.corrcoef(rx, ry)[0, 1])
    import math
    if math.factorial(n) <= max_perm:
        perms = (np.array(p) for p in permutations(range(n)))
        total = math.factorial(n)
    else:
        rng = np.random.default_rng(seed)
        perms = (rng.permutation(n) for _ in range(max_perm))
        total = max_perm
    hits = sum(abs(np.corrcoef(rx, ry[p])[0, 1]) >= abs(rho) - 1e-12 for p in perms)
    return {"rho": rho, "p_perm": hits / total, "n": n}
