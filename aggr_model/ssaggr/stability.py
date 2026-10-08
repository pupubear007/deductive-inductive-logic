"""Objective 1: aggressiveness and its consistency across hosts.

Joint regression of isolate aggressiveness on a host index (Finlay and Wilkinson 1963;
Eberhart and Russell 1966): for isolate i on host h,

    y_ih = mu_i + b_i * I_h + d_ih,      I_h = mean_i y_ih - grand mean,

b_i measures how strongly the isolate tracks host susceptibility and the deviation mean
square s2d_i = sum_h d_ih^2 / (H - 2) measures host-specific departures.

Note on scaling. If y is z-scored *within* each host, every host mean is 0 and I_h vanishes,
so the regression is undefined. The default therefore uses log(sAUDPC), which keeps host
differences in I_h and puts isolates on a ratio scale; the rank-based summaries used in CH4
(mean rank, rank range) are reported alongside.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def stability_table(pheno: pd.DataFrame) -> pd.DataFrame:
    """pheno: isolates x crops, mean sAUDPC. Returns one row per isolate."""
    y = np.log(pheno.clip(lower=1e-6))
    H = y.shape[1]
    I = (y.mean(0) - y.values.mean()).to_numpy()
    rows = []
    for iso, r in y.iterrows():
        v = r.to_numpy()
        mu = v.mean()
        b = float(np.dot(I, v - mu) / np.dot(I, I)) if np.dot(I, I) > 0 else float("nan")
        resid = v - mu - b * I
        s2d = float((resid ** 2).sum() / (H - 2)) if H > 2 else float("nan")
        rows.append({"isolate": iso, "mean_log_sAUDPC": mu, "slope_b": b, "deviation_ms": s2d})
    out = pd.DataFrame(rows).set_index("isolate")
    ranks = pheno.rank(axis=0)
    out["mean_rank"] = ranks.mean(axis=1)
    out["rank_range"] = ranks.max(axis=1) - ranks.min(axis=1)
    return out


def classify(st: pd.DataFrame, high_q: float = 0.5, dev_q: float = 0.5) -> pd.Series:
    """Provisional classes. The cut-offs are *decisions*, not results: set them before looking
    at expression data (proposal 5.8) and record them in the config."""
    hi = st["mean_log_sAUDPC"] >= st["mean_log_sAUDPC"].quantile(high_q)
    stable = st["deviation_ms"] <= st["deviation_ms"].quantile(dev_q)
    lab = np.where(hi & stable, "consistent_high", np.where(hi & ~stable, "host_variable_high",
                   np.where(~hi & stable, "consistent_low", "host_variable_low")))
    return pd.Series(lab, index=st.index, name="stability_class")


def interaction_share(lesions: pd.DataFrame) -> dict[str, float]:
    """Sequential sums of squares of log(sAUDPC + 1) on isolate, crop and isolate x crop from
    per-plant data (no block column available in the lesion files, so block is not fitted)."""
    d = lesions.dropna(subset=["sAUDPC"]).copy()
    if d.empty:
        return {}
    y = np.log1p(d["sAUDPC"].clip(lower=0).to_numpy())

    def ss_res(cols: list[str]) -> float:
        X = np.ones((len(d), 1))
        if cols:
            X = np.hstack([X, pd.get_dummies(d[cols].astype(str).agg("|".join, axis=1) if len(cols) > 1 else d[cols[0]].astype(str),
                                             drop_first=True).to_numpy(float)])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ beta) ** 2).sum())

    sst = float(((y - y.mean()) ** 2).sum())
    r0, r1, r2 = sst, ss_res(["isolate"]), None
    X_add = np.hstack([np.ones((len(d), 1)),
                       pd.get_dummies(d["isolate"].astype(str), drop_first=True).to_numpy(float),
                       pd.get_dummies(d["crop"].astype(str), drop_first=True).to_numpy(float)])
    beta, *_ = np.linalg.lstsq(X_add, y, rcond=None)
    r2 = float(((y - X_add @ beta) ** 2).sum())
    r3 = ss_res(["isolate", "crop"])  # full cell-means model = isolate x crop
    return {"isolate": (r0 - r1) / sst, "crop": (r1 - r2) / sst, "isolate_x_crop": (r2 - r3) / sst,
            "residual": r3 / sst, "n_plants": int(len(d))}
