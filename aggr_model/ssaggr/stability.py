"""Objective 1: aggressiveness and its consistency across hosts.

Joint regression of isolate aggressiveness on a host index (Finlay and Wilkinson 1963;
Eberhart and Russell 1966): for isolate i on host h,

    y_ih = mu_i + b_i * I_h + d_ih,      I_h = mean_i y_ih - grand mean,

b_i measures how strongly the isolate tracks host susceptibility (b = 1: average response) and
the deviation mean square s2d_i = sum_h d_ih^2 / (H - 2) measures host-specific departures.

Scale (``stability_scale``):
  raw  sAUDPC as measured (Eberhart and Russell's form). Default.
  log  log sAUDPC, closer to Finlay and Wilkinson's use of a log scale: a host that doubles every
       isolate's lesion then shifts every isolate equally.
z-scoring within host is not offered for this regression: it sets every host mean, and so I_h,
to zero.

Classes are *your hypothesis*: either listed explicitly (``hypothesis_classes``) or defined by
cut-offs on these statistics (``stability_cutoffs``), fixed before the expression analysis. The
code never invents cut-offs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def stability_table(pheno: pd.DataFrame, scale: str = "raw") -> pd.DataFrame:
    """pheno: isolates x crops, mean sAUDPC. Returns one row per isolate."""
    if scale == "log":
        y = np.log(pheno.clip(lower=1e-6))
    elif scale == "raw":
        y = pheno.astype(float)
    else:
        raise ValueError("stability_scale must be 'raw' or 'log'")
    H = y.shape[1]
    I = (y.mean(axis=0) - y.to_numpy().mean()).to_numpy()
    rows = []
    for iso, r in y.iterrows():
        v = r.to_numpy(float)
        mu = v.mean()
        b = float(np.dot(I, v - mu) / np.dot(I, I)) if np.dot(I, I) > 0 else float("nan")
        resid = v - mu - b * I
        s2d = float((resid ** 2).sum() / (H - 2)) if H > 2 else float("nan")
        rows.append({"isolate": iso, "mean_sAUDPC": float(pheno.loc[iso].mean()), "slope_b": b, "deviation_ms": s2d})
    out = pd.DataFrame(rows).set_index("isolate")
    ranks = pheno.rank(axis=0)
    out["mean_rank"] = ranks.mean(axis=1)
    out["rank_range"] = ranks.max(axis=1) - ranks.min(axis=1)
    out.attrs["host_index"] = dict(zip(pheno.columns, I))
    return out


def classes_from_cutoffs(st: pd.DataFrame, cut: dict) -> pd.Series | None:
    """Supported keys (any subset):
      mean_rank_min_high          isolates with mean_rank >= this are 'high'
      rank_range_max_consistent   isolates with rank_range <= this are 'consistent'
      deviation_ms_max_consistent isolates with deviation_ms <= this are 'consistent'
    """
    if not cut:
        return None
    lab = pd.Series("", index=st.index, dtype=object)
    if "mean_rank_min_high" in cut:
        lab += np.where(st["mean_rank"] >= cut["mean_rank_min_high"], "high", "low")
    cons = pd.Series(True, index=st.index)
    used = False
    if "rank_range_max_consistent" in cut:
        cons &= st["rank_range"] <= cut["rank_range_max_consistent"]
        used = True
    if "deviation_ms_max_consistent" in cut:
        cons &= st["deviation_ms"] <= cut["deviation_ms_max_consistent"]
        used = True
    if used:
        lab = lab.where(lab == "", lab + "_") + np.where(cons, "consistent", "host_variable")
    return lab.rename("class_from_cutoffs")


def interaction_share(lesions: pd.DataFrame) -> dict[str, float]:
    """Sequential sums of squares of per-plant sAUDPC on isolate, crop and isolate x crop
    (no block column in the per-plant file, so block is not fitted)."""
    d = lesions.dropna(subset=["sAUDPC"]).copy()
    if d.empty:
        return {}
    y = d["sAUDPC"].to_numpy(float)
    iso = pd.get_dummies(d["isolate"].astype(str), drop_first=True).to_numpy(float)
    crop = pd.get_dummies(d["crop"].astype(str), drop_first=True).to_numpy(float)
    cell = pd.get_dummies(d["isolate"].astype(str) + "|" + d["crop"].astype(str), drop_first=True).to_numpy(float)
    one = np.ones((len(d), 1))

    def rss(X):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ beta) ** 2).sum())

    sst = float(((y - y.mean()) ** 2).sum())
    r1, r2, r3 = rss(np.hstack([one, iso])), rss(np.hstack([one, iso, crop])), rss(np.hstack([one, cell]))
    return {"isolate": (sst - r1) / sst, "crop": (r1 - r2) / sst, "isolate_x_crop": (r2 - r3) / sst,
            "residual": r3 / sst, "n_plants": int(len(d))}
