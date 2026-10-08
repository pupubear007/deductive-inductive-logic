"""Count normalisation (DESeq2-style median-of-ratios with 'poscounts'), filtering and thinning.

The log2-normalised values approximate, but are not identical to, DESeq2's VST. If you prefer
the exact VST used in the chapters, export ``assay(vst(dds))`` from R and point the config at it
(see DATA_REQUIREMENTS.md); the rest of the pipeline is unchanged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def filter_genes(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def size_factors_poscounts(counts: pd.DataFrame) -> pd.Series:
    x = counts.to_numpy(float)
    n = x.shape[1]
    with np.errstate(divide="ignore"):
        logx = np.where(x > 0, np.log(x), 0.0)
    lgm = logx.sum(1) / n  # poscounts: zeros contribute 0 to the log sum, divide by n
    ok = np.isfinite(lgm) & (lgm > 0)
    sf = np.empty(n)
    for j in range(n):
        pos = ok & (x[:, j] > 0)
        sf[j] = np.exp(np.median(np.log(x[pos, j]) - lgm[pos])) if pos.any() else 1.0
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=counts.columns)


def log_normalize(counts: pd.DataFrame, by: pd.Series | None = None) -> pd.DataFrame:
    """log2(count / size factor + 1). Size factors are computed within each group of ``by``
    (e.g. host), because the two hosts were sequenced and quantified separately."""
    if by is None:
        sf = size_factors_poscounts(counts)
    else:
        sf = pd.concat([size_factors_poscounts(counts.loc[:, by.index[by == g]]) for g in by.unique()])
        sf = sf.reindex(counts.columns)
    return np.log2(counts / sf + 1.0)


def thin_to_depth(counts: pd.DataFrame, depth: int | None = None, seed: int = 0) -> pd.DataFrame:
    """Binomial thinning so every library has the same pathogen depth (proposal section 5.6),
    so that differences in fungal load cannot masquerade as expression differences."""
    rng = np.random.default_rng(seed)
    tot = counts.sum(axis=0)
    depth = int(tot.min()) if depth is None else depth
    out = counts.copy()
    for c in counts.columns:
        p = min(1.0, depth / tot[c]) if tot[c] > 0 else 0.0
        out[c] = rng.binomial(np.round(counts[c].to_numpy()).astype(np.int64), p)
    return out.astype(float)


def logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))
