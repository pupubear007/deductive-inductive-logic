"""Assemble the modelling panel: libraries x genes, metadata and isolate x host targets."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .io import load_expression, load_phenotype
from .model import Design
from .normalize import filter_genes, log_normalize, logit, thin_to_depth


@dataclass
class Panel:
    expr: pd.DataFrame        # genes x libraries, log2 normalised
    meta: pd.DataFrame        # libraries: hpi, rep, host, isolate, transcript_share
    target: pd.DataFrame      # isolates x hosts (centred log sAUDPC or within-host rank)
    hosts: list[str]
    isolates: list[str]

    def libs(self, isolates=None, hosts=None) -> pd.Index:
        m = self.meta
        keep = np.ones(len(m), bool)
        if isolates is not None:
            keep &= m["isolate"].isin(list(isolates)).to_numpy()
        if hosts is not None:
            keep &= m["host"].isin(list(hosts)).to_numpy()
        return m.index[keep]


def make_target(pheno: pd.DataFrame, isolates: list[str], hosts: list[str], host_to_crop: dict, scale: str) -> pd.DataFrame:
    missing = [i for i in isolates if i not in pheno.index]
    if missing:
        raise KeyError(f"isolates {missing} not in the phenotype table; add them or an alias in `isolate_aliases`")
    out = {}
    for h in hosts:
        crop = host_to_crop.get(h)
        if crop not in pheno.columns:
            raise KeyError(f"host {h!r} maps to crop {crop!r}, which is not a column of the phenotype table")
        v = pheno.loc[isolates, crop].astype(float)
        if scale == "z":      # standardised within host, as in the proposal
            out[h] = (v - v.mean()) / v.std(ddof=1)
        elif scale == "raw":  # centred sAUDPC
            out[h] = v - v.mean()
        elif scale == "log":
            out[h] = np.log(v) - np.log(v).mean()
        elif scale == "rank":
            out[h] = v.rank() - v.rank().mean()
        else:
            raise ValueError(f"target_scale {scale!r}: use z, raw, log or rank")
    return pd.DataFrame(out)


def build_panel(cfg: dict, counts: pd.DataFrame | None = None, meta: pd.DataFrame | None = None,
                pheno: pd.DataFrame | None = None, timepoints: list[str] | None = None) -> Panel:
    if counts is None:
        counts, meta = load_expression(cfg)
    if pheno is None:
        pheno = load_phenotype(cfg)
    tps = [str(t) for t in (timepoints or cfg["timepoints"])]
    m = meta[meta["inoculated"] & meta["hpi"].astype(str).isin(tps)].copy()
    c = counts[m.index]
    if cfg.get("thin_to_equal_depth"):
        c = pd.concat([thin_to_depth(c[m.index[m["host"] == h]], seed=cfg["seed"]) for h in m["host"].unique()], axis=1)[m.index]
    c = filter_genes(c, cfg["min_count"], cfg["min_samples"])
    expr = log_normalize(c, by=m["host"])
    if "transcript_share" not in m.columns:
        m["transcript_share"] = np.nan
    hosts = sorted(m["host"].unique())
    isolates = sorted(m["isolate"].unique())
    target = make_target(pheno, isolates, hosts, cfg["host_to_crop"], cfg["target_scale"])
    return Panel(expr, m, target, hosts, isolates)


# ------------------------------------------------------------------------------------------------
# Fold-wise preprocessing (fitted on training libraries only)
# ------------------------------------------------------------------------------------------------

@dataclass
class Prep:
    genes: list[str]
    mu: dict[str, np.ndarray]
    sd: dict[str, np.ndarray]
    cov_mu: dict[str, float]


def fit_prep(p: Panel, train: pd.Index, top_k: int | None) -> Prep:
    """Per-host centring/scaling and an unsupervised top-variance gene filter, learned on the
    training libraries only (so nothing about the held-out isolate leaks into the features)."""
    mt = p.meta.loc[train]
    var = []
    for h in mt["host"].unique():
        e = p.expr[mt.index[mt["host"] == h]]
        var.append(e.var(axis=1))
    v = pd.concat(var, axis=1).mean(axis=1)
    genes = list(v.sort_values(ascending=False).index[:top_k]) if top_k else list(v.index)
    mu, sd, cmu = {}, {}, {}
    for h in p.hosts:
        cols = mt.index[mt["host"] == h]
        if len(cols) == 0:  # host unseen in training (cross-host transfer): use its own libraries
            cols = p.meta.index[p.meta["host"] == h]
        e = p.expr.loc[genes, cols].to_numpy(float)
        mu[h] = e.mean(1)
        s = e.std(1)
        sd[h] = np.where(s < 1e-8, 1.0, s)
        ts = p.meta.loc[cols, "transcript_share"].to_numpy(float)
        cmu[h] = float(np.nanmean(logit(ts))) if np.isfinite(ts).any() else 0.0
    return Prep(genes, mu, sd, cmu)


def design(p: Panel, libs: pd.Index, prep: Prep, cfg: dict, target: pd.DataFrame | None = None) -> Design:
    target = p.target if target is None else target
    m = p.meta.loc[libs]
    X = np.zeros((len(libs), len(prep.genes)))
    for h in m["host"].unique():
        idx = np.where(m["host"].to_numpy() == h)[0]
        e = p.expr.loc[prep.genes, m.index[idx]].to_numpy(float).T
        X[idx] = (e - prep.mu[h]) / prep.sd[h]
    host_idx = np.array([p.hosts.index(h) for h in m["host"]])
    tps = sorted(p.meta["hpi"].astype(str).unique())
    group = np.array([p.hosts.index(h) * len(tps) + tps.index(str(t)) for h, t in zip(m["host"], m["hpi"])])
    if cfg.get("colonization_covariate") and np.isfinite(m["transcript_share"].to_numpy(float)).all():
        cov = np.array([logit(np.array([s]))[0] - prep.cov_mu[h] for s, h in zip(m["transcript_share"], m["host"])])
    else:
        cov = np.zeros(len(libs))
    y = np.array([target.loc[i, h] for i, h in zip(m["isolate"], m["host"])], float)
    cell = m["isolate"].astype(str) + "|" + m["host"].astype(str)
    w = 1.0 / cell.map(cell.value_counts()).to_numpy(float)
    w = w / w.sum()
    return Design(X, group, host_idx, cov, y, w, n_groups=len(p.hosts) * len(tps), n_hosts=len(p.hosts))
