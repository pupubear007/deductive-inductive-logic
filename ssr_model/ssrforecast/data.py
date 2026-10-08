"""Loading the four input tables into time-adapted training samples.

A *sample* is a pair (unit, issue day tau). Its inputs use only information available on or before
tau (adaptedness, Definition S3 in the mechanism notes): environment series over
[tau - lookback + 1, tau], the latest image of each modality taken on or before tau, and static
unit attributes. Its targets are assay readings in (tau, tau + horizon] for task = forecast, or
on day tau for task = nowcast (diagnosis).
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from PIL import Image
from torch.utils.data import Dataset

from .mechanism import COURT_INTERVENTIONS
from .schema import OUTCOME_OF_ASSAY, binary_assays, continuous_assays

DAY = np.timedelta64(1, "D")


def _s(v) -> str:
    """Category value as a string; missing values become the empty string."""
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA else str(v)


# ---------------------------------------------------------------------------------------------
# Table loading
# ---------------------------------------------------------------------------------------------

def read_tables(data_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Read whichever of the four CSV files exist. Missing optional tables become empty frames."""
    d = Path(data_dir)
    out: dict[str, pd.DataFrame] = {}
    for name in ("units", "observations", "timeseries", "images"):
        p = d / f"{name}.csv"
        out[name] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    for name, df in out.items():
        if "unit_id" in df.columns:
            df["unit_id"] = df["unit_id"].astype(str)
        if "station_id" in df.columns:
            df["station_id"] = df["station_id"].astype(str)
    return out


def _to_days(s: pd.Series) -> pd.Series:
    """Datetime-like series -> float days since 1970-01-01."""
    dt = pd.to_datetime(s, errors="coerce")
    return (dt - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)


def unit_start_days(units: pd.DataFrame) -> pd.Series:
    """Day 0 of each unit: planting date, else inoculation date, else NaN (filled later)."""
    start = pd.Series(np.nan, index=units["unit_id"].values, dtype=float)
    for col in ("planting_date", "inoculation_date"):
        if col in units.columns:
            v = _to_days(units[col]).values
            start = start.where(~start.isna(), pd.Series(v, index=start.index))
    return start


def add_time(df: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    """Give every row an absolute day ``t``. Rows may carry ``date`` (calendar) or ``dpi``
    (days post inoculation; requires ``inoculation_date`` in units.csv, else dpi is used as is)."""
    if df.empty:
        return df.assign(t=pd.Series(dtype=float))
    df = df.copy()
    t = pd.Series(np.nan, index=df.index, dtype=float)
    if "date" in df.columns:
        t = _to_days(df["date"])
    if "dpi" in df.columns:
        inoc = pd.Series(0.0, index=units["unit_id"].values)
        if "inoculation_date" in units.columns:
            inoc = pd.Series(_to_days(units["inoculation_date"]).fillna(0.0).values, index=units["unit_id"].values)
        key = df["unit_id"] if "unit_id" in df.columns else None
        base = key.map(inoc).fillna(0.0) if key is not None else 0.0
        t = t.where(~t.isna(), base + pd.to_numeric(df["dpi"], errors="coerce"))
    df["t"] = t.astype(float)
    return df.dropna(subset=["t"])


# ---------------------------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------------------------

@dataclass
class Meta:
    """Everything needed to rebuild tensors identically at evaluation time."""

    bin_names: list[str]
    cont_names: list[str]
    outcome_names: list[str]
    series_vars: list[str]
    num_cols: list[str]
    cat_cols: list[str]
    vocab: dict[str, dict[str, int]]
    num_mean: list[float]
    num_std: list[float]
    series_mean: list[float]
    series_std: list[float]
    modalities: dict[str, dict[str, int]]
    lookback: int
    map_size: int
    cont_scale: list[float]
    splits: dict[str, list[str]] = field(default_factory=dict)


def build_samples(tables: dict[str, pd.DataFrame], cfg: dict[str, Any]) -> pd.DataFrame:
    """One row per (unit, issue day) with lists of target readings."""
    units = tables["units"].copy()
    obs = add_time(tables["observations"], units)
    obs = obs[obs["assay"].isin(cfg["assays"].keys())]
    start = unit_start_days(units)
    first_obs = obs.groupby("unit_id")["t"].min()
    start = start.fillna(first_obs).fillna(0.0)

    H = float(cfg["horizon_days"])
    rows = []
    for uid, g in obs.groupby("unit_id"):
        days = np.sort(g["t"].unique())
        t0 = float(start.get(uid, days.min()))
        if cfg["task"] == "forecast":
            step = float(cfg.get("issue_every_days") or 0)
            grid = np.arange(t0, days.max(), step) if step > 0 else np.array([])
            issues = np.unique(np.concatenate([[t0], days, grid]))
            for tau in issues:
                win = g[(g["t"] > tau) & (g["t"] <= tau + H)]
                if len(win):
                    rows.append({"unit_id": uid, "tau": float(tau), "t0": t0, "targets": win})
        else:  # nowcast / diagnosis
            for tau in days:
                rows.append({"unit_id": uid, "tau": float(tau), "t0": t0, "targets": g[g["t"] == tau]})
    return pd.DataFrame(rows)


def split_units(units: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, list[str]]:
    """Group-wise split (default: by season), so the test set is an unobserved season
    (Theorem 5.4: the next season is an unobserved individual)."""
    sp = cfg["split"]
    key = sp.get("group_by", "season")
    if key == "site_season" and {"site", "season"} <= set(units.columns):
        groups = units["site"].map(_s) + "|" + units["season"].map(_s)
    elif key in units.columns:
        groups = units[key].map(_s)
    else:
        groups = units["unit_id"].astype(str)
    # Split within each experiment type, so growth-chamber, greenhouse and field units all
    # appear in every split (each type has its own held-out seasons / runs).
    strata = units["experiment_type"].map(_s) if "experiment_type" in units.columns else pd.Series("", index=units.index)
    rng = random.Random(cfg["seed"])
    assign: dict[tuple[str, str], str] = {}
    for st in sorted(strata.unique()):
        uniq = sorted(groups[strata == st].unique())
        rng.shuffle(uniq)
        n = len(uniq)
        n_test = max(1, round(n * sp["test_frac"])) if n >= 3 else 0
        n_val = max(1, round(n * sp["val_frac"])) if n >= 3 else 0
        for j, gk in enumerate(uniq):
            assign[(st, gk)] = "test" if j < n_test else "val" if j < n_test + n_val else "train"
    out = {"train": [], "val": [], "test": []}
    for uid, st, gk in zip(units["unit_id"].astype(str), strata, groups):
        out[assign[(st, gk)]].append(uid)
    if not out["train"]:  # tiny datasets: train on everything
        out["train"] = list(units["unit_id"].astype(str))
    return out


# ---------------------------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------------------------

def _load_image(path: Path, channels: int, size: int) -> torch.Tensor:
    """PNG/JPG/TIF via Pillow, or .npy arrays [C,H,W] / [H,W,C] for multispectral data."""
    if path.suffix.lower() == ".npy":
        a = np.load(path).astype(np.float32)
        if a.ndim == 2:
            a = a[None]
        if a.shape[-1] == channels and a.shape[0] != channels:
            a = np.moveaxis(a, -1, 0)
        x = torch.from_numpy(a)
        if x.max() > 1.5:
            x = x / 255.0
    else:
        im = Image.open(path)
        im = im.convert("L" if channels == 1 else "RGB")
        x = torch.from_numpy(np.asarray(im, dtype=np.float32) / 255.0)
        x = x[None] if x.ndim == 2 else x.permute(2, 0, 1)
    c = x.shape[0]
    if c < channels:
        x = torch.cat([x, torch.zeros(channels - c, *x.shape[1:])], 0)
    x = x[:channels]
    return torch.nn.functional.interpolate(x[None], size=(size, size), mode="bilinear", align_corners=False)[0]


class SSRDataset(Dataset):
    def __init__(self, samples: pd.DataFrame, tables: dict[str, pd.DataFrame], cfg: dict[str, Any],
                 meta: Meta, data_dir: str | Path, train: bool = False):
        self.cfg, self.meta, self.dir, self.train = cfg, meta, Path(data_dir), train
        self.samples = samples.reset_index(drop=True)
        units = tables["units"].set_index("unit_id")
        self.units = units
        self._cache: dict[tuple[str, int, int], torch.Tensor] = {}

        ts = add_time(tables["timeseries"], tables["units"]) if not tables["timeseries"].empty else tables["timeseries"]
        if not ts.empty and "unit_id" not in ts.columns and "station_id" in ts.columns and "station_id" in units.columns:
            st = units["station_id"].astype(str).reset_index()
            ts = ts.merge(st, on="station_id")
        self.series: dict[str, pd.DataFrame] = {}
        if not ts.empty:
            ts = ts[ts["variable"].isin(meta.series_vars)].copy()
            ts["day"] = np.floor(ts["t"]).astype(int)
            daily = ts.groupby(["unit_id", "day", "variable"])["value"].mean().unstack("variable")
            for uid, g in daily.groupby(level=0):
                self.series[uid] = g.droplevel(0).reindex(columns=meta.series_vars)

        im = add_time(tables["images"], tables["units"]) if not tables["images"].empty else tables["images"]
        self.images = {uid: g.sort_values("t") for uid, g in im.groupby("unit_id")} if not im.empty else {}
        self._targets = [self._compute_targets(r) for _, r in self.samples.iterrows()]  # once, not per epoch

    def _compute_targets(self, s: pd.Series):
        """Binary readings: any positive in the window (fractions: last reading > 0).
        Continuous readings: the last reading in the window, divided by the assay's scale."""
        tg: pd.DataFrame = s["targets"]
        order = np.argsort(tg["t"].values, kind="stable")
        assay = tg["assay"].values[order]
        vals = pd.to_numeric(tg["value"], errors="coerce").values[order]
        nb, nc = len(self.meta.bin_names), len(self.meta.cont_names)
        y_bin, m_bin = torch.zeros(nb), torch.zeros(nb)
        y_cont, m_cont = torch.zeros(nc), torch.zeros(nc)
        for j, a in enumerate(self.meta.bin_names):
            v = vals[(assay == a) & ~np.isnan(vals)]
            if len(v):
                kind = self.cfg["assays"][a]["kind"]
                y_bin[j] = float(v.max() > 0) if kind == "binary" else float(v[-1] > 0)
                m_bin[j] = 1.0
        for j, a in enumerate(self.meta.cont_names):
            v = vals[(assay == a) & ~np.isnan(vals)]
            if len(v):
                y_cont[j] = float(v[-1]) / self.meta.cont_scale[j]
                m_cont[j] = 1.0
        return y_bin, m_bin, y_cont, m_cont

    def __len__(self) -> int:
        return len(self.samples)

    # -- pieces ---------------------------------------------------------------------------
    def _series(self, uid: str, tau: float) -> tuple[torch.Tensor, torch.Tensor]:
        L, V = self.meta.lookback, len(self.meta.series_vars)
        x = np.zeros((L, V), np.float32)
        m = np.zeros((L, V), np.float32)
        df = self.series.get(uid)
        if df is not None and V:
            days = np.arange(int(np.floor(tau)) - L + 1, int(np.floor(tau)) + 1)
            w = df.reindex(days).values.astype(np.float32)
            m = (~np.isnan(w)).astype(np.float32)
            mu, sd = np.array(self.meta.series_mean, np.float32), np.array(self.meta.series_std, np.float32)
            x = np.nan_to_num((w - mu) / sd) * m
        return torch.from_numpy(x), torch.from_numpy(m)

    def _image(self, path: str, ch: int, size: int) -> torch.Tensor:
        key = (path, ch, size)
        if key not in self._cache:
            self._cache[key] = _load_image(self.dir / path, ch, size)
        return self._cache[key]

    def _tab(self, uid: str, tau: float, t0: float) -> tuple[torch.Tensor, torch.Tensor]:
        row = self.units.loc[uid]
        nums = []
        for c in self.meta.num_cols:
            v = tau - t0 if c == "days_since_start" else pd.to_numeric(row.get(c, np.nan), errors="coerce")
            nums.append(np.nan if v is None else float(v))
        nums = np.array(nums, np.float32)
        mu, sd = np.array(self.meta.num_mean, np.float32), np.array(self.meta.num_std, np.float32)
        nums = np.nan_to_num((nums - mu) / sd)
        cats = [self.meta.vocab[c].get(_s(row.get(c, "")), 0) for c in self.meta.cat_cols]
        return torch.from_numpy(nums), torch.tensor(cats, dtype=torch.long)

    # -- item -----------------------------------------------------------------------------
    def __getitem__(self, i: int) -> dict[str, Any]:
        s = self.samples.iloc[i]
        uid, tau, t0 = s["unit_id"], s["tau"], s["t0"]
        x_series, m_series = self._series(uid, tau)
        x_num, x_cat = self._tab(uid, tau, t0)

        imgs, present = {}, {}
        g = self.images.get(uid)
        for mod, spec in self.meta.modalities.items():
            imgs[mod] = torch.zeros(spec["channels"], spec["size"], spec["size"])
            present[mod] = torch.tensor(0.0)
            if g is not None:
                cand = g[(g["modality"] == mod) & (g["t"] <= tau) & (g["t"] > tau - self.meta.lookback)]
                if len(cand):
                    imgs[mod] = self._image(cand.iloc[-1]["path"], spec["channels"], spec["size"])
                    present[mod] = torch.tensor(1.0)

        y_bin, m_bin, y_cont, m_cont = self._targets[i]

        S = self.meta.map_size
        y_map, m_map = torch.zeros(1, S, S), torch.tensor(0.0)
        if g is not None and "target_map_path" in g.columns:
            if self.cfg["task"] == "forecast":
                win = g[(g["t"] > tau) & (g["t"] <= tau + self.cfg["horizon_days"])]
            else:
                win = g[g["t"] == tau]
            win = win[win["target_map_path"].notna() & (win["target_map_path"].astype(str) != "")]
            if len(win):
                y_map = self._image(str(win.iloc[-1]["target_map_path"]), 1, S).clamp(0, 1)
                m_map = torch.tensor(1.0)

        method = _s(self.units.loc[uid].get("inoculation_method", ""))
        return {
            "series": x_series, "series_mask": m_series, "num": x_num, "cat": x_cat,
            "images": imgs, "present": present,
            "do_court": torch.tensor(float(method in COURT_INTERVENTIONS)),
            "y_bin": y_bin, "m_bin": m_bin, "y_cont": y_cont, "m_cont": m_cont,
            "y_map": y_map, "m_map": m_map, "index": torch.tensor(i),
        }


def collate(batch: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k in batch[0]:
        if k in ("images", "present"):
            out[k] = {m: torch.stack([b[k][m] for b in batch]) for m in batch[0][k]}
        else:
            out[k] = torch.stack([b[k] for b in batch])
    return out


# ---------------------------------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------------------------------

def build_datasets(cfg: dict[str, Any], meta: Meta | None = None):
    """Return (train, val, test, meta). If ``meta`` is given (evaluation), reuse its
    vocabularies, normalisation and split instead of recomputing them."""
    tables = read_tables(cfg["data_dir"])
    if tables["units"].empty or tables["observations"].empty:
        raise FileNotFoundError(f"{cfg['data_dir']}: units.csv and observations.csv are required "
                                "(run `python -m ssrforecast.validate` for a full report)")
    units = tables["units"]
    samples = build_samples(tables, cfg)
    if samples.empty:
        raise ValueError("no (unit, issue day) pair has a target reading inside the horizon")

    if meta is None:
        splits = split_units(units, cfg)
        train_units = set(splits["train"])
        tr = samples[samples["unit_id"].isin(train_units)]

        num_cols = list(cfg["tabular"]["numeric"])
        num_vals = []
        for _, s in tr.iterrows():
            row = units.set_index("unit_id").loc[s["unit_id"]]
            num_vals.append([s["tau"] - s["t0"] if c == "days_since_start"
                             else pd.to_numeric(row.get(c, np.nan), errors="coerce") for c in num_cols])
        nv = np.array(num_vals, dtype=float) if num_vals else np.zeros((0, len(num_cols)))
        num_mean = np.nan_to_num(np.nanmean(nv, 0)) if len(nv) else np.zeros(len(num_cols))
        num_std = np.nan_to_num(np.nanstd(nv, 0)) if len(nv) else np.ones(len(num_cols))
        num_std[num_std < 1e-6] = 1.0

        cat_cols = [c for c in cfg["tabular"]["categorical"]]
        tru = units[units["unit_id"].isin(train_units)]
        vocab = {c: {v: i + 1 for i, v in enumerate(sorted(tru[c].map(_s).unique()))} if c in tru.columns else {}
                 for c in cat_cols}

        series_vars = list(cfg["series_variables"])
        ts = tables["timeseries"]
        s_mean, s_std = np.zeros(len(series_vars)), np.ones(len(series_vars))
        if not ts.empty:
            for j, v in enumerate(series_vars):
                vals = pd.to_numeric(ts.loc[ts["variable"] == v, "value"], errors="coerce").dropna()
                if len(vals):
                    s_mean[j], s_std[j] = vals.mean(), max(vals.std(), 1e-6) if len(vals) > 1 else 1.0

        cont = continuous_assays(cfg)
        meta = Meta(
            bin_names=binary_assays(cfg), cont_names=cont,
            outcome_names=[OUTCOME_OF_ASSAY.get(a, a) for a in cont],
            series_vars=series_vars, num_cols=num_cols, cat_cols=cat_cols, vocab=vocab,
            num_mean=num_mean.tolist(), num_std=num_std.tolist(),
            series_mean=s_mean.tolist(), series_std=s_std.tolist(),
            modalities={k: dict(v) for k, v in cfg["image_modalities"].items()},
            lookback=int(cfg["lookback_days"]), map_size=int(cfg["diffusion"]["map_size"]),
            cont_scale=[float(cfg["assays"][a].get("scale", 1.0)) for a in cont],
            splits=splits,
        )
    sets = []
    for name in ("train", "val", "test"):
        part = samples[samples["unit_id"].isin(set(meta.splits[name]))]
        sets.append(SSRDataset(part, tables, cfg, meta, cfg["data_dir"], train=(name == "train")))
    return sets[0], sets[1], sets[2], meta
