"""Reading count matrices, sample metadata, phenotypes and in-vitro assays.

Sample names follow the lab convention ``<hpi>_<rep>_<host>_<isolate>``, e.g. ``48_2_Gm_MNSS6``
and ``0_1_Gm_NC`` for non-inoculated controls (CH3 / CH4 pipeline, helpers/load_metadata.R).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import resolve


def canon(name: str, cfg: dict) -> str:
    """Canonical isolate name (aliases such as Ss1980 -> 1980)."""
    s = str(name).strip()
    return cfg["isolate_aliases"].get(s, s)


def parse_samples(samples: list[str], cfg: dict) -> pd.DataFrame:
    rows = []
    for s in samples:
        parts = s.split("_")
        if len(parts) < 4:
            raise ValueError(f"sample name {s!r} does not match <hpi>_<rep>_<host>_<isolate>")
        tp, rep, host, iso = parts[0], parts[1], parts[2], "_".join(parts[3:])
        rows.append({"sample": s, "hpi": tp, "rep": cfg["rep_remap"].get(rep, rep), "host": host,
                     "isolate": canon(iso, cfg), "inoculated": iso != "NC"})
    return pd.DataFrame(rows).set_index("sample")


def normalize_gene_id(g: str) -> str:
    """'Ss_jgi|Sclsc1|8988|SS1G_02403T0' -> 'SS1G_02403' (as in the R pipeline's normalize_gene_id);
    other identifiers are returned unchanged."""
    m = re.search(r"(SS1G_\d+)", str(g))
    return m.group(1) if m else str(g)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [normalize_gene_id(g) for g in df.index.astype(str)]
    df = df.fillna(0).astype(float)
    return df.groupby(level=0).sum()  # transcripts of one gene, if any, are summed


def load_expression(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All hosts' pathogen counts joined on shared genes (genes x samples) plus sample metadata."""
    mats = []
    for host, p in cfg["counts"].items():
        m = read_counts(resolve(cfg, p))
        mats.append(m)
    if not mats:
        raise ValueError("config: `counts` must map host codes to pathogen count matrices")
    shared = sorted(set.intersection(*[set(m.index) for m in mats]))
    counts = pd.concat([m.loc[shared] for m in mats], axis=1)
    meta = parse_samples(list(counts.columns), cfg)
    meta["pathogen_total"] = counts.sum(axis=0).values
    if cfg.get("library_totals"):
        lt = pd.read_csv(resolve(cfg, cfg["library_totals"]))
        tot = lt.set_index("sample")["host_total_counts"]
        meta["host_total"] = tot.reindex(meta.index).values
        meta["transcript_share"] = meta["pathogen_total"] / (meta["pathogen_total"] + meta["host_total"])
    return counts, meta


def load_phenotype(cfg: dict) -> pd.DataFrame:
    """Isolate x crop matrix of mean sAUDPC (rows: canonical isolate names). Accepts the long
    format written by CH4_multicrop_aggressiveness.Rmd (TableS1_long.csv: crop, isolate, mean, ...)
    or a wide isolate x crop table."""
    p = resolve(cfg, cfg["phenotype_table"])
    df = pd.read_csv(p)
    cols = {c.lower(): c for c in df.columns}
    if {"crop", "isolate", "mean"} <= set(cols):
        df = df.rename(columns={cols["crop"]: "crop", cols["isolate"]: "isolate", cols["mean"]: "mean"})
        df["isolate"] = df["isolate"].map(lambda x: canon(x, cfg))
        return df.pivot_table(index="isolate", columns="crop", values="mean", aggfunc="mean")
    df = df.set_index(df.columns[0])
    df.index = [canon(i, cfg) for i in df.index]
    crops = [c for c in df.columns if not c.lower().endswith("rank") and c.lower() not in ("mean rank", "rank range")]
    return df[crops].astype(float)


def sAUDPC(days: np.ndarray, values: np.ndarray) -> float:
    """Standardized AUDPC over the rated days, as in CH4: AUDPC / [D (n - 1) / n]."""
    ok = ~np.isnan(values)
    d, v = days[ok], values[ok]
    n = len(d)
    if n < 2:
        return float("nan")
    audpc = float(np.sum((v[1:] + v[:-1]) / 2 * np.diff(d)))
    D = float(d[-1] - d[0])
    return audpc / (D * (n - 1) / n) if D > 0 else float("nan")


def load_lesions(cfg: dict) -> pd.DataFrame:
    """Per-plant sAUDPC. Preferred: the per-plant table written by the CH4 pipeline
    (``per_plant_saudpc``, columns crop, isolate, sAUDPC). Fallback: recompute from raw lesion
    files (columns '<k>dpi'), which uses every rated day and may differ from the CH4 pipeline."""
    if cfg.get("per_plant_saudpc"):
        d = pd.read_csv(resolve(cfg, cfg["per_plant_saudpc"]))
        d["isolate"] = d["isolate"].map(lambda x: canon(x, cfg))
        return d[["crop", "isolate", "sAUDPC"]]
    rows = []
    for crop, p in (cfg.get("lesion_files") or {}).items():
        df = pd.read_csv(resolve(cfg, p), encoding="utf-8-sig")
        iso_col = "Ss_isolate" if "Ss_isolate" in df.columns else df.columns[1]
        day_cols = [c for c in df.columns if str(c).replace("dpi", "").strip().isdigit()]
        days = np.array([int(str(c).replace("dpi", "")) for c in day_cols], float)
        for _, r in df.iterrows():
            vals = pd.to_numeric(r[day_cols], errors="coerce").to_numpy(float)
            rows.append({"crop": crop, "isolate": canon(r[iso_col], cfg), "sAUDPC": sAUDPC(days, vals)})
    return pd.DataFrame(rows)


def load_in_vitro(cfg: dict) -> pd.DataFrame:
    """Isolate-level in vitro traits (oxalic acid, appressoria, radial growth ...) as columns."""
    out = {}
    for name, spec in (cfg.get("in_vitro") or {}).items():
        df = pd.read_csv(resolve(cfg, spec["file"]))
        for k, v in (spec.get("filter") or {}).items():
            df = df[df[k].astype(str) == str(v)]
        s = df.groupby(df[spec["isolate_col"]].map(lambda x: canon(x, cfg)))[spec["value_col"]].mean()
        out[name] = s
    return pd.DataFrame(out)


def load_annotation(cfg: dict) -> pd.DataFrame | None:
    if not cfg.get("annotation"):
        return None
    a = pd.read_csv(resolve(cfg, cfg["annotation"]), sep="\t")
    idc = cfg["annotation_id_col"]
    keep = [c for c in ("signal_peptide", "is_effector", "effectorp_pred", "cazyme", "cazyme_family", "GO_terms") if c in a.columns]
    return a.drop_duplicates(idc).set_index(idc)[keep]
