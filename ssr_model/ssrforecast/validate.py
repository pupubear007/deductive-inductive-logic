"""Check a data directory and say what is missing before training.

    python -m ssrforecast.validate --data DIR [--config CFG] [--report DIR/data_report.md]

The report lists, in order: files and columns, what each table contains, label audits against the
mechanism's hard constraints, which model components the data can train, and a prioritised
DATA REQUEST of what to collect next.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .data import add_time, read_tables
from .mechanism import COURT_INTERVENTIONS
from .schema import EXPERIMENT_TYPES, IMAGE_MODALITIES, REQUIRED_COLUMNS, load_config

LEVEL_ORDER = ["apothecia", "spore", "court", "infection"]
MIN_PER_CLASS = 20  # rough floor below which a binary head is not worth training


def validate(data_dir: str | Path, cfg: dict[str, Any]) -> tuple[list[str], dict[str, bool]]:
    d = Path(data_dir)
    T = read_tables(d)
    lines: list[str] = [f"# Data report for `{d}`", ""]
    request: list[str] = []
    ok: dict[str, bool] = {}

    def h(t: str) -> None:
        lines.extend(["", f"## {t}", ""])

    # 1. Files and columns ------------------------------------------------------------------
    h("1. Files and columns")
    for fname, cols in REQUIRED_COLUMNS.items():
        df = T[fname[:-4]]
        mandatory = fname in ("units.csv", "observations.csv")
        if df.empty:
            lines.append(f"- {'MISSING (required)' if mandatory else 'absent (optional)'}: `{fname}`")
            if mandatory:
                request.append(f"`{fname}` is required; see DATA_REQUIREMENTS.md for its columns.")
            continue
        need = list(cols)
        if fname == "timeseries.csv" and "unit_id" not in df.columns and "station_id" in df.columns:
            need = [c for c in need if c != "unit_id"] + ["station_id"]
        miss = [c for c in need if c not in df.columns]
        if fname != "units.csv" and not ({"date", "dpi"} & set(df.columns)):
            miss.append("date or dpi")
        lines.append(f"- `{fname}`: {len(df)} rows; " + (f"**missing columns {miss}**" if miss else "columns OK"))
        if miss:
            request.append(f"`{fname}` needs columns {miss}.")
    if T["units"].empty or T["observations"].empty:
        return lines + ["", "Cannot continue without units.csv and observations.csv."], ok

    units = T["units"]
    # 2. Units ----------------------------------------------------------------------------------
    h("2. Units")
    et = units.get("experiment_type", pd.Series(dtype=str)).astype(str)
    lines.append(f"- experiment types: {et.value_counts().to_dict()}")
    bad = sorted(set(et) - EXPERIMENT_TYPES)
    if bad:
        lines.append(f"- **unknown experiment_type values** {bad} (use {sorted(EXPERIMENT_TYPES)})")
    n_season = units["season"].nunique() if "season" in units.columns else 0
    lines.append(f"- seasons: {n_season}; sites: {units['site'].nunique() if 'site' in units.columns else 'n/a'}")
    if n_season < 3:
        request.append("At least 3 seasons (or site-seasons) are needed to hold out whole seasons for "
                       "testing; with fewer, test scores describe seen seasons only (Theorem 5.4).")
    for c in ("cultivar", "inoculation_method", "planting_date", "inoculation_date"):
        lines.append(f"- `{c}`: " + ("present" if c in units.columns else "absent"))
    if "inoculation_method" in units.columns:
        lines.append(f"- inoculation methods: {units['inoculation_method'].astype(str).value_counts().to_dict()} "
                     f"(court-level interventions: {sorted(COURT_INTERVENTIONS)})")

    # 3. Observations -------------------------------------------------------------------------
    h("3. Observations (assays)")
    obs = add_time(T["observations"], units)
    counts = obs["assay"].value_counts()
    known = [a for a in counts.index if a in cfg["assays"]]
    unknown = [a for a in counts.index if a not in cfg["assays"]]
    for a in known:
        spec = cfg["assays"][a]
        v = pd.to_numeric(obs.loc[obs["assay"] == a, "value"], errors="coerce").dropna()
        desc = f"n={len(v)}, range [{v.min():.3g}, {v.max():.3g}]" if len(v) else "no numeric values"
        note = ""
        if spec["kind"] == "binary":
            pos = int((v > 0).sum())
            desc += f", positives={pos}, negatives={len(v) - pos}"
            if not set(np.unique(v)) <= {0, 1}:
                note = " **values are not 0/1**"
        if spec["kind"] == "fraction" and len(v) and v.max() > 1:
            note = " **looks like percent: divide by 100 or set scale: 100**"
        lines.append(f"- `{a}` ({spec['kind']}{', level ' + spec['level'] if 'level' in spec else ''}): {desc}{note}")
    if unknown:
        lines.append(f"- **assays not in the config (ignored)**: {dict(counts[unknown])} -> add them under "
                     "`assays:` in the YAML with kind / level / scale")
    missing_assays = [a for a in cfg["assays"] if a not in counts.index]
    if missing_assays:
        lines.append(f"- configured but absent: {missing_assays}")

    # 4. Environment series ----------------------------------------------------------------
    h("4. Environment series")
    ts = T["timeseries"]
    if ts.empty:
        lines.append("- none")
        request.append("Environment series (`timeseries.csv`): daily air temperature, RH, rainfall and, "
                       "for the microclimate gate, canopy leaf wetness / RH and soil moisture.")
    else:
        found = set(ts["variable"].astype(str))
        want = set(cfg["series_variables"])
        lines.append(f"- variables found: {sorted(found)}")
        if want - found:
            lines.append(f"- **configured but absent**: {sorted(want - found)}")
        if not ({"leaf_wetness", "canopy_rh", "soil_moisture"} & found):
            request.append("No canopy-level moisture variable: the microclimate gate will be learned from "
                           "macro weather only (a weaker proxy).")
        key = "unit_id" if "unit_id" in ts.columns else "station_id"
        lines.append(f"- {ts[key].nunique()} {key}s with series")

    # 5. Images --------------------------------------------------------------------------------
    h("5. Images")
    im = T["images"]
    if im.empty:
        lines.append("- none")
    else:
        for m, g in im.groupby("modality"):
            missing_files = sum(not (d / str(p)).exists() for p in g["path"])
            flag = "" if m in IMAGE_MODALITIES else " **unknown modality**"
            conf = "" if m in cfg["image_modalities"] else " (not in config: ignored)"
            lines.append(f"- `{m}`: {len(g)} images, {g['unit_id'].nunique()} units, {missing_files} files missing{flag}{conf}")
            if missing_files:
                request.append(f"{missing_files} `{m}` image files listed in images.csv are not on disk.")
        n_maps = int(im["target_map_path"].notna().sum()) if "target_map_path" in im.columns else 0
        lines.append(f"- target disease maps (for map diffusion): {n_maps}")

    # 6. Label audit ----------------------------------------------------------------------------
    h("6. Label audit against hard constraints")
    lvl_of = {a: s.get("level") for a, s in cfg["assays"].items() if s["kind"] in ("binary", "fraction")}
    complete = {a for a, s in cfg["assays"].items() if s.get("se", 1.0) >= 1.0 and a in lvl_of}
    o = obs[obs["assay"].isin(lvl_of)].copy()
    o["pos"] = pd.to_numeric(o["value"], errors="coerce") > 0
    o["level"] = o["assay"].map(lvl_of)
    o["day"] = np.floor(o["t"])
    inter = set(units.loc[units.get("inoculation_method", pd.Series("", index=units.index)).astype(str)
                          .isin(COURT_INTERVENTIONS), "unit_id"].astype(str))
    n_bad = 0
    for (uid, day), g in o.groupby(["unit_id", "day"]):
        pos = set(g.loc[g["pos"], "level"])
        neg_complete = set(g.loc[~g["pos"] & g["assay"].isin(complete), "level"])
        # Inoculated units have the court set by intervention, so only apothecia -> spore applies.
        pairs = (("apothecia", "spore"),) if uid in inter else \
            (("infection", "court"), ("court", "spore"), ("apothecia", "spore"))
        n_bad += sum(hi in pos and lo in neg_complete for hi, lo in pairs)
    lines.append(f"- same-day readings where a higher level is positive but a complete (se = 1) assay "
                 f"reads the level it entails as negative: **{n_bad}**" +
                 (" -> either the assay is not complete (lower its `se`) or the label is wrong." if n_bad else ""))

    # 7. Trainable components ----------------------------------------------------------------
    h("7. What can be trained")
    v = pd.to_numeric(obs["value"], errors="coerce")
    def classes(a):
        x = v[obs["assay"] == a].dropna()
        return int((x > 0).sum()), int((x <= 0).sum())
    level_data = {}
    for a, lvl in lvl_of.items():
        p, n = classes(a)
        if p + n:
            level_data.setdefault(lvl, [0, 0])
            level_data[lvl][0] += p
            level_data[lvl][1] += n
    ok["chain_head"] = any(p >= MIN_PER_CLASS and n >= MIN_PER_CLASS for p, n in level_data.values())
    cont = [a for a, s in cfg["assays"].items() if s["kind"] in ("fraction", "continuous") and a in counts.index]
    ok["regression_head"] = bool(cont)
    ok["vector_diffusion"] = bool(cont) and bool(cfg["diffusion"].get("vector"))
    ok["map_diffusion"] = (not im.empty and "target_map_path" in im.columns and im["target_map_path"].notna().sum() >= 50)
    for k, val in ok.items():
        lines.append(f"- {k}: {'YES' if val else 'no'}")
    lines.append(f"- level coverage (positives, negatives): {level_data}")
    observed = set(level_data)
    unobserved = [l for l in LEVEL_ORDER if l not in observed]
    if unobserved:
        lines.append(f"- levels with no direct readings: {unobserved}. Their step kernels are only identified "
                     "through the product with the levels above them (Corollary 8.3 applied to parameters): "
                     "the model can fit P(infection) but cannot say which step failed.")
        if {"apothecia", "spore"} & set(unobserved):
            request.append("Inoculum-level readings (apothecia scouting and/or spore-trap qPCR) for at least "
                           "some field-seasons, so the lower steps of the chain are identifiable.")
        if "court" in unobserved:
            request.append("Petal-colonization readings (petal assay / qPCR at R1-R3) on a subset of fields.")
    if not ok["map_diffusion"]:
        request.append("For the spatial (map) diffusion head: >= 50 drone/satellite/field images with a "
                       "co-registered disease map (`target_map_path`, values 0-1 per pixel or grid cell).")

    h("DATA REQUEST")
    lines.extend(f"{i}. {r}" for i, r in enumerate(request, 1)) if request else lines.append("Nothing essential missing.")
    return lines, ok


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--report", default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    lines, _ = validate(a.data, cfg)
    text = "\n".join(lines) + "\n"
    print(text)
    Path(a.report or Path(a.data) / "data_report.md").write_text(text)


if __name__ == "__main__":
    main()
