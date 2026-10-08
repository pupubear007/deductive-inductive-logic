"""Configuration: defaults plus a YAML file. Every path is set in the YAML, never hard-coded."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

DEFAULTS: dict[str, Any] = {
    "seed": 0,
    "out_dir": "runs/default",
    # --- inputs ---------------------------------------------------------------------------
    "counts": {},            # host code -> pathogen-only count matrix (genes x samples, TSV)
    "library_totals": None,  # CSV host,sample,host_total_counts (for pathogen transcript share)
    "phenotype_table": None,  # mean sAUDPC: long (crop, isolate, mean) like CH4 TableS1_long.csv, or isolate x crop wide
    "per_plant_saudpc": None,  # per-plant sAUDPC (crop, isolate, sAUDPC), e.g. CH4 sAUDPC_per_plant.csv
    "host_to_crop": {"Gm": "Soybean", "Ha": "Sunflower"},  # count-matrix host -> phenotype column
    "lesion_files": {},      # crop -> per-plant lesion time course CSV (optional, for ANOVA)
    "in_vitro": {},          # trait name -> {file, isolate_col, value_col, filter (optional)}
    "annotation": None,      # gene annotation TSV (SS1G id, effector / CAZyme / signal peptide)
    "annotation_id_col": "SS1G",
    "isolate_aliases": {"Ss1980": "1980", "SsPotter": "SSPotter"},
    "classes": {             # a-priori classes on the RNA-seq panel (CH3 / CH4)
        "WISS47": "low", "JS659": "low", "MNSS6": "high", "Xtra7": "high",
        "SSPotter": "soy_adapted", "BN172": "sun_adapted",
    },
    "rep_remap": {"4": "1", "5": "2", "6": "3"},
    # --- preprocessing ----------------------------------------------------------------------
    "timepoints": ["24", "48", "96"],  # hpi in the panel
    # every set is analysed separately; a set with several time points is pooled (time-specific
    # intercepts, shared gene weights). Default: each time point alone, then all three pooled.
    "timepoint_sets": [["24"], ["48"], ["96"], ["24", "48", "96"]],
    "n_jobs": 2,             # parallel processes for the permutation null
    "min_count": 10,         # keep genes with >= min_count in >= min_samples libraries (CH4)
    "min_samples": 3,
    "top_variable_genes": 2000,  # unsupervised filter, recomputed inside every CV fold
    "colonization_covariate": True,  # include logit(pathogen transcript share) as a covariate
    "thin_to_equal_depth": False,    # binomial thinning of pathogen reads (proposal 5.6)
    # --- phenotype ----------------------------------------------------------------------------
    "target_scale": "z",     # z: sAUDPC standardised within host (proposal); raw: centred sAUDPC; rank; log
    "stability_scale": "raw",  # raw sAUDPC (Eberhart-Russell) or log (Finlay-Wilkinson's log scale)
    # Your hypothesis, stated before the expression analysis: either explicit classes ...
    "hypothesis_classes": {},  # e.g. {MNSS6: consistent, Xtra7: consistent, SSPotter: host_variable, ...}
    # ... or cut-offs on the stability statistics (any subset); isolates are then classified by them
    "stability_cutoffs": {},   # e.g. {mean_rank_min_high: 9, rank_range_max_consistent: 5}
    # --- model ----------------------------------------------------------------------------------
    "model": {
        "lambdas": [1.0, 0.5, 0.25, 0.1, 0.05],  # penalty path, relative to lambda_max
        "host_penalty_ratio": 2.0,  # host-specific deviations v are penalised this much more
        "l1_ratio": 0.5,          # elastic-net mixing: 1 = lasso, smaller = more grouping
        "ridge": 1e-3,
        "max_iter": 5000,         # coordinate-descent sweeps
        "tol": 1e-4,              # duality-gap tolerance (scikit-learn)
    },
    "n_permutations": 100,   # isolate-label permutations for the null distribution (720 = all for 6)
    "n_bootstrap": 50,       # bootstrap resamples of libraries for stability selection
    "selection_threshold": 0.6,
}


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict) and k not in ("counts", "lesion_files", "in_vitro", "classes"):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(path: str | Path | None = None, **over: Any) -> dict[str, Any]:
    cfg = copy.deepcopy(DEFAULTS)
    base = Path(".")
    if path is not None:
        with open(path) as f:
            cfg = _merge(cfg, yaml.safe_load(f) or {})
        base = Path(path).resolve().parent
    cfg = _merge(cfg, over)
    cfg["_config_dir"] = str(base)
    cfg["timepoints"] = [str(t) for t in cfg["timepoints"]]
    return cfg


def resolve(cfg: dict, p: str | None) -> Path | None:
    """Paths in the YAML are relative to the YAML file unless absolute."""
    if p is None:
        return None
    q = Path(p).expanduser()
    return q if q.is_absolute() else Path(cfg["_config_dir"]) / q
