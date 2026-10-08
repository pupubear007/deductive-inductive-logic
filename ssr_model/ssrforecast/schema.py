"""Input-data schema and configuration.

A dataset is a directory with up to four CSV files (only ``units.csv`` and ``observations.csv``
are mandatory) plus the image / map files they point to::

    units.csv         one row per experimental unit (pot, plant, plot, field-season)
    observations.csv  long format: one row per assay reading (disease assessment, qPCR, scouting ...)
    timeseries.csv    long format: environment readings (weather station, canopy sensors, chamber logs)
    images.csv        one row per image: growth chamber, greenhouse, field, drone, satellite

See DATA_REQUIREMENTS.md for the meaning of every column.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

REQUIRED_COLUMNS: dict[str, list[str]] = {
    "units.csv": ["unit_id", "experiment_type", "season"],
    "observations.csv": ["unit_id", "assay", "value"],  # plus 'date' or 'dpi'
    "timeseries.csv": ["unit_id", "variable", "value"],  # plus 'date' or 'dpi'
    "images.csv": ["unit_id", "modality", "path"],  # plus 'date' or 'dpi'; optional 'target_map_path'
}
OPTIONAL_UNIT_COLUMNS = [
    "site", "cultivar", "isolate", "inoculation_method", "inoculation_date", "planting_date",
    "lat", "lon",
]
EXPERIMENT_TYPES = {"growth_chamber", "greenhouse", "field"}
IMAGE_MODALITIES = {"gc_rgb", "greenhouse_rgb", "field_rgb", "drone_rgb", "drone_ms", "satellite"}

# Assay kinds:
#   binary      -> reading in {0,1}; observes a mechanism level through (se, sp)
#   fraction    -> reading in [0,1] (e.g. incidence); value > 0 observes `level` (hurdle model)
#   continuous  -> real reading; if `conditional_on` is set it is 0 whenever that level is 0
ASSAY_KINDS = {"binary", "fraction", "continuous"}

DEFAULT_CONFIG: dict[str, Any] = {
    "seed": 0,
    "data_dir": "data/synthetic",
    "out_dir": "runs/default",
    "task": "forecast",  # forecast: targets strictly after the issue date; nowcast: on/before it
    "horizon_days": 21,
    "issue_every_days": 7,  # also issue a forecast every k days (besides observation days)
    "lookback_days": 60,
    "split": {"group_by": "season", "val_frac": 0.15, "test_frac": 0.15},
    "series_variables": ["air_temp", "rh", "rain", "leaf_wetness", "soil_moisture"],
    "tabular": {
        "numeric": ["lat", "lon", "days_since_start"],
        "categorical": ["experiment_type", "cultivar", "isolate", "inoculation_method"],
    },
    "image_modalities": {
        "gc_rgb": {"channels": 3, "size": 64},
        "greenhouse_rgb": {"channels": 3, "size": 64},
        "field_rgb": {"channels": 3, "size": 64},
        "drone_rgb": {"channels": 3, "size": 64},
        "drone_ms": {"channels": 5, "size": 64},
        "satellite": {"channels": 4, "size": 32},
    },
    # Sensitivity / specificity of each assay for its level. 1.0 / 1.0 is the paper's idealised
    # (sound and complete) assay. Replace with measured values, or set learn: true.
    "assays": {
        "apothecia_scouting": {"kind": "binary", "level": "apothecia", "se": 1.0, "sp": 1.0},
        "spore_trap_qpcr": {"kind": "binary", "level": "spore", "se": 1.0, "sp": 1.0},
        "petal_qpcr": {"kind": "binary", "level": "court", "se": 1.0, "sp": 1.0},
        "stem_qpcr": {"kind": "binary", "level": "infection", "se": 1.0, "sp": 1.0},
        "incidence": {"kind": "fraction", "level": "infection", "se": 1.0, "sp": 1.0},
        "dsi": {"kind": "continuous", "conditional_on": "infection", "scale": 100.0},
        "lesion_length_mm": {"kind": "continuous", "conditional_on": "infection", "scale": 100.0},
        "yield_kg_ha": {"kind": "continuous", "scale": 5000.0},
    },
    "model": {"dim": 128, "series_hidden": 64, "image_width": 32, "modality_dropout": 0.2},
    "diffusion": {"vector": True, "map": True, "steps": 100, "map_size": 32, "hidden": 256},
    "loss_weights": {"chain": 1.0, "regression": 1.0, "vector_diffusion": 0.5, "map_diffusion": 0.5},
    "train": {"epochs": 30, "batch_size": 32, "lr": 1e-3, "weight_decay": 1e-4, "device": "auto",
              "num_workers": 0, "patience": 8, "diffusion_epochs": 30},
}

# Map continuous assays to the outcome names used by the generative (vector diffusion) head.
OUTCOME_OF_ASSAY = {"incidence": "incidence", "dsi": "dsi", "lesion_length_mm": "lesion_length", "yield_kg_ha": "yield"}


def _merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict) and k not in ("assays", "image_modalities"):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(path: str | Path | None = None, **overrides: Any) -> dict[str, Any]:
    """Load a YAML config on top of DEFAULT_CONFIG. ``assays`` and ``image_modalities`` replace
    the defaults wholesale when given, so a config lists exactly the assays it uses."""
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    if path is not None:
        with open(path) as f:
            cfg = _merge(cfg, yaml.safe_load(f) or {})
    cfg = _merge(cfg, overrides)
    for name, a in cfg["assays"].items():
        if a.get("kind") not in ASSAY_KINDS:
            raise ValueError(f"assay {name!r}: kind must be one of {sorted(ASSAY_KINDS)}")
    return cfg


def binary_assays(cfg: dict) -> list[str]:
    """Assays whose reading observes a mechanism level (binary, and the positive part of fractions)."""
    return [n for n, a in cfg["assays"].items() if a["kind"] in ("binary", "fraction")]


def continuous_assays(cfg: dict) -> list[str]:
    return [n for n, a in cfg["assays"].items() if a["kind"] in ("fraction", "continuous")]
