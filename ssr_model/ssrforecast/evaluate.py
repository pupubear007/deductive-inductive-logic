"""Evaluate a trained checkpoint on the held-out split.

    python -m ssrforecast.evaluate --run runs/default [--split test] [--samples 8]

Writes ``<run>/eval_<split>.json`` and ``<run>/predictions_<split>.csv``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from . import metrics as M
from .data import Meta, build_datasets, collate
from .mechanism import LEVELS, outcome_violations, project_outcomes
from .models.full import SSRModel


def to_device(b: dict[str, Any], dev: torch.device) -> dict[str, Any]:
    return {k: ({m: t.to(dev) for m, t in v.items()} if isinstance(v, dict) else v.to(dev)) for k, v in b.items()}


@torch.no_grad()
def evaluate(model: SSRModel, ds, cfg: dict[str, Any], dev: torch.device, n_samples: int = 0) -> tuple[dict, pd.DataFrame]:
    model.eval()
    meta = model.meta
    loader = DataLoader(ds, batch_size=cfg["train"]["batch_size"], shuffle=False, collate_fn=collate)
    P, YB, MB, EC, YC, MC, LV, idx = [], [], [], [], [], [], [], []
    viol_raw, viol_proj, map_mae = [], [], []
    for b in loader:
        b = to_device(b, dev)
        out = model(b)
        P.append(model.p_positive(out["chain"]).cpu()); YB.append(b["y_bin"].cpu()); MB.append(b["m_bin"].cpu())
        EC.append(model.expected_outcomes(out).cpu()); YC.append(b["y_cont"].cpu()); MC.append(b["m_cont"].cpu())
        LV.append(out["chain"].stack().cpu()); idx.append(b["index"].cpu())
        if n_samples and model.vdiff is not None:
            for _ in range(n_samples):
                y = model.vdiff.sample(out["z"])
                viol_raw.append(float(outcome_violations(y, meta.outcome_names)))
                viol_proj.append(float(outcome_violations(project_outcomes(y, meta.outcome_names), meta.outcome_names)))
        if n_samples and model.mdiff is not None and b["m_map"].sum() > 0:
            keep = b["m_map"] > 0
            s = model.mdiff.sample(out["z"][keep])
            map_mae.append(float((s - b["y_map"][keep]).abs().mean()))
    if not P:
        return {"n_samples": 0}, pd.DataFrame()
    P, YB, MB = torch.cat(P).numpy(), torch.cat(YB).numpy(), torch.cat(MB).numpy()
    EC, YC, MC = torch.cat(EC).numpy(), torch.cat(YC).numpy(), torch.cat(MC).numpy()
    LV, idx = torch.cat(LV).numpy(), torch.cat(idx).numpy()

    res: dict[str, Any] = {"n_samples": int(len(idx)), "assays": {}}
    for j, a in enumerate(meta.bin_names):
        m = MB[:, j] > 0
        p, y = P[m, j], YB[m, j]
        res["assays"][a] = {"n": int(m.sum()), "prevalence": float(y.mean()) if m.any() else None,
                            "brier": M.brier(p, y), "log_loss": M.log_loss(p, y), "auroc": M.auroc(p, y),
                            "ece": M.ece(p, y), **M.sens_exceeds_fpr(p, y)}
    for j, a in enumerate(meta.cont_names):
        m = MC[:, j] > 0
        sc = meta.cont_scale[j]
        res["assays"].setdefault(a, {})["mae_expected_value"] = M.mae(EC[m, j] * sc, YC[m, j] * sc)
    # Hard constraints hold by construction in the chain head; this is an audit.
    lv = dict(zip(LEVELS, LV.T))
    res["chain_monotone_violations"] = int(((lv["infection"] > lv["court"] + 1e-6)).sum())
    if viol_raw:
        res["diffusion_outcome_violation_rate_raw"] = float(np.mean(viol_raw))
        res["diffusion_outcome_violation_rate_projected"] = float(np.mean(viol_proj))
    if map_mae:
        res["map_diffusion_mae"] = float(np.mean(map_mae))

    s = ds.samples.iloc[idx]
    pred = pd.DataFrame({"unit_id": s["unit_id"].values, "issue_day": s["tau"].values,
                         "days_since_start": (s["tau"] - s["t0"]).values})
    for k, n in enumerate(LEVELS):
        pred[f"p_{n}"] = LV[:, k]
    for j, a in enumerate(meta.bin_names):
        pred[f"p_positive_{a}"] = P[:, j]
        pred[f"obs_{a}"] = np.where(MB[:, j] > 0, YB[:, j], np.nan)
    for j, a in enumerate(meta.cont_names):
        pred[f"expected_{a}"] = EC[:, j] * meta.cont_scale[j]
        pred[f"obs_{a}"] = np.where(MC[:, j] > 0, YC[:, j] * meta.cont_scale[j], np.nan)
    return res, pred


def load_run(run: str | Path, device: str = "cpu"):
    ck = torch.load(Path(run) / "model.pt", map_location=device, weights_only=False)
    meta = Meta(**ck["meta"])
    model = SSRModel(meta, ck["cfg"])
    model.load_state_dict(ck["state_dict"])
    return model.to(device), meta, ck["cfg"]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    ap.add_argument("--split", default="test", choices=["train", "val", "test"])
    ap.add_argument("--samples", type=int, default=4, help="diffusion draws per batch for the audit")
    ap.add_argument("--data", default=None, help="override data_dir")
    a = ap.parse_args(argv)
    dev = torch.device("cpu")
    model, meta, cfg = load_run(a.run)
    if a.data:
        cfg["data_dir"] = a.data
    tr, va, te, _ = build_datasets(cfg, meta)
    ds = {"train": tr, "val": va, "test": te}[a.split]
    res, pred = evaluate(model, ds, cfg, dev, a.samples)
    Path(a.run, f"eval_{a.split}.json").write_text(json.dumps(res, indent=2, default=float))
    pred.to_csv(Path(a.run, f"predictions_{a.split}.csv"), index=False)
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
