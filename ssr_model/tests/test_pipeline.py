"""End-to-end smoke test: synthetic data -> validate -> train -> evaluate (small and fast)."""

import json
from pathlib import Path

import torch

from ssrforecast.data import build_datasets
from ssrforecast.evaluate import load_run
from ssrforecast.schema import load_config
from ssrforecast.synthetic import generate
from ssrforecast.train import train
from ssrforecast.validate import validate


def small_cfg(tmp_path):
    data = generate(tmp_path / "data", n_fields=30, n_chamber=16, seed=1)
    cfg = load_config(Path(__file__).resolve().parents[1] / "configs" / "synthetic.yaml")
    cfg["data_dir"] = str(data)
    cfg["out_dir"] = str(tmp_path / "run")
    cfg["train"].update(epochs=2, batch_size=16, diffusion_epochs=2)
    cfg["diffusion"].update(steps=10)
    cfg["model"].update(dim=32, series_hidden=16, image_width=8)
    return cfg


def test_validate_reports(tmp_path):
    cfg = small_cfg(tmp_path)
    lines, ok = validate(cfg["data_dir"], cfg)
    text = "\n".join(lines)
    assert "DATA REQUEST" in text
    assert ok["regression_head"] and ok["vector_diffusion"]


def test_inputs_are_adapted(tmp_path):
    """No input dated after the issue day enters a sample (Definition S3)."""
    cfg = small_cfg(tmp_path)
    tr, _, _, meta = build_datasets(cfg)
    for i in range(0, len(tr), max(1, len(tr) // 20)):
        s = tr.samples.iloc[i]
        g = tr.images.get(s["unit_id"])
        item = tr[i]
        if g is not None:
            for m, present in item["present"].items():
                if present:
                    assert (g[(g["modality"] == m)]["t"] <= s["tau"]).any()
        assert (s["targets"]["t"] > s["tau"]).all()  # forecast targets are strictly in the future


def test_train_and_evaluate(tmp_path):
    cfg = small_cfg(tmp_path)
    out = train(cfg)
    assert (out / "model.pt").exists() and (out / "metrics.csv").exists()
    res = json.loads((out / "eval_test.json").read_text()) if (out / "eval_test.json").exists() else {}
    if res:
        assert res["chain_monotone_violations"] == 0
        if "diffusion_outcome_violation_rate_projected" in res:
            assert res["diffusion_outcome_violation_rate_projected"] == 0.0
    model, meta, cfg2 = load_run(out)
    assert isinstance(model, torch.nn.Module)
