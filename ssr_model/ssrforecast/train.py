"""Train the SSR model.

    python -m ssrforecast.train --config configs/synthetic.yaml
    python -m ssrforecast.train --config configs/default.yaml --data /path/to/your/data --out runs/real1

Before training on real data, run ``python -m ssrforecast.validate --data DIR --config CFG``;
it reports which inputs are missing and which model components can be trained.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from dataclasses import asdict
from pathlib import Path

import numpy as np
import torch
import yaml
from torch.utils.data import DataLoader

from .data import build_datasets, collate
from .evaluate import evaluate, to_device
from .models.full import SSRModel, compute_losses
from .schema import load_config


def set_seed(s: int) -> None:
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)


def pick_device(name: str) -> torch.device:
    if name == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(name)


@torch.no_grad()
def val_loss(model, ds, cfg, dev) -> float:
    if len(ds) == 0:
        return float("nan")
    model.eval()
    tot, n = 0.0, 0
    for b in DataLoader(ds, batch_size=cfg["train"]["batch_size"], collate_fn=collate):
        b = to_device(b, dev)
        out = model(b)
        L = compute_losses(model, b, out, cfg)
        # model selection uses the deterministic heads only (diffusion losses are noisy)
        v = sum(L[k] for k in ("chain", "regression") if k in L)
        tot += float(v) * len(b["index"])
        n += len(b["index"])
    return tot / max(n, 1)


def train_diffusion(model: SSRModel, ds, cfg: dict, dev: torch.device, out: Path) -> None:
    """Stage 2: fit the diffusion heads on the frozen embedding z of the selected stage-1 model.
    z is computed once (eval mode), so this stage is cheap."""
    heads = [h for h in (model.vdiff, model.mdiff) if h is not None]
    epochs = int(cfg["train"].get("diffusion_epochs", 0))
    if not heads or epochs <= 0 or len(ds) == 0:
        return
    model.eval()
    Z, YC, MC, YM, MM = [], [], [], [], []
    with torch.no_grad():
        for b in DataLoader(ds, batch_size=cfg["train"]["batch_size"], collate_fn=collate):
            b = to_device(b, dev)
            Z.append(model.embed(b)); YC.append(b["y_cont"]); MC.append(b["m_cont"])
            YM.append(b["y_map"]); MM.append(b["m_map"])
    Z, YC, MC, YM, MM = (torch.cat(x) for x in (Z, YC, MC, YM, MM))
    params = [p for h in heads for p in h.parameters()]
    opt = torch.optim.AdamW(params, lr=cfg["train"]["lr"], weight_decay=cfg["train"]["weight_decay"])
    bs, n = cfg["train"]["batch_size"], len(Z)
    w = cfg["loss_weights"]
    for h in heads:
        h.train()
    with open(out / "diffusion_metrics.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["epoch", "vector_diffusion", "map_diffusion"])
        for ep in range(1, epochs + 1):
            perm = torch.randperm(n, device=Z.device)
            tv, tm, k = 0.0, 0.0, 0
            for i in range(0, n, bs):
                j = perm[i:i + bs]
                lv = model.vdiff.loss(YC[j].clamp(-0.5, 1.5), MC[j], Z[j]) if model.vdiff is not None else Z.new_zeros(())
                lm = model.mdiff.loss(YM[j], MM[j], Z[j]) if model.mdiff is not None else Z.new_zeros(())
                loss = w.get("vector_diffusion", 1.0) * lv + w.get("map_diffusion", 1.0) * lm
                opt.zero_grad()
                loss.backward()
                opt.step()
                tv, tm, k = tv + float(lv.detach()), tm + float(lm.detach()), k + 1
            wr.writerow([ep, round(tv / k, 5), round(tm / k, 5)])
            if ep == 1 or ep % 10 == 0 or ep == epochs:
                print(f"diffusion epoch {ep:3d}  vector {tv / k:.4f}  map {tm / k:.4f}")
    model.eval()


def train(cfg: dict) -> Path:
    set_seed(cfg["seed"])
    out = Path(cfg["out_dir"])
    out.mkdir(parents=True, exist_ok=True)
    dev = pick_device(cfg["train"]["device"])
    tr, va, te, meta = build_datasets(cfg)
    print(f"samples: train={len(tr)} val={len(va)} test={len(te)} | device={dev}")
    print(f"binary assays: {meta.bin_names}\ncontinuous assays: {meta.cont_names}\nseries: {meta.series_vars}")

    model = SSRModel(meta, cfg).to(dev)
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=cfg["train"]["lr"], weight_decay=cfg["train"]["weight_decay"])
    loader = DataLoader(tr, batch_size=cfg["train"]["batch_size"], shuffle=True, collate_fn=collate,
                        num_workers=cfg["train"]["num_workers"])
    (out / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))

    best, bad = float("inf"), 0
    with open(out / "metrics.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["epoch", "train_total", "train_chain", "train_regression", "val_loss"])
        for ep in range(1, cfg["train"]["epochs"] + 1):
            model.train()
            agg: dict[str, float] = {}
            nb = 0
            for b in loader:
                b = to_device(b, dev)
                L = compute_losses(model, b, model(b), cfg)
                opt.zero_grad()
                L["total"].backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                opt.step()
                for k, v in L.items():
                    agg[k] = agg.get(k, 0.0) + float(v.detach())
                nb += 1
            agg = {k: v / max(nb, 1) for k, v in agg.items()}
            vl = val_loss(model, va, cfg, dev)
            wr.writerow([ep] + [round(agg.get(k, float("nan")), 5) for k in ("total", "chain", "regression")]
                        + [round(vl, 5)])
            f.flush()
            print(f"epoch {ep:3d}  train {agg.get('total', float('nan')):.4f}  val {vl:.4f}")
            score = vl if vl == vl else agg.get("total", 0.0)  # fall back to train loss if no val set
            if score < best - 1e-5:
                best, bad = score, 0
                torch.save({"state_dict": model.state_dict(), "meta": asdict(meta), "cfg": cfg}, out / "model.pt")
            else:
                bad += 1
                if bad >= cfg["train"]["patience"]:
                    print("early stopping")
                    break

    ck = torch.load(out / "model.pt", map_location=dev, weights_only=False)
    model.load_state_dict(ck["state_dict"])
    train_diffusion(model, tr, cfg, dev, out)
    torch.save({"state_dict": model.state_dict(), "meta": asdict(meta), "cfg": cfg}, out / "model.pt")
    for name, ds in (("val", va), ("test", te)):
        if len(ds):
            res, pred = evaluate(model, ds, cfg, dev, n_samples=2)
            (out / f"eval_{name}.json").write_text(json.dumps(res, indent=2, default=float))
            pred.to_csv(out / f"predictions_{name}.csv", index=False)
    print(f"done -> {out}")
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--data", default=None, help="override data_dir")
    ap.add_argument("--out", default=None, help="override out_dir")
    ap.add_argument("--epochs", type=int, default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    if a.data:
        cfg["data_dir"] = a.data
    if a.out:
        cfg["out_dir"] = a.out
    if a.epochs:
        cfg["train"]["epochs"] = a.epochs
    train(cfg)


if __name__ == "__main__":
    main()
