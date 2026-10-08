"""How many isolates does the fixed-program vs host-responsive test need? (simulation)

    python -m ssaggr.power --isolates 6 12 24 --seeds 5

For each panel size and each true regime, synthetic panels are simulated and analysed exactly as
the real data are: leave-one-isolate-out R2 per host, and cross-host transfer (fit on one host,
rank isolates on the other). Under a fixed program transfer should succeed; under host-responsive
regulation it should degrade. The output is the distribution of these statistics, i.e. how well a
panel of that size can tell the two regimes apart.

The effect size and noise levels are arbitrary simulation settings (see synthetic.py); vary
--effect to see how conclusions depend on them. This is a design aid, not an estimate.
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np
import pandas as pd

from . import cv
from .config import load_config
from .panel import build_panel
from .synthetic import write


def one(n_iso: int, regime: str, seed: int, effect: float, genes: int) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        write(tmp, n_iso=n_iso, n_genes=genes, regime=regime, seed=seed, effect=effect)
        cfg = load_config(f"{tmp}/config.yaml")
        cfg["model"]["lambdas"] = [0.5, 0.25, 0.1]
        p = build_panel(cfg)
        s = cv.summarize(cv.loio(p, p.hosts, "fixed", cfg))
        x = cv.cross_host(p, cfg)
    return {"isolates": n_iso, "true_regime": regime, "seed": seed,
            "loio_r2_mean": float(s["r2_loio"].mean()), "transfer_spearman_mean": float(x["spearman"].mean())}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--isolates", type=int, nargs="+", default=[6, 12, 24])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--effect", type=float, default=2.0)
    ap.add_argument("--genes", type=int, default=800)
    ap.add_argument("--out", default="runs/power.csv")
    a = ap.parse_args(argv)
    rows = [one(n, r, s, a.effect, a.genes) for n in a.isolates for r in ("fixed", "responsive") for s in range(a.seeds)]
    df = pd.DataFrame(rows)
    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    df.to_csv(a.out, index=False)
    summ = df.groupby(["isolates", "true_regime"])[["loio_r2_mean", "transfer_spearman_mean"]].agg(["mean", "std"])
    print(summ.round(3).to_string())
    # separation of the two regimes by the transfer statistic (AUROC over seeds)
    for n, g in df.groupby("isolates"):
        f = g.loc[g["true_regime"] == "fixed", "transfer_spearman_mean"].to_numpy()
        r = g.loc[g["true_regime"] == "responsive", "transfer_spearman_mean"].to_numpy()
        auc = np.mean([(fi > ri) + 0.5 * (fi == ri) for fi in f for ri in r]) if len(f) and len(r) else np.nan
        print(f"{n} isolates: P(transfer higher under fixed than responsive) = {auc:.2f}")


if __name__ == "__main__":
    main()
