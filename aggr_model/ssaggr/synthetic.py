"""Synthetic dual RNA-seq + phenotype data with a *known* answer, in the same file layout as the
real data (count TSVs named <hpi>_<rep>_<host>_<isolate>, library totals, phenotype table).

    python -m ssaggr.synthetic --out data/synthetic --regime responsive --isolates 6

Two regimes:
  fixed       one set of determinant genes tracks the isolate's host-independent aggressiveness
              on every host (a more aggressive isolate runs the same program harder);
  responsive  each host has its own determinant genes, and aggressiveness itself has a
              host-specific component (isolates reorder between hosts, as in CH4).

ALL PARAMETERS ARE ARBITRARY SIMULATION SETTINGS, NOT ESTIMATES FROM YOUR DATA. The generator
exists to check that the pipeline recovers a planted answer and to explore power (how many
isolates are needed), never to stand in for biology.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def simulate(n_iso: int = 6, n_genes: int = 1500, hosts=("Gm", "Ha"), reps: int = 3, hpi=("48",),
             regime: str = "responsive", n_det: int = 30, effect: float = 2.0, host_specific_sd: float = 0.8,
             seed: int = 0):
    rng = np.random.default_rng(seed)
    isos = [f"ISO{i + 1:02d}" for i in range(n_iso)]
    f = rng.normal(0, 1, n_iso)  # host-independent aggressiveness
    a = {h: f + (rng.normal(0, host_specific_sd, n_iso) if regime == "responsive" else 0.0) for h in hosts}
    mu = rng.normal(4, 1.5, n_genes)
    iso_noise = rng.normal(0, 0.3, (n_iso, n_genes))      # isolate identity unrelated to aggressiveness
    host_eff = {h: rng.normal(0, 0.5, n_genes) for h in hosts}
    det_shared = rng.choice(n_genes, n_det, replace=False)
    det = {}
    for h in hosts:
        det[h] = det_shared if regime == "fixed" else rng.choice(np.setdiff1d(np.arange(n_genes), det_shared), n_det, replace=False)
    sign = rng.choice([-1.0, 1.0], n_genes)
    counts, totals = {}, []
    for h in hosts:
        cols = {}
        for t in hpi:
            for r in range(1, reps + 1):
                for i, iso in enumerate(isos):
                    eta = mu + host_eff[h] + iso_noise[i]
                    driver = f[i] if regime == "fixed" else a[h][i]
                    eta[det[h]] += effect * sign[det[h]] * driver
                    eta += rng.normal(0, 0.25, n_genes)  # replicate noise
                    lib = rng.uniform(0.7, 1.3)
                    m = np.exp(eta) * lib
                    k = 10.0
                    cols[f"{t}_{r}_{h}_{iso}"] = rng.negative_binomial(k, k / (k + m))
                    share = 1 / (1 + np.exp(-(-2.5 - 0.3 * f[i] + rng.normal(0, 0.3))))
                    p_tot = cols[f"{t}_{r}_{h}_{iso}"].sum()
                    totals.append({"host": h, "sample": f"{t}_{r}_{h}_{iso}",
                                   "host_total_counts": float(p_tot * (1 - share) / share)})
        df = pd.DataFrame(cols, index=[f"SS1G_{g:05d}" for g in range(n_genes)])
        df.index.name = "gene_id"
        counts[h] = df
    crop = {"Gm": "Soybean", "Ha": "Sunflower"}
    pheno = pd.DataFrame({crop.get(h, h): np.exp(3 + 0.5 * a[h]) for h in hosts}, index=isos)
    truth = {"regime": regime, "aggressiveness": {h: dict(zip(isos, map(float, a[h]))) for h in hosts},
             "determinants": {h: [f"SS1G_{g:05d}" for g in det[h]] for h in hosts}}
    return counts, pd.DataFrame(totals), pheno, truth


def write(out: str | Path, **kw) -> Path:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    counts, totals, pheno, truth = simulate(**kw)
    for h, df in counts.items():
        df.to_csv(out / f"pathogen_counts_{h}.tsv", sep="\t")
    totals.to_csv(out / "library_totals.csv", index=False)
    pheno.to_csv(out / "phenotype.csv")
    (out / "truth.json").write_text(json.dumps(truth, indent=2))
    hosts = list(counts)
    crop = {"Gm": "Soybean", "Ha": "Sunflower"}
    cfg = {
        "out_dir": "run",
        "counts": {h: f"pathogen_counts_{h}.tsv" for h in hosts},
        "library_totals": "library_totals.csv",
        "phenotype_table": "phenotype.csv",
        "host_to_crop": {h: crop.get(h, h) for h in hosts},
        "classes": {},
        "timepoints": list(kw.get("hpi", ("48",))),
        "timepoint_sets": [list(kw.get("hpi", ("48",)))],
    }
    import yaml
    (out / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    (out / "README.txt").write_text("SYNTHETIC DATA from ssaggr.synthetic; parameters are arbitrary simulation settings.\n")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--regime", choices=["fixed", "responsive"], default="responsive")
    ap.add_argument("--isolates", type=int, default=6)
    ap.add_argument("--genes", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--effect", type=float, default=2.0, help="log-expression change per unit aggressiveness")
    a = ap.parse_args(argv)
    p = write(a.out, n_iso=a.isolates, n_genes=a.genes, regime=a.regime, seed=a.seed, effect=a.effect)
    print(f"synthetic data and config written to {p}")


if __name__ == "__main__":
    main()
