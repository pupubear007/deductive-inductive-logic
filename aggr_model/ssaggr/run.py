"""Run the whole analysis and write a report.

    python -m ssaggr.run --config configs/real_local.yaml [--quick]

Outputs in ``out_dir``: stability.csv, responsiveness.csv, loio_<model>.csv, permutation_*.csv,
cross_host.csv, determinants.csv, resolution.csv and report.md.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import cv
from .config import load_config
from .io import load_annotation, load_in_vitro, load_lesions, load_phenotype
from .panel import build_panel
from .resolution import resolution_report
from .responsiveness import exact_spearman, host_responsiveness, residualize
from .stability import classify, interaction_share, stability_table


def fmt(df: pd.DataFrame, floatfmt: str = "{:.3f}") -> str:
    if df is None or len(df) == 0:
        return "_(none)_\n"
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    cols = [d.index.name or ""] + list(d.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for idx, r in d.iterrows():
        lines.append("| " + " | ".join([str(idx)] + [str(v) for v in r.to_list()]) + " |")
    return "\n".join(lines) + "\n"


def run(cfg: dict, quick: bool = False, skip_perm: bool = False) -> Path:
    t0 = time.time()
    out = Path(cfg["out_dir"])
    if not out.is_absolute():
        out = Path(cfg["_config_dir"]) / out
    out.mkdir(parents=True, exist_ok=True)
    if quick:
        cfg["n_permutations"] = min(cfg["n_permutations"], 20)
        cfg["n_bootstrap"] = min(cfg["n_bootstrap"], 12)
    R: list[str] = ["# Aggressiveness determinants: analysis report", "",
                    f"Time points: {cfg['timepoints']} hpi. Target: {cfg['target_scale']} sAUDPC, centred within host.", ""]

    # 1. Objective 1: stability across crops -----------------------------------------------------
    pheno = load_phenotype(cfg)
    st = stability_table(pheno)
    st["stability_class"] = classify(st)
    st.to_csv(out / "stability.csv")
    R += ["## 1. Aggressiveness and its consistency across hosts (Objective 1)", "",
          "Joint regression on the host index, log sAUDPC. `stability_class` uses median cut-offs and "
          "is provisional: fix the cut-offs before looking at expression data.", "", fmt(st)]
    les = load_lesions(cfg)
    if len(les):
        share = interaction_share(les)
        R += ["Share of per-plant variation in log(sAUDPC+1) (sequential SS, no block term): " +
              ", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in share.items()), ""]

    # 2. Panel and Objective 2 index -------------------------------------------------------------
    p = build_panel(cfg)
    R += ["## 2. RNA-seq panel", "", f"{len(p.isolates)} isolates x {len(p.hosts)} hosts, "
          f"{p.meta.shape[0]} libraries, {p.expr.shape[0]} genes after filtering (shared between hosts).", "",
          "Targets (centred log sAUDPC):", "", fmt(p.target)]
    resp = host_responsiveness(residualize(p.expr, p.meta), p.meta)
    if len(resp):
        resp.to_csv(out / "responsiveness.csv")
        R += ["## 3. Host-responsiveness index (Objective 2)", "",
              "R2_host = share of the isolate's pathogen-expression variance explained by host, after removing "
              "colonization. Host and sequencing batch are confounded in the current data, so only differences "
              "between isolates are interpretable.", "", fmt(resp)]
        for metric in ("rank_range", "deviation_ms", "mean_log_sAUDPC"):
            if metric in st.columns:
                r = exact_spearman(resp["R2_host"], st[metric])
                R.append(f"- Spearman(R2_host, {metric}) = {r['rho']:.2f}, exact permutation p = {r['p_perm']:.3f} (n = {r['n']})")
        R.append("")

    # 3. Prediction: LOIO per host, joint fixed vs responsive, transfer ---------------------------
    R += ["## 4. Held-out isolate prediction (leave one isolate out)", ""]
    per_host = []
    for h in p.hosts:
        pred = cv.loio(p, [h], "fixed", cfg)
        pred.to_csv(out / f"loio_single_{h}.csv", index=False)
        per_host.append(pred)
    single = pd.concat(per_host)
    s_single = cv.summarize(single)
    R += ["Single-host sparse models:", "", fmt(s_single)]
    if not skip_perm and cfg["n_permutations"] > 0:
        perms = []
        for h in p.hosts:
            perms.append(cv.permutation_test(p, [h], "fixed", cfg, single[single["host"] == h], cfg["n_permutations"], cfg["seed"]))
        perm = pd.concat(perms)
        perm.to_csv(out / "permutation_single.csv")
        R += ["Permutation null (isolate labels shuffled, whole procedure rerun):", "", fmt(perm)]
    joint = {}
    if len(p.hosts) > 1:
        for kind in ("fixed", "responsive"):  # noqa: B007
            pr = cv.loio(p, p.hosts, kind, cfg)
            pr.to_csv(out / f"loio_joint_{kind}.csv", index=False)
            joint[kind] = cv.summarize(pr)
        cmp = pd.concat({k: v[["r2_loio", "spearman", "mse"]] for k, v in joint.items()}, axis=1)
        cmp.columns = [f"{a}_{b}" for a, b in cmp.columns]
        R += ["Joint model over both hosts: shared weights (fixed program) versus shared + host-specific "
              "weights (host-responsive). Lower held-out MSE favours that hypothesis; with 6 isolates treat "
              "this as descriptive.", "", fmt(cmp)]
        xh = cv.cross_host(p, cfg)
        xh.to_csv(out / "cross_host.csv", index=False)
        R += ["Cross-host transfer (fit on one host, rank isolates on the other):", "", fmt(xh.set_index("train_host"))]

    # 4. Determinants ----------------------------------------------------------------------------
    sel = [cv.stability_selection(p, h, cfg, cfg["seed"]) for h in p.hosts]
    det = cv.concordance(sel, p.hosts, cfg["selection_threshold"])
    ann = load_annotation(cfg)
    if ann is not None:
        det = det.join(ann, how="left")
    det.to_csv(out / "determinants.csv")
    counts = det["class"].value_counts().to_dict()
    R += ["## 5. Candidate determinants (stability selection)", "",
          f"Selection frequency >= {cfg['selection_threshold']} across LOIO training sets x bootstrap resamples. "
          f"Classes: {counts}. Full table: determinants.csv.", "",
          fmt(det[det["class"] != "not_selected"].head(25))]

    # 5. Resolution (Theorem 8.2) ----------------------------------------------------------------
    scores = load_in_vitro(cfg)
    ts = p.meta.groupby(["isolate", "host"])["transcript_share"].mean().unstack()
    for h in ts.columns:
        scores[f"transcript_share_{h}"] = ts[h]
    for h, g in single.groupby("host"):
        scores[f"loio_expression_score_{h}"] = g.set_index("isolate")["predicted"]
    labels = pd.DataFrame(index=sorted(set(scores.index) | set(p.isolates)))
    apriori = pd.Series(cfg.get("classes") or {}, dtype=object)
    if len(apriori):
        hl = apriori[apriori.isin(["high", "low"])]
        labels["apriori_high_vs_low"] = hl
    for h in p.hosts:
        t = p.target[h]
        labels[f"median_split_{h}"] = pd.Series(np.where(t > t.median(), "high", "low"), index=t.index)
    res = resolution_report(scores, labels)
    res.to_csv(out / "resolution.csv", index=False)
    R += ["## 6. Which assays resolve the aggressiveness class? (Theorem 8.2 / Corollary 8.3)", "",
          "An assay resolves the class if one threshold separates the classes with no exception; each "
          "witness pair is two isolates the assay cannot tell apart that differ in class. p_chance is the "
          "probability that a random score would also resolve it: with 4-6 isolates, a 'resolves' is weak "
          "evidence on its own.", "",
          fmt(res.set_index("assay")[["host", "resolves", "misclassified", "n", "p_chance", "witnesses"]])]

    R += ["## Caveats", "",
          "- With 6 isolates the isolate is the unit of evidence: n = 6, whatever the number of libraries.",
          "- Permutation p-values cannot fall below 1/720 with 6 isolates.",
          "- Soybean and sunflower libraries were sequenced and quantified separately; host and batch are confounded.",
          "- The phenotype table must refer to the host genotype that was sequenced (check `host_to_crop`).",
          f"- Runtime {time.time() - t0:.0f} s."]
    (out / "report.md").write_text("\n".join(R) + "\n")
    (out / "config_used.json").write_text(json.dumps({k: v for k, v in cfg.items() if not k.startswith("_")}, indent=2, default=str))
    print(f"report -> {out / 'report.md'}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--quick", action="store_true", help="few permutations / bootstraps, for a first look")
    ap.add_argument("--skip-perm", action="store_true")
    a = ap.parse_args(argv)
    run(load_config(a.config), quick=a.quick, skip_perm=a.skip_perm)


if __name__ == "__main__":
    main()
