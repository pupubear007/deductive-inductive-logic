"""Run the whole analysis and write a report.

    python -m ssaggr.run --config configs/real_local.yaml [--quick] [--skip-perm]

Objective 1 (stability across crops) is computed once. Everything that uses expression is run
for each time-point set in ``timepoint_sets`` (default: 24, 48 and 96 hpi separately, then all
three pooled), each in its own sub-folder, with a summary table across sets at the top of
``report.md``.
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
from .io import load_annotation, load_expression, load_in_vitro, load_lesions, load_phenotype
from .panel import build_panel
from .resolution import resolution_report
from .responsiveness import exact_spearman, host_responsiveness, residualize
from .stability import classes_from_cutoffs, interaction_share, stability_table


def fmt(df: pd.DataFrame | None, floatfmt: str = "{:.3f}") -> str:
    if df is None or len(df) == 0:
        return "_(none)_\n"
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    cols = [d.index.name or ""] + [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for idx, r in d.iterrows():
        lines.append("| " + " | ".join([str(idx)] + [str(v) for v in r.to_list()]) + " |")
    return "\n".join(lines) + "\n"


def hypothesis_labels(cfg: dict, st: pd.DataFrame) -> pd.Series | None:
    if cfg.get("hypothesis_classes"):
        return pd.Series(cfg["hypothesis_classes"], dtype=object).rename("hypothesis_class")
    lab = classes_from_cutoffs(st, cfg.get("stability_cutoffs") or {})
    return lab.rename("hypothesis_class") if lab is not None else None


def analyse_set(cfg: dict, tps: list[str], counts, meta, pheno, st, hyp, out: Path, skip_perm: bool) -> tuple[list[str], dict]:
    tag = "+".join(tps)
    p = build_panel(cfg, counts, meta, pheno, timepoints=tps)
    R = [f"## Time points: {tag} hpi", "",
         f"{len(p.isolates)} isolates x {len(p.hosts)} hosts, {p.meta.shape[0]} libraries, "
         f"{p.expr.shape[0]} shared genes after filtering. "
         + ("Pooled: intercepts per host x time point, gene weights shared across time points." if len(tps) > 1 else ""), ""]
    summary: dict = {"timepoints": tag}

    # host-responsiveness index (per time point, then averaged for pooled sets)
    resp_tabs = []
    for t in tps:
        sub = p.meta.index[p.meta["hpi"].astype(str) == t]
        r = host_responsiveness(residualize(p.expr[sub], p.meta.loc[sub]), p.meta.loc[sub])
        r["hpi"] = t
        resp_tabs.append(r)
    resp = pd.concat(resp_tabs)
    resp.to_csv(out / "responsiveness.csv")
    rmean = resp.groupby(level=0)["R2_host"].mean()
    R += ["**Host-responsiveness index** (R²_host within isolate, colonization removed; host and batch are "
          "confounded, so compare isolates, not absolute values):", "",
          fmt(resp.reset_index().pivot(index="isolate", columns="hpi", values="R2_host"))]
    for metric in ("rank_range", "deviation_ms", "slope_b"):
        r = exact_spearman(rmean, st[metric])
        R.append(f"- Spearman(R²_host, {metric}) = {r['rho']:.2f}, exact p = {r['p_perm']:.3f} (n = {r['n']})")
    if hyp is not None:
        R.append(f"- R²_host by hypothesis class: " + ", ".join(
            f"{c}: {rmean[rmean.index.isin(hyp[hyp == c].index)].mean():.3f}" for c in sorted(hyp.dropna().unique())))
    R.append("")

    # prediction
    single = pd.concat([cv.loio(p, [h], "fixed", cfg) for h in p.hosts])
    single.to_csv(out / "loio_single.csv", index=False)
    s_single = cv.summarize(single)
    R += ["**Leave-one-isolate-out prediction, one host at a time:**", "", fmt(s_single)]
    for h in p.hosts:
        summary[f"r2_loio_{h}"] = s_single.loc[h, "r2_loio"]
    if not skip_perm and cfg["n_permutations"] > 0:
        perm = pd.concat([cv.permutation_test(p, [h], "fixed", cfg, single[single["host"] == h],
                                              cfg["n_permutations"], cfg["seed"]) for h in p.hosts])
        perm.to_csv(out / "permutation_single.csv")
        R += ["Permutation null (isolates shuffled, whole procedure rerun):", "", fmt(perm)]
        for h in p.hosts:
            summary[f"p_perm_{h}"] = perm.loc[h, "p_perm"]
    if len(p.hosts) > 1:
        joint = {k: cv.summarize(cv.loio(p, p.hosts, k, cfg)) for k in ("fixed", "responsive")}
        cmp = pd.concat({k: v[["r2_loio", "spearman"]] for k, v in joint.items()}, axis=1)
        cmp.columns = [f"{a}_{b}" for a, b in cmp.columns]
        cmp.to_csv(out / "joint_fixed_vs_responsive.csv")
        xh = cv.cross_host(p, cfg)
        xh.to_csv(out / "cross_host.csv", index=False)
        R += ["**Fixed program vs host-responsive** (joint model over both hosts, held-out R²):", "", fmt(cmp),
              "**Cross-host transfer** (fit on one host, rank isolates on the other):", "",
              fmt(xh.set_index("train_host"))]
        for r in xh.itertuples():
            summary[f"transfer_{r.train_host}->{r.test_host}"] = r.spearman

    # determinants
    det = cv.concordance([cv.stability_selection(p, h, cfg, cfg["seed"]) for h in p.hosts], p.hosts,
                         cfg["selection_threshold"])
    ann = load_annotation(cfg)
    if ann is not None:
        det = det.join(ann, how="left")
    det.to_csv(out / "determinants.csv")
    vc = det["class"].value_counts().to_dict()
    summary.update({f"n_{k}": v for k, v in vc.items() if k != "not_selected"})
    keep = [c for c in det.columns if c.startswith(("freq_", "mean_weight_")) or c in ("class", "is_effector", "cazyme_family")]
    R += [f"**Candidate determinants** (selection frequency >= {cfg['selection_threshold']}): {vc}", "",
          fmt(det.loc[det["class"] != "not_selected", keep].head(20))]

    # resolution
    scores = load_in_vitro(cfg)
    ts = p.meta.groupby(["isolate", "host"])["transcript_share"].mean().unstack()
    for h in ts.columns:
        scores[f"transcript_share_{h}"] = ts[h]
    for h, g in single.groupby("host"):
        scores[f"expression_score_{h}"] = g.set_index("isolate")["predicted"]
    labels = pd.DataFrame(index=sorted(set(scores.index) | set(p.isolates)))
    ap = pd.Series(cfg.get("classes") or {}, dtype=object)
    if len(ap):
        labels["apriori_high_vs_low"] = ap[ap.isin(["high", "low"])]
    if hyp is not None and hyp.nunique() == 2:
        labels["hypothesis_class"] = hyp
    for h in p.hosts:
        t = p.target[h]
        labels[f"above_median_{h}"] = pd.Series(np.where(t > t.median(), "high", "low"), index=t.index)
    res = resolution_report(scores, labels)
    res.to_csv(out / "resolution.csv", index=False)
    R += ["**Which assays resolve the class (Theorem 8.2)?** Only rows that resolve are shown; full table in "
          "resolution.csv. p_chance = probability a random score would also resolve it.", "",
          fmt(res[res["resolves"] == True].set_index("assay")[["host", "n", "p_chance"]])]  # noqa: E712
    return R, summary


def run(cfg: dict, quick: bool = False, skip_perm: bool = False) -> Path:
    t0 = time.time()
    out = Path(cfg["out_dir"])
    if not out.is_absolute():
        out = Path(cfg["_config_dir"]) / out
    out.mkdir(parents=True, exist_ok=True)
    if quick:
        cfg["n_permutations"] = min(cfg["n_permutations"], 20)
        cfg["n_bootstrap"] = min(cfg["n_bootstrap"], 12)

    pheno = load_phenotype(cfg)
    st = stability_table(pheno, cfg["stability_scale"])
    hyp = hypothesis_labels(cfg, st)
    st_out = st.join(hyp) if hyp is not None else st
    st_out.to_csv(out / "stability.csv")
    head = ["# Aggressiveness determinants: analysis report", "",
            f"Phenotype: `{Path(str(cfg['phenotype_table'])).name}`. Expression target: sAUDPC, `{cfg['target_scale']}` within host. "
            f"Stability regression on `{cfg['stability_scale']}` sAUDPC.", ""]
    S = ["## Objective 1: aggressiveness and its consistency across crops", "",
         "Host index I_h (crop mean minus grand mean): " +
         ", ".join(f"{k} {v:.1f}" for k, v in st.attrs.get("host_index", {}).items()), "",
         fmt(st_out)]
    if hyp is None:
        S += ["_No hypothesis classes given. Set `hypothesis_classes` or `stability_cutoffs` in the config "
              "(before looking at expression results)._", ""]
    les = load_lesions(cfg)
    if len(les):
        sh = interaction_share(les)
        S += ["Per-plant sAUDPC, sequential SS shares: " + ", ".join(
            f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in sh.items()), ""]

    counts, meta = load_expression(cfg)
    body, summaries = [], []
    for tps in cfg["timepoint_sets"]:
        tps = [str(t) for t in tps]
        sub = out / ("tp_" + "_".join(tps))
        sub.mkdir(exist_ok=True)
        R, summ = analyse_set(cfg, tps, counts, meta, pheno, st, hyp, sub, skip_perm)
        body += R
        summaries.append(summ)
        print(f"done {summ['timepoints']} hpi ({time.time() - t0:.0f} s)")
    summ = pd.DataFrame(summaries).set_index("timepoints")
    summ.to_csv(out / "summary_by_timepoint.csv")
    head += ["## Summary across time points", "",
             "r2_loio > 0 means held-out isolates are predicted better than by the mean of the others; "
             "judge it by p_perm. Transfer = Spearman when a model fitted on one host ranks isolates on the other.", "",
             fmt(summ), ""]
    tail = ["## Caveats", "",
            "- The isolate is the unit of evidence: n = 6, whatever the number of libraries; the smallest possible "
            "permutation p is 1/720.",
            "- Several time points and hosts are tested; p-values are not adjusted across them.",
            "- Soybean and sunflower libraries were sequenced and quantified separately: host and batch are confounded.",
            "- Pooled time points share one set of gene weights; per-time-point sets show where the signal is.",
            f"- Runtime {time.time() - t0:.0f} s."]
    (out / "report.md").write_text("\n".join(head + S + body + tail) + "\n")
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
