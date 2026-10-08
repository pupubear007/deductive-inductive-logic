"""Figures from a finished run (PNG at 300 dpi + PDF), written to ``<run>/figures/``.

    python -m ssaggr.plots --run <out_dir of ssaggr.run>

Fig. 1  stability across crops (Objective 1): mean rank vs rank range, and Eberhart-Russell
        slope vs deviation; the six RNA-seq isolates highlighted.
Fig. 2  held-out (leave-one-isolate-out) R2 per time point and host, with permutation p.
Fig. 3  observed vs held-out predicted aggressiveness, one panel per host x time point.
Fig. 4  host-responsiveness index per isolate and time point, and against rank range.
Fig. 5  candidate determinants: selection frequency on soybean vs sunflower (pooled time points).
Fig. 6  which assays resolve the aggressiveness class (Theorem 8.2), with the chance rate.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Palette (validated reference palette, light mode): two hosts = categorical slots 1 and 2.
HOST_COLOR = {"Gm": "#2a78d6", "Ha": "#eb6834"}
HOST_NAME = {"Gm": "Soybean (Williams 82)", "Ha": "Sunflower (HA89)"}
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#a3a29c", "#e6e5e1", "#fcfcfb"
HILITE = "#2a78d6"


def _style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "figure.facecolor": "white", "axes.facecolor": "white", "legend.frameon": False,
        "axes.titleweight": "bold", "axes.titlelocation": "left", "text.color": INK,
    })


def _label_points(ax, xs, ys, texts, colors=None, fontsize=7):
    """Place point labels without overlapping each other or the points: for each label try a ring
    of offsets and keep the first whose box hits nothing already placed."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    frame = ax.get_window_extent(r)
    placed = []
    offsets = [(4, 3), (4, -9), (-4, 3), (-4, -9), (6, -3), (-6, -3), (0, 7), (0, -13),
               (10, 10), (-10, 10), (10, -16), (-10, -16), (14, 0), (-14, 0)]
    for k, (x, y, t) in enumerate(zip(xs, ys, texts)):
        col = colors[k] if colors is not None else INK2
        best = None
        for dx, dy in offsets:
            ha = "left" if dx >= 0 else "right"
            a = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=fontsize,
                            color=col, ha=ha)
            bb = a.get_window_extent(r).expanded(1.05, 1.1)
            inside = frame.x0 <= bb.x0 and bb.x1 <= frame.x1 and frame.y0 <= bb.y0 and bb.y1 <= frame.y1
            hit = (not inside) or any(bb.overlaps(o) for o in placed) or any(
                bb.contains(px, py) for j, (px, py) in enumerate(pts) if j != k)
            if not hit:
                best = a
                placed.append(bb)
                break
            a.remove()
        if best is None:  # fall back to the default offset
            a = ax.annotate(t, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=fontsize, color=col)
            placed.append(a.get_window_extent(r))


def _save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def _tp_dirs(run: Path) -> list[tuple[str, Path]]:
    dirs = sorted([d for d in run.glob("tp_*") if d.is_dir()], key=lambda d: (len(d.name), d.name))
    return [(d.name[3:].replace("_", "+"), d) for d in dirs]


def _hosts(df: pd.DataFrame) -> list[str]:
    return [h for h in ("Gm", "Ha") if h in set(df["host"])] + sorted(set(df["host"]) - {"Gm", "Ha"})


# ------------------------------------------------------------------------------------------------

def fig_stability(run: Path, out: Path, rnaseq: set[str]):
    st = pd.read_csv(run / "stability.csv", index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    for ax, (x, y, xl, yl) in zip(axes, [("mean_rank", "rank_range", "Mean rank across 4 crops (1 = least aggressive)",
                                          "Rank range across crops (0 = same rank everywhere)"),
                                         ("slope_b", "deviation_ms", "Eberhart–Russell slope b (tracks host susceptibility)",
                                          "Deviation mean square (host-specific departures)")]):
        on = st.index.isin(rnaseq)
        ax.scatter(st.loc[~on, x], st.loc[~on, y], s=36, color=MUTED, edgecolor="white", linewidth=1.5, zorder=2,
                   label="phenotyped only")
        ax.scatter(st.loc[on, x], st.loc[on, y], s=48, color=HILITE, edgecolor="white", linewidth=1.5, zorder=3,
                   label="in the RNA-seq panel")
        if x == "slope_b":
            ax.axvline(1.0, color=MUTED, lw=1, ls="--", zorder=1)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
    axes[0].set_title("A  Consistency of aggressiveness across crops")
    axes[1].set_title("B  Joint regression on the host index")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    for ax, (x, y) in zip(axes, [("mean_rank", "rank_range"), ("slope_b", "deviation_ms")]):
        _label_points(ax, st[x].to_numpy(), st[y].to_numpy(), list(st.index),
                      [INK if i in rnaseq else INK2 for i in st.index])
    _save(fig, out, "fig1_stability")


def fig_r2(run: Path, out: Path):
    s = pd.read_csv(run / "summary_by_timepoint.csv")
    s["timepoints"] = s["timepoints"].astype(str)
    hosts = [h for h in ("Gm", "Ha") if f"r2_loio_{h}" in s.columns]
    labels = ["pooled\n24+48+96" if "+" in t else f"{t} hpi" for t in s["timepoints"]]
    x = np.arange(len(s))
    w = 0.36
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for k, h in enumerate(hosts):
        v = s[f"r2_loio_{h}"].to_numpy(float)
        xs = x + (k - (len(hosts) - 1) / 2) * (w + 0.04)
        ax.bar(xs, v, width=w, color=HOST_COLOR.get(h, MUTED), label=HOST_NAME.get(h, h), zorder=2)
        pcol = f"p_perm_{h}"
        for xi, vi, pi in zip(xs, v, s[pcol] if pcol in s else [np.nan] * len(s)):
            if np.isfinite(pi):
                txt = f"p={pi:.2f}" if pi >= 0.01 else f"p={pi:.3f}"
                ax.annotate(txt, (xi, vi), xytext=(0, 3 if vi >= 0 else -11), textcoords="offset points",
                            ha="center", fontsize=7, color=INK if pi < 0.05 else INK2,
                            fontweight="bold" if pi < 0.05 else "normal")
    ax.axhline(0, color=INK2, lw=0.8, zorder=3)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.06 * (hi - lo))
    ax.set_xticks(x, labels)
    ax.set_ylabel("Held-out R² (leave one isolate out)")
    ax.set_title("Predicting a held-out isolate's aggressiveness from its pathogen expression")
    ax.legend(loc="lower right", fontsize=8)
    ax.text(0, -0.24, "R² > 0: better than predicting the mean of the other isolates.\np: isolate-label permutation "
            "null (whole procedure rerun); not adjusted across the 8 tests.", transform=ax.transAxes, fontsize=7,
            color=INK2, va="top")
    fig.tight_layout()
    _save(fig, out, "fig2_heldout_r2")


def fig_obs_pred(run: Path, out: Path):
    tps = _tp_dirs(run)
    data = [(t, pd.read_csv(d / "loio_single.csv")) for t, d in tps if (d / "loio_single.csv").exists()]
    if not data:
        return
    hosts = _hosts(data[0][1])
    fig, axes = plt.subplots(len(hosts), len(data), figsize=(2.55 * len(data), 2.6 * len(hosts)), squeeze=False)
    lim = max(max(df["observed"].abs().max(), df["predicted"].abs().max()) for _, df in data) * 1.15
    for j, (t, df) in enumerate(data):
        for i, h in enumerate(hosts):
            ax = axes[i, j]
            g = df[df["host"] == h]
            ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=1, ls="--", zorder=1)
            ax.scatter(g["observed"], g["predicted"], s=40, color=HOST_COLOR.get(h, MUTED), edgecolor="white",
                       linewidth=1.5, zorder=3)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            name = "pooled" if "+" in t else f"{t} hpi"
            ax.set_title(f"{HOST_NAME.get(h, h).split(' ')[0]}, {name}", fontsize=9)
            if i == len(hosts) - 1:
                ax.set_xlabel("Observed (z sAUDPC)")
            if j == 0:
                ax.set_ylabel("Held-out prediction")
    fig.suptitle("Observed vs held-out predicted aggressiveness (dashed: perfect prediction)", x=0.01, ha="left",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    for j, (t, df) in enumerate(data):
        for i, h in enumerate(hosts):
            g = df[df["host"] == h]
            _label_points(axes[i, j], g["observed"].to_numpy(), g["predicted"].to_numpy(), list(g["isolate"]),
                          fontsize=6.5)
    _save(fig, out, "fig3_observed_vs_predicted")


def fig_responsiveness(run: Path, out: Path):
    tabs = []
    for t, d in _tp_dirs(run):
        if "+" in t or not (d / "responsiveness.csv").exists():
            continue
        r = pd.read_csv(d / "responsiveness.csv")
        tabs.append(r)
    if not tabs:
        return
    r = pd.concat(tabs)
    r["hpi"] = r["hpi"].astype(str)
    piv = r.pivot(index="isolate", columns="hpi", values="R2_host")
    piv = piv[sorted(piv.columns, key=int)]
    st = pd.read_csv(run / "stability.csv", index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    order = piv.mean(axis=1).sort_values().index
    seq = ["#9ec5f4", "#3987e5", "#184f95"]  # one hue, light -> dark = early -> late
    for k, hpi in enumerate(piv.columns):
        ax.scatter(piv.loc[order, hpi], np.arange(len(order)), s=44, color=seq[k % 3], edgecolor="white",
                   linewidth=1.5, zorder=3, label=f"{hpi} hpi")
    for y, iso in enumerate(order):
        ax.plot([piv.loc[iso].min(), piv.loc[iso].max()], [y, y], color=GRID, lw=2, zorder=1)
    ax.set_yticks(np.arange(len(order)), order)
    ax.set_xlabel("R²_host: share of the isolate's pathogen-expression variance explained by host")
    ax.set_title("A  Host-responsiveness index", pad=18)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=8, handletextpad=0.2, borderaxespad=0.1)
    ax = axes[1]
    m = piv.mean(axis=1)
    rr = st.loc[m.index, "rank_range"]
    ax.scatter(rr, m, s=48, color=HILITE, edgecolor="white", linewidth=1.5, zorder=3)
    ax.set_xlabel("Rank range across 4 crops (low = consistent)")
    ax.set_ylabel("Mean R²_host over 24, 48, 96 hpi")
    ax.set_title("B  Responsiveness vs consistency", pad=18)
    ax.text(0, -0.24, "Host and sequencing batch are confounded:\ncompare isolates, not absolute values. n = 6.",
            transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    fig.tight_layout()
    _label_points(ax, rr.to_numpy(), m.to_numpy(), list(m.index))
    _save(fig, out, "fig4_host_responsiveness")


def fig_determinants(run: Path, out: Path, thr: float = 0.6):
    tps = _tp_dirs(run)
    pooled = [d for t, d in tps if "+" in t] or [tps[-1][1]]
    f = pooled[0] / "determinants.csv"
    if not f.exists():
        return
    d = pd.read_csv(f, index_col=0)
    if not {"freq_Gm", "freq_Ha"} <= set(d.columns):
        return
    eff = d.get("is_effector", pd.Series(False, index=d.index)).astype(str).str.lower().eq("true")
    caz = d.get("cazyme", pd.Series(False, index=d.index)).astype(str).str.lower().eq("true")
    sel = d[(d["freq_Gm"] >= thr) | (d["freq_Ha"] >= thr)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    jit = np.random.default_rng(0).uniform(-0.008, 0.008, (len(d), 2))
    x, y = d["freq_Gm"] + jit[:, 0], d["freq_Ha"] + jit[:, 1]
    plain = ~(eff | caz)
    ax.scatter(x[plain], y[plain], s=14, color=MUTED, alpha=0.6, edgecolor="none", zorder=2, label="other genes")
    ax.scatter(x[caz & ~eff], y[caz & ~eff], s=34, marker="s", color="#1baf7a", edgecolor="white", linewidth=1,
               zorder=3, label="CAZyme")
    ax.scatter(x[eff], y[eff], s=40, marker="D", color="#4a3aa7", edgecolor="white", linewidth=1, zorder=4,
               label="predicted effector")
    ax.scatter(sel["freq_Gm"], sel["freq_Ha"], s=70, facecolor="none", edgecolor=INK, linewidth=1, zorder=5,
               label=f"selected (≥ {thr})")
    ax.axvline(thr, color=HOST_COLOR["Gm"], lw=1, ls="--", zorder=1)
    ax.axhline(thr, color=HOST_COLOR["Ha"], lw=1, ls="--", zorder=1)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Selection frequency, soybean model")
    ax.set_ylabel("Selection frequency, sunflower model")
    ax.set_title("A  All genes (pooled 24+48+96 hpi)")
    ax.text(1.0, 0.30, "soybean-specific", color=HOST_COLOR["Gm"], fontsize=8, ha="right")
    ax.text(0.03, 0.97, "sunflower-specific", color=HOST_COLOR["Ha"], fontsize=8, va="top")
    ax.text(1.0, 0.97, "shared", color=INK2, fontsize=8, va="top", ha="right")
    ax.legend(loc="center right", fontsize=8)
    ax.text(0, -0.13, "Frequency = share of resamples (one isolate left out × replicate bootstrap)\nin which the "
            "gene was selected.", transform=ax.transAxes, fontsize=7, color=INK2, va="top")

    ax = axes[1]
    sel["host"] = np.where(sel["freq_Gm"] >= sel["freq_Ha"], "Gm", "Ha")
    sel["freq"] = sel[["freq_Gm", "freq_Ha"]].max(axis=1)
    sel["weight"] = np.where(sel["host"] == "Gm", sel.get("mean_weight_Gm", 0), sel.get("mean_weight_Ha", 0))
    sel = sel.sort_values(["host", "freq"], ascending=[True, True])
    yy = np.arange(len(sel))
    for k, (g, r) in enumerate(sel.iterrows()):
        col = HOST_COLOR[r["host"]]
        ax.plot([thr, r["freq"]], [k, k], color=col, lw=2, zorder=2, solid_capstyle="round")
        mk = "D" if eff.get(g, False) else ("s" if caz.get(g, False) else "o")
        ax.scatter(r["freq"], k, s=46, marker=mk, color=col, edgecolor="white", linewidth=1.2, zorder=3)
        ax.annotate("↑ higher in aggressive" if r["weight"] > 0 else "↓ lower in aggressive", (r["freq"], k),
                    xytext=(7, -3), textcoords="offset points", fontsize=6.5, color=INK2)
    ax.set_yticks(yy, list(sel.index), fontsize=7.5)
    ax.set_xlim(thr - 0.02, 1.25)
    ax.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Selection frequency on its host")
    ax.set_title("B  Selected genes (blue: soybean, orange: sunflower)")
    ax.text(0, -0.13, "Arrow: sign of the gene's weight (expression higher or lower in the more aggressive "
            "isolates).\nMarker: ◆ effector, ■ CAZyme, ● other.", transform=ax.transAxes, fontsize=7, color=INK2,
            va="top")
    fig.tight_layout()
    _save(fig, out, "fig5_determinants")


def fig_resolution(run: Path, out: Path):
    tps = _tp_dirs(run)
    rows = []
    for t, d in tps:
        f = d / "resolution.csv"
        if f.exists():
            r = pd.read_csv(f)
            r["tp"] = "pooled" if "+" in t else f"{t} hpi"
            rows.append(r)
    if not rows:
        return
    r = pd.concat(rows)
    r = r[~r["assay"].str.startswith("expression_score")]  # expression scores are already in Fig. 2
    r["label"] = (r["host"].str.replace("above_median_Gm", "above median on soybean")
                  .str.replace("above_median_Ha", "above median on sunflower")
                  .str.replace("apriori_high_vs_low", "a-priori high vs low"))
    r["assay"] = (r["assay"].str.replace("transcript_share_Gm", "transcript share (soybean)")
                  .str.replace("transcript_share_Ha", "transcript share (sunflower)")
                  .str.replace("oxalic_acid_per_mg", "oxalic acid per mg")
                  .str.replace("appressoria_per_field_96h", "appressoria per field, 96 h")
                  .str.replace("radial_growth_sAUGC", "radial growth (sAUGC)"))
    r["row"] = r["assay"] + "  |  " + r["label"]
    in_vitro = ~r["assay"].str.startswith("transcript share")
    # in vitro traits do not depend on time point: keep one copy
    r = pd.concat([r[in_vitro].drop_duplicates("row").assign(tp="any"), r[~in_vitro]])
    cols = ["any"] + [c for c in ["24 hpi", "48 hpi", "96 hpi", "pooled"] if c in set(r["tp"])]
    rws = list(dict.fromkeys(r["row"]))
    fig, ax = plt.subplots(figsize=(9.6, 0.34 * len(rws) + 1.2))
    for i, row in enumerate(rws):
        for j, c in enumerate(cols):
            g = r[(r["row"] == row) & (r["tp"] == c)]
            if g.empty:
                continue
            g = g.iloc[0]
            ok = bool(g["resolves"]) if str(g["resolves"]) != "nan" else False
            ax.scatter(j, i, s=110, marker="o", color=HILITE if ok else "white",
                       edgecolor=HILITE if ok else MUTED, linewidth=1.4, zorder=3)
            ax.annotate(f"{int(g['misclassified'])} wrong" if not ok else f"chance {g['p_chance']:.2f}",
                        (j, i), xytext=(9, -3), textcoords="offset points", fontsize=6.5,
                        color=INK if ok else INK2)
    ax.set_xticks(range(len(cols)), cols)
    ax.set_yticks(range(len(rws)), rws, fontsize=7.5)
    ax.set_xlim(-0.4, len(cols) - 0.35)
    ax.invert_yaxis()
    ax.grid(False)
    ax.set_title("Which assays resolve the aggressiveness class? (filled = one threshold separates all isolates)")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0)
    ax.spines["bottom"].set_visible(False)
    ax.text(0, -0.03, "chance: probability that a random score would separate the classes just as well "
            "(2 vs 2 isolates: 0.33; 3 vs 3: 0.10). 'any': in vitro traits do not depend on time point.",
            transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    fig.tight_layout()
    _save(fig, out, "fig6_assay_resolution")


def make_all(run: str | Path, rnaseq: list[str] | None = None) -> Path:
    _style()
    run = Path(run)
    out = run / "figures"
    out.mkdir(exist_ok=True)
    if rnaseq is None:
        tps = _tp_dirs(run)
        rnaseq = sorted(pd.read_csv(tps[0][1] / "loio_single.csv")["isolate"].unique()) if tps else []
    fig_stability(run, out, set(rnaseq))
    fig_r2(run, out)
    fig_obs_pred(run, out)
    fig_responsiveness(run, out)
    fig_determinants(run, out)
    fig_resolution(run, out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True, help="out_dir of a finished ssaggr.run")
    a = ap.parse_args(argv)
    out = make_all(a.run)
    print(f"figures -> {out}")


if __name__ == "__main__":
    main()
