"""Isolate-wise validation: the only honest test with a handful of isolates.

* Leave-one-isolate-out (LOIO): all libraries of the held-out isolate (on every host) are removed,
  the L1 weight is chosen by an inner LOIO over the remaining isolates, and the held-out isolate's
  aggressiveness is predicted from its own expression. With libraries as the unit, a flexible
  model would simply memorise isolate identity from replicates.
* Permutation null: isolate labels are permuted across the phenotype vector and the whole
  procedure (including the inner lambda choice) is rerun. With 6 isolates there are only 720
  orderings, so no p-value can be smaller than 1/720.
* Cross-host transfer: fit on one host, predict the isolates' ranking on the other. Under the
  fixed-program hypothesis this should work; under host-responsive regulation it should fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import fit, gene_weights, lambda_max, predict
from .panel import Panel, design, fit_prep


def _isolate_means(p: Panel, libs: pd.Index, pred: np.ndarray) -> pd.DataFrame:
    m = p.meta.loc[libs, ["isolate", "host"]].copy()
    m["pred"] = pred
    return m.groupby(["isolate", "host"])["pred"].mean().reset_index()


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def select_and_fit(p: Panel, train_iso: list[str], hosts: list[str], kind: str, cfg: dict,
                   target: pd.DataFrame | None = None):
    """Inner LOIO over ``train_iso`` to pick lambda, then refit on all of them."""
    mc = cfg["model"]
    target = p.target if target is None else target
    train = p.libs(train_iso, hosts)
    prep = fit_prep(p, train, cfg["top_variable_genes"])
    d_all = design(p, train, prep, cfg, target)
    lmax = lambda_max(d_all)
    lams = [r * lmax for r in mc["lambdas"]]
    best_lam, best_err = lams[0], np.inf
    if len(train_iso) >= 3:
        errs = np.zeros(len(lams))
        for iso in train_iso:
            inner = [i for i in train_iso if i != iso]
            tr, te = p.libs(inner, hosts), p.libs([iso], hosts)
            d_tr, d_te = design(p, tr, prep, cfg, target), design(p, te, prep, cfg, target)
            beta = None
            for k, lam in enumerate(lams):  # warm start down the path
                beta = fit(d_tr, lam, kind, mc["host_penalty_ratio"], mc["ridge"], mc["max_iter"], mc["tol"], beta, mc["l1_ratio"])
                im = _isolate_means(p, te, predict(d_te, beta, kind))
                obs = np.array([target.loc[r.isolate, r.host] for r in im.itertuples()])
                errs[k] += float(np.mean((im["pred"].to_numpy() - obs) ** 2))
        best_lam = lams[int(np.argmin(errs))]
    beta = fit(d_all, best_lam, kind, mc["host_penalty_ratio"], mc["ridge"], mc["max_iter"], mc["tol"], None, mc["l1_ratio"])
    return beta, prep, best_lam / lmax


def loio(p: Panel, hosts: list[str], kind: str, cfg: dict, target: pd.DataFrame | None = None) -> pd.DataFrame:
    """Held-out isolate predictions. Returns rows (isolate, host, observed, predicted, lambda)."""
    target = p.target if target is None else target
    rows = []
    for iso in p.isolates:
        train_iso = [i for i in p.isolates if i != iso]
        beta, prep, lam_rel = select_and_fit(p, train_iso, hosts, kind, cfg, target)
        te = p.libs([iso], hosts)
        im = _isolate_means(p, te, predict(design(p, te, prep, cfg, target), beta, kind))
        for r in im.itertuples():
            rows.append({"isolate": r.isolate, "host": r.host, "observed": float(target.loc[r.isolate, r.host]),
                         "predicted": float(r.pred), "lambda_rel": lam_rel})
    return pd.DataFrame(rows)


def summarize(pred: pd.DataFrame) -> pd.DataFrame:
    """Per host: Spearman and MSE of held-out predictions, and the MSE of the honest baseline that
    predicts each held-out isolate by the mean of the others. Under no signal, leave-one-out with an
    intercept gives *negative* correlations (the training mean moves away from the held-out value),
    so judge the model by R2_loio = 1 - mse / mse_baseline and by the permutation null, not by the
    sign of Spearman alone."""
    out = []
    for h, g in pred.groupby("host"):
        y = g["observed"].to_numpy()
        n = len(y)
        base = (y.sum() - y) / (n - 1) if n > 1 else np.zeros_like(y)
        mse = float(np.mean((y - g["predicted"].to_numpy()) ** 2))
        mse_b = float(np.mean((y - base) ** 2))
        out.append({"host": h, "spearman": _spearman(y, g["predicted"].to_numpy()), "mse": mse,
                    "mse_baseline": mse_b, "r2_loio": 1 - mse / mse_b if mse_b > 0 else np.nan, "n_isolates": int(n)})
    return pd.DataFrame(out).set_index("host")


_G: dict = {}


def _init(p, hosts, kind, cfg):
    _G.update(p=p, hosts=hosts, kind=kind, cfg=cfg)


def _perm_r2(perm: np.ndarray) -> dict:
    p, cfg = _G["p"], _G["cfg"]
    t = p.target.copy()
    t.index = [p.isolates[j] for j in perm]
    t = t.loc[p.isolates]
    return summarize(loio(p, _G["hosts"], _G["kind"], cfg, t))["r2_loio"].to_dict()


def permutation_test(p: Panel, hosts: list[str], kind: str, cfg: dict, observed: pd.DataFrame,
                     n_perm: int, seed: int = 0) -> pd.DataFrame:
    """Null distribution of the LOIO R2 per host under random isolate-to-phenotype mapping.
    For 6 isolates and n_perm >= 719 every non-identity ordering is used exactly once."""
    from itertools import permutations
    from math import factorial
    rng = np.random.default_rng(seed)
    n = len(p.isolates)
    if n_perm >= factorial(n) - 1:
        perms = [np.array(q) for q in permutations(range(n)) if list(q) != list(range(n))]
    else:
        perms = [rng.permutation(n) for _ in range(n_perm)]
    obs = summarize(observed)["r2_loio"]
    jobs = int(cfg.get("n_jobs") or 1)
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(jobs, initializer=_init, initargs=(p, hosts, kind, cfg)) as pool:
            res = pool.map(_perm_r2, perms, chunksize=max(1, len(perms) // (4 * jobs)))
    else:
        _init(p, hosts, kind, cfg)
        res = [_perm_r2(q) for q in perms]
    rows = []
    for h in obs.index:
        nv = np.array([r.get(h, np.nan) for r in res], float)
        nv = nv[np.isfinite(nv)]
        rows.append({"host": h, "r2_loio": obs[h], "p_perm": (1 + np.sum(nv >= obs[h] - 1e-12)) / (1 + len(nv)),
                     "n_perm": int(len(nv)), "null_mean": float(nv.mean()) if len(nv) else np.nan})
    return pd.DataFrame(rows).set_index("host")


def cross_host(p: Panel, cfg: dict) -> pd.DataFrame:
    """Train a fixed-program model on one host, predict isolate ranking on each other host."""
    rows = []
    for a in p.hosts:
        beta, prep, lam = select_and_fit(p, p.isolates, [a], "fixed", cfg)
        for b in p.hosts:
            if b == a:
                continue
            te = p.libs(p.isolates, [b])
            # the model's intercept groups belong to host a; ranking within host b is unaffected
            d = design(p, te, prep, cfg)
            d.group = np.zeros_like(d.group)
            d.cov = np.zeros_like(d.cov)
            im = _isolate_means(p, te, predict(d, beta, "fixed"))
            obs = np.array([p.target.loc[i, b] for i in im["isolate"]])
            rows.append({"train_host": a, "test_host": b, "spearman": _spearman(obs, im["pred"].to_numpy()),
                         "n_isolates": int(len(im))})
    return pd.DataFrame(rows)


def stability_selection(p: Panel, host: str, cfg: dict, seed: int = 0, lam_factors=(1.0, 0.5, 0.25)) -> pd.DataFrame:
    """Stability selection (Meinshausen and Buehlmann 2010), adapted to few isolates: subsamples are
    the LOIO training sets (all isolates but one) x bootstrap resamples of replicate libraries within
    isolate; each is fitted along a short lambda path below the LOIO-chosen lambda, and a gene's
    score is its maximum selection frequency over that path."""
    rng = np.random.default_rng(seed)
    mc = cfg["model"]
    beta, prep, lam_rel = select_and_fit(p, p.isolates, [host], "fixed", cfg)
    counts = {f: {} for f in lam_factors}
    wsum, nfit = {}, 0
    n_boot = max(1, cfg["n_bootstrap"] // max(1, len(p.isolates)))
    for iso in p.isolates:
        train_iso = [i for i in p.isolates if i != iso]
        base = p.libs(train_iso, [host])
        prep_f = fit_prep(p, base, cfg["top_variable_genes"])
        for _ in range(n_boot):
            m = p.meta.loc[base]
            pick = []
            for _, g in m.groupby("isolate"):
                pick.extend(rng.choice(g.index.to_numpy(), size=len(g), replace=True))
            d = design(p, pd.Index(pick), prep_f, cfg)
            lmax = lambda_max(d)
            b = None
            for f in lam_factors:
                b = fit(d, lam_rel * f * lmax, "fixed", mc["host_penalty_ratio"], mc["ridge"], mc["max_iter"], mc["tol"], b, mc["l1_ratio"])
                w = gene_weights(b, d, "fixed")["w"]
                for j in np.flatnonzero(w):
                    g_ = prep_f.genes[j]
                    counts[f][g_] = counts[f].get(g_, 0) + 1
                    wsum[g_] = wsum.get(g_, 0.0) + w[j]
            nfit += 1
    full = gene_weights(beta, design(p, p.libs(p.isolates, [host]), prep, cfg), "fixed")["w"]
    w_full = dict(zip(prep.genes, full))
    genes = set().union(*[set(c) for c in counts.values()])
    rows = []
    for g in genes:
        fr = max(counts[f].get(g, 0) for f in lam_factors) / nfit
        tot = sum(counts[f].get(g, 0) for f in lam_factors)
        rows.append({"gene": g, f"freq_{host}": fr, f"mean_weight_{host}": wsum[g] / tot,
                     f"full_weight_{host}": w_full.get(g, 0.0)})
    cols = [f"freq_{host}", f"mean_weight_{host}", f"full_weight_{host}"]
    return pd.DataFrame(rows).set_index("gene") if rows else pd.DataFrame(columns=cols)


def concordance(tables: list[pd.DataFrame], hosts: list[str], thr: float) -> pd.DataFrame:
    """Classify selected genes as shared-concordant, shared-opposite or host-specific (CH4 Fig. 10)."""
    df = pd.concat(tables, axis=1).fillna(0.0)
    sel = {h: df[f"freq_{h}"] >= thr for h in hosts}
    cls = []
    for g in df.index:
        on = [h for h in hosts if sel[h][g]]
        if len(on) == 0:
            cls.append("not_selected")
        elif len(on) == 1:
            cls.append(f"specific_{on[0]}")
        else:
            signs = {np.sign(df.loc[g, f"mean_weight_{h}"]) for h in on}
            cls.append("shared_concordant" if len(signs) == 1 else "shared_opposite")
    df["class"] = cls
    return df.sort_values([f"freq_{h}" for h in hosts], ascending=False)
