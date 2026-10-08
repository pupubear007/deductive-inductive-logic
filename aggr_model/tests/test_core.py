import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ssaggr import cv
from ssaggr.config import load_config
from ssaggr.io import parse_samples, sAUDPC
from ssaggr.model import Design, fit, gene_weights, lambda_max
from ssaggr.panel import build_panel, design, fit_prep
from ssaggr.resolution import resolves
from ssaggr.responsiveness import exact_spearman
from ssaggr.stability import stability_table
from ssaggr.synthetic import write


def test_parse_samples_and_rep_remap():
    cfg = load_config()
    m = parse_samples(["48_5_Ha_Xtra7", "0_1_Gm_NC", "24_2_Gm_Ss1980"], cfg)
    assert m.loc["48_5_Ha_Xtra7", "rep"] == "2"
    assert not m.loc["0_1_Gm_NC", "inoculated"]
    assert m.loc["24_2_Gm_Ss1980", "isolate"] == "1980"


def test_saudpc_formula():
    # constant lesion of 10 over 4 rating days 1..4: AUDPC = 30, D = 3, n = 4 -> 30 / (3*3/4) = 13.33
    assert sAUDPC(np.array([1, 2, 3, 4.0]), np.array([10, 10, 10, 10.0])) == pytest.approx(40 / 3)


def test_fista_recovers_sparse_signal():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 40))
    y = 2 * X[:, 3] - 1.5 * X[:, 7]
    d = Design(X, np.zeros(60, int), np.zeros(60, int), np.zeros(60), y, np.full(60, 1 / 60), 1, 1)
    b = fit(d, 0.05 * lambda_max(d), l1_ratio=1.0, max_iter=2000)
    w = gene_weights(b, d, "fixed")["w"]
    assert set(np.argsort(-np.abs(w))[:2]) == {3, 7}
    assert np.sum(w != 0) < 15


def test_stability_table_flags_reordering():
    pheno = pd.DataFrame({"A": [1, 2, 3.0], "B": [1, 2, 3.0], "C": [3, 2, 1.0]}, index=["i1", "i2", "i3"])
    st = stability_table(np.exp(pheno))
    assert st.loc["i2", "rank_range"] == 0 and st.loc["i1", "rank_range"] == 2


def test_resolution_witnesses():
    s = pd.Series({"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0})
    assert resolves(s, pd.Series({"a": "low", "b": "low", "c": "high", "d": "high"}))["resolves"]
    r = resolves(s, pd.Series({"a": "low", "b": "high", "c": "low", "d": "high"}))
    assert not r["resolves"] and r["witness_pairs"]


def test_exact_spearman_floor():
    x = pd.Series(range(6), index=list("abcdef"), dtype=float)
    r = exact_spearman(x, x)
    assert r["rho"] == pytest.approx(1.0) and r["p_perm"] == pytest.approx(2 / 720)  # +/-1 both extreme


@pytest.fixture(scope="module")
def fixed_panel(tmp_path_factory):
    d = write(tmp_path_factory.mktemp("fixed"), n_iso=10, n_genes=400, regime="fixed", seed=2, effect=3.0)
    cfg = load_config(Path(d) / "config.yaml")
    cfg["model"]["lambdas"] = [0.5, 0.25, 0.1]
    cfg["n_bootstrap"] = 10
    return build_panel(cfg), cfg, json.loads((Path(d) / "truth.json").read_text())


def test_preprocessing_uses_training_isolates_only(fixed_panel):
    p, cfg, _ = fixed_panel
    held = p.isolates[0]
    train = p.libs(p.isolates[1:], p.hosts)
    prep = fit_prep(p, train, 100)
    h = p.hosts[0]
    cols = p.meta.index[(p.meta["host"] == h) & (p.meta["isolate"] != held)]
    np.testing.assert_allclose(prep.mu[h], p.expr.loc[prep.genes, cols].mean(axis=1).to_numpy())


def test_planted_fixed_program_is_recovered_and_transfers(fixed_panel):
    p, cfg, truth = fixed_panel
    s = cv.summarize(cv.loio(p, p.hosts, "fixed", cfg))
    assert (s["r2_loio"] > 0.5).all()
    x = cv.cross_host(p, cfg)
    assert (x["spearman"] > 0.7).all()
    sel = cv.stability_selection(p, p.hosts[0], cfg).sort_values(f"freq_{p.hosts[0]}", ascending=False)
    top = set(sel.index[:20])
    assert len(top & set(truth["determinants"][p.hosts[0]])) >= 10


def test_run_end_to_end(tmp_path):
    from ssaggr.run import run
    d = write(tmp_path / "d", n_iso=6, n_genes=300, regime="responsive", seed=1)
    cfg = load_config(Path(d) / "config.yaml")
    cfg["model"]["lambdas"] = [0.5, 0.2]
    cfg["n_permutations"] = 3
    cfg["n_bootstrap"] = 6
    out = run(cfg)
    rep = (out / "report.md").read_text()
    assert "Held-out isolate prediction" in rep and (out / "determinants.csv").exists()
