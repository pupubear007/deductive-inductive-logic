"""The paper's assay-resolution theorem applied to isolates (Definition 8.1, Theorem 8.2, Cor. 8.3).

Worlds are isolates (on a given host). An *assay* a maps an isolate to a result: in vitro oxalic
acid, appressorium count, colony growth, pathogen transcript share, or an expression score. The
hypothesis H is the aggressiveness class (e.g. "high" versus "low"). A thresholded assay
[a >= theta] resolves H exactly when no two isolates on the same side of theta differ on H
(Theorem 8.2); any such pair is a witness of underdetermination (Corollary 8.3).

With few isolates a perfect split happens by chance often: for 2 vs 2 isolates a random score
separates them with probability 1/3, for 3 vs 3 with probability 1/10 (reported as p_chance).

Because H is host-indexed ("aggressive on soybean"), the same assay can resolve H on one host and
fail on another. That is the formal content of "aggressiveness is conditioned by the host".
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def resolves(score: pd.Series, label: pd.Series) -> dict:
    """Best threshold split of a binary label by a continuous isolate-level score."""
    d = pd.concat([score.rename("s"), label.rename("H")], axis=1).dropna()
    if d["H"].nunique() != 2:
        return {"resolves": None, "n": int(len(d)), "note": "label needs exactly two classes"}
    pos = sorted(d["H"].unique())[1]
    y = (d["H"] == pos).to_numpy()
    s = d["s"].to_numpy(float)
    best = None
    for direction in (1, -1):
        for th in np.unique(s):
            pred = direction * s >= direction * th
            err = int(np.sum(pred != y))
            if best is None or err < best[0]:
                best = (err, th, direction)
    err, th, direction = best
    side = direction * s >= direction * th
    witnesses = [(a, b) for i, a in enumerate(d.index) for j, b in enumerate(d.index)
                 if i < j and side[i] == side[j] and y[i] != y[j]]
    from math import comb
    k = int(y.sum())
    # chance that a random ordering of the isolates separates the classes with one threshold
    p_chance = min(1.0, 2 / comb(len(y), k)) if 0 < k < len(y) else float("nan")
    return {"resolves": err == 0, "misclassified": err, "p_chance": p_chance, "threshold": float(th),
            "direction": "high score = " + str(pos) if direction == 1 else "low score = " + str(pos),
            "witness_pairs": witnesses[:10], "n": int(len(d))}


def resolution_report(scores: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    """scores: isolates x assays; labels: isolates x hosts (binary classes per host)."""
    rows = []
    for a in scores.columns:
        for h in labels.columns:
            r = resolves(scores[a], labels[h])
            rows.append({"assay": a, "host": h, **{k: v for k, v in r.items() if k != "witness_pairs"},
                         "witnesses": "; ".join(f"{x}~{y}" for x, y in r.get("witness_pairs", []))})
    return pd.DataFrame(rows)
