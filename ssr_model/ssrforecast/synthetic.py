"""Mechanism-consistent synthetic data, in exactly the format real data must use.

    python -m ssrforecast.synthetic --out data/synthetic [--fields 240] [--chamber-units 48]

EVERY NUMBER IN THIS FILE IS AN ARBITRARY SIMULATION SETTING, NOT A BIOLOGICAL ESTIMATE.
The generator exists to (1) exercise the whole pipeline, (2) give a case where the true
mechanism is known, so one can check the model recovers it, and (3) show the file layout.

Worlds are generated along the mechanism (see mechanism.py):
    sclerotia -> apothecia (needs soil moisture) ; external spores ; spore = apothecia OR external
    court = spore AND flowering AND q_court ; infection = court AND microclimate AND q_inf
    incidence > 0 iff infection ; DSI, lesion length = 0 without infection
Growth-chamber units are inoculated by cut petiole, i.e. do(court = 1).
Hard constraints therefore hold in every generated world.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

SYNTHETIC_SE_SP = {  # used both to generate readings and in configs/synthetic.yaml
    "apothecia_scouting": (0.7, 1.0),
    "spore_trap_qpcr": (0.85, 0.97),
    "petal_qpcr": (0.8, 0.98),
    "stem_qpcr": (0.95, 0.99),
}


def _sig(x):
    return 1 / (1 + np.exp(-x))


def _blob_map(rng, incidence: float, size: int) -> np.ndarray:
    """Disease map with clustered foci whose mean is roughly the incidence."""
    if incidence <= 0:
        return np.zeros((size, size), np.float32)
    yy, xx = np.mgrid[0:size, 0:size]
    m = np.zeros((size, size))
    for _ in range(1 + int(rng.poisson(1 + 6 * incidence))):
        cy, cx = rng.uniform(0, size, 2)
        s = rng.uniform(1.5, size / 5)
        m += np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * s * s))
    m = m / (m.mean() + 1e-9) * incidence
    return np.clip(m, 0, 1).astype(np.float32)


def _canopy_image(rng, dmap: np.ndarray, size: int, channels: int = 3) -> np.ndarray:
    """Green canopy; diseased pixels pale / tan (RGB), or a lower NIR band (multispectral)."""
    d = np.kron(dmap, np.ones((size // dmap.shape[0], size // dmap.shape[1])))[:size, :size]
    base = rng.normal(0, 0.04, (size, size))
    if channels == 3:
        r = 0.25 + 0.5 * d + base
        g = 0.55 - 0.15 * d + base
        b = 0.2 + 0.45 * d + base
        img = np.stack([r, g, b], -1)
    else:
        bands = [0.08 + 0.1 * d, 0.12 + 0.05 * d, 0.1 + 0.12 * d, 0.6 - 0.35 * d, 0.3 - 0.1 * d][:channels]
        img = np.stack([x + base for x in bands], 0)
    return np.clip(img, 0, 1).astype(np.float32)


def _stem_image(rng, lesion_frac: float, size: int = 64) -> np.ndarray:
    """A growth-chamber plant: vertical green stem with a white/tan lesion of the given length."""
    img = np.full((size, size, 3), 0.85) + rng.normal(0, 0.02, (size, size, 3))
    x0 = size // 2 - 2
    img[4:size - 4, x0:x0 + 4] = [0.2, 0.5, 0.2]
    L = int(lesion_frac * (size - 8))
    if L > 0:
        top = size // 3
        img[top:min(top + L, size - 4), x0:x0 + 4] = [0.9, 0.85, 0.7]
    return np.clip(img, 0, 1).astype(np.float32)


def _save_png(arr: np.ndarray, path: Path) -> None:
    Image.fromarray((arr * 255).astype(np.uint8)).save(path)


def generate(out: str | Path, n_fields: int = 240, n_chamber: int = 48, seed: int = 0, map_size: int = 32) -> Path:
    rng = np.random.default_rng(seed)
    out = Path(out)
    (out / "images").mkdir(parents=True, exist_ok=True)
    (out / "maps").mkdir(exist_ok=True)
    units, obs, ts, ims = [], [], [], []
    cultivars = {"CV_A": 0.0, "CV_B": 0.4, "CV_C": 0.8, "CV_D": 1.2}  # synthetic resistance shift
    sites = {f"site{i}": (rng.uniform(42, 48), rng.uniform(-97, -91), rng.uniform(-2, 2)) for i in range(6)}
    seasons = [2019, 2020, 2021, 2022, 2023]

    # ---------------------------------------------------------------- field-seasons
    for k in range(n_fields):
        site = list(sites)[k % len(sites)]
        season = seasons[(k // len(sites)) % len(seasons)]
        lat, lon, toff = sites[site]
        uid = f"F{k:04d}"
        plant = pd.Timestamp(f"{season}-05-10") + pd.Timedelta(days=int(rng.integers(0, 15)))
        cv = rng.choice(list(cultivars))
        row = float(rng.choice([38.0, 76.0]))
        fung = rng.choice(["none", "R1", "R3"], p=[0.5, 0.25, 0.25])
        history = rng.random() < 0.6
        season_wet = rng.normal(0, 1) + 0.3 * (season % 3 - 1)  # shared within a season

        n_days = 130
        d = np.arange(n_days)
        temp = 17 + toff + 7 * np.sin(np.pi * d / n_days) + rng.normal(0, 2, n_days)
        rain = (rng.random(n_days) < _sig(-1.4 + 0.4 * season_wet)) * rng.gamma(2, 6, n_days)
        rh = np.clip(65 + 8 * season_wet + 0.8 * rain + rng.normal(0, 6, n_days), 30, 100)
        soil = np.zeros(n_days)
        for t in range(n_days):
            soil[t] = np.clip((soil[t - 1] if t else 25) * 0.93 + 0.9 * rain[t] + 1.2, 5, 60)
        closure = 1.0 if row == 38.0 else 0.6  # narrow rows close the canopy earlier
        lw = np.clip((rh - 70) / 3 + 0.3 * rain + 3 * closure * (d > 55) + rng.normal(0, 1, n_days), 0, 24)
        for t in range(n_days):
            day = (plant + pd.Timedelta(days=int(t))).strftime("%Y-%m-%d")
            for var, val in (("air_temp", temp[t]), ("rh", rh[t]), ("rain", rain[t]),
                             ("leaf_wetness", lw[t]), ("soil_moisture", soil[t])):
                ts.append((uid, day, var, round(float(val), 3)))

        r1 = int(rng.integers(48, 56))
        win = slice(r1, r1 + 21)  # R1-R3
        S = history or rng.random() < 0.15
        A = S and rng.random() < _sig((soil[r1 - 14:r1 + 7].mean() - 30) / 4)
        X = rng.random() < 0.12
        spore = A or X
        court = spore and rng.random() < 0.85
        micro = lw[win].mean() > 2.0
        q_inf = _sig(2.0 - cultivars[cv] - {"none": 0, "R1": 1.0, "R3": 0.6}[fung])
        inf = court and micro and rng.random() < q_inf
        inc = float(np.clip(rng.beta(2, 5) * (1.4 - 0.4 * cultivars[cv]) * (0.6 + 0.08 * lw[win].mean()), 0.01, 1)) if inf else 0.0
        dsi = float(np.clip(inc * 100 * rng.uniform(0.4, 0.9), 0, 100)) if inf else 0.0
        yld = float(4200 + 300 * rng.normal() - 2500 * inc)

        units.append({"unit_id": uid, "experiment_type": "field", "season": season, "site": site,
                      "cultivar": cv, "isolate": "", "inoculation_method": "natural",
                      "planting_date": plant.strftime("%Y-%m-%d"), "inoculation_date": "",
                      "lat": round(lat, 4), "lon": round(lon, 4), "row_spacing_cm": row,
                      "fungicide": fung, "sclerotia_history": int(history)})

        def date(t):
            return (plant + pd.Timedelta(days=int(t))).strftime("%Y-%m-%d")

        se, _ = SYNTHETIC_SE_SP["apothecia_scouting"]
        for t in range(r1 - 7, r1 + 21, 7):
            obs.append((uid, date(t), "apothecia_scouting", int(A and rng.random() < se)))
        se, sp = SYNTHETIC_SE_SP["spore_trap_qpcr"]
        if rng.random() < 0.6:  # not every field has a trap
            for t in range(r1, r1 + 21, 7):
                obs.append((uid, date(t), "spore_trap_qpcr", int(rng.random() < (se if spore else 1 - sp))))
        se, sp = SYNTHETIC_SE_SP["petal_qpcr"]
        if rng.random() < 0.5:
            obs.append((uid, date(r1 + 10), "petal_qpcr", int(rng.random() < (se if court else 1 - sp))))
        obs.append((uid, date(105), "incidence", round(inc, 4)))
        obs.append((uid, date(105), "dsi", round(dsi, 2)))
        obs.append((uid, date(128), "yield_kg_ha", round(yld, 1)))

        # Images: drone RGB at R3 and R6, satellite (4-band) at R6, with disease maps
        for t, mod, ch, size in ((r1 + 20, "drone_rgb", 3, 64), (100, "drone_rgb", 3, 64), (100, "satellite", 4, 32)):
            frac = inc * (0.3 if t < 90 else 1.0)
            dmap = _blob_map(rng, frac, map_size)
            img = _canopy_image(rng, dmap, size, ch)
            if ch == 3:
                p = f"images/{uid}_{mod}_{t}.png"
                _save_png(img, out / p)
            else:
                p = f"images/{uid}_{mod}_{t}.npy"
                np.save(out / p, img)
            mp = f"maps/{uid}_{t}.npy" if mod == "drone_rgb" else ""
            if mp:
                np.save(out / mp, dmap)
            ims.append((uid, date(t), mod, p, mp))

    # ---------------------------------------------------------------- growth chamber (cut petiole)
    isolates = {"ISO_aggr": 1.0, "ISO_low": 0.4}
    inoc0 = pd.Timestamp("2025-03-01")
    for k in range(n_chamber):
        uid = f"G{k:04d}"
        cv = list(cultivars)[k % 4]
        iso = list(isolates)[(k // 4) % 2]
        chamber_rh = float(rng.choice([70.0, 90.0]))
        inoc = inoc0 + pd.Timedelta(days=int(k // 16) * 30)
        units.append({"unit_id": uid, "experiment_type": "growth_chamber", "season": f"GC{k // 16}",
                      "site": "chamber1", "cultivar": cv, "isolate": iso, "inoculation_method": "cut_petiole",
                      "planting_date": (inoc - pd.Timedelta(days=28)).strftime("%Y-%m-%d"),
                      "inoculation_date": inoc.strftime("%Y-%m-%d"), "lat": "", "lon": "",
                      "row_spacing_cm": "", "fungicide": "none", "sclerotia_history": ""})
        for dpi in range(-7, 15):
            day = (inoc + pd.Timedelta(days=dpi)).strftime("%Y-%m-%d")
            ts.append((uid, day, "air_temp", round(float(21 + rng.normal(0, 0.5)), 3)))
            ts.append((uid, day, "rh", round(float(chamber_rh + rng.normal(0, 2)), 3)))
        micro = chamber_rh > 80 or rng.random() < 0.5
        inf = micro and rng.random() < _sig(2.5 - cultivars[cv] + isolates[iso])
        rate = (4.5 * isolates[iso] + 2) * (1 - 0.25 * cultivars[cv]) if inf else 0.0
        for dpi in (3, 5, 7, 9, 11, 14):
            L = max(0.0, rate * (dpi - 2) + rng.normal(0, 1.5)) if inf else 0.0
            obs.append((uid, None, "lesion_length_mm", round(L, 2), dpi))
            p = f"images/{uid}_gc_{dpi}.png"
            _save_png(_stem_image(rng, min(L / 100, 1.0)), out / p)
            ims.append((uid, None, "gc_rgb", p, "", dpi))
        se, sp = SYNTHETIC_SE_SP["stem_qpcr"]
        obs.append((uid, None, "stem_qpcr", int(rng.random() < (se if inf else 1 - sp)), 7))
        obs.append((uid, None, "dsi", round(min(100.0, rate * 6) if inf else 0.0, 2), 14))
        obs.append((uid, None, "incidence", 1.0 if inf else 0.0, 14))  # one plant: incidence is 0 or 1

    pd.DataFrame(units).to_csv(out / "units.csv", index=False)
    o = pd.DataFrame([r + (None,) * (5 - len(r)) for r in obs], columns=["unit_id", "date", "assay", "value", "dpi"])
    o.to_csv(out / "observations.csv", index=False)
    pd.DataFrame(ts, columns=["unit_id", "date", "variable", "value"]).to_csv(out / "timeseries.csv", index=False)
    i = pd.DataFrame([r + (None,) * (6 - len(r)) for r in ims],
                     columns=["unit_id", "date", "modality", "path", "target_map_path", "dpi"])
    i.to_csv(out / "images.csv", index=False)
    (out / "README.txt").write_text("SYNTHETIC DATA generated by ssrforecast.synthetic. All parameter values are "
                                    "arbitrary simulation settings, not biological estimates.\n")
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--fields", type=int, default=240)
    ap.add_argument("--chamber-units", type=int, default=48)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    p = generate(a.out, a.fields, a.chamber_units, a.seed)
    print(f"synthetic dataset written to {p}")


if __name__ == "__main__":
    main()
