# `ssaggr`: aggressiveness determinants of *Sclerotinia sclerotiorum*

`ssaggr` uses dual RNA-seq, growth-chamber lesion assays and in vitro assays to search for
pathogen genes that determine isolate aggressiveness. It also tests the two hypotheses in the
proposal, a **fixed program** versus **host-responsive regulation**, in a way that is honest
about how few isolates there are.

It reads the files your existing pipeline already produces:
- Salmon pathogen count matrices, sample names `<hpi>_<rep>_<host>_<isolate>`
- CH4 Table 2 sAUDPC means
- the per-plant lesion files
- oxalic acid, appressorium and radial growth tables
- the secretome annotation

## What it computes

| proposal | module | output |
|---|---|---|
| Objective 1: aggressiveness and its consistency across hosts | `stability.py` | Finlay–Wilkinson / Eberhart–Russell slope and deviation per isolate, mean rank, rank range, isolate × crop share of variance |
| Objective 2: host-responsiveness index | `responsiveness.py` | R²_host within each isolate after removing colonization, and its exact-permutation correlation with stability |
| Objective 2 as a prediction test | `model.py`, `cv.py` | sparse elastic-net models (FISTA): shared weights (fixed program) versus shared + host-specific weights (host-responsive); leave-one-isolate-out, permutation null, cross-host transfer |
| candidate determinants | `cv.stability_selection` | selection frequency per gene and host; class `shared_concordant`, `shared_opposite` or `specific_<host>` (compare CH4 Fig. 10); joined to effector / CAZyme annotation |
| paper, Theorem 8.2 | `resolution.py` | whether each assay (oxalic acid, appressoria, growth, transcript share, expression score) **resolves** the aggressiveness class on each host, the witness pairs, and the chance rate |
| design | `power.py` | simulated power: how well 6 / 12 / 24 isolates separate the two hypotheses |

## Why it is built this way

- **The isolate is the unit of evidence.** Six isolates is n = 6, whatever the number of
  libraries. A flexible model (deep net, random forest) would learn isolate identity from
  replicates and look perfect. Every reported number is therefore leave-one-isolate-out, with
  preprocessing fitted inside each fold. Every claim gets a permutation null that reruns the
  whole procedure. With 6 isolates there are only 720 orderings.
- **Colonization is a covariate, not a confounder.** Following CH3 and the proposal (5.6), logit
  pathogen transcript share enters every model unpenalised. Binomial thinning to equal depth is
  available (`thin_to_equal_depth`).
- **"Aggressive" is host-indexed.** In the paper's terms, the hypothesis H is "aggressive on host
  h", and an assay may resolve it on one host and not another (Theorem 8.2). The resolution table
  makes that explicit, which is the formal version of CH4's main claim.
- **No diffusion or deep generative model here.** With 6–12 isolates they cannot be trained
  honestly. The trial and phenotyping model is on branch `ssr-trial-phenotyping`.

## Quick start

```sh
cd aggr_model
pip install -e ".[dev]"
pytest -q

# synthetic data with a planted answer
python -m ssaggr.synthetic --out data/syn --regime responsive --isolates 12
python -m ssaggr.run --config data/syn/config.yaml --quick

# your data (paths in the config point at your RNAseq_paper folder)
python -m ssaggr.run --config configs/real_local.yaml            # ~30 min with 200 permutations
python -m ssaggr.run --config configs/real_local.yaml --quick    # first look, ~5 min

# design: isolates needed to separate fixed vs host-responsive
python -m ssaggr.power --isolates 6 12 24 --seeds 5
```

Outputs go to `out_dir`: `report.md`, `stability.csv`, `responsiveness.csv`, `loio_*.csv`,
`permutation_single.csv`, `cross_host.csv`, `determinants.csv`, `resolution.csv`.

## Before trusting a result

1. **Phenotype and RNA-seq hosts must match.** CH4 Table 2 soybean values come from
   Dwight / 52-82B; the RNA-seq host is Williams 82. Point `phenotype_table` at the matching
   phenotype if one exists.
2. **Decide the stability cut-offs** (`stability.classify`) before looking at expression results.
3. **Host and batch are confounded.** Soybean and sunflower were sequenced and quantified
   separately, so R²_host and cross-host transfer include batch. A shared control library across
   runs (proposal 5.4) is what fixes this.
4. **Read `p_chance`** next to every "resolves": with 4–6 isolates, perfect separations happen by
   chance.
