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
python -m ssaggr.run --config configs/real_local.yaml --quick    # first look, ~10 min
python -m ssaggr.run --config configs/real_local.yaml            # 24, 48, 96 hpi and pooled; ~1 h with 200 permutations

# design: isolates needed to separate fixed vs host-responsive
python -m ssaggr.power --isolates 6 12 24 --seeds 5
```

Outputs go to `out_dir`: `report.md` (with a summary table across time points),
`summary_by_timepoint.csv` and `stability.csv`, plus one folder per time-point set (`tp_24`,
`tp_48`, `tp_96`, `tp_24_48_96`) holding `responsiveness.csv`, `loio_single.csv`,
`permutation_single.csv`, `joint_fixed_vs_responsive.csv`, `cross_host.csv`, `determinants.csv`
and `resolution.csv`.

## Scales

- **Expression target:** sAUDPC standardized within host (z-score), as in the proposal
  (`target_scale: z`). `raw`, `rank` and `log` are options.
- **Stability regression:** sAUDPC as measured (`stability_scale: raw`, Eberhart–Russell). A
  within-host z-score cannot be used here, because it sets every host mean, and so the host
  index, to zero. `log` is an option.

## Time points

Each time point (24, 48, 96 hpi) is analysed on its own, then all three are pooled (intercepts
per host × time point, shared gene weights). The per-time-point runs show *when* a signal
appears. The pooled run uses all libraries.

## Before trusting a result

1. **Phenotype source.** The default is the final Table S1 (`outputs_fig2/TableS1_long.csv`).
   Its soybean column is the same dataset as the CH3 Fig. 1 sAUDPC.
2. **Your hypothesis goes in the config** (`hypothesis_classes` or `stability_cutoffs`) before
   looking at expression results. The code never invents class cut-offs.
3. **Host and batch are confounded.** Soybean and sunflower were sequenced and quantified
   separately, so R²_host and cross-host transfer include batch. A shared control library across
   runs (proposal 5.4) is what fixes this.
4. **Read `p_chance`** next to every "resolves": with 4–6 isolates, perfect separations happen by
   chance.
