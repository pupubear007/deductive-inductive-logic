# Prompt for Claude Code

Paste everything below the line into Claude Code, opened in a clone of
`pupubear007/deductive-inductive-logic`. Edit the bracketed parts first.

---

You are working in my repository `deductive-inductive-logic`, on branch `ssr-forecast-model`
(run `git checkout ssr-forecast-model`). The folder `ssr_model/` holds `ssrforecast`, a PyTorch
model for Sclerotinia stem rot (SSR) of soybean. Its structure follows my paper (`paper.pdf`,
Section 8, plus the Lean files in `lean/WangLogic/`): hard gates are entailments, step kernels
are defeasible rules, level probabilities are the Bayesian layer, and two conditional diffusion
heads act as a learned prior over worlds.

**Before changing anything**
1. Read `ssr_model/CLAUDE.md` (the invariants you must keep), `ssr_model/README.md`,
   `ssr_model/DATA_REQUIREMENTS.md` and `ssr_model/ssrforecast/mechanism.py`.
2. `cd ssr_model && pip install -e ".[dev]" && pytest -q`. Report the result.

**Task: get my data into the model**
3. My data are in `[PATH TO MY DATA FOLDER]`. They are [describe: e.g. growth-chamber lesion
   lengths at 3-14 dpi for 4 checklines x 7 isolates in an Excel sheet; greenhouse photos; chamber
   temperature/RH logs; field incidence/DSI from colleagues as CSV; drone images].
   Inspect the files and tell me how each column maps onto `units.csv`, `observations.csv`,
   `timeseries.csv` and `images.csv`. **Ask me** about anything ambiguous (units of measure,
   percent vs fraction, dpi vs date, which column is the unit id, band order of multispectral
   images) before converting. Do not guess.
4. Write a converter `ssr_model/scripts/convert_<source>.py` that produces the four tables in
   `data/real/`. Never modify my original files.
5. Copy `configs/default.yaml` to `configs/real.yaml`, keeping only the assays, series
   variables and image modalities I actually have. Leave se/sp at 1.0 unless I give measured
   values.
6. Run `python -m ssrforecast.validate --data data/real --config configs/real.yaml`. Show me
   the DATA REQUEST section and **ask me** which missing items I can supply. Stop here until I
   answer.

**Task: train and report**
7. Train with `python -m ssrforecast.train --config configs/real.yaml`. Use a GPU if
   available; otherwise reduce `model.dim` and `diffusion.steps` and tell me.
8. Report test-season results from `eval_test.json` as a table: per assay n, prevalence,
   Brier, AUROC, ECE, sensitivity, FPR and `positive_raises_risk`; MAE of continuous outcomes;
   the constraint audits. Also compare against two baselines you add: prevalence-only, and
   logistic regression on the tabular plus weather-summary features. Say plainly if the model does
   not beat them.
9. Run the **null control**: shuffle the targets across units (or use my file
   `synthetic_Field_trial_data.csv`, whose columns are independent random draws). Performance
   must drop to chance; if it does not, find the leak.

**Rules**
- Keep every invariant in `ssr_model/CLAUDE.md`. In particular, never replace
  `chain_probabilities` with an unconstrained classifier, never use data after the issue day as
  input, and never split rows at random for reported metrics.
- Do not invent biological parameter values. If something needs one, ask me.
- Commit on `ssr-forecast-model` in small commits with clear messages. Do not push to `main`.
- When you finish a step, tell me in a few lines what changed and what you need from me next.

**Later extensions, not now:** a pretrained image backbone for drone and satellite data; an
explicit between-field dispersal kernel at the spore level; do-operators for fungicide timing
(counterfactual loss avoided); Lean statements for the chained-forecast and forecast-horizon
results in `lean/WangLogic/`.
