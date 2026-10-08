# `ssrforecast`: a trainable Sclerotinia stem rot forecast

A PyTorch model for Sclerotinia stem rot (SSR) of soybean whose structure follows the paper's
framework. The learned parts are fitted inside a fixed logical skeleton.

| paper | here |
|---|---|
| worlds W, hard constraints (Def. 3.1, S1) | outputs can only describe admissible worlds (`mechanism.py`) |
| entailment / modus tollens (Thm 4.2, Prop 4.4) | **gates**: closed gate ⇒ next level has probability 0 by construction |
| levels (Def. 6.3) | chain spore ≤ court ≤ infection, with apothecia and external spores as sources |
| defeasible rules (Def. 5.12) | learned **step kernels** `q_court`, `q_inf` (normally progresses, can fail) |
| Bayesian confirmation (Def. 5.6) | calibrated level probabilities and assay likelihoods |
| assays (Def. 8.1) | observation channels with sensitivity / specificity |
| forecast horizon (S3/S4) | inputs strictly up to the issue day; targets strictly after it |
| next season unobserved (Thm 5.4) | test split = held-out whole seasons |
| prior over worlds (S12) | **conditional diffusion** over outcome vectors and over 2-D disease maps |

```
images (chamber, greenhouse, field, drone, satellite) ─┐
environment series (≤ issue day) ──────────────────────┼─► fused embedding z ─┬─► gated chain head ─► P(apothecia, spore, court, infection)
cultivar / isolate / management / inoculation ─────────┘                      ├─► outcome head ─────► E[incidence, DSI, lesion, yield]
                                                                              ├─► vector diffusion ─► samples of outcomes (projected onto W_K)
                                                                              └─► map diffusion ────► samples of disease maps
```

## Quick start

```sh
cd ssr_model
pip install -e ".[dev]"                                   # torch, numpy, pandas, pyyaml, pillow
python -m ssrforecast.synthetic --out data/synthetic      # example data in the required format
python -m ssrforecast.validate  --data data/synthetic --config configs/synthetic.yaml
python -m ssrforecast.train     --config configs/synthetic.yaml
python -m ssrforecast.evaluate  --run runs/synthetic --samples 8
pytest -q
```

On your own data: copy `configs/default.yaml`, point `data_dir` at your folder, run `validate`,
fix what it reports, then train. See **DATA_REQUIREMENTS.md** for every column.

## Outputs (`runs/<name>/`)

* `model.pt`: weights plus the normalisation, vocabularies and split, so evaluation reuses them
* `metrics.csv`: per-epoch losses
* `eval_val.json`, `eval_test.json`: per assay Brier score, log loss, AUROC, ECE, sensitivity and
  FPR (with the Prop. 8.6 check `positive_raises_risk`); MAE of expected outcomes; a hard-constraint
  audit; the diffusion constraint-violation rate before and after projection
* `predictions_*.csv`: per sample, every level probability and every expected reading next to the
  observed one

## Status and limits

* Version 0.1 is a scaffold. On synthetic data it shows the pipeline runs and the constraints
  hold. It is not evidence about real SSR.
* The chain assumes the gate and the lower level are conditionally independent given the inputs.
* Interventions are limited to inoculation, do(court = 1). Management (fungicide, row spacing)
  enters as covariates, not as do-operations, so the model does not answer counterfactual
  questions such as "loss avoided by spraying".
* The map head learns spatial patterns but is not yet an explicit dispersal kernel between fields.
* Every number in `synthetic.py` is an arbitrary simulation setting.
