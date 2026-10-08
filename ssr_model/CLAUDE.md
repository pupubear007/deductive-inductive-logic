# Instructions for Claude Code in `ssr_model/`

This package is the trainable layer of the paper in `../paper/` (Lean proofs in `../lean/`).
Keep these invariants. If a change would break one, stop and say so instead of working around it.

1. **Hard constraints by construction.** Level probabilities come only from
   `mechanism.chain_probabilities`. Never predict P(infection) with an unconstrained head.
   `tests/test_mechanism.py` must pass.
2. **Adaptedness.** A sample's inputs use data dated on or before its issue day; targets are
   strictly after it (forecast) or on it (nowcast). Never add a feature computed from the target
   window.
3. **Observations are assay readings, not facts.** Binary readings enter the loss through the
   (se, sp) channel in `SSRModel.p_positive`. Do not feed a reading in as a ground-truth level.
4. **Interventions are not observations.** Inoculation methods in `COURT_INTERVENTIONS` set
   do(court = 1). Do not let such units train the spore or court steps.
5. **Split by season** (or chamber run), stratified by experiment type. Do not use random
   row-level splits for reported metrics.
6. **No invented biology.** Parameters in `synthetic.py` and the se/sp defaults are placeholders.
   Never present them as estimates, and never hard-code a biological value without a citation the
   user has confirmed.
7. Generated outcomes must pass through `mechanism.project_outcomes` before use, and evaluation
   reports the violation rate before projection.

Run `pytest -q` (in `ssr_model/`) before committing. For a quick check use
`python -m ssrforecast.train --config configs/synthetic.yaml --epochs 2`.
