# Research roadmap

Hsuan Fu Wang, Department of Plant Pathology, University of Minnesota

**Goal.** Use one framework for reasoning to make plant disease diagnosis, forecasting and
surveillance more transparent and better designed. The framework has deduction ("this hypothesis
predicts that observation"), induction ("these data leave these hypotheses open") and assay
resolution ("can this measurement tell the hypotheses apart at all?"). Each stage below is a
separate paper. Each one builds on the previous stage and adds real data, a domain question, and a
formal result checked in Lean.

## Stage 1: Foundation (this repository)

*Deduction and Induction in One Diagram*: the semantics, the Bayesian and defeasible extensions,
the assay-resolution theorem (Thm 8.2), and a quantum model. Every theorem is verified in Lean 4.

- **Status:** complete; to be posted as an arXiv preprint after advisor review.
- **Role:** gives later papers definitions and theorems to cite. It is not the final word on any
  application.

## Stage 2: Forecasting and its evaluation

**Question.** For a weather-based disease forecast, when does a positive forecast justify action,
and which additional observations would make it more reliable?

- **Pathosystem:** to choose with the advisor, based on available data. Candidates are
  Sclerotinia stem rot of soybean, Fusarium head blight of wheat, and potato late blight.
- **Model, in three layers:**
  1. disease-triangle constraints as hard gates (no susceptible stage or no inoculum means no
     risk);
  2. defeasible expert rules, with their ranking estimated from field records;
  3. a calibrated probability model, P(disease | weather, field, scouting).
- **Evaluation:**
  - discrimination: ROC, Youden's J (Prop 8.6);
  - calibration: reliability diagrams, Brier score;
  - decision value: cost–loss threshold, net benefit compared with always or never spraying;
  - transfer: leave-one-site-year-out validation.
- **Use of Thm 8.2:** find field-seasons that share all predictor values but differ in outcome.
  Rank candidate extra assays (spore traps, apothecia scouting, canopy sensors, qPCR) by how many
  of these ambiguities each resolves.
- **New formal results:** the cost–loss decision threshold, and resolution under a noisy assay
  (sensitivity < 1).
- **Collaborators needed:** statistics or epidemiology.

## Stage 3: Dispersal and surveillance design

**Question.** Where should spore traps or sensors be placed so that the source and timing of
inoculum can be identified?

- **Model:** airborne spore transport. Start with advection–diffusion driven by measured or
  modelled wind; use more detailed flow models (canopy airflow, Navier–Stokes-based) only with
  collaborators who have that expertise.
- **Use of Thm 8.2:** two inoculum scenarios that give the same trap readings cannot be told
  apart, however long the traps run. This turns trap placement into a resolution problem.
- **Deliverable:** a method for comparing surveillance designs by what they can resolve, tested
  on simulated and, if available, real trap data.
- **Collaborators needed:** atmospheric or fluid-dynamics modelling.

## Skills to build along the way

- Be able to explain every proof in the Stage 1 paper without notes. Most are short.
- Bayesian statistics and forecast verification for Stage 2.
- Basics of transport (advection–diffusion) models for Stage 3.

## Questions for the advisor

1. Which pathosystem and datasets are realistic for Stage 2?
2. Is the Stage 1 preprint appropriate to post now, and under which affiliation and
   acknowledgements?
3. Who in the department or in the wider university could collaborate on statistics and on
   dispersal modelling?
