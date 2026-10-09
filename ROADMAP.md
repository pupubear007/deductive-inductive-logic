# Research roadmap

Hsuan Fu Wang, Department of Plant Pathology, University of Minnesota

**Goal.** Make plant disease reasoning more transparent and better designed at every scale, from
field to cell. The guiding question at every scale is *what can a measurement tell apart?* In the
paper this is assay resolution (Thm 8.2); in modelling it is called identifiability.

Each stage below is a separate paper. Each one builds on the previous ones and adds real data, a
domain question, and a formal result checked in Lean.

| Stage | Scale | Question | Priority |
|---|---|---|---|
| 1 | Logic | Foundation: deduction, induction, assay resolution | Done |
| 2 | Field | Forecasting and its evaluation | **Core** |
| 3 | Cell | Cell models from sequencing data | **Next** (data in hand) |
| 4 | Device and landscape | Transport: microfluidics, spore dispersal, AFM mechanics | Later (after coursework and access) |

## Stage 1: Foundation (this repository)

*Deduction and Induction in One Diagram*: the semantics, the Bayesian and defeasible extensions,
the assay-resolution theorem (Thm 8.2), and a quantum model. Every theorem is verified in Lean 4.

- **Status:** complete; to be posted as an arXiv preprint after advisor review.
- **Role:** gives later papers definitions and theorems to cite. It is not the final word on any
  application.

## Stage 2 (core): Forecasting and its evaluation

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

## Stage 3 (next): Cell models from sequencing data

**Question.** Which parts of a plant–pathogen interaction model are actually determined by the
sequencing data, and which further experiment would determine the rest?

- **Data:** existing sequencing data. Record the organism, bulk or single-cell, time course or
  not, and the conditions.
- **Model:** a small mechanistic model of one pathway or response, chosen from the data.
  - Time-course expression data: a gene-regulatory model, as differential equations or a Boolean
    network.
  - A genome and expression profiles of a pathogen: a constraint-based metabolic model.
  - Start with a model small enough to understand completely.
- **Use of Thm 8.2:** test identifiability. Look for parameter settings, or network structures,
  that give the same predicted data but differ in a biological conclusion. Then rank candidate
  experiments (extra time points, mutants, conditions, a protein-level assay) by how many of
  these ambiguities each resolves. This turns "the model has too many parameters" into a design
  result: which experiment to run next.
- **New formal results:** identifiability for a simple class of models (for example linear
  differential equations), stated as assay resolution.
- **Collaborators needed:** systems biology or bioinformatics.

## Stage 4 (later): Transport and mechanics, from device to landscape

The same transport physics appears at two very different scales, so this stage treats them
together.

**Microfluidics and biochips.**
- In microchannels the Reynolds number is small, so the Navier–Stokes equations reduce to the
  linear Stokes equations. Flow is laminar, and mixing happens mainly by diffusion (Péclet
  number).
- Question: which channel design and readouts resolve a pathogen behaviour of interest, such as
  spore germination, hyphal guidance, root colonization, or on-chip pathogen detection?
- Simple exact flow solutions, such as pressure-driven flow in a channel, are small enough to
  formalize in Lean.

**Spore dispersal and surveillance design.**
- Model airborne spore transport as advection–diffusion driven by measured or modelled wind.
  Use more detailed flow models (canopy airflow, Navier–Stokes-based) only with collaborators who
  have that expertise.
- Use of Thm 8.2: two inoculum scenarios that give the same trap readings cannot be told apart,
  however long the traps run. Trap placement becomes a resolution problem.

**Atomic force microscopy.**
- Measure stiffness, adhesion and turgor-related mechanics of fungal cell walls, spores and
  infection structures, and fit contact-mechanics models.
- Question: can force data separate wall stiffness from turgor pressure, and which additional
  measurement would?

**Prerequisites and collaborators.**
- Coursework first: transport phenomena or fluid mechanics, and biophysics or quantitative cell
  biology.
- Access to microfabrication and AFM facilities.
- Collaborators in bioengineering, and in atmospheric or fluid modelling.

## Skills to build along the way

- Be able to explain every proof in the Stage 1 paper without notes. Most are short.
- Stage 2: Bayesian statistics and forecast verification.
- Stage 3: dynamical models of gene regulation or metabolism, parameter estimation, and
  identifiability analysis.
- Stage 4: transport phenomena (Stokes flow, advection–diffusion) and contact mechanics.

## Questions for the advisor

1. Which pathosystem and datasets are realistic for Stage 2?
2. Which part of the existing sequencing data could support a Stage 3 model, and which pathway or
   response is most interesting?
3. Is the Stage 1 preprint appropriate to post now, and under which affiliation and
   acknowledgements?
4. Which courses fit Stage 4 prerequisites, and are microfluidics or AFM facilities available?
5. Who in the department or in the wider university could collaborate on statistics, systems
   biology, and transport modelling?
