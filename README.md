# Deduction and Induction in One Diagram: A Machine-Checked Semantics with Applications to Quantum Information and Plant Disease Diagnosis

[![CI](https://github.com/pupubear007/deductive-inductive-logic/actions/workflows/ci.yml/badge.svg)](https://github.com/pupubear007/deductive-inductive-logic/actions/workflows/ci.yml)

**Hsuan Fu Wang**, Department of Plant Pathology, University of Minnesota ·
[wan00965@umn.edu](mailto:wan00965@umn.edu)

**Paper:** [`paper.pdf`](paper.pdf) (submission version, October 2026) ·
**Lean formalization:** [`lean/`](lean/) · **Submission package:** [`submission/`](submission/)

> We give a diagrammatic representation of deductive and inductive reasoning in which deduction is a
> certain passage from thought to existence, and induction is a tentative, iterated passage from the
> study of existence to the theory of thought. We make the diagram precise in a possible-worlds
> semantics: thoughts are propositions about worlds, deduction is semantic entailment, polarity
> distinguishes affirmation from modus tollens, and an inductive theory is supported by a study when
> some world is consistent with both. Within this setting we prove that deduction is sound and
> transitive, that iterated observation can only narrow the supported theories, that a single
> counterexample falsifies a universal theory, and that no finite study entails a universal theory
> over an unobserved individual. The last statement is a formal version of Hume's problem of
> induction. As an application we realize the framework in quantum information theory: thoughts
> become subspaces, deduction becomes certain measurement and unitary evolution, induction becomes
> state tomography, and the classical setting is recovered as the case of commuting observables. The
> quantum model reproduces refinement, falsification and underdetermination; measurements in a
> single basis never determine a quantum state. As a second application we read plant disease
> diagnosis and forecasting in the framework: diagnosis is elimination and refinement, Koch's
> postulates are hypothetico-deductive confirmation, forecasting rules are defeasible, and a new
> theorem characterizes when an assay can resolve a diagnosis, of which the quantum basis limitation
> is an instance. We further give a Bayesian refinement of support, under which every genuine test a
> theory predicts confirms it, a ranked semantics for defeasible inference, and an order on levels
> under which iteration moves upward. Every theorem of the paper is verified in the Lean 4 proof
> assistant with Mathlib.

## Repository layout

| Path | Contents |
|---|---|
| [`paper.pdf`](paper.pdf) | The paper |
| [`paper/`](paper/) | LaTeX source |
| [`lean/`](lean/) | Lean 4 formalization of every theorem (package `WangLogic`) |
| [`submission/`](submission/) | arXiv source bundle, arXiv metadata, pre-submission checklist |
| [`aggr_model/`](aggr_model/) | Aggressiveness determinants of *S. sclerotiorum* from dual RNA-seq and growth-chamber assays: isolate-wise sparse models, fixed-program vs host-responsive test, assay resolution (Thm 8.2) (branch `ss-aggressiveness-determinants`, work in progress) |
| [`archive/`](archive/) | Earlier versions: the February 2025 manuscripts with the list of changes, and version 1 of this paper without the plant-pathology application |

## Building

**Lean** (needs [elan](https://github.com/leanprover/elan); the toolchain `v4.34.1` is installed
automatically):

```sh
cd lean
lake exe cache get     # download prebuilt Mathlib
lake build
lake env lean AxiomCheck.lean   # axioms used by every cited theorem
```

`Basic.lean` (Sections 3–6) uses only the Lean core library; the other files use Mathlib, pinned
to commit `d13f23b723b8a846827a245b89c10fc7d3f11612`. No file contains `sorry`, and every theorem
depends at most on the axioms `propext`, `Quot.sound` and `Classical.choice`. CI checks both on
every push.

**Paper** (needs a TeX Live installation with `latexmk`, or run the four commands):

```sh
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

## Theorems and their Lean declarations

All declarations are in the namespace `WangLogic`. Appendix A of the paper has the same table.

| File | Paper |
|---|---|
| [`Basic.lean`](lean/WangLogic/Basic.lean) | Sections 3–6 |
| [`Confirmation.lean`](lean/WangLogic/Confirmation.lean) | Remark 3.5; §5.3 Bayesian confirmation |
| [`Defeasible.lean`](lean/WangLogic/Defeasible.lean) | §5.4 defeasible inference |
| [`Levels.lean`](lean/WangLogic/Levels.lean) | §6.1 ordered levels |
| [`Quantum.lean`](lean/WangLogic/Quantum.lean) | Props. 7.1, 7.3, 7.5 |
| [`QuantumExtensions.lean`](lean/WangLogic/QuantumExtensions.lean) | Props. 7.2, 7.4, 7.6; entropy invariance |
| [`Plant.lean`](lean/WangLogic/Plant.lean) | §8 plant disease diagnosis and forecasting |

| Paper, Sections 3–6 | Lean |
|---|---|
| Def. 3.2, entailment `φ ⊨ ψ` | `Entails` |
| Def. 3.4, support `s ⟳ t` | `Supported` |
| Prop. 4.1, `(p → q) ⇔ (¬p ∨ q)` | `imp_iff_not_or` |
| Thm. 4.2, deduction is certain | `deduction_sound` |
| Prop. 4.3, chains of deduction | `Entails.trans` |
| Prop. 4.4, negative polarity is modus tollens | `modus_tollens` |
| Example 4.5, deductive swan | `swan_syllogism` |
| Prop. 5.1, Thm. 5.2, refinement is monotone | `Supported.of_cons`, `Supported.mono` |
| Prop. 5.3, falsification | `falsification` |
| Thm. 5.4, induction is not deduction | `induction_not_deduction` |
| Prop. 5.5, inductive swan | `allWhite_supported` |
| Def. 6.1, Prop. 6.2, foundation map | `m`, `m_spec` |

| Draft, Sections 3–6 extensions | Lean |
|---|---|
| Remark 3.5, `P(t ∧ s) > 0 ⇔ s ⟳ t` | `Prior.prob_and_pos_iff_supported` |
| Def. 5.6, confirmation `P(t \| s) > P(t)` (`tϕ↑`) and disconfirmation (`tϕ↓`) | `Prior`, `Confirms`, `Disconfirms` |
| Prop. 5.7, confirmation implies support | `Confirms.supported` |
| Thm. 5.8, hypothetico-deductive confirmation | `confirms_of_entails` |
| Prop. 5.9, iteration raises confirmation | `condProb_mono` |
| Prop. 5.10, falsification, quantitatively | `Disconfirms.of_refutes` |
| Cor. 5.11, Bayesian swan | `allWhite_confirmed`, `allWhite_condProb_mono` |
| Def. 5.12–Thm. 5.13, ranked defeasible entailment and rules of system **P** | `Normally`, `NormallyFrom`, `Normally.of_entails`, `.and`, `.or`, `.cut`, `.cautious_mono` |
| Prop. 5.14, defeasible inference is nonmonotone | `normally_nonmonotone` |
| Def. 6.3–Thm. 6.4, ordered levels; iteration moves upward | `Hierarchy`, `Hierarchy.Established.mono`, `.of_le`, `.tenable`, `.iterate`, `Hierarchy.Tenable.anti` |

| Draft, Section 7 | Lean |
|---|---|
| Prop. 7.1, quantum deductive certainty | `quantum_deduction_sound` (with `starProjection_comp_eq_iff`, `starProjection_orthogonal_eq_sub`) |
| Prop. 7.2, unitary deduction | `unitary_deduction`, `born_unitary_conj` |
| Prop. 7.3(1), `P_φ ≤ P_ψ ⇔ φ ⊨ ψ`, `P_¬φ = I − P_φ` | `classicalProj_le_iff`, `classicalProj_not`, `classicalProj_eq_sum` |
| Prop. 7.3(2), `Tr(ρ P_t P_o₁ ⋯ P_oₖ) > 0 ⇔ s ⟳ t` | `classical_trace_pos_iff` |
| Prop. 7.4, quantum logic is not distributive | `quantum_not_distributive` |
| Prop. 7.5(1), refinement | `QSupported.of_concat`, `likelihood_concat_le` |
| Prop. 7.5(2), underdetermination and maximum likelihood | `spinState_supported`, `spinUp_maxLikelihood`, `likelihood_spinState_up` |
| Prop. 7.5(3), falsification | `spinUp_refuted` |
| Prop. 7.5(4), the phase is invisible | `likelihood_spinState_phase` |
| §7.6, `S(UρU†) = S(ρ)` | `vonNeumannEntropy_unitary_conj`, `vonNeumannEntropy_eq_sum_eigenvalues` |
| Prop. 7.6, deduction and falsification under depolarizing noise | `born_depolarize_ge`, `noisy_spinUp_not_refuted` |

| Paper, Section 8 | Lean |
|---|---|
| Def. 8.1, assay and resolution | `AssayResult`, `Resolves` |
| Thm. 8.2, assay resolution | `resolves_iff` |
| Cor. 8.3, underdetermination by an assay | `undetermined`, `Resolves.pair_left` |
| Example 8.4, formae speciales of *Fusarium oxysporum* | `morphology_not_resolves`, `hostTest_resolves` |
| Prop. 8.5, Koch's postulates | `confirms_of_entails`, `modus_tollens` |
| Prop. 8.6, when a positive result raises the risk | `positive_raises_risk_iff` |
| Prop. 8.7, a defeasible forecasting rule | `sclerotinia_nonmonotone`, `apothecia_normally_infected`, `noApothecia_normally_not_infected` |

In Prop. 7.1 a projector is Mathlib's orthogonal projection `K.starProjection` onto a subspace of a
finite-dimensional complex inner product space, and `P ≤ Q` is subspace inclusion. Elsewhere in
Section 7 operators are complex matrices: a projector is a Hermitian idempotent matrix, a density
operator is positive semidefinite with trace 1, and `P ≤ Q` is written `Q * P = P`. In §5.3 the set
of worlds is finite and a prior is a strictly positive weight. Not formalized: Landauer's
principle (a physical law), `H(T | S) ≤ H(T)`, and the implementations and experimental proposals
of Section 7. The Discussion of the paper lists the open problems, with proposed approaches.

[`lean/ComparatorChallenges/InductionNotDeduction.json`](lean/ComparatorChallenges/InductionNotDeduction.json)
states Theorem 5.4 for independent checking with
[Comparator](https://github.com/leanprover/comparator), and
[`lean/formalization.yaml`](lean/formalization.yaml) describes the project in the
[formalization.yaml](https://github.com/mathlib-initiative/formalization.yaml) format.

## Citation

```bibtex
@unpublished{Wang2026DeductiveInductive,
  author = {Hsuan Fu Wang},
  title  = {Deduction and Induction in One Diagram: A Machine-Checked Semantics with
            Applications to Quantum Information and Plant Disease Diagnosis},
  note   = {Submitted. Lean formalization at
            \url{https://github.com/pupubear007/deductive-inductive-logic}},
  year   = {2026}
}
```

GitHub's "Cite this repository" button reads [`CITATION.cff`](CITATION.cff).

## License

The Lean code is licensed under the [Apache License 2.0](LICENSE). The paper text, its LaTeX
source and the PDFs are licensed under [CC BY 4.0](LICENSE-paper).
