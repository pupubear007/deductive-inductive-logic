# [A Philosophical-Mathematical Representation of Deductive and Inductive Logic: A Diagrammatic and Formal Analysis](paper.pdf)

**Author:** Hsuan Fu Wang  
**Date:** Preprint draft, October 7, 2026 (revised from the [February 2025 version](original.pdf))

Section 7 is an application to quantum information theory, adapted from the companion manuscript
[*A Quantum Information Theoretic Model for the Philosophical-Mathematical Representation of
Deductive and Inductive Logic*](original-quantum.pdf) (February 2025). Its three propositions are
formalized in Lean with Mathlib; see below.

## Citation

```bibtex
@misc{Wang:Deductive-Inductive-Logic-2026,
  author = {Hsuan Fu Wang},
  title = {{A Philosophical-Mathematical Representation of Deductive and Inductive Logic:
            A Diagrammatic and Formal Analysis}},
  note = {Preprint draft},
  year = {2026}
}
```

## Building the PDF

The LaTeX source is in [`build/`](build/) and follows the layout of the preprints in
[`../../preprints`](../../preprints).

```sh
cd build
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

## Lean formalization

[`lean/WangLogic/Basic.lean`](lean/WangLogic/Basic.lean) formalizes Sections 3–6 in plain Lean 4
(toolchain `v4.34.1`, the same one as [`../../lean`](../../lean)), without Mathlib.
[`lean/WangLogic/Quantum.lean`](lean/WangLogic/Quantum.lean) formalizes Propositions 7.1–7.3 with
Mathlib, pinned to commit `d13f23b723b8a846827a245b89c10fc7d3f11612` (the same pin as
[`../../lean/lakefile.lean`](../../lean/lakefile.lean)). Neither file has a `sorry`, and every
declaration uses at most the axioms `propext`, `Quot.sound` and `Classical.choice`.

```sh
cd lean
lake exe cache get
lake build
```

| 2025 version | Lean |
|---|---|
| §3.1 `Φ`, `∃` | `Thought W := W → Prop` over a type of possible worlds `W` |
| §3.3 `ϕ → ∃` as `¬ϕ ∨ ∃` | `imp_iff_not_or` |
| §3.3 `Φ_ϕ ⊆ f⁻¹(∃_∃)` | `Entails`, `deduction_sound`, `Entails.trans` |
| §2.2 deductive swan example | `swan_syllogism` |
| `↑`, `↓`; `ϕ → ∃↓` | `Polarity`, `Polarity.apply`, `modus_tollens` |
| §3.4 `P(T_ϕ \| S_∃) > 0` | `Supported` (some world fits the study and the theory) |
| §3.4 iterative refinement `i` | `Supported.of_cons`, `Supported.mono` |
| §4 Popper's falsifiability | `falsification` |
| §4 Hume's problem of induction | `induction_not_deduction` |
| §2.2 inductive swan example | `allWhite_supported` |
| §3.5 `m : I × H × L × {0,1} → Φ` | `m`, `m_spec` |

| Draft, Section 7 | Lean (`WangLogic/Quantum.lean`) |
|---|---|
| Prop. 7.1, quantum deductive certainty | `quantum_deduction_sound` (with `starProjection_comp_eq_iff`, `starProjection_orthogonal_eq_sub`) |
| Prop. 7.2(1), `P_φ ≤ P_ψ ⇔ φ ⊨ ψ`, `P_¬φ = I − P_φ` | `classicalProj_le_iff`, `classicalProj_not`, `classicalProj_eq_sum` |
| Prop. 7.2(2), `Tr(ρ P_t P_o₁ ⋯ P_oₖ) > 0 ⇔ s ⟳ t` | `classical_trace_pos_iff` |
| Prop. 7.3(1), refinement | `QSupported.of_concat`, `likelihood_concat_le` |
| Prop. 7.3(2), underdetermination and maximum likelihood | `spinState_supported`, `spinUp_maxLikelihood`, `likelihood_spinState_up` |
| Prop. 7.3(3), falsification | `spinUp_refuted` |
| Prop. 7.3(4), the phase is invisible | `likelihood_spinState_phase` |

In Prop. 7.1 a projector is Mathlib's orthogonal projection `K.starProjection` onto a subspace of a
finite-dimensional complex inner product space, and `P ≤ Q` is subspace inclusion. In Props. 7.2
and 7.3 operators are complex matrices: a projector is a Hermitian idempotent matrix, a density
operator is positive semidefinite with trace 1, and `P ≤ Q` is written `Q * P = P`. The informal
parts of Section 7 (unitary evolution, failure of distributivity, entropy, implementations) are not
formalized.

[`lean/ComparatorChallenges/InductionNotDeduction.json`](lean/ComparatorChallenges/InductionNotDeduction.json)
states the headline theorem in the repository's
[Comparator](../../lean/ComparatorChallenges/README.md) format, and
[`lean/formalization.yaml`](lean/formalization.yaml) catalogues it in the format of
[`../../lean/formalization.yaml`](../../lean/formalization.yaml).

## Revisions to the manuscript

The draft in [`paper.pdf`](paper.pdf) makes each of these changes to the [2025 version](original.pdf).
Each item is a place where the 2025 text could not be stated in Lean as written; the Lean file
shows the replacement.

1. **Rename the set of existences.** `∃` is used for the existential quantifier, the set of
   outcomes, and an element of that set. Use `W` (possible worlds) for the set and `w` for an
   element.
2. **Make Axiom 1 a relation, not a function.** `f : Φ → ∃` must send each thought to exactly one
   existence, but a thought entails many. Define entailment `ϕ ⊨ ψ` as
   `{w | ϕ(w)} ⊆ {w | ψ(w)}` (`Entails`) and drop `f`.
3. **Replace "Thus, `ϕ → ∃` is valid" with theorems that have content.** §3.3–3.5 assume their
   conclusions in the axioms. State instead what the framework proves: soundness
   (`deduction_sound`), modus tollens for `↓` (`modus_tollens`), falsification
   (`falsification`), and that no finite study entails a universal theory
   (`induction_not_deduction`). The last is the paper's strongest formal point and supports the
   Hume discussion in §4.
4. **Define the inductive criterion precisely.** `P(T_ϕ | S_∃) > 0` is undefined when
   `P(S_∃) = 0` and is met by every theory consistent with the data. Present it as
   "the theory is consistent with the study" (`Supported`) and state monotonicity under new
   observations (`Supported.mono`) as the meaning of iteration `i`.
5. **Define the inductive polarity.** "Ascending probability" for `tϕ↑` needs a definition. A
   natural one is Bayesian confirmation, `P(t | s) > P(t)`; formalizing that needs Mathlib's
   probability library.
6. **Give observations a type.** `s∃ ∈ S × ∃` makes an observation a pair, and `⊨` between pairs
   is undefined. Treat a study as a list of thoughts known to hold (`List (Thought W)`).
7. **Fix the codomain of `m`.** `m(i, h, l, 1) = ϕ or ∃` mixes thoughts and existences. Let `m`
   return a thought, with state `0` giving its negation; `m_spec` then proves the binary state is
   the truth value.
8. **Retitle §3** from "Mathematical Proof" to "Formalization" and cite the Lean file, or add a
   short appendix listing the theorem names above.
