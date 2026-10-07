# [A Philosophical-Mathematical Representation of Deductive and Inductive Logic: A Diagrammatic and Formal Analysis](paper.pdf)

**Author:** Hsuan Fu Wang  
**Date:** February 2025

## Citation

```bibtex
@misc{Wang:Deductive-Inductive-Logic-2025,
  author = {Hsuan Fu Wang},
  title = {{A Philosophical-Mathematical Representation of Deductive and Inductive Logic:
            A Diagrammatic and Formal Analysis}},
  year = {2025}
}
```

## Lean formalization

[`lean/WangLogic/Basic.lean`](lean/WangLogic/Basic.lean) formalizes Section 3 in plain Lean 4
(toolchain `v4.34.1`, the same one as [`../../lean`](../../lean)). It needs no Mathlib, has no
`sorry`, and uses only the axioms `propext`, `Quot.sound` and `Classical.choice`.

```sh
cd lean
lake build
```

| Paper | Lean |
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

[`lean/ComparatorChallenges/InductionNotDeduction.json`](lean/ComparatorChallenges/InductionNotDeduction.json)
states the headline theorem in the repository's
[Comparator](../../lean/ComparatorChallenges/README.md) format, and
[`lean/formalization.yaml`](lean/formalization.yaml) catalogues it in the format of
[`../../lean/formalization.yaml`](../../lean/formalization.yaml).

## Revisions to the manuscript

Each item below is a place where the current text cannot be stated in Lean as written; the Lean
file shows the replacement.

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
