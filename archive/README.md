# Archive

| File | |
|---|---|
| [`wang-2025-deductive-inductive-logic.pdf`](wang-2025-deductive-inductive-logic.pdf) | *A Philosophical-Mathematical Representation of Deductive and Inductive Logic: A Diagrammatic and Formal Analysis*, February 2025 |
| [`wang-2025-quantum-model.pdf`](wang-2025-quantum-model.pdf) | *A Quantum Information Theoretic Model for the Philosophical-Mathematical Representation of Deductive and Inductive Logic*, February 2025; adapted as Section 7 of the paper |

## Changes from the February 2025 manuscripts

The paper in [`../paper.pdf`](../paper.pdf) makes each of these changes to the
[February 2025 version](wang-2025-deductive-inductive-logic.pdf).
Each item is a place where the 2025 text could not be stated in Lean as written; the Lean files
in [`../lean`](../lean) show the replacement.

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
   natural one is Bayesian confirmation, `P(t | s) > P(t)`. The draft now does this in §5.3
   (`Confirms`), and proves that white sightings raise the probability of "all swans are white"
   at every step (`allWhite_confirmed`, `allWhite_condProb_mono`).
6. **Give observations a type.** `s∃ ∈ S × ∃` makes an observation a pair, and `⊨` between pairs
   is undefined. Treat a study as a list of thoughts known to hold (`List (Thought W)`).
7. **Fix the codomain of `m`.** `m(i, h, l, 1) = ϕ or ∃` mixes thoughts and existences. Let `m`
   return a thought, with state `0` giving its negation; `m_spec` then proves the binary state is
   the truth value.
8. **Retitle §3** from "Mathematical Proof" to "Formalization" and cite the Lean file, or add a
   short appendix listing the theorem names below.

How the 2025 sections map to Lean:

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
