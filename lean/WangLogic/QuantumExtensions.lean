import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.LinearAlgebra.Matrix.Charpoly.Basic
import WangLogic.Quantum

/-!
# Section 7: unitary evolution, non-distributivity, entropy and noise

The informal claims of Section 7, stated and proved in the matrix model of `WangLogic.Quantum`.

* `unitary_deduction`: if `U P_φ U† ≤ P_ψ`, unitary evolution carries every state with property
  `φ` to one with property `ψ`; `born_unitary_conj`: evolving state and observable together
  leaves Born probabilities unchanged.
* `quantum_not_distributive`: the subspaces of a qubit do not form a distributive lattice.
* `vonNeumannEntropy_unitary_conj`: unitary evolution preserves von Neumann entropy.
* `born_depolarize`, `born_depolarize_ge`, `noisy_spinUp_not_refuted`: under depolarizing noise
  deductive certainty degrades by at most the noise rate, and a single outcome no longer refutes
  a theory outright.
-/

open Matrix
open scoped ComplexOrder

namespace WangLogic

variable {n : Type*} [Fintype n] [DecidableEq n]

/-! ## Unitary evolution is deductive -/

/-- **Proposition 7.2 (unitary deduction).** If `U P_φ U† ≤ P_ψ` (as `P_ψ (U P_φ U†) = U P_φ U†`)
and `P_φ v = v`, then `P_ψ (U v) = U v`. -/
theorem unitary_deduction {U Pφ Pψ : Matrix n n ℂ} (hU : Uᴴ * U = 1)
    (hle : Pψ * (U * Pφ * Uᴴ) = U * Pφ * Uᴴ) {v : n → ℂ} (hv : Pφ *ᵥ v = v) :
    Pψ *ᵥ (U *ᵥ v) = U *ᵥ v := by
  have key : (U * Pφ * Uᴴ) *ᵥ (U *ᵥ v) = U *ᵥ v := by
    rw [mulVec_mulVec, Matrix.mul_assoc (U * Pφ), hU, Matrix.mul_one, ← mulVec_mulVec, hv]
  rw [← key, mulVec_mulVec, hle]

/-- Evolving the state and the observable by the same unitary leaves Born probabilities
unchanged. -/
theorem born_unitary_conj {U ρ P : Matrix n n ℂ} (hU : Uᴴ * U = 1) :
    born (U * ρ * Uᴴ) (U * P * Uᴴ) = born ρ P := by
  unfold born
  have h : U * ρ * Uᴴ * (U * P * Uᴴ) = U * (ρ * P) * Uᴴ := by
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc Uᴴ U, hU, Matrix.one_mul]
  rw [h, Matrix.trace_mul_comm, ← Matrix.mul_assoc, hU, Matrix.one_mul]

/-! ## Quantum logic is not distributive -/

/-- The state space of a qubit. -/
abbrev Qubit := Fin 2 → ℂ

/-- The subspace spanned by `|↑⟩` (spin up along `z`). -/
noncomputable def zUp : Submodule ℂ Qubit := ℂ ∙ ![1, 0]

/-- The subspace spanned by `|↓⟩` (spin down along `z`). -/
noncomputable def zDown : Submodule ℂ Qubit := ℂ ∙ ![0, 1]

/-- The subspace spanned by `(|↑⟩ + |↓⟩)/√2` (spin up along `x`). -/
noncomputable def xPlus : Submodule ℂ Qubit := ℂ ∙ ![1, 1]

theorem zUp_sup_zDown : zUp ⊔ zDown = ⊤ := by
  refine eq_top_iff.2 fun v _ => Submodule.mem_sup.2 ?_
  refine ⟨v 0 • ![1, 0], Submodule.mem_span_singleton.2 ⟨v 0, rfl⟩,
    v 1 • ![0, 1], Submodule.mem_span_singleton.2 ⟨v 1, rfl⟩, ?_⟩
  funext i
  fin_cases i <;> simp

theorem xPlus_inf_eq_bot {a b : ℂ} (hab : a ≠ b) : xPlus ⊓ (ℂ ∙ ![a, b]) = ⊥ := by
  refine (Submodule.eq_bot_iff _).2 fun v hv => ?_
  obtain ⟨c, hc⟩ := Submodule.mem_span_singleton.1 hv.1
  obtain ⟨d, hd⟩ := Submodule.mem_span_singleton.1 hv.2
  have h0 := congrFun (hc.trans hd.symm) 0
  have h1 := congrFun (hc.trans hd.symm) 1
  simp at h0 h1
  have hc0 : c = 0 := by
    by_contra hc0
    apply hab
    have : d * a = d * b := by rw [← h0, ← h1]
    have hd0 : d ≠ 0 := by
      rintro rfl
      simp at h0
      exact hc0 h0
    exact mul_left_cancel₀ hd0 this
  rw [← hc, hc0, zero_smul]

/-- **Proposition 7.4 (failure of distributivity).** `P_{x+} ∧ (P_{z+} ∨ P_{z-}) = P_{x+}`, while
`(P_{x+} ∧ P_{z+}) ∨ (P_{x+} ∧ P_{z-}) = 0`. -/
theorem quantum_not_distributive :
    xPlus ⊓ (zUp ⊔ zDown) = xPlus ∧ (xPlus ⊓ zUp) ⊔ (xPlus ⊓ zDown) = ⊥ ∧ xPlus ≠ ⊥ := by
  refine ⟨by rw [zUp_sup_zDown, inf_top_eq], ?_, ?_⟩
  · rw [zUp, zDown, xPlus_inf_eq_bot one_ne_zero, xPlus_inf_eq_bot zero_ne_one, bot_sup_eq]
  · intro h
    have : (![1, 1] : Qubit) ∈ xPlus := Submodule.mem_span_singleton_self _
    rw [h, Submodule.mem_bot] at this
    have := congrFun this 0
    simp at this

/-! ## Von Neumann entropy -/

/-- The von Neumann entropy `S(ρ) = -Tr ρ ln ρ = -∑ λ ln λ`, summed over the roots of the
characteristic polynomial (the eigenvalues, with multiplicity). -/
noncomputable def vonNeumannEntropy (ρ : Matrix n n ℂ) : ℝ :=
  (ρ.charpoly.roots.map fun z => Real.negMulLog z.re).sum

/-- For a Hermitian `ρ` the definition is `∑ᵢ -λᵢ ln λᵢ` over its eigenvalues. -/
theorem vonNeumannEntropy_eq_sum_eigenvalues {ρ : Matrix n n ℂ} (hρ : ρ.IsHermitian) :
    vonNeumannEntropy ρ = ∑ i, Real.negMulLog (hρ.eigenvalues i) := by
  rw [vonNeumannEntropy, hρ.roots_charpoly_eq_eigenvalues, Multiset.map_map]
  rfl

/-- **Deduction preserves entropy.** `S(U ρ U†) = S(ρ)` for unitary `U`. -/
theorem vonNeumannEntropy_unitary_conj {U : Matrix n n ℂ} (hU : Uᴴ * U = 1)
    (ρ : Matrix n n ℂ) : vonNeumannEntropy (U * ρ * Uᴴ) = vonNeumannEntropy ρ := by
  unfold vonNeumannEntropy
  rw [charpoly_mul_comm, ← Matrix.mul_assoc, hU, Matrix.one_mul]

/-! ## Depolarizing noise -/

/-- The depolarizing channel `ρ ↦ (1 - ε) ρ + ε I / d`. -/
noncomputable def depolarize (ε : ℝ) (ρ : Matrix n n ℂ) : Matrix n n ℂ :=
  ((1 - ε : ℝ) : ℂ) • ρ + ((ε / Fintype.card n : ℝ) : ℂ) • (1 : Matrix n n ℂ)

theorem born_depolarize (ε : ℝ) (ρ P : Matrix n n ℂ) :
    born (depolarize ε ρ) P = (1 - ε) * born ρ P + ε / Fintype.card n * P.trace.re := by
  simp only [born, depolarize, Matrix.add_mul, Matrix.smul_mul, Matrix.one_mul,
    Matrix.trace_add, Matrix.trace_smul, smul_eq_mul, Complex.add_re, Complex.re_ofReal_mul]

/-- **Proposition 7.6, first part (noisy deductive certainty).** If `ρ` passes the test `P` with
certainty, its depolarized version passes with probability at least `1 - ε`. -/
theorem born_depolarize_ge {ε : ℝ} (hε : 0 ≤ ε) {ρ P : Matrix n n ℂ} (hP : IsProjector P)
    (h : born ρ P = 1) : 1 - ε ≤ born (depolarize ε ρ) P := by
  rw [born_depolarize, h, mul_one]
  have htr : 0 ≤ P.trace.re := by
    simpa [born] using born_nonneg (Matrix.PosSemidef.one (R := ℂ)) hP
  have : 0 ≤ ε / Fintype.card n * P.trace.re :=
    mul_nonneg (div_nonneg hε (Nat.cast_nonneg _)) htr
  linarith

/-- **Proposition 7.6, second part (falsification becomes statistical under noise).** After
depolarizing noise of rate `ε > 0`, the state `|↑⟩` gives the outcome `↓` with probability `ε / 2 >
0`, so one outcome `↓` no longer refutes it. -/
theorem noisy_spinUp_not_refuted {ε : ℝ} (hε : 0 < ε) :
    born (depolarize ε (pureState spinUp)) (basisProj 1) = ε / 2 ∧
      0 < born (depolarize ε (pureState spinUp)) (basisProj 1) := by
  have h : born (depolarize ε (pureState spinUp)) (basisProj 1) = ε / 2 := by
    rw [born_depolarize, born_pureState_basisProj]
    simp [spinUp, basisProj, Matrix.trace_diagonal]
  exact ⟨h, h ▸ half_pos hε⟩

end WangLogic
