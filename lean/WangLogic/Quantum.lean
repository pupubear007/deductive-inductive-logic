import Mathlib.Analysis.InnerProductSpace.Projection.Basic
import Mathlib.Analysis.InnerProductSpace.Projection.Submodule
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.LinearAlgebra.Matrix.PosDef
import WangLogic.Basic

/-!
# Section 7: a quantum information model

Companion to Section 7 of the preprint draft. Two models of the same quantum mechanics are used.

* **Projections** (Proposition 7.1). A thought is a subspace `K` of a finite-dimensional complex
  inner product space `E`, realized by its orthogonal projection `K.starProjection`. The paper's
  order `P ≤ Q` on projectors is inclusion of ranges, `K ≤ L`, equivalently `Q ∘ P = P`
  (`starProjection_comp_eq_iff`), and `I - P` is the projection onto `Kᗮ`.
* **Matrices** (Propositions 7.2 and 7.3). Operators on `ℂ^W` are matrices `Matrix W W ℂ`, a
  density operator is a positive semidefinite matrix of trace `1`, and the Born probability of
  outcome `P` in state `ρ` is `Tr(ρ P)`. For the classical case the order `P ≤ Q` is taken in its
  algebraic form `Q * P = P`.
-/

open Matrix
open scoped InnerProductSpace ComplexOrder

namespace WangLogic

/-! ## Proposition 7.1: deduction as certain measurement -/

section Projection

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E] [FiniteDimensional ℂ E]

/-- The paper's two readings of `P ≤ Q` for projectors agree: `ran P ⊆ ran Q` iff `Q P = P`. -/
theorem starProjection_comp_eq_iff (K L : Submodule ℂ E) :
    L.starProjection ∘L K.starProjection = K.starProjection ↔ K ≤ L := by
  constructor
  · intro h v hv
    rw [← Submodule.starProjection_eq_self_iff, ← (Submodule.starProjection_eq_self_iff).2 hv]
    exact congrArg (fun T => T v) h
  · intro h
    ext v
    exact (Submodule.starProjection_eq_self_iff).2 (h (K.starProjection_apply_mem v))

/-- Negative polarity `P↓ = I - P` is the projection onto the orthogonal complement. -/
theorem starProjection_orthogonal_eq_sub (K : Submodule ℂ E) :
    Kᗮ.starProjection = ContinuousLinearMap.id ℂ E - K.starProjection := by
  ext v
  have := K.starProjection_add_starProjection_orthogonal v
  simp only [_root_.sub_apply, ContinuousLinearMap.id_apply]
  rw [eq_sub_iff_add_eq, add_comm, this]

/-- **Proposition 7.1 (quantum deductive certainty).** If `P_φ ≤ P_ψ` and the unit vector `v`
lies in the range of `P_φ`, then measuring `{P_ψ, I - P_ψ}` on `v` gives `P_ψ` with probability
`⟨v, P_ψ v⟩ = 1`; moreover `I - P_ψ ≤ I - P_φ`. -/
theorem quantum_deduction_sound {K L : Submodule ℂ E} (hKL : K ≤ L) {v : E} (hv : v ∈ K)
    (hnorm : ‖v‖ = 1) : ⟪v, L.starProjection v⟫_ℂ = 1 ∧ Lᗮ ≤ Kᗮ := by
  refine ⟨?_, Submodule.orthogonal_le hKL⟩
  rw [(Submodule.starProjection_eq_self_iff).2 (hKL hv), inner_self_eq_norm_sq_to_K, hnorm]
  simp

end Projection

/-! ## Born probabilities for matrices -/

section Born

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- A projector: a Hermitian idempotent matrix. -/
def IsProjector (P : Matrix n n ℂ) : Prop := P.IsHermitian ∧ P * P = P

/-- A density operator: positive semidefinite with trace `1`. -/
def IsDensity (ρ : Matrix n n ℂ) : Prop := ρ.PosSemidef ∧ ρ.trace = 1

/-- The Born probability `Tr(ρ P)` of the outcome `P` in the state `ρ`. -/
noncomputable def born (ρ P : Matrix n n ℂ) : ℝ := (ρ * P).trace.re

/-- The likelihood `L(ρ) = ∏ₖ Tr(ρ P_{n_k})` of a record of outcomes. -/
noncomputable def likelihood (ρ : Matrix n n ℂ) (outcomes : List (Matrix n n ℂ)) : ℝ :=
  (outcomes.map (born ρ)).prod

theorem IsProjector.one_sub {P : Matrix n n ℂ} (hP : IsProjector P) : IsProjector (1 - P) := by
  refine ⟨?_, ?_⟩
  · simp only [IsHermitian, conjTranspose_sub, conjTranspose_one, hP.1.eq]
  · rw [sub_mul, mul_sub, mul_sub, hP.2]
    simp

omit [DecidableEq n] in
/-- `Tr(ρ P) ≥ 0` for a positive semidefinite `ρ` and a projector `P`. -/
theorem trace_mul_nonneg {ρ P : Matrix n n ℂ} (hρ : ρ.PosSemidef) (hP : IsProjector P) :
    0 ≤ (ρ * P).trace := by
  have h : (ρ * P).trace = (Pᴴ * ρ * P).trace := by
    rw [hP.1.eq, Matrix.mul_assoc, Matrix.trace_mul_comm P, Matrix.mul_assoc, hP.2]
  rw [h]
  exact (hρ.conjTranspose_mul_mul_same P).trace_nonneg

omit [DecidableEq n] in
theorem born_nonneg {ρ P : Matrix n n ℂ} (hρ : ρ.PosSemidef) (hP : IsProjector P) :
    0 ≤ born ρ P :=
  (Complex.nonneg_iff.1 (trace_mul_nonneg hρ hP)).1

/-- Born probabilities are at most `1`. -/
theorem born_le_one {ρ P : Matrix n n ℂ} (hρ : IsDensity ρ) (hP : IsProjector P) :
    born ρ P ≤ 1 := by
  have h := born_nonneg hρ.1 hP.one_sub
  simp only [born, Matrix.mul_sub, Matrix.mul_one, Matrix.trace_sub, hρ.2,
    Complex.sub_re, Complex.one_re] at h
  unfold born
  linarith

omit [DecidableEq n] in
theorem likelihood_nonneg {ρ : Matrix n n ℂ} (hρ : ρ.PosSemidef) {outcomes : List (Matrix n n ℂ)}
    (h : ∀ P ∈ outcomes, IsProjector P) : 0 ≤ likelihood ρ outcomes := by
  unfold likelihood
  apply List.prod_nonneg
  intro x hx
  obtain ⟨P, hP, rfl⟩ := List.mem_map.1 hx
  exact born_nonneg hρ (h P hP)

theorem likelihood_le_one {ρ : Matrix n n ℂ} (hρ : IsDensity ρ) {outcomes : List (Matrix n n ℂ)}
    (h : ∀ P ∈ outcomes, IsProjector P) : likelihood ρ outcomes ≤ 1 := by
  induction outcomes with
  | nil => simp [likelihood]
  | cons P Ps ih =>
    have hP := h P (by simp)
    have hPs : ∀ Q ∈ Ps, IsProjector Q := fun Q hQ => h Q (by simp [hQ])
    simp only [likelihood, List.map_cons, List.prod_cons] at ih ⊢
    calc born ρ P * likelihood ρ Ps ≤ 1 * 1 :=
          mul_le_mul (born_le_one hρ hP) (ih hPs) (likelihood_nonneg hρ.1 hPs) zero_le_one
      _ = 1 := one_mul 1

/-- A state `ρ` is *supported* by a record of outcomes when its likelihood is positive. -/
def QSupported (ρ : Matrix n n ℂ) (outcomes : List (Matrix n n ℂ)) : Prop :=
  0 < likelihood ρ outcomes

end Born

/-! ## Proposition 7.2: the classical framework is the commuting case -/

section Classical

variable {W : Type} [Fintype W] [DecidableEq W]

open Classical in
/-- `P_φ = ∑_{n : φ n} |n⟩⟨n|`: the diagonal projector onto the worlds where `φ` holds. -/
noncomputable def classicalProj (φ : Thought W) : Matrix W W ℂ :=
  diagonal fun w => if φ w then 1 else 0

open Classical in
/-- `classicalProj φ` is the sum of the rank-one projectors `|n⟩⟨n|` over the worlds `n` where
`φ` holds. -/
theorem classicalProj_eq_sum (φ : Thought W) :
    classicalProj φ = ∑ w with φ w, Matrix.single w w (1 : ℂ) := by
  ext i j
  rw [Matrix.sum_apply]
  by_cases hij : i = j
  · subst hij
    by_cases hi : φ i <;> simp [classicalProj, Matrix.single_apply, hi]
  · rw [classicalProj, diagonal_apply_ne _ hij, Finset.sum_eq_zero]
    intro w _
    rw [Matrix.single_apply, ite_eq_right_of_eq_false _ _ (eq_false ?_)]
    rintro ⟨rfl, rfl⟩
    exact hij rfl

theorem classicalProj_mul (φ ψ : Thought W) :
    classicalProj φ * classicalProj ψ = classicalProj fun w => φ w ∧ ψ w := by
  classical
  simp only [classicalProj, diagonal_mul_diagonal]
  congr 1
  funext w
  by_cases hφ : φ w <;> by_cases hψ : ψ w <;> simp [hφ, hψ]

theorem classicalProj_isProjector (φ : Thought W) : IsProjector (classicalProj φ) := by
  classical
  refine ⟨?_, ?_⟩
  · simp only [IsHermitian, classicalProj, diagonal_conjTranspose]
    congr 1
    funext w
    by_cases h : φ w <;> simp [h]
  · rw [classicalProj_mul]
    simp

/-- **Proposition 7.2 (1).** `P_φ ≤ P_ψ` (in the form `P_ψ P_φ = P_φ`) iff `φ ⊨ ψ`. -/
theorem classicalProj_le_iff (φ ψ : Thought W) :
    classicalProj ψ * classicalProj φ = classicalProj φ ↔ Entails φ ψ := by
  classical
  rw [classicalProj_mul]
  simp only [classicalProj]
  constructor
  · intro h w hφ
    have := congrFun (diagonal_injective h) w
    by_contra hψ
    simp [hφ, hψ] at this
  · intro h
    congr 1
    funext w
    by_cases hφ : φ w
    · simp [hφ, h w hφ]
    · simp [hφ]

omit [Fintype W] in
/-- **Proposition 7.2 (1), negation.** `P_{¬φ} = I - P_φ`. -/
theorem classicalProj_not (φ : Thought W) :
    classicalProj (fun w => ¬φ w) = 1 - classicalProj φ := by
  classical
  simp only [classicalProj, ← diagonal_one, diagonal_sub]
  congr 1
  funext w
  by_cases h : φ w <;> simp [h]

theorem classicalProj_prod (obs : List (Thought W)) :
    (obs.map classicalProj).prod = classicalProj fun w => ∀ o ∈ obs, o w := by
  classical
  induction obs with
  | nil => simp [classicalProj]
  | cons o obs ih =>
    rw [List.map_cons, List.prod_cons, ih, classicalProj_mul]
    simp

/-- **Proposition 7.2 (2).** For a diagonal density `ρ = ∑ pₙ |n⟩⟨n|` with every `pₙ > 0`,
`Tr(ρ P_t P_{o₁} ⋯ P_{o_k}) > 0` iff the study `o₁, …, o_k` supports `t`. -/
theorem classical_trace_pos_iff (p : W → ℝ) (hp : ∀ w, 0 < p w) (obs : List (Thought W))
    (t : Thought W) :
    0 < (diagonal (fun w => (p w : ℂ)) * classicalProj t * (obs.map classicalProj).prod).trace ↔
      Supported obs t := by
  classical
  rw [classicalProj_prod, Matrix.mul_assoc, classicalProj_mul]
  have htr : (diagonal (fun w => (p w : ℂ)) *
      classicalProj (fun w => t w ∧ ∀ o ∈ obs, o w)).trace =
      ((∑ w, if t w ∧ ∀ o ∈ obs, o w then p w else 0 : ℝ) : ℂ) := by
    simp only [classicalProj, diagonal_mul_diagonal, trace_diagonal]
    push_cast
    refine Finset.sum_congr rfl fun w _ => ?_
    split_ifs <;> simp
  rw [htr, Complex.zero_lt_real]
  constructor
  · intro h
    by_contra hns
    have : ∀ w, (if t w ∧ ∀ o ∈ obs, o w then p w else 0) = 0 := fun w => by
      exact ite_eq_right_of_eq_false _ _ (eq_false fun hw => hns ⟨w, hw⟩)
    simp [this] at h
  · rintro ⟨w, hw⟩
    apply Finset.sum_pos' (fun x _ => by split_ifs <;> simp [(hp x).le])
    exact ⟨w, Finset.mem_univ w, by rw [ite_eq_left_of_eq_true _ _ (eq_true hw)]; exact hp w⟩

end Classical

/-! ## Proposition 7.3: induction as state tomography -/

section Tomography

variable {n : Type*} [Fintype n] [DecidableEq n]

omit [DecidableEq n] in
/-- **Proposition 7.3 (1).** A further outcome can only shrink the set of supported states:
`L_{N+1}(ρ) > 0` implies `L_N(ρ) > 0`. -/
theorem QSupported.of_concat {ρ : Matrix n n ℂ} (hρ : ρ.PosSemidef)
    {outcomes : List (Matrix n n ℂ)} {P : Matrix n n ℂ} (hP : IsProjector P)
    (h : QSupported ρ (outcomes ++ [P])) : QSupported ρ outcomes := by
  unfold QSupported likelihood at *
  rw [List.map_append, List.prod_append] at h
  simp only [List.map_cons, List.map_nil, List.prod_cons, List.prod_nil, mul_one] at h
  exact pos_of_mul_pos_left h (born_nonneg hρ hP)

/-- The likelihood of a density operator never increases with a further outcome. -/
theorem likelihood_concat_le {ρ : Matrix n n ℂ} (hρ : IsDensity ρ)
    {outcomes : List (Matrix n n ℂ)} (h : ∀ Q ∈ outcomes, IsProjector Q) {P : Matrix n n ℂ}
    (hP : IsProjector P) : likelihood ρ (outcomes ++ [P]) ≤ likelihood ρ outcomes := by
  have hL := likelihood_nonneg hρ.1 h
  unfold likelihood at *
  rw [List.map_append, List.prod_append]
  simp only [List.map_cons, List.map_nil, List.prod_cons, List.prod_nil, mul_one]
  exact mul_le_of_le_one_right hL (born_le_one hρ hP)

/-- The pure state `|v⟩⟨v|`. -/
noncomputable def pureState (v : n → ℂ) : Matrix n n ℂ := vecMulVec v (star v)

/-- The projector `|k⟩⟨k|` onto a basis vector. -/
def basisProj (k : n) : Matrix n n ℂ := diagonal (Pi.single k 1)

/-- The Born rule in a basis: outcome `k` has probability `|v_k|²`. -/
theorem born_pureState_basisProj (v : n → ℂ) (k : n) :
    born (pureState v) (basisProj k) = Complex.normSq (v k) := by
  simp only [born, pureState, basisProj, trace, diag, mul_diagonal, vecMulVec_apply,
    Pi.single_apply, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true,
    Pi.star_apply, Complex.star_def, Complex.mul_conj, Complex.ofReal_re]

theorem basisProj_isProjector (k : n) : IsProjector (basisProj k) := by
  refine ⟨?_, ?_⟩
  · simp only [IsHermitian, basisProj, diagonal_conjTranspose]
    congr 1
    funext i
    by_cases h : i = k <;> simp [h]
  · simp only [basisProj, diagonal_mul_diagonal]
    congr 1
    funext i
    by_cases h : i = k <;> simp [h]

theorem likelihood_pureState (v : n → ℂ) (outcomes : List n) :
    likelihood (pureState v) (outcomes.map basisProj) =
      (outcomes.map fun k => Complex.normSq (v k)).prod := by
  simp only [likelihood, List.map_map]
  congr 1
  exact List.map_congr_left fun k _ => born_pureState_basisProj v k

omit [DecidableEq n] in
theorem pureState_isDensity (v : n → ℂ) (hv : ∑ k, Complex.normSq (v k) = 1) :
    IsDensity (pureState v) := by
  refine ⟨posSemidef_vecMulVec_self_star v, ?_⟩
  simp only [pureState, trace, diag, vecMulVec_apply, Pi.star_apply, Complex.star_def,
    Complex.mul_conj]
  exact_mod_cast hv

/-- The spin state `|θ, ϕ⟩ = cos θ |↑⟩ + e^{iϕ} sin θ |↓⟩`, with `|↑⟩ = |0⟩` and `|↓⟩ = |1⟩`. -/
noncomputable def spinState (θ ϕ : ℝ) : Fin 2 → ℂ :=
  ![Real.cos θ, Complex.exp (ϕ * Complex.I) * Real.sin θ]

theorem normSq_spinState (θ ϕ : ℝ) (k : Fin 2) :
    Complex.normSq (spinState θ ϕ k) = if k = 0 then Real.cos θ ^ 2 else Real.sin θ ^ 2 := by
  fin_cases k
  · show Complex.normSq (Real.cos θ : ℂ) = if (0 : Fin 2) = 0 then _ else _
    rw [Complex.normSq_ofReal, ite_eq_left_of_eq_true _ _ (eq_true rfl), sq]
  · show Complex.normSq (Complex.exp (ϕ * Complex.I) * (Real.sin θ : ℂ)) =
      if (1 : Fin 2) = 0 then _ else _
    rw [Complex.normSq_mul, Complex.normSq_ofReal, Complex.normSq_eq_norm_sq,
      Complex.norm_exp_ofReal_mul_I, ite_eq_right_of_eq_false _ _ (eq_false (by decide)), one_pow,
      one_mul, sq]

theorem spinState_isDensity (θ ϕ : ℝ) : IsDensity (pureState (spinState θ ϕ)) := by
  apply pureState_isDensity
  simp [Fin.sum_univ_two, normSq_spinState, Real.cos_sq_add_sin_sq]

/-- The basis state `|↑⟩`. -/
def spinUp : Fin 2 → ℂ := Pi.single 0 1

theorem spinUp_eq (ϕ : ℝ) : spinUp = spinState 0 ϕ := by
  funext k
  fin_cases k <;> simp [spinUp, spinState]

/-- **Proposition 7.3 (2), likelihood.** After `N` outcomes `↑`, the state `|θ, ϕ⟩` has
likelihood `cos^{2N} θ`. -/
theorem likelihood_spinState_up (θ ϕ : ℝ) (N : ℕ) :
    likelihood (pureState (spinState θ ϕ)) (List.replicate N (basisProj 0)) =
      Real.cos θ ^ (2 * N) := by
  rw [← List.map_replicate, likelihood_pureState, List.map_replicate, List.prod_replicate,
    normSq_spinState, ite_eq_left_of_eq_true _ _ (eq_true rfl), pow_mul]

/-- **Proposition 7.3 (2), underdetermination.** Every `|θ, ϕ⟩` with `0 ≤ θ < π/2` is supported
by any number of outcomes `↑`. -/
theorem spinState_supported {θ : ℝ} (h₀ : 0 ≤ θ) (h₁ : θ < Real.pi / 2) (ϕ : ℝ) (N : ℕ) :
    QSupported (pureState (spinState θ ϕ)) (List.replicate N (basisProj 0)) := by
  unfold QSupported
  rw [likelihood_spinState_up]
  exact pow_pos (Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], h₁⟩) _

/-- **Proposition 7.3 (2), maximum likelihood.** `|↑⟩` maximizes the likelihood of `N` outcomes
`↑` over all density operators. -/
theorem spinUp_maxLikelihood (N : ℕ) :
    IsDensity (pureState spinUp) ∧
      ∀ ρ : Matrix (Fin 2) (Fin 2) ℂ, IsDensity ρ →
        likelihood ρ (List.replicate N (basisProj 0)) ≤
          likelihood (pureState spinUp) (List.replicate N (basisProj 0)) := by
  refine ⟨spinUp_eq 0 ▸ spinState_isDensity 0 0, fun ρ hρ => ?_⟩
  rw [spinUp_eq 0, likelihood_spinState_up, Real.cos_zero, one_pow]
  exact likelihood_le_one hρ fun P hP => (List.eq_of_mem_replicate hP) ▸ basisProj_isProjector 0

/-- **Proposition 7.3 (3), falsification.** A single outcome `↓` refutes `|↑⟩`. -/
theorem spinUp_refuted (outcomes : List (Fin 2)) (h : 1 ∈ outcomes) :
    likelihood (pureState spinUp) (outcomes.map basisProj) = 0 := by
  rw [likelihood_pureState]
  apply List.prod_eq_zero
  exact List.mem_map.2 ⟨1, h, by simp [spinUp]⟩

/-- **Proposition 7.3 (4), the phase is invisible.** No record of outcomes in the `↑, ↓` basis
distinguishes `|θ, ϕ⟩` from `|θ, ϕ'⟩`. -/
theorem likelihood_spinState_phase (θ ϕ ϕ' : ℝ) (outcomes : List (Fin 2)) :
    likelihood (pureState (spinState θ ϕ)) (outcomes.map basisProj) =
      likelihood (pureState (spinState θ ϕ')) (outcomes.map basisProj) := by
  simp only [likelihood_pureState, normSq_spinState]

end Tomography

end WangLogic
