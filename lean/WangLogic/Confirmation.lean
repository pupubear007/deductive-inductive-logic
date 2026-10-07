import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Fintype.Pi
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Tactic.Linarith
import WangLogic.Basic

/-!
# Quantitative support: Bayesian confirmation

Over a finite type of worlds `W`, a *prior* gives every world a strictly positive weight. The
probability of a thought is the normalized weight of the worlds where it holds, and conditional
probability is a ratio of weights. A study `s` *confirms* a theory `t` when `P(t | s) > P(t)`; this
is the quantitative inductive polarity `tϕ↑`, and `P(t | s) < P(t)` is `tϕ↓`.

* `prob_and_pos_iff_supported` proves Remark 3.5: under a prior with full support,
  `P(t ∧ s) > 0` holds exactly when `s ⟳ t`.
* `Confirms.supported`: confirmation implies (qualitative) support.
* `confirms_of_entails`: hypothetico-deductive confirmation. If the theory entails the study
  and the study could have failed, the study confirms the theory.
* `condProb_mono`: iteration raises confirmation. For such a theory, more observations give a
  higher posterior: the "ascending probability" of the 2025 version.
* `condProb_eq_zero_of_refutes`, `Disconfirms.of_refutes`: a refuting study sends the posterior
  to `0`.
* `allWhite_confirmed`, `allWhite_condProb_mono`: the swan example.
-/

open Classical

namespace WangLogic

/-- A prior on a finite set of worlds: a strictly positive weight on every world. -/
structure Prior (W : Type) [Fintype W] where
  weight : W → ℝ
  pos : ∀ w, 0 < weight w

namespace Prior

variable {W : Type} [Fintype W] (P : Prior W)

/-- The (unnormalized) weight of the worlds where `φ` holds. -/
noncomputable def mass (φ : Thought W) : ℝ := ∑ w, if φ w then P.weight w else 0

/-- The probability `P(φ)`. -/
noncomputable def prob (φ : Thought W) : ℝ := P.mass φ / P.mass fun _ => True

/-- The conditional probability `P(t | s) = P(t ∧ s) / P(s)`. -/
noncomputable def condProb (t s : Thought W) : ℝ := P.mass (fun w => t w ∧ s w) / P.mass s

variable {P}

theorem mass_nonneg (φ : Thought W) : 0 ≤ P.mass φ :=
  Finset.sum_nonneg fun w _ => by split_ifs <;> simp [(P.pos w).le]

theorem mass_pos_iff (φ : Thought W) : 0 < P.mass φ ↔ ∃ w, φ w := by
  constructor
  · intro h
    by_contra hn
    have : P.mass φ = 0 :=
      Finset.sum_eq_zero fun w _ => ite_eq_right_of_eq_false _ _ (eq_false fun hw => hn ⟨w, hw⟩)
    linarith
  · rintro ⟨w, hw⟩
    exact Finset.sum_pos' (fun x _ => by split_ifs <;> simp [(P.pos x).le])
      ⟨w, Finset.mem_univ w, by simp [hw, P.pos w]⟩

/-- Weights add over a case split on `ψ`. -/
theorem mass_split (φ ψ : Thought W) :
    P.mass φ = P.mass (fun w => φ w ∧ ψ w) + P.mass (fun w => φ w ∧ ¬ψ w) := by
  unfold mass
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun w _ => ?_
  by_cases hφ : φ w <;> by_cases hψ : ψ w <;> simp [hφ, hψ]

theorem mass_congr {φ ψ : Thought W} (h : ∀ w, φ w ↔ ψ w) : P.mass φ = P.mass ψ := by
  unfold mass
  exact Finset.sum_congr rfl fun w _ => by simp only [h w]

theorem mass_mono {φ ψ : Thought W} (h : Entails φ ψ) : P.mass φ ≤ P.mass ψ := by
  unfold mass
  refine Finset.sum_le_sum fun w _ => ?_
  by_cases hφ : φ w
  · simp [hφ, h w hφ]
  · by_cases hψ : ψ w <;> simp [hφ, hψ, (P.pos w).le]

theorem prob_nonneg (φ : Thought W) : 0 ≤ P.prob φ :=
  div_nonneg (mass_nonneg φ) (mass_nonneg _)

/-- **Remark 3.5.** Under a prior with full support, `P(t ∧ s) > 0` iff `s ⟳ t`. -/
theorem prob_and_pos_iff_supported (obs : List (Thought W)) (t : Thought W) :
    0 < P.prob (fun w => t w ∧ ∀ o ∈ obs, o w) ↔ Supported obs t := by
  unfold prob Supported
  constructor
  · intro h
    by_contra hn
    have h0 : P.mass (fun w => t w ∧ ∀ o ∈ obs, o w) = 0 := by
      have := (mass_pos_iff (P := P) (fun w => t w ∧ ∀ o ∈ obs, o w)).not.2 hn
      linarith [mass_nonneg (P := P) (fun w => t w ∧ ∀ o ∈ obs, o w)]
    rw [h0, zero_div] at h
    exact lt_irrefl 0 h
  · rintro ⟨w, hw⟩
    exact div_pos ((mass_pos_iff _).2 ⟨w, hw⟩) ((mass_pos_iff _).2 ⟨w, trivial⟩)

end Prior

variable {W : Type} [Fintype W]

/-- The conjunction of the observations of a study. -/
def StudyHolds (obs : List (Thought W)) : Thought W := fun w => ∀ o ∈ obs, o w

/-- **Definition 5.6.** Inductive polarity `↑`: the study confirms the theory, `P(t | s) > P(t)`. -/
def Confirms (P : Prior W) (obs : List (Thought W)) (t : Thought W) : Prop :=
  P.prob t < P.condProb t (StudyHolds obs)

/-- Inductive polarity `↓`: the study disconfirms the theory, `P(t | s) < P(t)`. -/
def Disconfirms (P : Prior W) (obs : List (Thought W)) (t : Thought W) : Prop :=
  P.condProb t (StudyHolds obs) < P.prob t

open Prior

/-- **Proposition 5.7.** Confirmation implies support: quantitative `↑` refines qualitative `⟳`. -/
theorem Confirms.supported {P : Prior W} {obs : List (Thought W)} {t : Thought W}
    (h : Confirms P obs t) : Supported obs t := by
  have hpos : 0 < P.condProb t (StudyHolds obs) := lt_of_le_of_lt (prob_nonneg t) h
  unfold condProb at hpos
  have hm : 0 < P.mass (fun w => t w ∧ StudyHolds obs w) := by
    by_contra hn
    have h0 : P.mass (fun w => t w ∧ StudyHolds obs w) = 0 :=
      le_antisymm (not_lt.1 hn) (mass_nonneg _)
    rw [h0, zero_div] at hpos
    exact lt_irrefl 0 hpos
  obtain ⟨w, ht, hs⟩ := (mass_pos_iff _).1 hm
  exact ⟨w, ht, hs⟩

/-- When `t ⊨ s`, the posterior is `P(t | s) = mass t / mass s`. -/
theorem condProb_eq_of_entails {P : Prior W} {t s : Thought W} (h : Entails t s) :
    P.condProb t s = P.mass t / P.mass s := by
  unfold condProb
  rw [mass_congr fun w => ⟨fun hw => hw.1, fun hw => ⟨hw, h w hw⟩⟩]

/-- **Theorem 5.8 (hypothetico-deductive confirmation).** If the theory is possible, entails the
study, and the study could have failed, then the study confirms the theory. -/
theorem confirms_of_entails (P : Prior W) {obs : List (Thought W)} {t : Thought W}
    (hent : Entails t (StudyHolds obs)) (ht : ∃ w, t w) (hfail : ∃ w, ¬StudyHolds obs w) :
    Confirms P obs t := by
  unfold Confirms prob
  rw [condProb_eq_of_entails hent]
  have hmt : 0 < P.mass t := (mass_pos_iff t).2 ht
  have hms : 0 < P.mass (StudyHolds obs) := by
    obtain ⟨w, hw⟩ := ht
    exact (mass_pos_iff _).2 ⟨w, hent w hw⟩
  have hsplit := mass_split (P := P) (fun _ => True) (StudyHolds obs)
  have hfailpos : 0 < P.mass (fun w => True ∧ ¬StudyHolds obs w) := by
    obtain ⟨w, hw⟩ := hfail
    exact (mass_pos_iff _).2 ⟨w, trivial, hw⟩
  have hsT : P.mass (fun w => True ∧ StudyHolds obs w) = P.mass (StudyHolds obs) :=
    mass_congr fun w => ⟨fun hw => hw.2, fun hw => ⟨trivial, hw⟩⟩
  apply div_lt_div_of_pos_left hmt hms
  linarith

/-- **Proposition 5.9 (iteration raises confirmation).** If the theory entails the larger study
`s'`, which entails the smaller study `s`, then `P(t | s) ≤ P(t | s')`. -/
theorem condProb_mono {P : Prior W} {t s s' : Thought W} (ht : ∃ w, t w) (h₁ : Entails t s')
    (h₂ : Entails s' s) : P.condProb t s ≤ P.condProb t s' := by
  rw [condProb_eq_of_entails (Entails.trans h₁ h₂), condProb_eq_of_entails h₁]
  have hmt : 0 < P.mass t := (mass_pos_iff t).2 ht
  obtain ⟨w, hw⟩ := ht
  have hms' : 0 < P.mass s' := (mass_pos_iff _).2 ⟨w, h₁ w hw⟩
  exact div_le_div_of_nonneg_left hmt.le hms' (mass_mono h₂)

/-- A study that refutes the theory sends its posterior to `0`. -/
theorem condProb_eq_zero_of_refutes {P : Prior W} {t s : Thought W} (h : Entails s fun w => ¬t w) :
    P.condProb t s = 0 := by
  unfold condProb
  have : P.mass (fun w => t w ∧ s w) = 0 := by
    have hn : ¬∃ w, t w ∧ s w := fun ⟨w, ht, hs⟩ => h w hs ht
    have := (mass_pos_iff (P := P) (fun w => t w ∧ s w)).not.2 hn
    linarith [mass_nonneg (P := P) (fun w => t w ∧ s w)]
  rw [this, zero_div]

/-- **Proposition 5.10 (falsification, quantitatively).** A study refuting a possible theory
disconfirms it. -/
theorem Disconfirms.of_refutes (P : Prior W) {obs : List (Thought W)} {t : Thought W}
    (h : Entails (StudyHolds obs) fun w => ¬t w) (ht : ∃ w, t w) : Disconfirms P obs t := by
  unfold Disconfirms prob
  rw [condProb_eq_zero_of_refutes h]
  obtain ⟨w, hw⟩ := ht
  exact div_pos ((mass_pos_iff _).2 ⟨w, hw⟩) ((mass_pos_iff _).2 ⟨w, trivial⟩)

/-! ## The swan example -/

variable {Bird : Type} [Fintype Bird] [DecidableEq Bird]

omit [Fintype Bird] [DecidableEq Bird] in
theorem allWhite_entails_seenWhite (seen : List Bird) :
    Entails (AllWhite : Thought (Coloring Bird)) (StudyHolds [SeenWhite seen]) := by
  intro c hc o ho
  simp only [List.mem_singleton] at ho
  subst ho
  exact fun x _ => hc x

/-- **Corollary 5.11 (Bayesian swan).** White sightings confirm "all swans are white".** Under any
prior with full support on the colourings, a non-empty study of white birds raises the probability
that all birds are white. -/
theorem allWhite_confirmed (P : Prior (Coloring Bird)) {seen : List Bird} (hseen : seen ≠ []) :
    Confirms P [SeenWhite seen] AllWhite := by
  apply confirms_of_entails P (allWhite_entails_seenWhite seen) ⟨fun _ => true, fun _ => rfl⟩
  obtain ⟨b, hb⟩ := List.exists_mem_of_ne_nil seen hseen
  refine ⟨fun _ => false, fun h => ?_⟩
  have := h (SeenWhite seen) (List.mem_singleton_self _) b hb
  simp at this

/-- **Corollary 5.11, second part.** More white sightings, more confirmation: If `seen ⊆ seen'`, the
posterior of "all birds are white" after `seen'` is at least the posterior after `seen`. -/
theorem allWhite_condProb_mono (P : Prior (Coloring Bird)) {seen seen' : List Bird}
    (hsub : ∀ x ∈ seen, x ∈ seen') :
    P.condProb AllWhite (StudyHolds [SeenWhite seen]) ≤
      P.condProb AllWhite (StudyHolds [SeenWhite seen']) := by
  apply condProb_mono ⟨fun _ => true, fun _ => rfl⟩ (allWhite_entails_seenWhite seen')
  intro c hc o ho
  simp only [List.mem_singleton] at ho
  subst ho
  exact fun x hx => hc (SeenWhite seen') (List.mem_singleton_self _) x (hsub x hx)

end WangLogic
