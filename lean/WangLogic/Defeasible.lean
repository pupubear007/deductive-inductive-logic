import Mathlib.Order.Basic
import Mathlib.Data.Nat.Find
import WangLogic.Basic

/-!
# Multi-premise defeasible reasoning

Defeasible ("normally") inference is given the ranked semantics of Kraus, Lehmann and Magidor:
a rank function `r : W → ℕ` orders worlds by normality (lower is more normal), and premises `Γ`
defeasibly entail `ψ` when `ψ` holds in the most normal worlds satisfying all of `Γ`.

* `Normally.of_entails`: every deduction is a defeasible inference (supraclassicality).
* `Normally.right_weakening`, `Normally.and`, `Normally.or`, `Normally.cut`,
  `Normally.cautious_mono`: the rules of the preferential system **P**.
* `normally_nonmonotone`: unlike `Entails`, adding a premise can retract a conclusion
  (birds normally fly; penguins normally do not).
-/

namespace WangLogic

variable {W : Type}

/-- **Definition 5.12.** `φ |~ ψ` in the ranked model `r`: the `r`-minimal `φ`-worlds satisfy `ψ`.
-/
def Normally (r : W → ℕ) (φ ψ : Thought W) : Prop :=
  ∀ w, φ w → (∀ v, φ v → r w ≤ r v) → ψ w

/-- Multi-premise defeasible entailment `Γ |~ ψ`: the premises are conjoined. -/
def NormallyFrom (r : W → ℕ) (Γ : List (Thought W)) (ψ : Thought W) : Prop :=
  Normally r (fun w => ∀ γ ∈ Γ, γ w) ψ

/-- Every satisfiable thought has an `r`-minimal world. -/
theorem exists_minimal (r : W → ℕ) {φ : Thought W} (h : ∃ w, φ w) :
    ∃ w, φ w ∧ ∀ v, φ v → r w ≤ r v := by
  classical
  obtain ⟨w₀, hw₀⟩ := h
  have hex : ∃ n, ∃ w, φ w ∧ r w = n := ⟨r w₀, w₀, hw₀, rfl⟩
  obtain ⟨w, hw, hr⟩ := Nat.find_spec hex
  exact ⟨w, hw, fun v hv => hr ▸ Nat.find_min' hex ⟨v, hv, rfl⟩⟩

variable {r : W → ℕ} {φ ψ χ : Thought W}

/-- Supraclassicality: deduction is a special case of defeasible inference. -/
theorem Normally.of_entails (h : Entails φ ψ) : Normally r φ ψ :=
  fun w hw _ => h w hw

theorem Normally.refl : Normally r φ φ := fun _ hw _ => hw

theorem Normally.right_weakening (h : Normally r φ ψ) (hψχ : Entails ψ χ) : Normally r φ χ :=
  fun w hw hmin => hψχ w (h w hw hmin)

theorem Normally.and (h₁ : Normally r φ ψ) (h₂ : Normally r φ χ) :
    Normally r φ fun w => ψ w ∧ χ w :=
  fun w hw hmin => ⟨h₁ w hw hmin, h₂ w hw hmin⟩

theorem Normally.or (h₁ : Normally r φ χ) (h₂ : Normally r ψ χ) :
    Normally r (fun w => φ w ∨ ψ w) χ := by
  intro w hw hmin
  rcases hw with hφ | hψ
  · exact h₁ w hφ fun v hv => hmin v (Or.inl hv)
  · exact h₂ w hψ fun v hv => hmin v (Or.inr hv)

theorem Normally.cut (h₁ : Normally r (fun w => φ w ∧ ψ w) χ) (h₂ : Normally r φ ψ) :
    Normally r φ χ :=
  fun w hw hmin => h₁ w ⟨hw, h₂ w hw hmin⟩ fun v hv => hmin v hv.1

/-- **Theorem 5.13**, the non-trivial rule. Cautious monotonicity: adding a premise that was itself
a defeasible conclusion keeps every other conclusion. -/
theorem Normally.cautious_mono (h₁ : Normally r φ ψ) (h₂ : Normally r φ χ) :
    Normally r (fun w => φ w ∧ ψ w) χ := by
  intro w hw hmin
  apply h₂ w hw.1
  intro v hv
  obtain ⟨m, hm, hmmin⟩ := exists_minimal r ⟨v, hv⟩
  exact le_trans (hmin m ⟨hm, h₁ m hm hmmin⟩) (hmmin v hv)

/-- Deductive entailment from premises is monotone: more premises, no fewer conclusions. -/
theorem entails_premises_mono {Γ Γ' : List (Thought W)} (hsub : ∀ γ ∈ Γ, γ ∈ Γ')
    (h : Entails (fun w => ∀ γ ∈ Γ, γ w) ψ) : Entails (fun w => ∀ γ ∈ Γ', γ w) ψ :=
  fun w hw => h w fun γ hγ => hw γ (hsub γ hγ)

/-! ## Tweety: defeasible inference is nonmonotone -/

/-- Worlds `(penguin, flies)` for a single bird. -/
abbrev BirdWorld := Bool × Bool

/-- Normal birds are flying non-penguins; flightless penguins are the next most normal. -/
def birdRank : BirdWorld → ℕ
  | (false, true) => 0
  | (true, false) => 1
  | _ => 2

def Penguin : Thought BirdWorld := fun w => w.1 = true
def Flies : Thought BirdWorld := fun w => w.2 = true

theorem bird_normally_flies : NormallyFrom birdRank [] Flies := by
  intro w _ hmin
  have := hmin (false, true) (by simp)
  rcases w with ⟨_ | _, _ | _⟩ <;> simp_all [birdRank, Flies]

theorem penguin_normally_not_flies : NormallyFrom birdRank [Penguin] fun w => ¬Flies w := by
  intro w hw hmin
  have := hmin (true, false) (by simp [Penguin])
  rcases w with ⟨_ | _, _ | _⟩ <;> simp_all [birdRank, Flies, Penguin]

/-- **Proposition 5.14 (defeasible inference is nonmonotone).** From no premises a bird normally
flies, but from the premise "penguin" it does not. -/
theorem normally_nonmonotone :
    NormallyFrom birdRank [] Flies ∧ ¬NormallyFrom birdRank [Penguin] Flies := by
  refine ⟨bird_normally_flies, fun h => ?_⟩
  have h₁ := h (true, false) (by simp [Penguin]) fun v hv => by
    rcases v with ⟨_ | _, _ | _⟩ <;> simp_all [birdRank, Penguin]
  simp [Flies] at h₁

end WangLogic
