/-!
# Deductive and inductive logic, formalized

Companion to *A Philosophical-Mathematical Representation of Deductive and Inductive Logic:
A Diagrammatic and Formal Analysis* (H. F. Wang, February 2025), Section 3.

Reading of the paper's notation used throughout:

* `W` is a type of possible worlds (the paper's `∃`, the set of possible existences).
* A thought `ϕ` is a proposition about worlds, `Thought W := W → Prop` (the paper's `Φ`).
* `Entails ϕ ψ` is semantic entailment: every world where `ϕ` holds also has `ψ`. This is the
  paper's set inclusion `Φ_ϕ ⊆ f⁻¹(∃_∃)` with the worlds as the underlying set.
* A study `s∃` is a list of observations, each a thought known to hold in the actual world.
* A theory `tϕ` is *supported* by a study when some world satisfies the study and the theory.
  This is the qualitative content of the paper's condition `P(T_ϕ | S_∃) > 0`.
-/

namespace WangLogic

/-- A thought `ϕ ∈ Φ`: a proposition about possible worlds. -/
abbrev Thought (W : Type) := W → Prop

/-- Semantic entailment `ϕ ⊨ ψ`. -/
def Entails {W : Type} (φ ψ : Thought W) : Prop := ∀ w, φ w → ψ w

/-! ## §3.3 Deduction `ϕ → ∃` -/

/-- The paper's definition of `ϕ → ∃` as `¬ϕ ∨ ∃` agrees with implication (classically). -/
theorem imp_iff_not_or (p q : Prop) : (p → q) ↔ (¬p ∨ q) := by
  constructor
  · intro h
    cases Classical.em p with
    | inl hp => exact Or.inr (h hp)
    | inr hp => exact Or.inl hp
  · intro h hp
    cases h with
    | inl hnp => exact absurd hp hnp
    | inr hq => exact hq

/-- Deduction is certain: if `ϕ ⊨ ∃` and `ϕ` holds in a world, `∃` holds there. -/
theorem deduction_sound {W : Type} {φ ψ : Thought W} (h : Entails φ ψ) {w : W} (hw : φ w) :
    ψ w :=
  h w hw

/-- Deductions chain: `ϕ ⊨ ψ` and `ψ ⊨ χ` give `ϕ ⊨ χ`. -/
theorem Entails.trans {W : Type} {φ ψ χ : Thought W} (h₁ : Entails φ ψ) (h₂ : Entails ψ χ) :
    Entails φ χ :=
  fun w hw => h₂ w (h₁ w hw)

/-- §2.2, deductive example: all swans are white; this is a swan; hence it is white. -/
theorem swan_syllogism {Bird : Type} (Swan White : Bird → Prop)
    (all_white : ∀ b, Swan b → White b) (a : Bird) (ha : Swan a) : White a :=
  all_white a ha

/-! ## Polarity `↑`, `↓` -/

/-- Polarity: `up` affirms a thought, `down` negates it. -/
inductive Polarity where
  | up
  | down
  deriving DecidableEq, Repr

/-- Apply a polarity to a proposition: `p↑ = p`, `p↓ = ¬p`. -/
def Polarity.apply : Polarity → Prop → Prop
  | .up, p => p
  | .down, p => ¬p

/-- Negative flow `ϕ → ∃↓` is modus tollens: if `ϕ ⊨ ψ` and `ψ` fails, so does `ϕ`. -/
theorem modus_tollens {W : Type} {φ ψ : Thought W} (h : Entails φ ψ) {w : W}
    (hψ : Polarity.down.apply (ψ w)) : Polarity.down.apply (φ w) :=
  fun hφ => hψ (h w hφ)

/-! ## §3.4 Induction `s∃ ⟳ tϕ` -/

/-- A theory `t` is supported by a study `obs` when some world satisfies every observation
and the theory. This is the paper's `P(T_ϕ | S_∃) > 0`. -/
def Supported {W : Type} (obs : List (Thought W)) (t : Thought W) : Prop :=
  ∃ w, t w ∧ ∀ o ∈ obs, o w

/-- Iteration `i`: dropping an observation never removes support. Equivalently, each new
observation can only shrink the set of supported theories. -/
theorem Supported.of_cons {W : Type} {o : Thought W} {obs : List (Thought W)} {t : Thought W}
    (h : Supported (o :: obs) t) : Supported obs t := by
  obtain ⟨w, ht, hobs⟩ := h
  exact ⟨w, ht, fun o' ho' => hobs o' (List.mem_cons_of_mem o ho')⟩

/-- Iterative refinement is monotone: a larger study supports fewer theories. -/
theorem Supported.mono {W : Type} {obs obs' : List (Thought W)} {t : Thought W}
    (hsub : ∀ o ∈ obs, o ∈ obs') (h : Supported obs' t) : Supported obs t := by
  obtain ⟨w, ht, hobs⟩ := h
  exact ⟨w, ht, fun o ho => hobs o (hsub o ho)⟩

/-- Popper's falsification: one non-white swan refutes "all swans are white" deductively. -/
theorem falsification {Bird : Type} (Swan White : Bird → Prop) (b : Bird) (hs : Swan b)
    (hw : ¬White b) : ¬∀ x, Swan x → White x :=
  fun h => hw (h b hs)

/-- A world for the swan examples: which birds are white. -/
abbrev Coloring (Bird : Type) := Bird → Bool

/-- The study "each bird in `seen` was observed white". -/
def SeenWhite {Bird : Type} (seen : List Bird) : Thought (Coloring Bird) :=
  fun c => ∀ x ∈ seen, c x = true

/-- The theory "all swans are white". -/
def AllWhite {Bird : Type} : Thought (Coloring Bird) :=
  fun c => ∀ x, c x = true

/-- Hume's problem: if some bird is unobserved, no number of white sightings entails that all
birds are white. Induction is not deduction. -/
theorem induction_not_deduction {Bird : Type} [DecidableEq Bird] (seen : List Bird) (b : Bird)
    (hb : b ∉ seen) : ¬Entails (SeenWhite seen) AllWhite := by
  intro h
  let c : Coloring Bird := fun x => if x = b then false else true
  have hseen : SeenWhite seen c := by
    intro x hx
    have hxb : x ≠ b := fun e => hb (e ▸ hx)
    simp [c, hxb]
  have := h c hseen b
  simp [c] at this

/-- §2.2, inductive example: the theory "all swans are white" stays supported by any number of
white sightings, so it may be held tentatively (`⟳`) though never deduced. -/
theorem allWhite_supported {Bird : Type} (seen : List Bird) :
    Supported [SeenWhite seen] AllWhite :=
  ⟨fun _ => true, fun _ => rfl, fun o ho => by
    simp only [List.mem_singleton] at ho
    subst ho
    intro _ _
    rfl⟩

/-! ## §3.5 Foundation `m : I × H × L × {0,1} → Φ` -/

/-- The foundation map: binary state `1` (`true`) gives the thought, `0` gives its negation. -/
def m {I H L W : Type} (base : I → H → L → Thought W) (i : I) (h : H) (l : L) :
    Bool → Thought W
  | true => base i h l
  | false => fun w => ¬base i h l w

/-- The binary state of `m` records the truth value: `m … b` holds exactly when the base
thought's truth agrees with `b`. -/
theorem m_spec {I H L W : Type} (base : I → H → L → Thought W) (i : I) (h : H) (l : L)
    (b : Bool) (w : W) : m base i h l b w ↔ (base i h l w ↔ b = true) := by
  cases b with
  | true => simp [m]
  | false => simp [m]

end WangLogic
