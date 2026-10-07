namespace WangLogic

abbrev Thought (W : Type) := W → Prop

def Entails {W : Type} (φ ψ : Thought W) : Prop := ∀ w, φ w → ψ w

abbrev Coloring (Bird : Type) := Bird → Bool

def SeenWhite {Bird : Type} (seen : List Bird) : Thought (Coloring Bird) :=
  fun c => ∀ x ∈ seen, c x = true

def AllWhite {Bird : Type} : Thought (Coloring Bird) :=
  fun c => ∀ x, c x = true

theorem induction_not_deduction {Bird : Type} [DecidableEq Bird] (seen : List Bird) (b : Bird)
    (hb : b ∉ seen) : ¬Entails (SeenWhite seen) AllWhite := by
  sorry

end WangLogic
