import Mathlib.Order.Basic
import WangLogic.Basic

/-!
# Ordered levels: iteration moves upward

The foundation `m : I × H × L × {0,1} → Φ` leaves the set of levels `L` uninterpreted. Here `L`
carries a preorder, read as "observational below theoretical", and a *hierarchy* assigns a thought
to each level so that a higher-level thought entails every lower-level one (a theory entails its
observable consequences).

Given a study, a level is

* *established* when the study entails its thought (it is reached by the straight arrow), and
* *tenable* when the study supports its thought (it may still be held along the circular arrow).

Both sets of levels are down-sets. Under iteration (adding observations) the established levels
can only grow and the tenable levels can only shrink: iteration pushes the floor of established
knowledge upward and the ceiling of tenable theory downward, and for a consistent study the
first stays below the second.
-/

namespace WangLogic

/-- **Definition 6.3.** A hierarchy of thoughts indexed by ordered levels: higher levels entail
lower ones. -/
structure Hierarchy (L W : Type) [Preorder L] where
  base : L → Thought W
  entails_of_le : ∀ {l l' : L}, l ≤ l' → Entails (base l') (base l)

variable {L W : Type} [Preorder L] (Hy : Hierarchy L W)

/-- Level `l` is established by the study: the study entails the level's thought. -/
def Hierarchy.Established (obs : List (Thought W)) (l : L) : Prop :=
  Entails (fun w => ∀ o ∈ obs, o w) (Hy.base l)

/-- Level `l` is tenable given the study: the study supports the level's thought. -/
def Hierarchy.Tenable (obs : List (Thought W)) (l : L) : Prop :=
  Supported obs (Hy.base l)

variable {Hy}

/-- Established levels form a down-set. -/
theorem Hierarchy.Established.of_le {obs : List (Thought W)} {l l' : L} (h : l ≤ l')
    (hl' : Hy.Established obs l') : Hy.Established obs l :=
  Entails.trans hl' (Hy.entails_of_le h)

/-- Tenable levels form a down-set. -/
theorem Hierarchy.Tenable.of_le {obs : List (Thought W)} {l l' : L} (h : l ≤ l')
    (hl' : Hy.Tenable obs l') : Hy.Tenable obs l := by
  obtain ⟨w, hw, hobs⟩ := hl'
  exact ⟨w, Hy.entails_of_le h w hw, hobs⟩

/-- **Theorem 6.4 (iteration moves upward).** Adding observations never loses an established level.
-/
theorem Hierarchy.Established.mono {obs obs' : List (Thought W)}
    (hsub : ∀ o ∈ obs, o ∈ obs') {l : L} (h : Hy.Established obs l) :
    Hy.Established obs' l :=
  fun w hw => h w fun o ho => hw o (hsub o ho)

/-- Adding observations never creates a tenable level (`Supported.mono`). -/
theorem Hierarchy.Tenable.anti {obs obs' : List (Thought W)}
    (hsub : ∀ o ∈ obs, o ∈ obs') {l : L} (h : Hy.Tenable obs' l) : Hy.Tenable obs l :=
  Supported.mono hsub h

/-- For a consistent study, every established level is tenable. -/
theorem Hierarchy.Established.tenable {obs : List (Thought W)}
    (hcons : ∃ w, ∀ o ∈ obs, o w) {l : L} (h : Hy.Established obs l) : Hy.Tenable obs l := by
  obtain ⟨w, hw⟩ := hcons
  exact ⟨w, h w hw, hw⟩

/-- Along an iteration `i ↦ obs i` that only adds observations, the established levels at step
`i` remain established at every later step `j`. -/
theorem Hierarchy.Established.iterate (obs : ℕ → List (Thought W))
    (hgrow : ∀ i, ∀ o ∈ obs i, o ∈ obs (i + 1)) {i j : ℕ} (hij : i ≤ j) {l : L}
    (h : Hy.Established (obs i) l) : Hy.Established (obs j) l := by
  induction hij with
  | refl => exact h
  | step _ ih => exact Hierarchy.Established.mono (hgrow _) ih

end WangLogic
