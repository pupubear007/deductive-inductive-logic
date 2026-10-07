import WangLogic.Basic
import WangLogic.Defeasible

/-!
# Application: plant disease diagnosis and forecasting (Section 8)

Worlds are candidate states of a plant or field (causal agent, host, conditions). A diagnosis is a
thought; an *assay* is a map from worlds to results. Elimination and refinement are already
`modus_tollens` and `Supported.mono`, and Koch's postulates are an instance of
`confirms_of_entails`; this file adds:

* `resolves_iff`: **Theorem 8.2 (assay resolution).** An assay resolves a hypothesis exactly
  when each of its results either entails the hypothesis or rules it out.
* `undetermined`, `Resolves.pair_left`: **Corollary 8.3.** A result shared by two worlds that
  disagree on the hypothesis settles nothing; adding an assay never loses resolution.
* `morphology_not_resolves`, `hostTest_resolves`: **Example 8.4**, formae speciales of
  *Fusarium oxysporum*.
* `positive_raises_risk_iff`: **Proposition 8.6.** In frequency form, a positive test or
  forecast raises the risk of disease exactly when sensitivity exceeds the false-positive rate.
* `apothecia_normally_infected`, `noApothecia_normally_not_infected`,
  `sclerotinia_nonmonotone`: **Proposition 8.7**, a defeasible forecasting rule for
  Sclerotinia stem rot.

Everything here uses only the core library besides `Basic` and `Defeasible`.
-/

namespace WangLogic

/-! ## Assays and resolution -/

/-- The observation "assay `a` gave result `r`". -/
def AssayResult {W R : Type} (a : W → R) (r : R) : Thought W := fun w => a w = r

/-- **Definition 8.1.** Assay `a` *resolves* `H` when worlds with the same result agree on `H`. -/
def Resolves {W R : Type} (a : W → R) (H : Thought W) : Prop :=
  ∀ w w', a w = a w' → (H w ↔ H w')

/-- An observation *decides* `H` when it entails `H` or entails `¬H`. -/
def Decides {W : Type} (o H : Thought W) : Prop :=
  Entails o H ∨ Entails o fun w => ¬H w

/-- **Theorem 8.2 (assay resolution).** An assay resolves `H` exactly when every result it can
give decides `H`. -/
theorem resolves_iff {W R : Type} (a : W → R) (H : Thought W) :
    Resolves a H ↔ ∀ w, Decides (AssayResult a (a w)) H := by
  constructor
  · intro hres w
    cases Classical.em (H w) with
    | inl hH => exact Or.inl fun v hv => (hres v w hv).mpr hH
    | inr hH => exact Or.inr fun v hv hHv => hH ((hres v w hv).mp hHv)
  · intro hdec w w' hww'
    cases hdec w with
    | inl hE => exact ⟨fun _ => hE w' hww'.symm, fun _ => hE w rfl⟩
    | inr hE =>
      exact ⟨fun hHw => absurd hHw (hE w rfl), fun hHw' => absurd hHw' (hE w' hww'.symm)⟩

/-- **Corollary 8.3 (underdetermination by an assay).** If two worlds give the same result but
disagree on `H`, that result neither entails `H` nor rules it out. -/
theorem undetermined {W R : Type} (a : W → R) {H : Thought W} {w w' : W}
    (hsame : a w = a w') (hH : H w) (hnH : ¬H w') :
    ¬Entails (AssayResult a (a w)) H ∧ ¬Entails (AssayResult a (a w)) fun v => ¬H v :=
  ⟨fun hE => hnH (hE w' hsame.symm), fun hE => hE w rfl hH⟩

/-- **Corollary 8.3, second part.** Running a further assay `b` alongside `a` never loses
resolution. -/
theorem Resolves.pair_left {W R S : Type} {a : W → R} (b : W → S) {H : Thought W}
    (h : Resolves a H) : Resolves (fun w => (a w, b w)) H :=
  fun w w' hab => h w w' (congrArg Prod.fst hab)

/-! ### Example 8.4: formae speciales of *Fusarium oxysporum* -/

/-- Two formae speciales, as the worlds of a toy model. -/
inductive FormaSpecialis
  | cubense
  | lycopersici
  deriving DecidableEq

/-- The host on which an isolate causes wilt. -/
inductive Host
  | banana
  | tomato
  deriving DecidableEq

/-- In the model, culture morphology gives the same result for both formae speciales. -/
def morphology : FormaSpecialis → Unit := fun _ => ()

/-- A host test records the host on which the isolate causes wilt. -/
def hostTest : FormaSpecialis → Host
  | .cubense => .banana
  | .lycopersici => .tomato

/-- The diagnosis "the isolate is f. sp. *cubense*". -/
def IsCubense : Thought FormaSpecialis := fun f => f = .cubense

theorem morphology_not_resolves : ¬Resolves morphology IsCubense := by
  intro h
  have := (h .cubense .lycopersici rfl).mp rfl
  cases this

theorem hostTest_resolves : Resolves hostTest IsCubense := by
  intro w w' h
  cases w <;> cases w' <;> simp_all [hostTest, IsCubense]

/-! ## Tests and forecasts -/

/-- **Proposition 8.6 (when a positive result raises risk).** Of `D` diseased and `N` healthy
units, `tp` diseased and `fp` healthy units test (or are forecast) positive. The share of
diseased units among the positives, `tp / (tp + fp)`, exceeds the prevalence `D / (D + N)`
exactly when the sensitivity `tp / D` exceeds the false-positive rate `fp / N`. Both sides are
stated cross-multiplied, so no positivity assumption is needed. -/
theorem positive_raises_risk_iff (D N tp fp : Nat) :
    D * (tp + fp) < tp * (D + N) ↔ D * fp < tp * N := by
  rw [Nat.mul_add, Nat.mul_add, Nat.mul_comm D tp]
  exact Nat.add_lt_add_iff_left

/-! ### Proposition 8.7: a defeasible forecasting rule for Sclerotinia stem rot

Worlds `(noApothecia, infected)` for a soybean field at flowering under weather favourable to
disease. Most normal: apothecia present and the crop infected. Next: no apothecia and no
infection. Other combinations are least normal, not impossible. -/

/-- Worlds `(noApothecia, infected)`. -/
abbrev FieldWorld := Bool × Bool

def fieldRank : FieldWorld → Nat
  | (false, true) => 0
  | (true, false) => 1
  | _ => 2

def NoApothecia : Thought FieldWorld := fun w => w.1 = true
def Infected : Thought FieldWorld := fun w => w.2 = true

theorem apothecia_normally_infected : NormallyFrom fieldRank [] Infected := by
  intro w _ hmin
  have := hmin (false, true) (by simp)
  rcases w with ⟨_ | _, _ | _⟩ <;> simp_all [fieldRank, Infected]

theorem noApothecia_normally_not_infected :
    NormallyFrom fieldRank [NoApothecia] fun w => ¬Infected w := by
  intro w hw hmin
  have := hmin (true, false) (by simp [NoApothecia])
  rcases w with ⟨_ | _, _ | _⟩ <;> simp_all [fieldRank, Infected, NoApothecia]

/-- **Proposition 8.7.** Under favourable weather a field normally becomes infected, but not once
it is known that no apothecia formed: the forecasting rule is defeasible. -/
theorem sclerotinia_nonmonotone :
    NormallyFrom fieldRank [] Infected ∧ ¬NormallyFrom fieldRank [NoApothecia] Infected := by
  refine ⟨apothecia_normally_infected, fun h => ?_⟩
  have h₁ := h (true, false) (by simp [NoApothecia]) fun v hv => by
    rcases v with ⟨_ | _, _ | _⟩ <;> simp_all [fieldRank, NoApothecia]
  simp [Infected] at h₁

end WangLogic
