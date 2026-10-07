import WangLogic

/-!
Prints the axioms used by every theorem cited in the paper. CI fails if any output line mentions
`sorryAx` or an axiom other than `propext`, `Classical.choice` and `Quot.sound`.

Run with `lake env lean AxiomCheck.lean`.
-/

open WangLogic

-- Sections 3–6 (`Basic.lean`)
#print axioms imp_iff_not_or
#print axioms deduction_sound
#print axioms Entails.trans
#print axioms swan_syllogism
#print axioms modus_tollens
#print axioms Supported.of_cons
#print axioms Supported.mono
#print axioms falsification
#print axioms induction_not_deduction
#print axioms allWhite_supported
#print axioms m_spec

-- Remark 3.5 and Section 5.3 (`Confirmation.lean`)
#print axioms Prior.prob_and_pos_iff_supported
#print axioms Confirms.supported
#print axioms confirms_of_entails
#print axioms condProb_mono
#print axioms Disconfirms.of_refutes
#print axioms allWhite_confirmed
#print axioms allWhite_condProb_mono

-- Section 5.4 (`Defeasible.lean`)
#print axioms Normally.of_entails
#print axioms Normally.right_weakening
#print axioms Normally.and
#print axioms Normally.or
#print axioms Normally.cut
#print axioms Normally.cautious_mono
#print axioms entails_premises_mono
#print axioms normally_nonmonotone

-- Section 6.1 (`Levels.lean`)
#print axioms Hierarchy.Established.of_le
#print axioms Hierarchy.Tenable.of_le
#print axioms Hierarchy.Established.mono
#print axioms Hierarchy.Tenable.anti
#print axioms Hierarchy.Established.tenable
#print axioms Hierarchy.Established.iterate

-- Section 7 (`Quantum.lean`, `QuantumExtensions.lean`)
#print axioms starProjection_comp_eq_iff
#print axioms starProjection_orthogonal_eq_sub
#print axioms quantum_deduction_sound
#print axioms unitary_deduction
#print axioms born_unitary_conj
#print axioms classicalProj_eq_sum
#print axioms classicalProj_le_iff
#print axioms classicalProj_not
#print axioms classical_trace_pos_iff
#print axioms quantum_not_distributive
#print axioms QSupported.of_concat
#print axioms likelihood_concat_le
#print axioms likelihood_spinState_up
#print axioms spinState_supported
#print axioms spinUp_maxLikelihood
#print axioms spinUp_refuted
#print axioms likelihood_spinState_phase
#print axioms vonNeumannEntropy_eq_sum_eigenvalues
#print axioms vonNeumannEntropy_unitary_conj
#print axioms born_depolarize_ge
#print axioms noisy_spinUp_not_refuted
