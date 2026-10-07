import Lake

open Lake DSL

package WangLogic where
  leanOptions := #[⟨`autoImplicit, false⟩]

@[default_target] lean_lib WangLogic
lean_lib ComparatorChallenges where globs := #[`ComparatorChallenges.+]
