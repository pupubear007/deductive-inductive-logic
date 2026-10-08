# Instructions for Claude Code in `aggr_model/`

Unpublished RNA-seq and phenotype data. Keep these rules; if a request would break one, say so.

1. **Never commit data or results.** `data/` and `runs/` are git-ignored; configs point at the
   user's folders by path. Do not copy count matrices, phenotype tables or reports into tracked files.
2. **The isolate is the unit.** Report only leave-one-isolate-out results, with preprocessing
   (gene filter, scaling, covariate centring) fitted inside each fold. Never split libraries at
   random. Never report training-set fit as evidence.
3. **Every claim gets its null.** Permutation tests shuffle isolates, not libraries, and rerun the
   whole procedure, including the lambda choice. State the floor (1/720 for 6 isolates).
4. **Colonization stays in.** Keep pathogen transcript share as an unpenalised covariate, or use
   thinning. Do not interpret expression differences that colonization could explain.
5. **No model that can memorise isolates** (deep nets, boosted trees) unless the number of
   isolates makes it testable, and then only under rules 2–3.
6. **Decisions before data.** Stability cut-offs, time points and thresholds are set in the
   config before looking at expression results. Sensitivity checks are reported, not cherry-picked.
7. **No invented biology.** Synthetic parameters are simulation settings. Do not state biological
   values without a source the user has confirmed.

Run `pytest -q` before committing.
