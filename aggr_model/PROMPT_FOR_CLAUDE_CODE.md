# Prompt for Claude Code

Paste everything below the line into Claude Code on your Mac, in a clone of
`pupubear007/deductive-inductive-logic` (it needs to reach `/Users/andrew20445/RNAseq_paper`).

---

Work on branch `ss-aggressiveness-determinants` (`git checkout ss-aggressiveness-determinants`).
The folder `aggr_model/` contains `ssaggr`, which tests aggressiveness determinants of
*Sclerotinia sclerotiorum* from my dual RNA-seq (soybean Williams 82 and sunflower HA89; 6
isolates; 24/48/96 hpi) and my growth-chamber cut-petiole assays (17 isolates, 4 crops). It
follows my NSF proposal: Objective 1 is stability across hosts, and Objective 2 is fixed program
versus host-responsive regulation.

**Before changing anything**
1. Read `aggr_model/CLAUDE.md` (rules you must keep), `README.md` and `DATA_REQUIREMENTS.md`.
2. `cd aggr_model && pip install -e ".[dev]" && pytest -q`. Report the result.

**Run on my data**
3. Check that every path in `configs/real_local.yaml` exists. **Ask me** before changing any of
   them. Ask me for my hypothesis classes (or cut-offs) and put them in `hypothesis_classes` /
   `stability_cutoffs` before running.
4. `python -m ssaggr.run --config configs/real_local.yaml` (runs 24, 48 and 96 hpi separately,
   then pooled).
5. Summarise `summary_by_timepoint.csv` for me: LOIO R² and permutation p per host and time
   point, fixed versus responsive, cross-host transfer, the number of determinants per class, and
   which assays resolve the class (with `p_chance`). State plainly what is not significant.
6. Cross-check the top determinants against my CH3 candidate tables (Table 5 hubs, Table 6
   effectors / CAZymes, Table S1.1) in `/Users/andrew20445/RNAseq_paper`. List overlaps; do not
   edit those files.

**Then, only if I agree**
7. Export the exact DESeq2 VST from `Sclerotinia_dualRNAseq/results/01_counts_prepared/*_dds.rds`
   (R: `assay(vst(dds))`) and add an option to load it, so the features match the chapters.
8. Run `python -m ssaggr.power --isolates 6 12 24 --seeds 10` and give me one paragraph for the
   proposal's "Potential problems" section on how many isolates Objective 2 needs. Say that the
   effect sizes are simulation settings.

**Rules**
- Keep every rule in `aggr_model/CLAUDE.md`. Never commit data or results.
- Commit on `ss-aggressiveness-determinants` only, in small commits. Do not push to `main`.
- After each step, tell me in a few lines what changed and what you need from me.
