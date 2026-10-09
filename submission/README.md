# Submission package

| File | Use |
|---|---|
| [`../paper.pdf`](../paper.pdf) | Submission PDF (21 pages). Rename it on upload if the venue asks, e.g. `wang-deductive-inductive-logic.pdf`. |
| [`arxiv-source.tar.gz`](arxiv-source.tar.gz) | arXiv source: `paper.tex`, `preamble.tex`, `sections/`, and the pre-built `paper.bbl`. Compiles with `pdflatex` alone. |

To rebuild both from [`../paper`](../paper):

```sh
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
tar czf ../submission/arxiv-source.tar.gz paper.tex preamble.tex paper.bbl sections
cp paper.pdf ../paper.pdf
```

## arXiv metadata

- **Title:** Deduction and Induction in One Diagram: A Machine-Checked Semantics with Applications to Quantum Information and Plant Disease Diagnosis
- **Authors:** Hsuan Fu Wang
- **Primary category:** math.LO (Logic). **Cross-lists:** cs.LO (Logic in Computer Science),
  quant-ph (Quantum Physics).
- **MSC class:** 03A05, 03B48, 68T27, 68V20, 81P10
- **Comments:** 21 pages, 1 figure, 3 tables. All results are formalized in Lean 4 with Mathlib;
  the code is at https://github.com/pupubear007/deductive-inductive-logic
- **License:** CC BY 4.0, matching [`../LICENSE-paper`](../LICENSE-paper).

**Abstract** (plain text for the submission form; 1451 characters, within arXiv's 1,920 limit):

> We give a diagrammatic representation of deductive and inductive reasoning: deduction is a certain
> passage from thought to existence, and induction a tentative, iterated passage from the study of
> existence to the theory of thought. We make the diagram precise in a possible-worlds semantics, in
> which deduction is semantic entailment and a theory is supported by a study when some world is
> consistent with both. We prove that deduction is sound and transitive, that iterated observation
> can only narrow the supported theories, that a single counterexample falsifies a universal theory,
> and that no finite study entails a universal theory over an unobserved individual, which is the
> logical part of Hume's problem of induction. A Bayesian refinement shows that every genuine test a
> theory predicts confirms it; a ranked semantics treats defeasible inference; and ordered levels
> show iteration moving upward. We apply the framework to plant disease diagnosis and forecasting:
> diagnosis is elimination and refinement, Koch's postulates are hypothetico-deductive confirmation,
> forecasting rules are defeasible, and a new theorem characterizes when an assay can resolve a
> diagnosis. A quantum information model realizes the same structure, with the classical setting as
> the case of commuting observables and the single-basis limit of state tomography as an instance of
> assay resolution. Every theorem is verified in the Lean 4 proof assistant with Mathlib.

## Before submitting

- [x] Code in a dedicated public repository, cited in the paper's Code availability statement.
- [x] Licences: Apache-2.0 for the Lean code, CC BY 4.0 for the paper.
- [x] CI builds the Lean project, checks axioms and `sorry`, and builds the PDF on every push.
- [ ] **Archive a release on Zenodo.** Sign in at zenodo.org with GitHub, enable this repository
      under *GitHub* settings, then create a GitHub release (e.g. tag `v1.0-submitted`). Zenodo
      mints a DOI from `CITATION.cff`. Add the DOI to the Code availability statement in
      `paper/sections/08-discussion.tex` and rebuild.
- [ ] **Check the AI-use statement** in "Statements and declarations" against your venue's policy.
- [ ] **Spot-check the references** (author spellings, editions, DOIs), especially the editions
      chosen for Aristotle, Descartes, Hume and Kant.
- [ ] **arXiv endorsement.** First-time submitters to math.LO may need an endorsement from an
      established author in that category.
- [ ] **Journal version.** The source uses the standard `article` class with numeric references.
      For a journal, move the sections into its template; the `\lean{...}` names and the
      appendix table carry over unchanged.
