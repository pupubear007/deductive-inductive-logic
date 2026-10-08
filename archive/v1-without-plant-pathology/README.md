# Version 1: without the plant-pathology application

*A Philosophical-Mathematical Representation of Deductive and Inductive Logic: A Diagrammatic and
Formal Analysis*, submission version of October 2026 (17 pages). Sections 1–7 cover the
framework, its extensions and the quantum information model. It does not contain Section 8 (plant
disease diagnosis and forecasting), which the current paper in [`../../paper.pdf`](../../paper.pdf)
adds.

| File | |
|---|---|
| [`paper.pdf`](paper.pdf) | The paper |
| [`paper/`](paper/) | LaTeX source |
| [`arxiv-source.tar.gz`](arxiv-source.tar.gz) | arXiv source bundle (compiles with `pdflatex` alone) |

The author block uses the current affiliation (Department of Plant Pathology, University of
Minnesota). Otherwise the text is identical to commit `6328f85`.

Every theorem of this version is formalized in [`../../lean`](../../lean). The only Lean file this
version does not use is `WangLogic/Plant.lean`. The formalization exactly as it stood for this
version is at commit `6328f85`.

To rebuild:

```sh
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```
