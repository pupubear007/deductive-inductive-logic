# Version 2: with the plant-pathology application, before the revision pass

*Deduction and Induction in One Diagram: A Machine-Checked Semantics with Applications to Quantum
Information and Plant Disease Diagnosis*, October 2026 (22 pages), exactly as merged in pull
request #1 (commit `210ae35`).

| File | |
|---|---|
| [`paper.pdf`](paper.pdf) | The paper |
| [`paper/`](paper/) | LaTeX source |
| [`arxiv-source.tar.gz`](arxiv-source.tar.gz) | arXiv source bundle (compiles with `pdflatex` alone) |

The current paper in [`../../paper.pdf`](../../paper.pdf) revises this version. It narrows the
Hume and Descartes claims, discusses the tacking problem, replaces the quantum analogies with a
short factual paragraph, rewrites the closing sentence, and shortens the abstract. The theorems and
their Lean proofs are unchanged.

To rebuild:

```sh
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```
