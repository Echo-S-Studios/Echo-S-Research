# Figures -- *The Generative Content of a Conserved Emptiness*

Figures for `papers/2026-06-generative-emptiness/generative_emptiness.tex`.

The paper contains exactly one `tikzpicture`, reconstructed here by
`code/2026-06-generative-emptiness/make_figures.py`, which extracts the TikZ
source, wraps it in the paper's own preamble (so the color macros `growc`,
`capc`, `openc` and the LaTeX definitions all resolve), and compiles it with
XeLaTeX via the `preview` package.

| File | Paper figure | Caption / label |
|---|---|---|
| `figure1.png`, `figure1.pdf` | **Figure 1** (`\label{fig:charge}`, Sec 3) | "The Z/4Z angle charge." The four directions `{0, pi/2, pi, 3pi/2}` are the roots of the content polynomial `x^4-1`. Realised on-circle: `+-1` (charges 0,2, green). The `+-i` sector (charges 1,3) is realised only *off* the circle, as K's place `+-i*beta` (blue). Every irrational angle (red) is the empty Salem sector. |

Reconstruction: 1 of 1 figures rendered; no figure was unreconstructable.

Regenerate with:

```
py code/2026-06-generative-emptiness/make_figures.py
```
