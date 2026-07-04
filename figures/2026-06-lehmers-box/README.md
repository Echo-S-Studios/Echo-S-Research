# Figures — *Lehmer's Box*

Reconstructed from the inline TikZ in
**`papers/2026-06-lehmers-box/lehmers_box.tex`**. Built by
`code/2026-06-lehmers-box/make_figures.py`, which extracts each
`tikzpicture` block, wraps it in the paper's own preamble (so every custom
macro and colour resolves identically), compiles with XeLaTeX under the
`preview` package, and rasterises to PNG + single-page PDF.

The paper contains **1** `tikzpicture` environment; **1** figure was
reconstructed.

| File | Paper figure | `\label` | Caption (abridged) |
|---|---|---|---|
| `figure1.png`, `figure1.pdf` | Figure 1 | `fig:box` | Lehmer's Box (shaded columns): the green floor wall at `log phi`; the four blue lattice walls at the on-circle directions `(pi/2)Z`; the forbidden strip `(1,phi)`; the golden seed `phi` and the minimal degree-four Salem `beta_4>phi` on the `theta=0` post; Lehmer's number `L` outside the box on both walls (irrational on-circle angle and measure below `phi`). |

Regenerate:

```sh
py code/2026-06-lehmers-box/make_figures.py
```

All reconstructable: nothing in this paper's figures was left unrendered.
