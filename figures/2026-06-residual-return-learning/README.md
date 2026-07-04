# Figures — Residual Return: Exact Learning Dynamics and Language over the Vector Substrate

Figures for `papers/2026-06-residual-return-learning/residual_return_learning.tex`.

The paper contains three inline `tikzpicture` environments. They are reconstructed by
`code/2026-06-residual-return-learning/make_figures.py`, which reuses the paper's own
preamble (so every macro and colour resolves) and compiles each picture under the
`preview` package with XeLaTeX, then rasterises to PNG (200 dpi) and single-page PDF.

Run: `py code/2026-06-residual-return-learning/make_figures.py`

| File | Paper figure | `\label` | Caption (short) |
|---|---|---|---|
| `figure1.png` / `figure1.pdf` | Figure 1 | `fig:twofaces` | The two faces of residual return, glued by the golden law x²−x−1: learning (number field K) on the left, language (Cl(2,0)≅M₂(R)) on the right. |
| `figure2.png` / `figure2.pdf` | Figure 2 | `fig:nondisjoint` | The non-disjoint witness Q(√2)(√2+√3)=Q(√2,√3): true compositum degree 4 (e′=2), not the tensor degree 8 that does not exist. |
| `figure3.png` / `figure3.pdf` | Figure 3 | `fig:loop` | The automatic detector-driven loop: detect → auto-gain → propose → gate (effective degree) → confirm → re-home → r=0. |

**Figures in paper: 3. Reconstructed: 3. Unreconstructable: none.**
