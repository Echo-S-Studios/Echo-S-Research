# Figures

Figures used in the papers, one folder per paper (`figures/<shortname>/` mirrors `papers/<shortname>/`).

The papers contain **no external image files** — every figure is an inline TikZ/pgfplots picture. Each was reconstructed by reusing that paper's *own* preamble plus the `preview` package to compile the `tikzpicture` standalone with XeLaTeX, then rendering to PNG (200 dpi) and per-figure PDF with `pymupdf`; the generator is `code/<shortname>/make_figures.py`. Each folder's `README.md` maps every `figureN.*` to its paper Figure number. Only **5 of the 14 papers** contain figures; the other 9 are pure-algebra/table papers (their `figures/<shortname>/README.md` records that no figures exist, and their tabular data is emitted under `data/`).

## Index by paper

| Paper (`figures/<shortname>/`) | Figures | Contents |
|---|--:|---|
| `2026-06-generative-emptiness` | 1 | Figure 1 — the Z/4Z angle-charge diagram (`fig:charge`) |
| `2026-06-lehmers-box` | 1 | Figure 1 — the height–angle box (`fig:box`) |
| `2026-06-residual-return-learning` | 3 | Figures 1–3 — two-faces, non-disjoint compositum, return loop |
| `2026-06-salem-slot` | 2 | Figures 1–2 — the trace-line trifurcation (`fig:trace`) and the √-edge branch point (`fig:edge`) |
| `2026-06-vector-substrate` | 4 | Figures 1–4 (TikZ/pgfplots); + `figure3_regen`/`figure4_regen` regenerated from data with matplotlib |
| _the other 9 papers_ | 0 | no `tikzpicture` in the source (README notes this per folder) |

**Totals: 11 figures reconstructed across 5 papers (0 unreconstructable); 24 image files (PNG + per-figure PDF, plus 2 matplotlib regenerations).**
