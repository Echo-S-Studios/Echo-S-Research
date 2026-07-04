# Figures — *The Occupant of the Salem Slot*

Rendered from the two inline `tikzpicture` environments of
**`papers/2026-06-salem-slot/salem_slot.tex`** by
`code/2026-06-salem-slot/make_figures.py`, which reuses the paper's own
preamble verbatim (so every macro/color/library — `\trace`, `growc`, `capc`,
`openc`, `arrows.meta`, … — resolves) and previews each picture on a tight page
with XeLaTeX + PyMuPDF.

The paper numbers figures with the article-class **flat** counter (only its
*theorems* are section-numbered), so the two figures are **Figure 1** and
**Figure 2** in source order.

| File | Paper figure | Label | Section | Caption (abridged) |
|------|--------------|-------|---------|--------------------|
| `figure1.png` / `figure1.pdf` | **Figure 1** | `fig:trace` | §3 (The redirection map) | The trifurcation lives on the trace line, cut by the flip `t=±2`. A would-be Salem is redirected by `trace` to `τ₀=β+β⁻¹` in **grow**; its conjugates become **captured** roots in `(−2,2)`. Shown: Lehmer's five trace roots; as `β→φ` the redirection `→√5=φ+φ⁻¹`. |
| `figure2.png` / `figure2.pdf` | **Figure 2** | `fig:edge` | §4 (The action just past the flip) | The grow root `β(t)=½(t+√(t²−4))` on the grow side of the flip. The square-root edge at `t=2` (vertical tangent) is the branch point: `β−1∼√(t−2)`. The interval `(2,√5]` corresponds to real pairs `{β,β⁻¹}` with `β∈(1,φ]`; endpoint `t=√5` lifts to `{φ,φ⁻¹}`. Lehmer's `τ₀=2.026` sits just past the flip. |

## Reconstruction status

Both inline figures reconstructed exactly (2 of 2). Each is a schematic
`tikzpicture` drawn directly in the source with literal coordinates, so the
render is a faithful copy of the paper's own figure. There are no external
image includes, plots-from-data, or unreconstructable figures.

## Regenerate

```
py code/2026-06-salem-slot/make_figures.py
```

Produces `figure1.{png,pdf}` and `figure2.{png,pdf}` at 200 dpi; intermediate
`_figbuild.*` files are cleaned up automatically.
