# Figures — *The Vector Substrate: Number Fields as Exact Learning Geometry*

Source paper: `papers/2026-06-vector-substrate/vector_substrate.tex`.

The paper has **4 inline figures** (all `tikzpicture`, the last also `pgfplots`).
`code/2026-06-vector-substrate/make_figures.py` extracts each `tikzpicture` in
document order and typesets it with the paper's own preamble (so every macro and
color resolves), producing `figureN.png` / `figureN.pdf`. All 4 are
reconstructed. The two data-driven figures are additionally regenerated from
computed data with matplotlib (`figureN_regen.png`).

Run: `py code/2026-06-vector-substrate/make_figures.py`

| Output | Paper figure | Label | Caption (abridged) | Type |
|---|---|---|---|---|
| `figure1.png` / `.pdf` | Figure 1 (§1) | `fig:triptych` | Three views of K=Q(θ) — Q-vector space, matrix algebra ρ, Euclidean lattice — glued by the single trace-form Gram `G=MᵀM`. | Conceptual diagram (not a data plot) |
| `figure2.png` / `.pdf` | Figure 2 (§3) | `fig:proj` | Projector and residual (Ex. 3.4/7.15) in the Minkowski plane of Q(√5): the line R·σ(φ), σ(√5) projecting to (5/3)σ(φ), residual r with ‖r‖²=5/3. | Geometric diagram (exact coordinates) |
| `figure3.png` / `.pdf` | Figure 3 (§5) | `fig:mahler` | Mahler picture for the plastic number ψ (root of x³−x−1): eigenvalues of ρ(ψ)=C(x³−x−1) with the unit circle; only λ₁=μ_S=1.3247… lies outside. | Data plot (companion spectrum) |
| `figure4.png` / `.pdf` | Figure 4 (§8) | `fig:threshold` | Growth decision in the (cost, gain) plane (log-log): 2√6 (96 vs 6.356) and √7 (56 vs 3.892) GROW; the lattice-aligned tiny residual (0.1) STOPs below the 0.562 floor. | Data plot (pgfplots axis) |

## Regenerated data plots (matplotlib)

| Output | Reconstructs | Regenerated from |
|---|---|---|
| `figure3_regen.png` | Figure 3 (`fig:mahler`) | Eigenvalues of `numpy.roots([1,0,-1,-1])` (companion spectrum of x³−x−1); Mahler = ∏_{|λ|>1}|λ| = 1.324718 = μ_S, one root outside the unit circle (red), two inside. |
| `figure4_regen.png` | Figure 4 (`fig:threshold`) | `data/2026-06-vector-substrate/threshold_decisions.csv` (the GROW/STOP points), with the gain=cost boundary and the 0.56 / 4.5 floors. |

Figures 1 and 2 are conceptual/geometric diagrams (no underlying dataset), so
they are reconstructed only from the paper's TikZ source, not regenerated from
data. All four paper figures are reconstructed; none are unreconstructable.
