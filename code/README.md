# Code

Scripts and programs that produce the results, data, and figures for the papers — the computational products of the research.

`code/<shortname>/` mirrors `papers/<shortname>/`: each folder holds clean **producer** scripts that recompute that paper's results from its stated premises and **emit** machine-readable output to `data/<shortname>/` (and figures to `figures/<shortname>/`). These are the source-of-truth producers; the independent **verifiers** live in [`tests/`](../tests) and are kept separate — the code does not import the tests and vice versa, so the two are genuinely independent derivations. Each folder has its own `README.md` with the exact run command per script.

Run any producer with the `py` launcher from the repo root, e.g.:

```bash
py code/2026-06-salem-slot/make_golden_rate.py   # -> data/2026-06-salem-slot/geometric_rate.csv
py code/2026-06-salem-slot/make_figures.py       # -> figures/2026-06-salem-slot/figure*.png/pdf
```

Dependencies (installed for `py`): `sympy`, `mpmath`, `numpy`, `matplotlib`, `pymupdf`; figure scripts also use MiKTeX XeLaTeX. Python cache/venv files are ignored by [`.gitignore`](../.gitignore). All artifacts use the corrected errata values (degree-12 census 256/46, z5 reciprocal 4, salem-slot `−φ⁻³`, lehmers-box 27 subfields, `R₁∼−R`).

## Index by paper

| Paper (`code/<shortname>/`) | Scripts | Structure → what it produces |
|---|--:|---|
| `2026-06-charge-measure-coupling` | 10 | `cmc_core`+`cmc_io` + 8 producers → Appendix-A ledger, group structure, parity floor, corrected tensor law, floor bounds, Salem commutator |
| `2026-06-emission-gap` | 10 | `emgap_core` + 9 producers → catalog, angle confinement, degree-2 Mahler gap, Salem=flip-straddle, signatures, self-action |
| `2026-06-generative-emptiness` | 8 | `ge_core` + 6 producers + `make_figures` → the five objects, empty (1,φ) gap, orbit measures, Figure 1 |
| `2026-06-lambda-2c` | 10 | `lambda2c_common` + 9 producers → λ=2c identity, trifurcation, gate ladder, keystone powers, det G=4D flip, Kuramoto |
| `2026-06-lehmers-box` | 9 | `box_core`+`box_io` + 6 producers + `make_figures` → constants, empty strip, 27-subfield census, lattice walls, Figure 1 |
| `2026-06-operator-algebra` | 9 | `opalg_core` + 8 producers → semiring/Adams/plethysm laws, floor monoid, `S_k` rates, fixed points |
| `2026-06-residual-return-learning` | 8 | `rrl_core` + 6 producers + `make_figures` → Gram/projector, witness digest, compositum, Fisher, R₁∼−R, Figures 1–3 |
| `2026-06-salem-slot` | 10 | `salem_core`+`salem_io` + 7 producers + `make_figures` → benchmarks, corrected Prop 6.4 Taylor, rate table, Figures 1–2 |
| `2026-06-vector-substrate` | 10 | `vsub_core` + 9 producers + `make_figures` → regular rep, invariant factors, Grams, Fisher, threshold, Figures 1–4 |
| `2026-06-z5-no-salem-dichotomy` | 8 | `z5_core`+`z5_io` + 6 producers → pentagon minimizer, pure-power μ(5)=2, corrected 4/13 window |
| `2026-07-emission-algebra-primer` | 12 | `eap_core`+`eap_io` + 10 producers → power law, `sl2` bracket table, Trace-Form Duality, partner K |
| `2026-07-helix-orthogonal-partner` | 9 | `helix_core`+`helix_io` + 7 producers → `ad_R` spectrum, K-formation, D₄ registry quartics, harness |
| `2026-07-pisot-residue` | 11 | `pisot_lib` + 10 producers → Rat_p, 3125-quintic census, degree-20 residue, corrected 256/46 Salem census |
| `2026-07-relational-charge` | 15 | `relcharge_core`+`relcharge_io` + 13 producers → contact signatures A–X, descent, parity floor, corrected 256/46 census |

**Totals: 139 producer/library scripts across 14 papers.**
