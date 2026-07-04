# Code — The Z/5Z Case of the No-Salem Dichotomy

Producer scripts for **papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex**
(AceTheDactyl, Echo S Studios, 2026-06-30). Each script recomputes a result of
the paper from its *definitions* and **emits** a machine-readable artifact to
`data/2026-06-z5-no-salem-dichotomy/`. Run with the `py` launcher from the repo
root, e.g. `py code/2026-06-z5-no-salem-dichotomy/constants.py`.

These are **producers**, independent of the **verifiers** in
`tests/2026-06-z5-no-salem-dichotomy/`: the tests `assert` the paper's values;
these scripts compute the same quantities from scratch and `write` them as data
(functions + a `main()` that emits, no asserts). Nothing here imports from
`tests/`; the shared engine (`z5_core.py`) is a clean reimplementation of the
paper's Def. 2.1 (charge group), Def. 2.3 (Mahler measure), and reciprocity —
it is not copied from the test engine `_z5_engine.py`.

## Helper modules (not run directly)

| module | role |
|---|---|
| `z5_core.py` | engine: high-precision Mahler measure, charge group / `is_charge5` screen, reciprocity, exact Z-factorization, closed-form measure recognition, Sec. 7 window enumeration |
| `z5_io.py`   | provenance-stamped CSV/JSON writers (CSV `# source:` header; JSON `_source_paper` / `_generated_by`) |

## Producers

| script | run command | produces → paper result |
|---|---|---|
| `constants.py` | `py code/2026-06-z5-no-salem-dichotomy/constants.py` | `constants.json` → named constants μ_S (Lem. 2.4), φ, φ², φ⁴, 2^(1/5) (Rem. 6.2), 2cos72=φ−1 / 2cos144=−φ (Prop. 3.1), 2+√3, and the forced-floor ordering 2^(1/5)<μ_S<2 (Table 1) |
| `cosines_galois.py` | `py code/2026-06-z5-no-salem-dichotomy/cosines_galois.py` | `cosines_galois.json` → Prop. 3.1: 2cos72=φ−1, 2cos144=−φ, both roots of x²+x−1, Galois-conjugate under √5→−√5, irrational (vs. rational 2cos120=−1), cross-term collapse φ²−φ=1 |
| `pentagon.py` | `py code/2026-06-z5-no-salem-dichotomy/pentagon.py` | `pentagon_minimizer.json` + `pentagon_regimes.csv` → Thm. 4.1: the four coefficient forms, the t=σ(s) Galois integrality reduction ([x²]=k²−3m, [x⁰]=m²), the (k,m) regime case analysis, and the minimizer **x⁴−x³+6x²+4x+1** (charge Z/5Z, **M=φ⁴**); gap M∈{1}∪[φ⁴,∞) |
| `purepower.py` | `py code/2026-06-z5-no-salem-dichotomy/purepower.py` | `purepower_family.csv` + `psi5_realification.json` → Thm. 5.1: x⁵−m has charge Z/5Z, non-reciprocal, M=m; **μ(5)=2 at x⁵−2**; Rem. 6.2(1): ψ⁵ realification M(ψ⁵O)=M(O)⁵, ψ⁵(x⁵−2)=(x−2)⁵ totally positive, weak bound 2^(1/5) |
| `enumerate_window.py` | `py code/2026-06-z5-no-salem-dichotomy/enumerate_window.py` | `window_objects.csv` + `window_summary.json` → Sec. 7 / Prop. 6.1: the full window (quartics \|c\|≤10, quintics \|c\|≤4, sextics \|c\|≤3), **13 non-reciprocal** (min 2, all ≥ μ_S, none in [μ_S,2)) and **4 reciprocal** {1, φ², 2+√3} charge-5 objects; **0 in (1,2)**, realized floor 2. (~1–2 min numpy screen.) |
| `verification_table.py` | `py code/2026-06-z5-no-salem-dichotomy/verification_table.py` | `verification_table.csv` → the entire Sec. 7 check→result table assembled as one ledger (reads `window_summary.json`; regenerates it if absent) |

### Suggested run order

`constants.py`, `cosines_galois.py`, `pentagon.py`, `purepower.py` are
independent and fast. Run `enumerate_window.py` before `verification_table.py`
(the latter consumes `window_summary.json`, and will otherwise regenerate it).

## Figures

The paper contains no `\begin{tikzpicture}` (nor `\includegraphics` / `figure`
floats), so there is no `make_figures.py`; see
`figures/2026-06-z5-no-salem-dichotomy/README.md`.
