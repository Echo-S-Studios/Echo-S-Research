# Code — *The Occupant of the Salem Slot*

Producer scripts for **`papers/2026-06-salem-slot/salem_slot.tex`**
("The Occupant of the Salem Slot: the Trace Redirection, the Grow Channel, the
√5 Limit at the Floor, and Its Rate", AceTheDactyl / Echo S Studios).

Each `make_*.py` script computes a result **from the paper's own premises** and
writes machine-readable output to `data/2026-06-salem-slot/`. These are
*producers* (a `main()` that emits data); the suite under
`tests/2026-06-salem-slot/` are the independent *verifiers* (they `assert` the
same numbers). The producers do **not** import the tests, and vice versa.
The single shared map is the trace substitution `trace(θ)=θ+1/θ`.

Run everything from the repo root with the `py` launcher.

## Shared library (not run directly)

| File | Role |
|------|------|
| `salem_core.py` | The trace map, its inverse (the lift `L`), the exact **trace-down** `T` of a reciprocal `P` (power-sum recurrence), Mahler measure, the canonical Salem family `S_n=x^nP−P*`, and the angle charge `A`. |
| `salem_io.py` | Provenance-stamped CSV/JSON writers (`# source:` comment / `_source_paper`+`_generated_by` fields). |

## Producers

| Script | Run command | Paper result produced | Output |
|--------|-------------|-----------------------|--------|
| `make_trace_regions.py` | `py code/2026-06-salem-slot/make_trace_regions.py` | §2 downstairs structure: flip `D=t²−4` (Def 2.1), Lehmer trace-down + round-trip lift, three regions (Lem 2.2), no-fourth-channel (Cor 2.3), lattice lift `{−2,0,2}`→roots of unity, interval endpoints `t=2→x=1`, `t=√5→{φ,φ⁻¹}` (Lem 4.4) | `trace_regions.json` |
| `make_benchmarks.py` | `py code/2026-06-salem-slot/make_benchmarks.py` | §5 benchmark table: Lehmer / β₄ / deg-6 Salems with `τ₀`, trace-down `T`, `Mah(T)`, entropy-trade logs (+ Lehmer numeric §4) | `benchmarks.csv`, `benchmarks.json` |
| `make_redirection.py` | `py code/2026-06-salem-slot/make_redirection.py` | §3/5/6/7 exact identities: AM–GM (Lem 3.1), occupant bound `τ₀>2>φ` (Thm 3.3), entropy trade (Cor 5.1), golden trace `φ+φ⁻¹=√5` (Lem 6.1), **forced** self-action = trace-down (Prop 7.1) | `redirection_identities.json` |
| `make_branch_point.py` | `py code/2026-06-salem-slot/make_branch_point.py` | §4 branch point: √-edge series `β=1+s+s²/2+s³/8`, `s=√(t−2)` (Lem 4.1), quadratic redirection `τ₀−2=(β−1)²/β` (Prop 4.2), monodromy swap `β↔β⁻¹` | `branch_point.json` |
| `make_golden_rate.py` | `py code/2026-06-salem-slot/make_golden_rate.py` | §6 golden limit + rate: `β_n→φ`, `τ₀→√5` (Thm 6.2), linear slope `φ⁻¹` + curvature `√5−2=φ⁻³`, **corrected** Prop 6.4 Taylor (quadratic term `−(√5−2)=−φ⁻³`), geometric-rate table `n=9..27` (Thm 6.6) | `prop64_taylor.json`, `geometric_rate.csv` |
| `make_superselection.py` | `py code/2026-06-salem-slot/make_superselection.py` | §"entry" superselection: angle charge `A` of generators (Prop charge), tensor reducibility `φ⊗φ`, `φ⁴⊗φ⁴`, `(π/2)ℤ` lattice closure (Thm superselect), β₄ off-lattice lift, the six-route entry table (Thm rejects) | `angle_charges.csv`, `tensor_products.json`, `entry_routes.csv` |
| `make_pisot_factorization.py` | `py code/2026-06-salem-slot/make_pisot_factorization.py` | §8 Pisot-trace accumulation (golden→√5, plastic `μ_P`→2.079596) + §9 `S_n` cyclotomic×grow factorizations (`n=6,10,12`), Salem-factor Mahler climbing 1.5061→1.6134 | `pisot_traces.csv`, `sn_factorizations.csv` |
| `make_figures.py` | `py code/2026-06-salem-slot/make_figures.py` | The two inline TikZ figures (renders using the paper's own preamble) | `figures/2026-06-salem-slot/figure{1,2}.{png,pdf}` |

## Reproduce all

```
py code/2026-06-salem-slot/make_trace_regions.py
py code/2026-06-salem-slot/make_benchmarks.py
py code/2026-06-salem-slot/make_redirection.py
py code/2026-06-salem-slot/make_branch_point.py
py code/2026-06-salem-slot/make_golden_rate.py
py code/2026-06-salem-slot/make_superselection.py
py code/2026-06-salem-slot/make_pisot_factorization.py
py code/2026-06-salem-slot/make_figures.py
```

## Note on the Prop 6.4 sign (errata 2026-07-04)

The quadratic coefficient in the expansion `√5 − τ₀ = φ⁻¹u − (√5−2)u² + O(u³)`
(with `u = φ − β`) is **negative**: `−(√5−2) = −φ⁻³`. `make_golden_rate.py`
computes it as `2 − √5` and records the correction. The linear slope `φ⁻¹` and
the curvature *magnitude* `√5−2 = φ⁻³` are unchanged.
