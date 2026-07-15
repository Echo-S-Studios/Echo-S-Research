# code/2026-07-complex-rung-generator/

The verification **pipes** for **papers/2026-07-complex-rung-generator/**: self-contained
exact-arithmetic softcheck harnesses. Each decides every claim over ℚ / ℚ(√5) (extended by
i, K decided at the squared level), with mpmath confined to separately-counted interval
guards or numeric displays, and exits 0 iff all its checks pass. They are *verifiers*, not
data producers — there is no `data/` tree; the drift check simply runs each and requires
exit 0 (`tests/2026-07-complex-rung-generator/test_harnesses.py`).

Every claim in these harnesses passed a two-lane adversarial verification (harness-discipline
audit + independent mathematical refutation) during authoring.

| Pipe | Layer | Checks (exact + guards) |
|------|-------|-------------------------|
| `nx_softcheck.py` | v1.5 third harness: replication, relational type, sphere | 24 + 3 displays |
| `qx_softcheck.py` | v1.6: K-family, quarter twist, rail bounds, Γ-transfer | 62 + 3 |
| `tx_softcheck.py` | v1.7: transcendence pack (κ transcendental given G–S; W[open]4 closed) | 30 + 6 |
| `mx_softcheck.py` | v1.7: multiplier torus, systolic 7th characterization, geodesic ledger | 39 + 3 |
| `rx_softcheck.py` | v1.7: RAD-0 parity no-go, area-transfer law, RAD-4 kills | 25 + 4 |
| `ux_softcheck.py` | v1.7: universal seed ladder (Vein 6), Theorem U, F1 falsifier | 35 + 1 |
| `ex_softcheck.py` | v1.7: incomplete-elliptic rail closed form + non-elementarity | 17 + 5 |
| `cx_softcheck.py` | v1.7: χ-selection chain, pentagon-by-image, monomial no-reach | 34 + 2 |
| `lx_softcheck.py` | v1.7: lens scoping lemma + bounded-negative scan | 26 + 2 (+7 scan) |
| `dx_softcheck.py` | v1.7: k-invariance, axiom H, ledger counter-models, R6 seed | 28 + 6 |
| `rd_softcheck.py` | v1.8: OP-RADIUS apex-ramification certificate + bridge obstruction | 22 + 2 |
| `cl_softcheck.py` | v1.8: OP-RATE R6 classification no-go (κ declared-up-to-ℤ) | 20 + 5 |
| `ch_softcheck.py` | v1.8: χ doctrine — on-circle image of S = μ₂, so ±i ∉ S | 26 + 2 |
| `rn_softcheck.py` | v1.8: relational note — P4 recovered/reduced, odd floor bracketed, P1′ scan | 20 + 4 |
| `pm_softcheck.py` | v1.8: period-relation map (τ₀ algebraic; κ-vs-τ₀ negative; L_res period) | 27 + 7 |
| `ra_softcheck.py` | v1.9: OP-RADIUS reduced to the declared atom D4 + differential strengthening | 23 + 4 |
| `qt_softcheck.py` | v1.9: D2 relocation to the K-seed rotation axis (advance, g1-disjoint) | 19 + 3 |
| `of_softcheck.py` | v1.9: odd relational floor bracketed + reduced to the lemma (ODD-2) | 20 + 6 |
| `pf_softcheck.py` | v1.9: period frontier — forced ledger, Baker/Schanuel map, PSLQ | 13 + 5 (+13 PSLQ) |
| `decimal_audit.py` | register decimal literals under the 0.51-ulp / truncation matcher | 4 literals |

Session totals across the 17 tx…pf harnesses: **424 exact checks + 67 certified guards + 13 PSLQ corroborations**.
Two harnesses the paper cites but does not ship (cite-not-certify): `ozy_softcheck.py` (91),
`relational_softcheck.py` (27). Run one pipe directly: `py code/2026-07-complex-rung-generator/tx_softcheck.py`.
Pinned stack: Python 3.12, sympy 1.14.0, mpmath 1.3.0.
