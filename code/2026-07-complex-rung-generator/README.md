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
| `bl_softcheck.py` | v2.0→v2.1: the active-blockers ledger — one forced fact per blocker class (ODD-2 retired) | 16 + 6 labelled walls |
| `nd_softcheck.py` | v2.1: the ODD-2 closure certificate — relative-norm descent, elementary TP-2 M(N)≥2^r (reads `pop.json`) | 28 · 565 objects |
| `pop.json` | fixture: the 565 forced on-ray objects `nd_softcheck.py` certifies (deterministic) | — |
| `xs_softcheck.py` | v2.2: the cross-shell no-tie certificate — for a multi-shell object, cross-shell coherent ⟺ charge-admissible (conjugation-closure), the P1′ ν-criterion resolved; validated fail-first vs `thm:reltype`, corroborated on the same 855-object window (0 exceptions) | 10 + 3 |
| `ev_softcheck.py` | v2.3–v2.5: the EVEN relational floor M(α) ≥ φ (`thm:evenfloor`) — the even mirror of the odd descent with **Schinzel**'s totally-real bound in place of TP-2; EV-S the Schinzel battery (0 violations, tight at x²−x−1), EV-D the descent certificate, EV-F the parity-split floor scan (even-min φ, odd-min 2), EV-A the p₁=q₂ identity, **EV-P** (v2.5) the N² peel localizing Schinzel to r=1 (r≥2 elementary via TP-2 on N²). Six fail-first falsifiers | 16 + 6 fal. |
| `ax_softcheck.py` | v2.4: C=4/9 is not an emission eigenvalue (`rem:lens49`) — √3 is the unique catalog 3-carrier (w₁=+2), emission monomials keep w₁≥0, and 4/9=2²·3⁻² has w₁=−8<0; the exact prime-3 mirror of `prop:unimon`'s 5-valuation | 4 |
| `decimal_audit.py` | register decimal literals under the 0.51-ulp / truncation matcher | 4 literals |

Session totals across the 18 tx…bl harnesses: **440 exact checks + 67 certified guards + 13 PSLQ corroborations + 6 labelled walls**, plus the v2.1 closure `nd_softcheck.py` (28 checks, 565 objects) that proves ODD-2 and the v2.2 closure `xs_softcheck.py` (10 exact + 3 guards, 855-object window) that proves the P1′ cross-shell no-tie. The ODD-2 closure sources (generator + notes) live in the source workspace's `odd2-closure/`.
Two harnesses the paper cites but does not ship (cite-not-certify): `ozy_softcheck.py` (91),
`relational_softcheck.py` (27). Run one pipe directly: `py code/2026-07-complex-rung-generator/tx_softcheck.py`.
Pinned stack: Python 3.12, sympy 1.14.0, mpmath 1.3.0, numpy 2.4.6 (`nd_softcheck.py` only).
