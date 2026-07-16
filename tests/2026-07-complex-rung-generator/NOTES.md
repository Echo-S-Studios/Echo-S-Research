# NOTES — papers/2026-07-complex-rung-generator/

`test_harnesses.py` runs every `*.py` in `code/2026-07-complex-rung-generator/` as an
isolated subprocess and asserts exit 0. Each is an exact-arithmetic softcheck harness
(the paper's verification "pipes"); exit 0 means all of that harness's checks passed.
Wired to the manifest as `verify: pytest tests/2026-07-complex-rung-generator`, so the
`drift` workflow re-certifies the whole verification layer on every push. The harnesses
are verifiers, not producers — no `data/` outputs, so there is nothing to diff/drift.

| Claim group | How tested | Verified / xfail / untestable | Notes |
|-------|-----------|-------------------------------|-------|
| Generator, chart, K-family, sphere, rail (v1.0–v1.6) | nx/qx (+ cited ozy/relational) | Verified (exit 0) | exact over ℚ(√5)[i], K at K² |
| Rate arithmetic, systole, parity no-go, seed ladder, rail closed form (v1.7) | tx mx rx ux ex cx lx dx | Verified (exit 0) | 234 exact + 29 guards |
| OP-RATE classification, radius apex-ramification, χ doctrine, relational note, period map (v1.8) | rd cl ch rn pm | Verified (exit 0) | 115 exact + 20 guards |
| ODD-2 closure (v2.1) + P1′ cross-shell no-tie (v2.2) | nd_softcheck (565 objects) · xs_softcheck (855-object window) | Verified (exit 0) | forced theorems; 0 exceptions |
| Register decimal literals | decimal_audit | Verified (exit 0) | 0.51-ulp / truncation matcher |

- **Verified:** every shipped harness exits 0 under the pinned stack (sympy 1.14.0, mpmath 1.3.0).
- **Flagged (xfail):** none.
- **Untestable (by design):** the transcendence results are FORCED-GIVEN Gelfond–Schneider (an
  external theorem the harness cannot check — only its algebraic side conditions and branch
  enumeration are certified); OP-RATE/OP-RADIUS/W[open]2 remain open/partial as derivations and are
  logged in `../../papers/2026-07-complex-rung-generator/` companions and the corpus open-problem
  register (not in this repo). Two cited harnesses (ozy_softcheck, relational_softcheck) are
  cite-not-certify and are NOT shipped here.
