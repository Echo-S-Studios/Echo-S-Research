# Verification notes — *The Occupant of the Salem Slot*

**Paper:** `papers/2026-06-salem-slot/salem_slot.tex` — "The Occupant of the Salem
Slot: A Positive Characterization: the Trace Redirection, the Grow Channel, the
√5 Limit at the Floor, and Its Rate" (AceTheDactyl / Echo S Studios).

**Suite:** 47 tests across 7 files, all independent re-derivations
(sympy exact symbolic + mpmath ≥40 dps + numpy). Final run: **46 passed,
1 xfailed** (the xfail is a genuine sign discrepancy, see below).

The paper is a chain of small, exact results built on one map — the trace
substitution `trace(θ)=θ+1/θ`. Almost every numeric/symbolic claim reproduces
*exactly* from the paper's own premises. The single defect found is a
second-order **sign typo** in the rate expansion (Prop 6.4); its headline
consequences (linear slope φ⁻¹, geometric rate) are nonetheless correct.

---

## Claims table

| # | Claim (paper location) | How tested | Status |
|---|------------------------|------------|--------|
| 1 | Flip `D(t)=t²−4` is the discriminant of `x²−tx+1` (Def 2.1) | `sympy.discriminant` | verified |
| 2 | Trace-down `P(x)=xᵐ T(x+1/x)`, Lehmer `T=t⁵+t⁴−5t³−5t²+4t+3` (Def 2.1) | substitute T, expand, compare to Lehmer poly | verified |
| 3 | Three regions: on-circle→[−2,2], θ>1→t>2, θ<−1→t<−2 (Lem 2.2) | sample θ=e^{iψ}, θ>1, θ<−1; check t=2cosψ | verified |
| 4 | No fourth channel: `(−2,2)∩(2,∞)=∅` (Cor 2.3) | sympy Interval intersection | verified |
| 5 | Salem trace-down totally real: 1 root >2, m−1 in (−2,2) (Cor 2.3) | nroots of Lehmer T, count/locate | verified |
| 6 | Lift on lattice {−2,0,2} = roots of unity `(x∓1)²`, `x²+1` (Sec entry) | build `L(T)=x^{deg}T(x+1/x)`, compare | verified |
| 7 | `t=2↦x=1`; `t=√5↦{φ,φ⁻¹}`; `(2,√5]=trace((1,φ])` (Lem 4.4) | solve lifts, check sum √5, product 1 | verified |
| 8 | Redirection injective (trace strictly increasing on β>1) | `d/dβ(β+1/β)=1−1/β²>0` | verified |
| 9 | AM–GM: `β+1/β−2=(√β−1/√β)²` (Lem 3.1) | symbolic simplify | verified |
| 10 | `β+1/β>2`, →2⁺ as β→1⁺ (Lem 3.1) | sympy limit + samples | verified |
| 11 | `τ₀−2=(β−1)²/β` exactly (Prop 4.2) | symbolic simplify | verified |
| 12 | Near floor `τ₀=2+δ²+O(δ³)` (Prop 4.2) | series in δ | verified |
| 13 | Occupant `τ₀>2>φ` for all Salem (Thm 3.3) | symbolic `2−φ>0`, `(β−1)²/β>0` | verified |
| 14 | `Mah(T)≥τ₀` (Thm 3.3) | build Mahler from Lehmer T roots | verified |
| 15 | Entropy trade `log τ₀>log β` (Cor 5.1) | difference `=1/β>0` + numeric | verified |
| 16 | `φ+φ⁻¹=√5` (Lem 6.1) | symbolic | verified |
| 17 | Forced: `φₐ−ψₐ=φₐ+1/φₐ=√(a²+4)`, a=1→√5 (Prop 7.1) | solve `x²−ax−1`, symbolic identities | verified |
| 18 | Square-root edge series `β=1+s+s²/2+s³/8`, `s=√(t−2)` (Lem 4.1) | sympy series | verified |
| 19 | `β−1∼√(t−2)`, `log β∼β−1` (Lem 4.1) | mpmath limits t→2⁺ | verified |
| 20 | Monodromy round t=2 swaps β↔1/β (Sec entry) | two branches at t=2.3, product=1 | verified |
| 21 | Lehmer: β=1.1762808, τ₀=2.026418, Mah(T)=5.615601, logs 0.162→0.706 | roots of Lehmer poly + trace-down | verified |
| 22 | β₄: T=t²−t−3, τ₀=2.302776, β=1.7220838, Mah=3.000000, logs 0.544→0.834 | root of T, lift, cross-check `x⁴−x³−x²−x+1` | verified |
| 23 | deg-6: T=t³−t²−3t+1, τ₀=2.170086, β=1.5061357, Mah=3.214320, logs 0.410→0.775 | root of T, lift, cross-check `x⁶−x⁵−x³−x+1` | verified |
| 24 | Proof-sketch: `P*=1−x−x²`, `P*(φ)=−2φ`, `P′(φ)=2φ−1=√5` (Thm 6.5) | symbolic | verified |
| 25 | Golden limit: βₙ→φ, τ₀→√5 (Thm 6.2) | β₃₀ via S₃₀ root | verified |
| 26 | Linear rate slope `trace′(φ)=1−φ⁻²=(√5−1)/2=φ⁻¹` (Prop 6.4) | symbolic | verified |
| 27 | Curvature magnitude `½ trace″(φ)=√5−2=φ⁻³` (Prop 6.4) | symbolic | verified |
| 28 | **Quadratic term of `√5−τ₀` written `+(√5−2)(φ−β)²`** (Prop 6.4) | independent series coeff = `−(√5−2)` | **xfail (sign error)** |
| 29 | Geometric-rate table (gap, gap·φⁿ, ratio) n=9…27 (Thm 6.6) | βₙ roots of Sₙ, high precision | verified |
| 30 | Asymptotics `gap·φⁿ→2/√5`, ratio→φ (Thm 6.6) | limits from βₙ family | verified |
| 31 | Charge A: `K=x⁴+5x²−5` & real seeds have all roots in (π/2)ℤ (Prop charge) | nroots args vs lattice | verified |
| 32 | `φ⊗φ=(x+1)²(x²−3x+1)`, on-circle root −1 (Sec entry) | pairwise products of `x²−x−1` roots | verified |
| 33 | `φ⁴⊗φ⁴=(x−1)²(x²−47x+1)`, on-circle root +1 (Sec entry) | products of `x²−7x+1` roots (min poly of φ⁴) | verified |
| 34 | Operators preserve (π/2)ℤ: doubling & addition closed (Thm superselect) | lattice arithmetic mod 4 | verified |
| 35 | β₄ off-lattice trace root `(1−√13)/2` lifts to conjugate at 0.726π (Sec entry) | lift, |·|=1, arg/π≈0.7258, off (π/2)ℤ | verified |
| 36 | Plastic `μ_P` (root `x³−x−1`)=1.324718, `μ_P+1/μ_P=2.079596` (Prop 8.1) | mpmath root | verified |
| 37 | `Sₙ=xⁿP−P*` factorizations n=6,10,12 (Sec 9) | sympy factor vs stated factors | verified |
| 38 | Sₙ Salem factor Mahler climbs 1.5061→1.6054→1.6134<φ (Sec 9) | largest root of each Salem factor | verified |
| 39 | Lattice-forcing ⇒ reducible (cyclotomic×grow) — sanity (Prop 9.1) | `is_irreducible` on the split examples | verified (instances) |
| U1 | Salem's theorem (every Pisot a 2-sided limit of Salems) (Prop 8.1) | — | untestable (external theorem) |
| U2 | φ is the least limit point of the Pisot set (Prop 8.1) | — | untestable (structural) |
| U3 | Mahler spectrum discrete/finitely-generated, no accumulation at 1 (Sec entry) | — | untestable (companion paper) |
| U4 | General: on-circle roots at 4th-roots-of-unity ⇒ cyclotomic factor (Prop 9.1) | instances only | untestable (in full generality) |
| U5 | "Slot rejects filling" as a global proof (Thm rejects/superselect) | closure core verified | untestable (argument-level) |

---

## VERIFIED (38 claims)

Claims 1–27 and 29–39 above all reproduce exactly from the paper's own
premises. Highlights of the independence:
- **Trace-down** verified by *reconstructing* Lehmer's degree-10 polynomial from
  the stated degree-5 trace-down via `x⁵T(x+1/x)`.
- **Benchmarks** (Lehmer, β₄, deg-6): τ₀, Mah(T), and logs rebuilt from roots;
  each β additionally cross-checked as a root of its *upstairs* Salem polynomial
  (`x⁴−x³−x²−x+1`, `x⁶−x⁵−x³−x+1`).
- **Geometric-rate table**: every entry (gap, gap·φⁿ, and the consecutive-ratio
  column, which we identified as `gap(n−1)/gap(n)`) matches to ≥5 significant
  figures, from independently computed Salem roots of `Sₙ`.
- **⊗ products**: `φ⊗φ` and `φ⁴⊗φ⁴` rebuilt as the multiset of pairwise products
  of the seed's roots (the min poly of φ⁴ derived via `minimal_polynomial`),
  giving exactly the stated factorizations `(x+1)²(x²−3x+1)` and
  `(x−1)²(x²−47x+1)`.
- **Superselection core**: the (π/2)ℤ closure under angle-doubling (squaring) and
  angle-addition (⊗) verified as lattice arithmetic.

## FLAGGED FOR HUMAN REVIEW (1 claim — xfail)

**Prop 6.4 (Proposition "Linear rate, slope φ⁻¹"), second-order term sign.**
The paper displays
> √5 − τ₀(β) = φ⁻¹(φ−β) **+ (√5−2)**(φ−β)² + O((φ−β)³).

An independent Taylor expansion of `√5 − (β+1/β)` at β=φ in powers of `u=φ−β`
gives coefficient of `u¹` = **+φ⁻¹** (matches the paper) but coefficient of
`u²` = **−φ⁻³ = −(√5−2)** — i.e. the quadratic term should be **negative**.
The stated *curvature magnitude* `½ trace″(φ)=√5−2=φ⁻³` is itself correct (that
is the curvature of `trace`); the sign error arises only when inserting it into
the expansion of `√5−τ₀ = −(τ₀−√5)`, which flips the even-order term.

This is a typo-level slip with **no downstream impact**: the headline claims —
linear approach with slope φ⁻¹, geometric rate `√5−τ₀ ∼ (2/√5)φ⁻ⁿ`, and
consecutive ratio →φ — are all correct and independently verified (tests 26, 27,
29, 30, and `test_linear_rate_expansion_leading_term`). Captured as
`test_golden_rate.py::test_linear_rate_expansion_quadratic_sign_as_written`
(xfail), which asserts the paper's `+` form and therefore fails as expected.

## UNTESTABLE (5 claims)

- **U1** Salem's theorem (every Pisot is a two-sided limit of Salems) — a deep
  external theorem. We verify only the concrete `βₙ→φ` instance.
- **U2** "φ is the least limit point of the Pisot set" — a structural fact about
  the Pisot set, not a finite computation.
- **U3** The Mahler spectrum being a discrete, finitely-generated sub-semigroup of
  [φ,∞) with no accumulation at 1 — deferred to the companion papers.
- **U4** Prop 9.1 in full generality (all-on-circle-roots-at-4th-roots-of-unity ⇒
  divisible by a cyclotomic) — a general algebraic theorem; only specific
  instances (the `Sₙ` splits and the ⊗ reducibility) are checked.
- **U5** Theorem "The slot rejects filling" as a *global* proof — the conceptual
  superselection argument. Its computational core (charge on generators + lattice
  closure of the operators) is verified (tests 31–35); the "no finite word of
  operators reaches A=0" universal statement is argument-level.

---

## Reproduce

```
py -m pytest tests\2026-06-salem-slot -v -p no:cacheprovider
```
Expected: 46 passed, 1 xfailed.

## Correction applied 2026-07-04

Prop 6.4 quadratic term sign corrected +(sqrt5-2)(phi-beta)^2 -> -(sqrt5-2)(phi-beta)^2 (= -phi^-3); test test_linear_rate_expansion_quadratic_sign now asserts the minus coefficient. Linear slope phi^-1 and curvature magnitude sqrt5-2 unchanged.
