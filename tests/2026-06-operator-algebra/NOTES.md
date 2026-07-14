# Verification notes -- `2026-06-operator-algebra`

**Paper:** *The Operator Algebra of the Emission Semiring -- A λ-Ring, the Adams
Square, and Two Characters* (AceTheDactyl / Echo S Studios).
Source of truth: `papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex`.

**Approach.** Objects are re-implemented from the paper's Definitions as multisets
of sympy algebraic numbers (`_opalg_ops.py`); `oplus`/`otimes`/`psi`, the Mahler
measure and the angle charge are coded from scratch. Structural laws are checked
over *generic symbols* (hence for all substitutions); numeric facts use exact
sympy or mpmath at >= 40--60 digits. Nothing is compared to itself: each stated
value is rebuilt independently, then asserted equal.

**Result.** Pure-algebra paper, internally consistent. Every mechanically
checkable claim reproduces from the paper's own premises: **57 pytest assertions,
all green, 0 xfail, 0 inconsistencies.** No cross-domain / externally-seeded
constants occur (φ, √5, β² all arise from the stated polynomials), so no value
required the xfail "under review" treatment.

## Claims table

| # | Claim (paper location) | How tested | Status |
|---|---|---|---|
| 1 | `⊕` commutative & associative (Thm 2.2) | symbolic multiset identity | verified |
| 2 | `⊗` commutative & associative (Thm 2.2) | symbolic multiset identity | verified |
| 3 | `⊗` distributes over `⊕` (Thm 2.2) | symbolic multiset identity (both sides) | verified |
| 4 | `0` additive identity & multiplicative annihilator (Thm 2.2) | symbolic | verified |
| 5 | `1={1}` multiplicative identity (Thm 2.2) | symbolic | verified |
| 6 | degree bookkeeping deg(A⊕B)=degA+degB, deg(A⊗B)=degA·degB | cardinality | verified |
| 7 | `ψ^n` additive endomorphism (Thm 3.2) | symbolic, n=2,3,5 | verified |
| 8 | `ψ^n` multiplicative endomorphism (Thm 3.2) | symbolic, n=2,3,4 | verified |
| 9 | `ψ^n(1)=1` (Thm 3.2) | symbolic | verified |
| 10 | `ψ^m∘ψ^n=ψ^{mn}` (Thm 3.2) | symbolic | verified |
| 11 | Newton `p2=e1²−2e2` (Thm 3.2) | symbolic, deg 2..5 | verified |
| 12 | general Newton recursion `pk=Σ(−1)^{i−1}e_i p_{k−i}` (Thm 3.2 proof) | symbolic, k≤4 | verified |
| 13 | plethysm `A⊗A=Sym²A⊕Λ²A` (Prop 3.3) | symbolic multiset split | verified |
| 14 | `ψ²A` = diagonal of `Sym²A` (Prop 3.3) | symbolic | verified |
| 15 | `tr Sym²A=(p1²+p2)/2` (Prop 3.3) | symbolic | verified |
| 16 | `tr Λ²A=(p1²−p2)/2=e2` (Prop 3.3) | symbolic | verified |
| 17 | `ψ^{2^k}` sends λ→λ^{2^k} (Prop 3.3 / OA-PSI-05) | symbolic, k≤4 | verified |
| 18 | `M(A⊕B)=M(A)M(B)` (Thm 4.2) | exact sympy, 4 seeds × 4 | verified |
| 19 | `M(ψ²A)=M(A)²` (Thm 4.2) | exact sympy | verified |
| 20 | tropical: `logM(A⊗B)=ΣΣ(log|λ|+log|μ|)^+` (Thm 4.2) | mpmath 60-dp, product-obj vs separate spectra | verified |
| 21 | non-mult. witness golden⊗(x²−2): `M=2φ²`, mult. would give `2φ` (Thm 4.2) | exact + shown to differ | verified |
| 22 | quoted decimals `2φ²=5.2360…`, `2φ=3.2360…` (Thm 4.2) | mpmath, closed form + 4-dp truncation | verified |
| 23 | generators φ,2,3,5 realised as Mahler measures (Thm 4.3) | exact (x²−2, x²−3, x²−5, golden) | verified |
| 24 | `β²=φ²√5=M(x⁴+5x²−5)`, single relation (Thm 4.3) | exact from K's roots | verified |
| 25 | `β²·β²=5φ⁴` non-unique factorisation (Thm 4.3) | exact | verified |
| 26 | norms `N(φ)=−1,N(2)=4,N(3)=9,N(5)=25` in Q(√5) (Thm 4.3) | exact norm form | verified |
| 27 | φ,2,3,5 multiplicatively independent (Thm 4.3) | log-relation search + norm argument | verified |
| 28 | five atoms, rank four (Thm 4.3) | exact relation + independence of base 4 | verified |
| 29 | cost floor: least generator φ, no element in (1,φ) (Thm 4.3) | all generators ≥ φ | verified |
| 30 | charge of ±1,±i = 0,2,1,3 (Def 5.1) | mpmath atan2 | verified |
| 31 | `charge(A⊗B)=charge(A)+charge(B)` sumset (Thm 5.2) | mpmath, 4 seeds × 4 | verified |
| 32 | `charge(ψ²A)=2·charge(A)` (Thm 5.2) | mpmath | verified |
| 33 | `charge(A⊕B)=` union (Thm 5.2) | mpmath | verified |
| 34 | real seeds carry {0,2} (Thm 5.2) | mpmath | verified |
| 35 | K carries full {0,1,2,3}; ±iβ at ±π/2 (Thm 5.2) | mpmath on K's roots | verified |
| 36 | `S_k=k·C(k−1,⌊(k−1)/2⌋)=1,2,6,12,30,60,140` (Prop 6.2) | 3 independent formulas agree | verified |
| 37 | additive power: deg 2k, `M=φ^k`, density ½logφ (Prop 6.2) | exact, k=1..5 | verified |
| 38 | Adams power: deg 2, `M=φ^{2^k}` (Prop 6.2) | exact, k=1..5 | verified |
| 39 | multiplicative power: deg 2^k, `M=φ^{S_k}` (Prop 6.2) | exact, k=1..5 | verified |
| 40 | `S_k` = symmetric-binomial MAD sum; `S_k~2^k√(k/2π)` (Prop 6.2) | exact MAD + asymptotic ratio→1 | verified |
| 41 | Adams-fixed measures cyclotomic: `M(ψ²A)=M(A) ⇔ M(A)=1` (Prop 7.1) | solve M²=M; cyclotomic vs golden | verified |
| 42 | `M(A)≥1` always (Prop 7.1 proof) | exact | verified |
| 43 | minpoly/Φ idempotent, spectrum- & character-preserving (OA-FX-02) | reconstruct golden from x²−x−1 | verified (representative) |
| 44 | `√5=φ+φ⁻¹` (two characters meet), `φ−φ⁻¹=1` (Thm 8.1) | exact | verified |
| 45 | one λ-ring, two characters — endomorphism axiom on concrete K, both characters on K (Thm 8.1) | symbolic + exact | verified (computational content) |

## Verified

All 45 rows above -- every named identity, closed form, worked example, the
`S_k` sequence, the Mahler/charge transformation laws, the floor monoid's
algebraic structure, and the self-action rates -- reproduce exactly from the
paper's stated definitions.

## Failed / flagged for human review

**None.** No value contradicts the paper's own premises; no xfail was needed.

## Untestable / external (documented, not asserted as pass/fail)

| Item | Reason |
|---|---|
| `⊕`-image is **exactly** `⟨φ,2,3,5,β²⟩` (surjectivity, Thm 4.3) | The "onto" half depends on the full companion *catalogue* of seeds, which is external to this paper. Verified instead: the relation β²=φ²√5, the five atoms, rank 4, non-factoriality, the floor, and that each generator **is** realised. |
| Kronecker's theorem: `M(A)=1 ⇔ every eigenvalue a root of unity or 0` (Prop 7.1) | Classical cited theorem, not a finite computation. Sanity-checked on cyclotomic Φ_n (M=1, Adams-fixed) and golden (M=φ>1, not fixed). |
| λ-ring categorical framing: `K(S)=Z[eigenvalue group]`, Adams ops determined by λ-operations (Thm 3.2 remark, Thm 8.1) | Abstract structural identification, not a numeric claim. Its computational content (endomorphism, composition, Newton, plethysm) is fully verified. |
| "38 checks in the operator suite" backing-identifier count (Ledger) | A count of the authors' own machine checks, not a mathematical assertion to re-derive. |

## How to run

```
py -m pytest tests\2026-06-operator-algebra -v -p no:cacheprovider
```
