# Verification notes — *The Dissolved Helix and Its Orthogonal Partner*

**Paper:** `papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex`
**Subtitle:** *Terrain, rotation, and the completion of Z/4Z*
**Author:** AceTheDactyl (Echo S Studios Research Developments)

**Final suite status:** 61 passed, 0 failed, 0 xfailed, 0 skipped, 0 errors
(6 test files; 45 distinct paper claims covered). Run:
`py -m pytest tests/2026-07-helix-orthogonal-partner -v -p no:cacheprovider`.

All derivations below were **re-built independently** from the paper's stated
premises (the Fibonacci companion `R = [[0,1],[1,1]]`; `phi, psi` the roots of
`x^2-x-1`; the quartic `Kf = x^4+5x^2-5`; the registry quartics `cons`, `res`)
using exact `sympy` symbolic algebra, `mpmath` (dps 50-60) for high-precision
numeric cross-checks, and `sympy.galois_group` / `field_isomorphism` for the
number-theoretic claims. No paper number was ever compared to itself; each value
is constructed from scratch and then matched against the paper's stated value.

The paper's own discipline ("no float crosses a decision boundary") makes almost
every displayed constant an **exact** algebraic quantity over Q, Q(sqrt5) or
Q(5^{1/4}). Consequently everything reproduces exactly and **no test needed an
xfail**: there are no cross-domain "seeded" constants that fail to follow from the
paper's premises. The single `[open]`-tagged item (the *provenance* of the labels
`cons`/`res`) is not a mathematical claim and is recorded as untestable; every
mathematical claim *about* those polynomials is derived and verified.

## Claims table

| # | Claim (paper location) | How tested | Verdict |
|---|------------------------|-----------|---------|
| 1 | charpoly(R)=λ²-λ-1, Tr=1, det=-1, spec={φ,ψ} (Prop. 2.2) | build R, symbolic charpoly/trace/det/eigenvals | verified |
| 2 | R²=R+I (Prop. 2.2, Cayley–Hamilton) | matrix algebra | verified |
| 3 | Vieta: φ+ψ=1, φψ=-1 (Prop. 2.2) | symbolic | verified |
| 4 | eq.(1) φ-ψ=√5=φ+φ⁻¹ (Prop. 2.2) | symbolic simplify | verified |
| 5 | ψ=1-φ=-φ⁻¹ (Def. 2.1) | symbolic | verified |
| 6 | eq.(2) spec(ad_R)={0,0,+√5,-√5} (Prop. 2.3) | build 4×4 ad_R, eigenvals + charpoly λ²(λ²-5) | verified |
| 7 | coupling √5=φ+φ⁻¹=φ-ψ irrational (Prop. 2.3) | symbolic | verified |
| 8 | Mah(A_φ)=φ (|ψ|=φ⁻¹<1) (Def. 2.4) | product of max(1,|λ|) | verified |
| 9 | χ(A_φ)={0,2} (Def. 2.4) | exact arg + round(2arg/π) mod 4 | verified |
| 10 | operator laws: ψ² squares Mah, doubles χ mod4 (Def. 2.4) | symbolic on A_φ | verified |
| 11 | eq.(3) Φ_R(R)=R²-R-I=0 (Prop. 3.1) | matrix algebra | verified |
| 12 | [R,R]=0 (Prop. 3.1) | matrix algebra | verified |
| 13 | Mah(A⊕A)=Mah(A)², Mah(ψ²A)=Mah(A)² (Lem. 3.2) | symbolic | verified |
| 14 | S10: Mah(ψ²A)=Mah(A) ⟺ Mah(A)=1 (Lem. 3.2) | two instances (φ≠1; {i,-i}=1) | verified |
| 15 | self-composites Mah∈{φ²,φ⁴,…}, never 1 (Lem. 3.2) | iterate squaring, closed form φ^(2^k) | verified |
| 16 | {0,2} closed under ∪,+mod4,×2; {1,3} unreachable (Prop. 3.3) | exhaustive Z/4Z | verified |
| 17 | on |θ|=1, θ+θ⁻¹=2cosα∈[-2,2] (Prop. 3.4) | symbolic identity + range | verified |
| 18 | θ+θ⁻¹-2=(θ-1)²/θ, >0 for θ>1 (Prop. 3.4) | symbolic identity | verified |
| 19 | eq.(4) τ₀=β+β⁻¹>2>φ, τ₀-2=(β-1)²/β (Prop. 3.4) | symbolic identity + mpmath ordering | verified |
| 20 | Def. 4.1 D=1+4C for x²+x-C (Def. 4.1) | discriminant | verified |
| 21 | roots (-1±√D)/2; double root -1/2 at D=0 (Thm. 4.2) | solve | verified |
| 22 | eq.(5) Gram [[2,-1],[-1,1+2C]], det G=1+4C=D (Thm. 4.2) | field trace over Q(θ), det | verified |
| 23 | signature flip: pos-def D>0, indefinite (1,1) D<0 (Thm. 4.2) | leading minor + eigenvalue signs both sides | verified |
| 24 | rotation channel {-i√|D|,0,+i√|D|} for D<0 (Thm. 4.2) | ad_M spectrum at C=-1 (D=-3) | verified |
| 25 | flip at C=-1/4; golden face C=1 gives D=5=(√5)² (Thm. 4.2) | symbolic | verified |
| 26 | eq.(6) y±=(-5±3√5)/2, straddle 0 (Prop. 5.1) | solve y²+5y-5, signs | verified |
| 27 | eq.(7) K=√y₊=5^{1/4}/φ≈0.9242 (Prop. 5.1) | rationalisation + minimal_polynomial=Kf | verified |
| 28 | β=√\|y₋\|≈2.4195, β>1 (Prop. 5.1) | symbolic + mpmath | verified |
| 29 | eq.(8) Mah(Kf)=β²=(5+3√5)/2=φ²√5≈5.8541 (Prop. 5.1) | product over exact roots + closed form | verified |
| 30 | K<1<β; Kf has 2 real + 2 imaginary roots (Prop. 5.1) | all_roots, re/im split | verified |
| 31 | χ(Kf)={0,1,2,3}; imag roots Re=0 (Prop. 6.1) | exact arg of ±K, ±iβ | verified |
| 32 | Kf even; cons,res non-even (Prop. 6.2) | p(-x)=p(x) test | verified |
| 33 | cons factorisation over Q(√5) (Prop. 6.2) | expand stated product = cons | verified |
| 34 | res factorisation over Q(√5) (Prop. 6.2) | expand stated product = res | verified |
| 35 | cons complex-block Re part = φ²=(3+√5)/2 (Prop. 6.2) | block discriminant sign + b/2 | verified |
| 36 | res complex-block Re part = -φ (Prop. 6.2) | block discriminant sign + -a₁/2 | verified |
| 37 | Kf,cons,res all generate Q(5^{1/4}) (Prop. 6.2) | 5^{1/4}∈Q(K); field_isomorphism for cons,res | verified |
| 38 | signature (2,1) for all three (Prop. 6.2) | count real roots = 2 | verified |
| 39 | Galois group D_4 order 8 for all three (Prop. 6.2) | sympy galois_group order + transitive | verified |
| 40 | none reciprocal ⟹ none Salem (Prop. 6.2) | palindrome test on coeffs | verified |
| 41 | Mah(cons)=2φ⁵=11+5√5 (Prop. 6.2) | mpmath Mahler + symbolic identity | verified |
| 42 | Mah(res)=12+19φ=(43+19√5)/2 (Prop. 6.2) | mpmath Mahler + symbolic identity | verified |
| 43 | filter asymmetry: 1 counterexample refutes, finite range never promotes (Prop. 8.4) | Z/4Z witness + Euler n²+n+41 demotion | verified (sanity) |
| 44 | generator inert on closed sublattice: C_⊣=∅ (Prop. 8.6) | {0,2} closure ∩ {1,3}=∅ | verified (sanity) |
| 45 | kernel identity √5=φ+φ⁻¹ (Rem. 8.9) | symbolic | verified |

## VERIFIED (independent re-derivation matches the paper)

Claims 1–45 above. Highlights of genuinely independent reconstruction:
- **ad_R spectrum** (eq. 2) built as an explicit 4×4 commutator matrix on M₂≅R⁴;
  eigenvalues and characteristic polynomial λ²(λ²-5) recovered from scratch.
- **Gram/flip** (eq. 5): the field trace over Q(θ) computed as a sum over the two
  conjugate roots (-1±√D)/2, giving det G = 1+4C = D symbolically; the signature
  toggle confirmed by eigenvalue signs on both the D>0 (C=1) and D<0 (C=-1) sides.
- **K = 5^{1/4}/φ** confirmed two independent ways: (i) `(5^{1/4}/φ)² = y₊` by
  rationalising `√5/φ² = (3√5-5)/2`, and (ii) `minimal_polynomial(5^{1/4}/φ) = Kf`.
- **Field claim**: `5^{1/4} = K(K²+4)/3 ∈ Q(K)` shown exactly, so Kf generates
  Q(5^{1/4}); `field_isomorphism` confirms cons, res generate isomorphic fields.
- **Galois D_4**: `sympy.galois_group` returns an order-8 transitive group for all
  three quartics (the unique transitive order-8 subgroup of S₄ is D_4).
- **Mahler measures**: recomputed at 60-digit precision as ∏max(1,|root|) over the
  numerically found roots, matched to the exact golden-field closed forms
  φ²√5, 11+5√5, (43+19√5)/2; the identities 2φ⁵=11+5√5 and 12+19φ=(43+19√5)/2
  verified symbolically.

## FLAGGED / xfail

None. Every mathematical claim reproduces exactly from the paper's own premises.

## UNTESTABLE (documented, not mechanically checkable here)

- **Provenance of the labels `cons`/`res`** (Prop. 6.2, `[open]`): the paper itself
  marks the *origin* of these registry polynomials as open ("consolidation /
  resonance quartic" labels in an external `seeds.py`; absent from the emission-
  algebra audit). This is a repository-provenance question, not a derivation; the
  mathematical facts about the polynomials are all verified (claims 32–42).
- **Repository cross-references** (`test_p1_05_keystone.py`, `test_p2_07_uniform.py`,
  `KL_DTA.py`, `polytope.py`, `hull_certificate.py`, and the "suites pass 10/10 and
  7/7" in App. A): assertions about external files not present in this archive.
- **Metalogical framing of the harness** (Defs 8.1–8.3, 8.7; Props 8.4–8.6 as
  general statements; Remarks 8.8–8.10; the ledger Table 2 verdict language): the
  epistemology of ⊢/⊣/? adjudication and the "survivors/casualties/residue"
  partition are interpretive/logical, not arithmetic. Their concrete computational
  kernels (filter asymmetry, generator inertia, √5=φ+φ⁻¹) are checked as sanity
  claims 43–45; the surrounding narrative is a reading, not a computation.
- **The `[posited]` helix overlay** (Table 1; Rems. 5.2, 7.3): the "helix / terrain /
  rotation" language is explicitly an interpretive overlay on the `[forced]`
  substrate. Its algebraic substrate (real vs. imaginary roots, the D=1+4C flip,
  χ∈{0,2} vs. {1,3}) is verified; the metaphor itself is not a checkable statement.
