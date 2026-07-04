# Verification notes — 2026-06-emission-gap

**Paper:** *The Emission–Gap Theorem: No Salem Number in the Spectral Image, and
the Cost Floor It Forces* — A Companion to *The Exchange Rate λ = 2c*
(AceTheDactyl / Echo S Studios). Source: `emission_gap_paper.tex`.

**Verdict:** every mechanically-checkable claim reproduced from the paper's own
premises via independent re-derivation. 47 tests, all green, no xfail, no
inconsistency flagged. The paper's "machine-verified in exact arithmetic"
banner holds up under independent re-computation.

Methods: `sympy` (exact factorization, Sturm root-counting, minimal
polynomials, Newton power sums), `mpmath` at 50 dps (roots, Mahler measures,
transcendental constants), `numpy` (eigenvalues, Kronecker/self-action
operators, trace-form signatures). Helpers live in `emgap_util.py`; every
re-derivation builds the quantity from a defining polynomial/expression *before*
comparing to the paper's printed value.

## Claims table

| # | Claim (paper location) | How tested | Status |
|---|------------------------|------------|--------|
| 1 | μ_S = 1.3247179572 = root of x³−x−1; log μ_S = 0.281200 (Sec.1, Cor.6.1, App.A) | findroot on x³−x−1, mpmath log | verified |
| 2 | φ = 1.618034, log φ = 0.481212; φ²=2.618, φ⁴=(7+3√5)/2 (Cor.6.1, App.A) | mpmath from √5 | verified |
| 3 | Lehmer number 1.1762808 = largest root & Mahler of Lehmer poly, < μ_S (Sec.1, Lem.5.1) | mp_roots + Mahler of degree-10 poly | verified |
| 4 | β₄ = 1.722084 = root of x⁴−x³−x²−x+1, > φ (Cor.10.4) | findroot, compare to φ | verified |
| 5 | Catalog min-polys; seed values are their roots (Def.2.1) | substitute each seed into stated poly | verified |
| 6 | Catalog eigenvalue arguments ∈ (π/2)ℤ (Lem.4.1, App.A) | roots' args mod 90° ≈ 0 | verified |
| 7 | K-seed x⁴+5x²−5 spectrum {±K, ±iβ}, β=2.4195; |i·5^{1/4}|=1.4953 (Lem.4.1, Prop.9.1) | root geometry | verified |
| 8 | gap-seed x²−7x+1 realises Mahler φ⁴ (App.A) | dominant root & Mahler = φ⁴ | verified |
| 9 | Catalog Mahler measures {φ,φ,2,3,5,φ⁴,β²}, min φ, none in (1,μ_S) (Cor.6.1) | Mahler of each seed | verified |
| 10 | (π/2)ℤ closed under + (kron) and doubling (squaring) (Lem.4.2) | ℤ/4 group check | verified |
| 11 | φ⊗φ has on-circle eigenvalue −1 (=φ·ψ) at 180° (App.A) | eigvals of C⊗C | verified |
| 12 | On-circle eigenvalues of spectral products ∈ {1,i,−1,−i} (Lem.4.3) | kron/dsum/square of catalog | verified |
| 13 | Lehmer irreducible w/ root>1; 8 on-circle conj. at 62.81…297.19°, none ⋅90° (Lem.5.1, App.A) | factor_list + root args | verified |
| 14 | Degree-2 Mahler gap (1,φ) empty; min-above-1 = φ over 625 quadratics (Lem.6.1) | brute scan |b|,|c|≤12 | verified |
| 15 | x²−bx+1: Mahler 1 for |b|≤2, ≥φ² for |b|≥3, =φ² at |b|=3 (Lem.6.1) | Mahler by |b| | verified |
| 16 | M(p⊕q)=M(p)M(q); M(φ⊕√2)=2φ=3.236 (Lem.6.2) | Mahler of product poly | verified |
| 17 | Squaring squares Mahler; M(x²−3x+1)=φ²=M(φ)² (Lem.6.2) | Mahler compare | verified |
| 18 | Kronecker sample: no value in (1,φ), min φ² (Lem.6.2) | scan catalog kron spectra | verified |
| 19 | Cyclotomic quadratics x²−x+1, x²+1, x−1 have Mahler 1 (App.A) | Mahler | verified |
| 20 | Lehmer trace-down roots {2.0264,0.9137,−0.5847,−1.4689,−1.8866}: flip-straddle (Lem.8.2) | trace_down + Sturm | verified |
| 21 | β₄ trace-down t²−t−3, flip-straddle; Salem, >φ (Lem.8.2, Cor.10.4) | trace_down exact | verified |
| 22 | Emission delta fires on Salem (Lehmer, β₄), not on Φ₁₀ / non-reciprocal (Thm.8.3) | salem_factors detector | verified |
| 23 | 4th-roots {1,i,−1,−i} trace-downs {2,0,−2}, never interior (Prop.8.4) | ρ(z)=z+1/z | verified |
| 24 | Spectral on-circle trace-downs land in {2,0,−2} (Prop.8.4) | eigvals of kron/square | verified |
| 25 | Traceless quartic x⁴+bx²+1: trace-down t²+(b−2), symmetric roots, never straddle (Prop.11.1) | trace_down for range of b | verified |
| 26 | Every degree-4 Salem poly has trace a≠0 (≠ commutator) (Prop.11.1) | scan Salem quartics | verified |
| 27 | Salem field sig (2,m−1), trace form (m+1,m−1); β₄ (3,1), Lehmer (6,4) (Prop.9.1) | real-root count + Newton trace form | verified |
| 28 | K-seed sig (2,1), trace form (3,1), complex place off-circle 2.4195 (Prop.9.1) | signature + trace form | verified |
| 29 | Totally-real seeds sig (2,0), definite trace form (Prop.9.1) | signature + trace form | verified |
| 30 | Subfield signatures Q(5^{1/4})=(2,1), Q(√2,5^{1/4})=(4,2), Q(√2,√3,√5)=(8,0), K=(8,4) (Thm.10.3) | minimal_polynomial + real-root count | verified |
| 31 | Every non-tot-real subfield of K has shape (2k,k); Salem shape ⇒ deg 4 (Thm.10.3) | shape check on computed sigs | verified |
| 32 | Q(ζ_n) totally complex (0,φ(n)/2), disjoint from Salem sig (Prop.12.3) | cyclotomic poly signature | verified |
| 33 | disc(x²−x−1)=5 (Lem.10.1) | b²−4ac | verified |
| 34 | ad_R spectrum {−√5,0,√5}, off-circle (Lem.10.1, App.A) | kron self-action operator eigvals | verified |
| 35 | ad_M spectrum = difference set {μ_i−μ_j} (Lem.10.2) | generic matrix vs ad operator | verified |
| 36 | 2K min-poly x⁴+20x²−80, sig (2,1), complex mod 2β≈4.84, not Salem (Lem.10.2) | substitution + signature | verified |
| 37 | Size-2 commutator charpoly x²+det, Mahler image {1}∪[2,∞), no (1,2) (Prop.8.5) | scan traceless 2×2 | verified |
| 38 | Circulants (n=4,5,6) emit no Salem (Prop.12.3) | scan 400 random integer circulants | verified |
| 39 | Commutators over catalog & kron/dsum products: 0 Salem factors (Prop.8.6) | build [A,B], factor, detect | verified |
| 40 | Cartan A_n eigenvalues 2−2cos(kπ/(n+1)) ∈ [0,4], totally real, no Salem (Prop.12.4) | eigvalsh of A₃,A₅,A₈ | verified |
| 41 | Entropy h = log M; exp maps entropy gap to height gap (Sec.7 table, LSW) | exp(log φ)=φ | verified |
| 42 | Three distinct gap ends √5=2.23607, φ=1.61803, log φ=0.48121 (Rem.7.2) | mpmath, pairwise separation | verified |

## VERIFIED (all 42 checkable claims above)

Highlights of genuinely independent re-derivations:
- **Trace-down machinery** rebuilt from the Chebyshev recurrence
  (cₖ: c₀=2, c₁=t, c_{k+1}=t·cₖ−c_{k−1}), not from the paper's printed T; it
  reproduces the Lehmer T = t⁵+t⁴−5t³−5t²+4t+3 and its five roots to 3–4 dp.
- **Salem detection** implemented two independent ways (numeric root geometry
  *and* exact palindromy + Sturm flip-straddle) and cross-checked on Lehmer/β₄.
- **Trace-form signatures** computed from exact Newton power sums (Gram matrix
  G_{ij}=Tr(θ^{i+j})), a computation independent of the real-embedding count
  used for the field signature — the two agree with (r₁+r₂, r₂)=(m+1,m−1).
- **Self-action** ad_M built as the flattened operator M⊗I−I⊗Mᵀ; its spectrum
  reproduces the eigenvalue difference set and, for R, {−√5,0,√5}.
- **Degree-16 field** K=Q(√2,√3,5^{1/4}) signature (8,4) obtained from the exact
  minimal polynomial of a primitive element (Sturm root count).

## FAILED / flagged for human review

None. No claim failed to reproduce; no xfail was required.

## UNTESTABLE (documented, not mechanically checkable)

| Item | Reason |
|------|--------|
| Geometry row "det G = 4D" in the Sec.7 non-local table | The Gram matrix G and D are inherited from the companion paper's metric and not defined in this paper; the relation cannot be pinned to a specific object here. The other three rows of the table (height, entropy, channel) are covered by tests 41–42, 34. |
| Uniform bound "at every size" (Thm.10.5), floor forced for all dimensions | An induction/closure statement over unbounded dimension. Its three load-bearing computational facts (difference spectrum, field-confinement to K, signature lattice ⇒ deg-4 ⇒ Mahler ≥ β₄ > φ) are each verified (tests 34–36, 30–31, 4); the "all sizes" quantifier itself is not finitely checkable. |
| Exact scan counts: "400 circulants → 1126 factors", "16 Kronecker products", Prop.8.6 "4 reciprocal factors of deg ≥4" | Sample/enumeration-dependent tallies; the paper does not fully pin the generating slice. The *structural* conclusion each supports (no Salem factor / no value in the gap) is reproduced with independent scans (tests 38, 39, 18). |
| General Lehmer problem (Rem.5.2, ledger last row) | Explicitly stated by the paper as OPEN and outside the system; nothing to verify. |
| Proof-level logical steps: "conjugates travel together" closure (Sec.2), the four-forms equivalence as a single structural identity (Thm.3.2), Shoda's theorem, Kronecker/Smyth/Dobrowolski/LSW citations | External theorems / prose arguments. Sanity-checked where they carry computational content (e.g. Smyth bound μ_S value, LSW h=log M, Salem⇒indefinite signature), but the arguments themselves are not machine claims. |

## Note on cross-domain / seeded constants

The catalog seeds (φ, √2, √3, √5, τ, gap, K) are *stated free parameters* given
by their minimal polynomials (Def.2.1). Per the verification discipline these
were tested for **procedural consistency** — that each named value is a root of
its stated polynomial and yields the stated Mahler measure / argument /
signature — not audited as arithmetic to be "corrected". All such consistency
checks passed, so none needed xfail.
