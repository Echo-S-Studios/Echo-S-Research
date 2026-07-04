# Verification notes — *The Z/5Z Case of the No-Salem Dichotomy*

**Paper:** `papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex`
(AceTheDactyl, Echo S Studios, 2026-06-30).
**Subtitle:** "A charge-coupled golden floor that confines the pentagon sector, the
Smyth bound that sharpens the general floor, and why μ(5)=2 stays *computed*."

Every test re-derives a value from the paper's **definitions** (Mahler measure
Def. 2.3, charge group Def. 2.1, reciprocal test) and only then compares to the
paper's stated result. The window enumeration is an independent reproduction of
the Sec. 7 verification table, using a from-scratch charge/measure engine
(`_z5_engine.py`).

Final run: **54 passed, 1 xfailed** (`py -m pytest <folder> -p no:cacheprovider`).

## Claims table

| # | Claim (paper location) | How tested | Status | Notes |
|---|---|---|---|---|
| 1 | μ_S = 1.3247179… is the real root of x³−x−1 (Lem. 2.4) | mpmath root of x³−x−1, compare digits | verified | |
| 2 | μ_S is the plastic number, a Pisot number (Lem. 2.4) | conjugates have \|·\|<1, real root>1 | verified | "smallest" Pisot = Siegel, cited not re-proved |
| 3 | realification value 2^(1/5)=1.1487 (Rem. 6.2) | mpmath 5th root, (·)^5=2 | verified | |
| 4 | φ=(1+√5)/2, φ²−φ−1=0 (Prop. 3.1) | sympy simplify | verified | |
| 5 | φ²=(3+√5)/2≈2.618, φ⁴=(7+3√5)/2=3φ+2≈6.854102 | sympy + mpmath | verified | |
| 6 | forced floor rises 2^(1/5)→μ_S, both <2 (Table 1) | ordering 2^(1/5)<μ_S<2 | verified | |
| 7 | 2cos72°=φ−1 (Prop. 3.1) | exact sympy cos | verified | |
| 8 | 2cos144°=−φ (Prop. 3.1) | exact sympy cos | verified | |
| 9 | {2cos72°,2cos144°} are the roots of x²+x−1 (Prop. 3.1) | sympy solve + annihilation | verified | |
| 10 | the two cosines are Galois-conjugate under √5→−√5 (Prop. 3.1) | σ swaps them; sum,prod rational | verified | |
| 11 | both cosines irrational (Prop. 3.1); 2cos120°=−1 rational (Z/3Z) | sympy is_rational | verified | motivates why Z/3Z proof breaks |
| 12 | cross-term collapse φ²−φ=1 (Thm. 4.1 proof) | sympy | verified | |
| 13 | pentagon lattice {0,72,144,216,288}° omits 180° (Lem. 2.2) | fractions of a turn; 1/2 absent | verified | |
| 14 | on-unit charge-5 roots are the 5th roots of unity (Lem. 2.2) | \|z\|=1, 5θ∈ℤ; Φ₅ charge 5, M=1 | verified | |
| 15 | negative real root ⇒ even charge, not 5 (Lem. 2.2) | (x+1)Φ₅ has charge 10; x⁵−2 real root>0 | verified | |
| 16 | reciprocal unit trace 2 ⇒ r=1 (Lem. 2.5) | sympy solve r+1/r=2 | verified | |
| 17 | trace 3 ⇒ r=(3+√5)/2=φ²; trace 4 ⇒ 2+√3 (Lem. 2.5) | quadratic roots; monotone in k | verified | |
| 18 | real reciprocal unit ⇒ M∈{1}∪[φ²,∞) (Lem. 2.5) | scan integer traces, none in (1,φ²) | verified | |
| 19 | pentagon expansion [x³],[x²],[x¹],[x⁰] (Thm. 4.1) | sympy expand of the factored product | verified | 4 separate tests |
| 20 | Galois conjugacy t=σ(s) clears √5 from all coeffs (Thm. 4.1) | σ-invariance of each coeff | verified | |
| 21 | [x²]=k²−3m, [x⁰]=m² with k=s+t,m=st (Thm. 4.1) | substitute roots of y²−ky+m | verified | |
| 22 | minimizer (k,m)=(3,1)⇒x⁴−x³+6x²+4x+1 (Thm. 4.1/Sec.7) | build from s=φ²,t=φ⁻²; match coeffs | verified | |
| 23 | minimizer has charge Z/5Z and M=φ⁴ (Sec. 7) | independent charge + hp Mahler | verified | irreducible confirmed |
| 24 | M=max(1,s)²max(1,t)² (Thm. 4.1) | =φ⁴ for the minimizer | verified | |
| 25 | regime s,t>1 minimum M=16 at s=t=2,(k,m)=(4,4) (Thm. 4.1) | direct | verified | above global min |
| 26 | regime s>1>t,m=1 minimum M=φ⁴ at k=3 (Thm. 4.1) | s=(k+√(k²−4))/2, monotone | verified | |
| 27 | degree-4: M∈{1}∪[φ⁴,∞), so ∉(1,2) (Thm. 4.1) | φ⁴>2; global min φ⁴ | verified | |
| 28 | x⁵−m (m≥2): charge Z/5Z, non-reciprocal, M=m (Thm. 5.1) | charge + hp Mahler + palindrome test | verified | m=2,3,4,5,7,10 |
| 29 | roots of x⁵−m are m^(1/5)·(5th roots of unity) (Thm. 5.1) | mpmath roots, equal moduli | verified | |
| 30 | μ(5)=2 realized at x⁵−2 (Thm. 5.1) | hp Mahler of x⁵−2 = 2 | verified | |
| 31 | family realizes M∈{1}∪[2,∞) (Thm. 5.1) | m=1→1, m≥2→≥2, none in (1,2) | verified | |
| 32 | ψ⁵(x⁵−2): fifth powers all=2, totally positive, M(O)=32^(1/5)=2 (Rem. 6.2) | mpmath 5th powers; M(ψ⁵O)=32 | verified | |
| 33 | M(ψ⁵O)=M(O)⁵ generally (Rem. 6.2) | minimizer ψ⁵ totally positive, M=φ²⁰ | verified | |
| 34 | realification bound 2^(1/5)<2, strictly weaker (Rem. 6.2) | ordering | verified | |
| 35 | charge-5 quartics \|c\|≤10: distinct M={1,φ⁴} (Sec. 7) | independent enumeration | verified | |
| 36 | charge-5 quartics with M∈(1,2): 0 (Sec. 7) | independent enumeration | verified | |
| 37 | 13 non-reciprocal charge-5 objects in window (Sec. 7) | independent enumeration | verified | **exact count reproduced** |
| 38 | non-reciprocal realized min M=2 (x⁵−2 and kin) (Sec. 7) | enumeration min | verified | |
| 39 | all non-reciprocal charge-5 obey M≥μ_S (Sec. 7) | enumeration | verified | |
| 40 | non-reciprocal charge-5 in [μ_S,2): 0 (Sec. 7/Prop. 6.1) | enumeration | verified | the residual window is empty |
| 41 | reciprocal charge-5 measures = {1,φ²,2+√3} (Sec. 7) | independent enumeration | verified | distinct measures reproduce |
| 42 | reciprocal charge-5 in (1,2): 0 (Sec. 7) | enumeration | verified | |
| 43 | overall: no charge-5 object with M∈(1,2), floor 2 (Sec. 7/Prop. 6.1) | enumeration | verified | the μ(5)=2 window result |
| 44 | every reciprocal charge-5 object = Φ₅ × real reciprocal unit (Thm. 3.1b) | Φ₅ divides each; M in {1}∪[φ²,∞) | verified | |
| 45 | **5 reciprocal charge-5 objects in window (Sec. 7)** | independent enumeration | **xfail** | reproduced **4**, not 5 (see below) |
| 46 | Mahler measure is multiplicative M(pq)=M(p)M(q) (Def. 2.3) | hp Mahler on products | verified | |
| 47 | Kronecker: M=1 exactly on cyclotomics/powers of x (Def. 2.3) | Φ_n all have M=1; non-cyclotomic >1 | verified | |
| 48 | Smyth: non-reciprocal, M>1 ⇒ M≥μ_S (Lem. 2.4) | μ_S=M(x³−x−1); example scan | verified | sanity, not full proof |
| 49 | Z/3Z collapses to x³−2, M=2 (Sec. 1.2) | 2cos120°=−1; charge 3; M=2 | verified | |
| 50 | totally-positive (1,2) gap, floor 2 (Rem. 6.2 premise) | scan totally-positive deg≤3 | verified | small-case sanity of a cited result |

## VERIFIED (independently re-derived, then matched)

Claims 1–44 and 46–50 above (54 passing tests). Highlights of genuinely
independent re-derivations that reproduced the paper exactly:

- **Prop. 3.1** cosine identities, the quadratic x²+x−1, Galois conjugacy, and
  irrationality — exact symbolic.
- **Thm. 4.1** the full pentagon expansion (all four coefficients), the Galois
  integrality mechanism, and the explicit minimizer **x⁴−x³+6x²+4x+1 with charge
  Z/5Z and M=φ⁴** — built from the two moduli s=φ², t=φ⁻² with a from-scratch
  Mahler engine.
- **Thm. 5.1** x⁵−m ⇒ M=m and μ(5)=2 at x⁵−2; the ψ⁵ realification identity
  M(ψ⁵O)=M(O)⁵ with 2^(1/5)=1.1487.
- **Sec. 7 window** (quartics |c|≤10, quintics |c|≤4, sextics |c|≤3): an
  independent brute-force enumeration reproduced **no charge-5 object with M in
  (1,2)** (realized floor 2), **13 non-reciprocal objects** (exact), the quartic
  measure set **{1, φ⁴}**, and the reciprocal measure set **{1, φ², 2+√3}**.

## FLAGGED for human review

- **Claim 45 — "5 reciprocal charge-5 objects in the window" (Sec. 7 table).**
  Marked `@pytest.mark.xfail`. An independent enumeration of the paper's exact
  window finds **exactly 4** reciprocal charge-Z/5Z objects:
  `Φ₅`, `x⁵−1`, `Φ₅·(x²−3x+1)`, `Φ₅·(x²−4x+1)` — with measures `1, 1, φ², 2+√3`.
  The paper's **three distinct measures {1, φ², 2+√3} DO reproduce**; only the
  *object count* is off by one (paper 5 vs reproduced 4). This does not affect
  any mathematical conclusion of the paper (all reciprocal objects are
  Φ₅ × real-reciprocal-unit with M∈{1}∪[φ²,∞), and none lands in (1,2)). Most
  likely a benign counting-convention difference in the paper's engine (e.g.
  also counting the power-of-x object x·Φ₅, whose root at 0 has undefined angle
  and is excluded here, or double-counting a reflection). Notably the
  *non-reciprocal* count (13) reproduced **exactly** under the same conventions,
  which is what makes the reciprocal off-by-one worth a human glance.

## UNTESTABLE (documented, not mechanically decidable here)

- **The universal no-Salem dichotomy / μ(5)=2 as a *forced* claim.** The paper
  itself holds this at `[computed]` (Prop. 6.1, Rem. 6.2): the residual
  "no non-reciprocal charge-Z/5Z object has M∈[μ_S,2)" is only verified over the
  finite window (which we reproduced empty), not proven universally. The
  universal statement is an open proof obligation (Schur–Siegel–Smyth
  trace-problem territory) and is not mechanically decidable.
- **"μ_S is the *smallest* Pisot number."** We verify the Pisot property and the
  value; minimality is Siegel's theorem (1944), cited, not re-derivable in-session.
- **The general "totally-positive (1,2) gap"** used by the realification. We
  sanity-checked degree ≤3; the general theorem is a cited external result.
- **Corpus cross-references** — *Emission-Gap* Cor. 5.x, *vector_substrate*
  Thm. smythfloor + Rem. reciprocal, *Lehmer's Box*, *Occupant of the Salem
  Slot*. These live in other papers and are out of scope for this folder.
- **Provenance sha256 / epistemic-tag bookkeeping** — not mathematical content.

## Correction applied 2026-07-04

Sec. 7 reciprocal charge-Z/5Z object count corrected 5 -> 4; test test_reciprocal_object_count_is_4 now asserts 4. Distinct measures {1, phi^2, 2+sqrt3} and the non-reciprocal count 13 unchanged.
