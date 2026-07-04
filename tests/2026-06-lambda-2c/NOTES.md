# Verification notes — `2026-06-lambda-2c`

**Paper:** *The Exchange Rate λ = 2c: A Conformal Identity, Its Gate, and Its Flip*
(AceTheDactyl / Echo S Studios Research Developments).
Source of truth: `papers/2026-06-lambda-2c/lambda_2c_paper.tex`.

**Method.** Every mechanically-checkable claim is INDEPENDENTLY re-derived from the
paper's stated premises (sympy exact symbolic; mpmath at dps 40–50 for transcendental
values; numpy for enumerations) and then compared to the paper's stated value. No test
restates the paper's final number and compares it to itself — values are built up
independently first (e.g. Mahler measures from computed roots, R⁴ from R, min-Mahler by
enumeration, KL from its integral, the ad-matrix spectrum from `[R,·]`).

**Result:** 72 tests, **72 passed, 0 failed, 0 xfailed, 0 skipped**. Every load-bearing
computational identity in the paper reproduced exactly. No inconsistencies found — this
paper's abstract claim that "every load-bearing identity is machine-verified in exact
arithmetic" holds up under independent re-derivation.

---

## Claims table

| # | Claim (paper ref) | How tested | Verdict | Notes |
|---|---|---|---|---|
| 1 | 2nd-order KL = ½·Fisher form = (1/2c)‖r‖²_G (eq:kl2c) | Gaussian-location KL closed form vs ½rᵀΣ⁻¹r; exp-family KL Hessian = Fisher | verified | Gaussian is exactly quadratic; exp-family confirms O(‖r‖³) statement |
| 2 | λ = 2c (thm:lambda2c) | solve balance KL=logM for the gain coefficient | verified | coefficient of logM is 2c |
| 3 | c=1→λ=2, c=n→λ=2n (thm:lambda2c pf) | substitute | verified | |
| 4 | σ = 1/(2c) = 1/λ (rem:sigma) | symbolic | verified | |
| 5 | trace form = conjugate Gram, G=BᵀB (lem:fisher substrate) | Q(√2), Q(√5): G_ij=Tr(θ^{i+j})=(BᵀB)_ij | verified | full "G = n·I_exp" statistical reading is model-specific → row U2 |
| 6 | disc(x²+x−C) = 1+4C (eq:gengate) | sympy discriminant | verified | |
| 7 | R_C=[[0,C],[1,−1]]: charpoly x²+x−C, eigs (−1±√(1+4C))/2, gap √(1+4C) (eq:companion) | eigenvals of R_C | verified | |
| 8 | gate table C∈{¼,½,1}→D∈{2,3,5}, seeds √D (Sec 5) | per-row D and companion gap | verified | |
| 9 | golden gate roots {1/φ, −φ} (Sec 5) | roots of x²+x−1 | verified | |
| 10 | disc(x²+x−1)=disc(x²−x−1)=5 (Sec 5) | discriminant | verified | |
| 11 | spec(ad_{R_C})={−√(1+4C),0,+√(1+4C)}, charpoly t²(t²−(1+4C)), 0 mult 2 (prop:trifurcation) | build 4×4 ad matrix, eigenvals | verified | golden gate {−√5,0,√5} checked |
| 12 | 0-eigenspace = span{I, R_C} (prop:trifurcation) | [R,I]=[R,R]=0, ad rank 2 | verified | |
| 13 | gate balance λC=r(R_C)=√(1+4C) at C=1 (eq:gatebalance) | √5·1=√5 | verified | |
| 14 | Mah(x²−D)=D; r(R_C)=√D=√Mah (rem:mahler) | Mahler from roots | verified | |
| 15 | frame-shift c=√(1+4C)/(2C) (eq:frameshift) | solve 2cC=√(1+4C) | verified | |
| 16 | at C=1: c=√5/2, λ=√5=φ−ψ (def:frameshift) | symbolic | verified | |
| 17 | three canonicalizations c∈{1,n,√(1+4C)/(2C)} (canon table) | selectors + cost floors | verified | |
| 18 | ternary ⇔ degree 2: d²−d+1=3 iff d=2 (prop:ternary) | solve over ℤ⁺; Golomb-ruler channel count | verified | d=3→7 channels |
| 19 | min Mahler>1 over integer quadratics = φ, only at disc 5 (lem:mincost/lem:tie) | enumerate x²+bx+c | verified | attained only (b,c)∈{(−1,−1),(1,−1)} |
| 20 | plastic μ_S=1.32472 (x³−x−1) < φ; min over cubics is plastic (rem:cubic) | mpmath root + enumerate cubics | verified | |
| 21 | squaring firewall x⁴+5x²−5→y²+5y−5, disc 45=9·5 (prop:firewallimage) | substitution y=x², discriminant | verified | |
| 22 | firewall Φ=companion∘charpoly fixes companions; φ⊕φ derogatory (Sec 15, def:firewall) | Φ(R_φ)=R_φ; φ⊕φ minpoly deg 2 < charpoly deg 4 | verified | |
| 23 | keystone R²=R+I unique to golden; τ:R²=I−R; radicand:R²=DI (lem:keystone) | matrix identities | verified | R²=R+I forces eigs {φ,ψ} |
| 24 | forced value C=1, c=√5/2, λ=√5=φ−ψ (thm:gateforced) | symbolic | verified | |
| 25 | Perron breaks tie: x²−x−1 dominant +φ (R²=R+I); τ dominant −φ (R²=I−R); detR=−1 (lem:perron) | dominant-eig sign, relations, det | verified | |
| 26 | φ = smallest Perron root of 2×2 primitive nonneg int matrix (Fibonacci) (thm:keystonederived) | enumerate 0..3 entries, primitivity, Perron | verified | argmin = Fibonacci matrix/transpose |
| 27 | swap [[0,1],[1,0]] Perron root 1, excluded by growth (thm:keystonederived pf) | eigenvalues {1,−1} | verified | |
| 28 | Rⁿ=F_n R+F_{n−1}I; eigs φⁿ,ψⁿ; tr=L_n∈ℤ; det=(−1)ⁿ (lem:keypowers) | n=1..10 | verified | |
| 29 | Pell L_n²−5F_n²=4(−1)ⁿ; charpoly x²−L_n x+(−1)ⁿ, disc 5F_n² (lem:pell) | n=1..10 | verified | |
| 30 | R⁴=3R+2I=[[2,3],[3,5]], charpoly x²−7x+1, entries {2,3,5}={F₃,F₄,F₅}, L₄=7 (prop:L4forced) | compute R⁴; Pell at n=4 → L₄=7 | verified | forced three ways |
| 31 | roots of x²−7x+1 are φ⁴,φ⁻⁴; z_c=√(L₄−4)/2=√3/2 (prop:L4forced) | solve + symbolic | verified | |
| 32 | drop-one witnesses: growth→cyclotomic Mah(x²+1)=1; integrality→no floor (App A) | Mahler; Perron→1⁺ as ε→0 | verified | degree-drop→plastic tested in #20 |
| 33 | App A identity (iv): {2,3,5}={F₃,F₄,F₅}, 7²=5·3²+4=49 (App A) | integer arithmetic | verified | |
| 34 | flip table: C∈{1,0,−¼,−1,−2}→D∈{5,1,0,−3,−7}; roots real/double/complex (thm:flip) | discriminant + roots per row | verified | |
| 35 | double root −½ at C=−¼; roots {−1,0} at C=0 (thm:flip) | roots dict | verified | |
| 36 | C=−1: x²+x+1 roots e^{±2πi/3} on unit circle, D=−3 (thm:flip/prop:meet) | roots, r³=1, \|r\|=1, disc | verified | |
| 37 | det G = 4D in gap basis; G=diag(2,2D) (prop:detG) | Tr(1)=2,Tr(√D)=0,Tr(D)=2D | verified | |
| 38 | signature flip: PD iff D>0, Lorentzian iff D<0, degenerate at 0 (prop:detG) | eigen-signs of diag(2,2D) | verified | |
| 39 | Q(i) basis {1,i}: diag(2,−2), det=−4=4·(−1) (prop:detG) | conjugate Gram | verified | |
| 40 | two flips meet at C=−1: Mah(x²+x+1)=1 (prop:meet) | Mahler from roots | verified | |
| 41 | K-formation y-roots (−5±3√5)/2 straddle 0 (prop:kform) | solve y²+5y−5 | verified | 0.8541, −5.8541 |
| 42 | K=√(y₊)=5^{1/4}/φ≈0.9242 (inside circle) (prop:kform) | symbolic equality + numeric | verified | |
| 43 | β=√\|y₋\|≈2.4195; Mah(x⁴+5x²−5)=β²=(5+3√5)/2 (prop:kform) | independent Mahler from roots | verified | |
| 44 | complex place off circle: \|iβ\|=2.4195, \|5^{1/4}i\|=1.4953 (ssec:kformface) | mpmath | verified | |
| 45 | Lehmer number L=1.17628 (root of Lehmer poly) (rem:lehmer) | mpmath findroot | verified | |
| 46 | Smyth floor μ_S=1.32472; band ordering 1<L<μ_S<φ (Smyth) | mpmath | verified | |
| 47 | catalog roots have argument in (π/2)ℤ (ssec:angle) | np.angle of 7 catalog minpolys | verified | |
| 48 | (π/2)ℤ closed under add/double; on-circle emit ⇒ 4th root of unity (ssec:angle) | subgroup closure | verified | |
| 49 | Salem conjugate has irrational angle (not in (π/2)ℤ) (ssec:angle) | angle of deg-4 Salem's on-circle conjugate | verified | mechanism of the no-Salem argument |
| 50 | spec(ad_R)={0,±(φ−ψ)}={0,±√5} difference set (ssec:uniform) | difference set of {φ,ψ} | verified | |
| 51 | smallest deg-4 Salem β₄=1.72208 > φ, a genuine Salem (ssec:uniform) | roots of x⁴−x³−x²−x+1, Salem signature | verified value | "smallest" is literature/[EG] → row U8 |
| 52 | engine table: Mah(x²−D)=D=1+4C; det G=4D (sec:engine) | Mahler + det | verified | |
| 53 | engine cost floors: Mah(x²−24)=24 (2√6), Mah(x²−7)=7 (√7); 2log24, 2log7 (sec:engine) | Mahler; 2√6=√24 | verified | |
| 54 | Clifford unity τ=x²+x−1 = golden gate (C=1) (sec:engine) | gate poly at C=1 | verified | |
| 55 | Kuramoto Lorentzian threshold K_c=2/(πg(0))=2γ (sec:kuramoto) | g(0)=1/(πγ); density integral = 1 | verified | |
| 56 | mean-field order-parameter exponent ½: r~√(K−K_c) (sec:kuramoto) | Landau normal-form steady state | verified | full self-consistency integral not derived → U7 |
| 57 | z_c=√3/2 ↔ C=½ gate; D=(2z)²=4z², C=z²−¼ (prop:square) | algebra | verified | "negentropy-peak" linkage → U7 |

**Distinct checkable claims with an assertion test: 57** (across 72 test functions).

---

## VERIFIED (independently reproduced) — all 57 rows above

Highlights where independence matters most:
- **λ = 2c** re-derived by solving the MDL balance for the gain coefficient, not by
  assertion (rows 1–3).
- **min Mahler = φ over integer quadratics**, attained *only* at disc 5 — by brute
  enumeration + Mahler-from-roots, matching the paper's uniqueness claim (row 19).
- **φ = smallest Perron root of a 2×2 primitive non-negative integer matrix** — by
  enumeration; argmin is exactly the Fibonacci matrix (row 26).
- **R⁴ = [[2,3],[3,5]], L₄ = 7** computed from R = [[0,1],[1,1]] and cross-checked
  against the Pell relation at n=4 — the "forced three ways" claim holds (row 30).
- **det G = 4D** and the **Q(i) → diag(2,−2)** signature flip re-derived from field
  traces (rows 37–39).
- **No-Salem mechanism**: catalog arguments confined to (π/2)ℤ, closure under
  doubling/addition, and the deg-4 Salem's on-circle conjugate landing at an
  irrational angle — the computational core of the emission-gap argument (rows 47–51).

## FAILED / FLAGGED for human review

**None.** No computational claim failed to reproduce; no xfail was warranted. Four of my
own tests initially failed on transcription/re-derivation bugs (a `Matrix == 0` idiom, a
non-generic eigenvalue set with colliding differences, an `exp(2πi/3)` vs rectangular
comparison, and a fragile `minimal_polynomial` call) — all fixed in the tests, none
implicating the paper.

## UNTESTABLE (documented, not mechanically checkable here)

- **U1 — Čencov's theorem (thm:cencov).** A cited uniqueness/impossibility meta-theorem
  ("Fisher is the unique Markov-invariant metric up to scale"); not a finite computation.
- **U2 — lem:fisher full statistical statement.** "G = n·I_exp = n × empirical Fisher of
  the location family over the n conjugates" under trace-zero and 1∈col(B). The paper's
  location-family model isn't fully specified here; I verified its load-bearing algebraic
  substrate (G_ij = Tr(θ^{i+j}) = (BᵀB)_ij, trace form = conjugate covariance, rows 5),
  but not the exact "n×" normalization/conditions.
- **U3 — two-part MDL principle & Northcott justification** that log M(θ) is the right
  code length (eq:balance rationale). Decision-theoretic justification, not a computation;
  the *algebra* of the balance IS tested (row 2).
- **U4 — thm:collapse** "the ternary arity and the keystone are the same datum." A logical
  identification; its two computational components (deg-2 realizes ternary, row 18;
  min positive-growth companion at deg 2 is R²=R+I, rows 25–26) are verified, the
  "sameness" itself is conceptual.
- **U5 — thm:floorresolved / emission no-Salem theorem in full.** "𝒮 emits no Salem
  number at any matrix size", the signature-lattice of ℚ(√2,√3,5^{1/4}), and the uniform
  bound — proofs cite the companion paper [EG]. I verified the computational mechanisms
  (angle confinement, difference-set spectrum, β₄>φ, K-formation complex place off the
  circle; rows 47–51, 44), not the full theorem.
- **U6 — prop:square empirical `\computed` claim** "verified across sub-/super-critical
  Kuramoto / Stuart–Landau and gate-seeking gradient flows." Specific flow simulations
  are not reproducible from the paper; the algebraic pinning D=4z²≥0 IS tested (row 57).
- **U7 — Kuramoto z_c=√3/2 "where the negentropy peaks."** The negentropy-peak criterion
  that pins the critical coherence to √3/2 is not fully specified in this paper. Tested:
  K_c=2γ (row 55), the ½ exponent (row 56), and z=√3/2 ↔ C=½ (row 57); NOT the derivation
  that Kuramoto's critical coherence is specifically √3/2.
- **U8 — "smallest degree-4 Salem number" = β₄=1.72208.** The *value* reproduces (verified
  it is a genuine Salem number, root of x⁴−x³−x²−x+1, exceeding φ; row 51), but the
  *minimality* over all degree-4 Salem numbers is a literature/[EG] enumeration claim not
  independently proven here.

Untestable items: **8** (several with their computational core nonetheless verified).
