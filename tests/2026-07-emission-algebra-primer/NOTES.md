# Verification notes -- *The Emission Algebra A* (2026-07 primer)

**Paper:** `papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex`
**Title:** *The Emission Algebra A -- A self-contained primer: one matrix, one relation, three layers*
**Author:** AceTheDactyl (@AceTheDactyl), Echo S Studios Research Developments

## Method

Every displayed constant, identity, worked example and exercise was
**re-derived independently** from the paper's premises and only then compared
to the paper's stated value. Independence discipline:

* Fibonacci / Lucas numbers come from the **bare recurrences** (`_helpers.fib`,
  `_helpers.luc`), not from any closed form the paper asserts.
* The power law `R^n = F_n R + F_{n-1} I` is checked by computing `R**n` by
  **actual matrix exponentiation** (sympy) and comparing to the independently
  built `F_n R + F_{n-1} I` -- the two routes never share code.
* The `ad_R` spectrum is obtained by assembling the **4x4 matrix** of the
  commutator map on `M_2` and diagonalising it, not by assuming the eigenvalues.
* Mahler measures of the golden-field objects are formed **symbolically** in
  `Q(sqrt5)` / `Q(5^{1/4})`; Lehmer's number is computed from its degree-10
  polynomial with `mpmath.polyroots` at `dps=50`.
* All decision-boundary equalities are settled with `sympy.simplify(...) == 0`
  (exact); floats appear only to confirm the paper's illustrative decimals.

Result: **every** mechanically-checkable claim reproduces exactly from the
paper's own premises. No arithmetic discrepancies were found -- consistent with
the paper's "every displayed constant is verified in exact arithmetic" claim.
Nothing needed to be flagged `xfail`.

## Final pytest run

```
104 passed  (0 failed, 0 xfailed, 0 skipped, 0 errors)
```

Files: `test_core.py` (19), `test_lie_layer.py` (20), `test_semiring.py` (16),
`test_reps.py` (9), `test_traceform.py` (14), `test_questions.py` (15),
`test_floor.py` (11).

## Claims table

| # | Claim (paper location) | How tested | Status |
|---|------------------------|------------|--------|
| 1 | `det(xI-R)=x^2-x-1`, keystone `R^2=R+I`, `Tr R=1`, `det R=-1`, eigenvalues phi,psi (Def 1.1 / Prop 1.2) | char poly + direct matrix mult; eigenvals | VERIFIED |
| 2 | Vieta `phi+psi=1, phi*psi=-1`; `phi-psi=sqrt5=phi+1/phi` (eq. 2) | exact symbolic | VERIFIED |
| 3 | Low powers `R^2..R^5 = {1,1},{2,1},{3,2},{5,3}` in (R,I) (Sec 1.3) | matrix power vs coeff form | VERIFIED |
| 4 | Power law `R^n=F_nR+F_{n-1}I=[[F_{n-1},F_n],[F_n,F_{n+1}]]`, all n (Thm 1.3) | matrix exp vs Fibonacci, n in [-16,16] | VERIFIED |
| 5 | `R^5=[[3,5],[5,8]]`, `R^{-1}=[[-1,1],[1,0]]`, `R R^{-1}=I` (Ex 1.4) | exact | VERIFIED |
| 6 | Lucas trace `Tr(R^n)=L_n=phi^n+psi^n` (Cor 1.5a) | trace of power vs Lucas recurrence & eigen-sum | VERIFIED |
| 7 | `det(R^n)=(-1)^n` (Cor 1.5b) | exact, n in [-14,14] | VERIFIED |
| 8 | Cassini `F_{n-1}F_{n+1}-F_n^2=(-1)^n` (Cor 1.5c) | exact | VERIFIED |
| 9 | Binet `phi^n=F_n phi+F_{n-1}`, `F_n=(phi^n-psi^n)/sqrt5` (Cor 1.5d) | exact | VERIFIED |
| 10 | Ex 1.6 (n=5): Tr=11=L5, det=-1, Cassini=-1, phi^5=5phi+3 (~11.09) | exact + decimal | VERIFIED |
| 11 | Ex 1.9: `R^{-2}=[[2,-1],[-1,1]]` two ways, det=1, Tr=3=L_{-2} | exact, both routes | VERIFIED |
| 12 | Z[R] free rank 2 (I,R independent) (Thm 1.7) | linear solve | VERIFIED |
| 13 | Maximal order: `disc(x^2-x-1)=5=disc(Q(sqrt5))`, index 1 (Thm 1.7) | discriminant + fundamental-disc rule | VERIFIED |
| 14 | Clumsy seed 2phi: `disc=20=2^2*5`, index 2 (Intuition after 1.7) | minimal_polynomial + discriminant | VERIFIED |
| 15 | `spec(ad_R)={0,0,+/-sqrt5}` on M_2 (Prop 2.2) | 4x4 commutator matrix eigenvalues | VERIFIED |
| 16 | Root vectors `N_+,N_-` explicit; eigen-eqns `[R,N_+/-]=+/-sqrt5 N`; 0-space=span{I,R}; v_phi.v_psi=0 (Thm 2.5) | exact matrix algebra + centraliser solve | VERIFIED |
| 17 | sl2 triple `[H,N_+]=2sqrt5 N_+`, `[H,N_-]=-2sqrt5 N_-`, `[N_+,N_-]=sqrt5 H`; standard normalisation (Thm 2.6) | exact commutators | VERIFIED |
| 18 | Signature (2,1): `H^2=S^2=5I`, `J^2=-I`, `N_+^2=N_-^2=0`; eigenvalue types (Prop 2.7) | exact squares + eigenvalues | VERIFIED |
| 19 | Rational structure `[H,S]=10J`,`[H,J]=2S`,`[S,J]=-2H`; `HS`,`SH`; adH on {S,J}=[[0,2],[10,0]], charpoly lambda^2-20 (Prop 2.9) | exact | VERIFIED |
| 20 | Null-frame `N_+/-=(S+/-sqrt5 J)/2`; transition matrices inverse & `det=-2/sqrt5`; realises basis change (Thm 2.10) | exact matrix identities | VERIFIED |
| 21 | Conjugator `V=[[1,1],[phi,psi]]`: `V^-1 H V=diag(sqrt5,-sqrt5)`, `V^-1 N_+ V=[[0,(5-sqrt5)/2],[0,0]]`, `V^-1 N_- V=[[0,0],[(5+sqrt5)/2,0]]` (Thm 2.10) | exact conjugation | VERIFIED |
| 22 | Ex 2.12: N_+N_-, N_-N_+ entrywise, `[N_+,N_-]=sqrt5 H`, `phi^2-psi^2=sqrt5` | exact | VERIFIED |
| 23 | Golden object `M(A_phi)=phi`, `chi={0,2}` (Ex 3.3) | symbolic Mahler + charge | VERIFIED |
| 24 | Roots of unity weightless (M=1); mu_4 charges 0,1,2,3; cube root -> charge 1 (Ex 3.4) | Mahler + charge rounding | VERIFIED |
| 25 | Grading laws: `M(A+B)=M(A)M(B)`, chi union; `chi(AxB)=chi(A)+chi(B) mod4`; `M(psi^n A)=M(A)^n`, `chi=n chi mod4` (Prop 3.6) | symbolic Mahler + charge sumsets | VERIFIED |
| 26 | Tropical coupling: `AxB={phi^3,-phi,phi^{-1},-phi^{-3}}`, `M=phi^4` (not naive phi^3); `+` does factor to phi^3 (Rem 3.7 / Ex 3.8) | raw products vs forms; branch-decided product | VERIFIED |
| 27 | Ex 3.9 Adams psi^3 on golden object; psi^6 tower | symbolic | VERIFIED |
| 28 | Adams pullback: Binet, `psi^2` is keystone, multiplicative-not-additive (differ by 2phi) (Prop 3.10) | exact | VERIFIED |
| 29 | sl2 irreps: `dim V_m=m+1`, weights {m..-m}, Casimir `m(m+2)/2` (Thm 4.1) | formula + weight strings | VERIFIED |
| 30 | Sym^n eigenvalue tower `{phi^{n-k}psi^k}`; H/sqrt5 integer weights; `M_2=V_2+V_0` (Prop 4.3) | symbolic tower + weight arithmetic | VERIFIED |
| 31 | Clebsch-Gordan dimension identity `(a+1)(b+1)=sum(j+1)`; V1xV1, V2xV2, V3xV3 (Thm 4.4 / Ex 4.5) | decomposition + dim bookkeeping | VERIFIED |
| 32 | Casimir on coupled targets: Cas(V3)=15/2, Cas(V4)=12 (Prop 5.10) | CG + Casimir formula | VERIFIED |
| 33 | Sym^2 eigenvalues {phi^2,-1,psi^2}, weights {2,0,-2} (Ex 4.7) | symbolic | VERIFIED |
| 34 | Trace form diag(10,10,-2), signature (2,1), invariance (Prop 5.2) | Gram matrix + trace cyclicity | VERIFIED |
| 35 | Deviation `X_n=2R^n-L_n I=F_n H` (integer, traceless); X_1,X_3,X_5 (Prop 5.4 / Ex 5.5) | matrix build vs F_n H | VERIFIED |
| 36 | Trace-Form Duality `(1/2)Tr(X_n^2)=5F_n^2=L_n^2-4(-1)^n=(phi^n-psi^n)^2` (Thm 5.6) | exact, n in [-8,8] | VERIFIED |
| 37 | psi^n dilation: gap F_n sqrt5; `ad_{R^n}=F_n ad_R`; index multiplicative; seed lift (Thm 5.7) | exact on test matrices | VERIFIED |
| 38 | Deviation ladder `(1/2)Tr_{V_m}(X_n^2)=5F_n^2 C(m+2,3)=(5F_n^2/3)dim Cas`; V3->50F_n^2, V4->100F_n^2 (Thm 5.11) | weight-sum-of-squares, both forms | VERIFIED |
| 39 | Ex 5.12: X_4=3H, (1/2)Tr=45, ladder into V_2 = 180 | exact | VERIFIED |
| 40 | Obstruction 1: `Phi_R(R)=R^2-R-I=0`, `[R,R]=0` (Prop 6.1) | exact | VERIFIED |
| 41 | Obstruction 2: {0,2} closed under union/+/double, never reaches {1,3} (Prop 6.2) | Z/4Z arithmetic | VERIFIED |
| 42 | Obstruction 3: `theta+1/theta` in [-2,2] on circle, >2 for real theta>1 (excess (theta-1)^2/theta); golden combos stay real (Prop 6.3) | symbolic + samples | VERIFIED |
| 43 | The flip: `D=1+4C`, Gram det `=4C+1=D`, flip at C=-1/4; C=1 real golden roots (D=5), C=-1 cube roots (D=-3) (Prop 6.4 / Ex 6.5) | discriminant + field-trace Gram | VERIFIED |
| 44 | Partner `K=x^4+5x^2-5`: y-roots straddle 0; real roots `+/-K=5^{1/4}/phi`; imaginary `+/-i beta`; `M(K)=beta^2=phi^2 sqrt5`; `chi={0,1,2,3}` (Thm 6.6) | quartic solve + symbolic Mahler/charge | VERIFIED |
| 45 | Parity criterion: even quartic w/ negative x^2-root -> imaginary complex roots (charges {1,3}); non-even -> off-axis (Prop 6.7) | several even quartics + one non-even | VERIFIED (see notes) |
| 46 | Ex 6.8: `x^4-1` roots {1,-1,i,-i}, charges Z/4Z, M=1 (floor charge completion) | exact | VERIFIED |
| 47 | Kronecker: cyclotomic <-> M=1; non-cyclotomic M>1 (Thm 7.1) | illustrated on cyclotomic & non-cyclotomic examples | VERIFIED (instances) |
| 48 | Boundary conditions B1-B5 (floor/identity, sub-semiring closure, collapse to charge group mu_4=Z/4Z, golden gap, one-way) (Thm 7.2) | Mahler + charge-group homomorphism checks | VERIFIED |
| 49 | Golden gap `phi-1~0.618`; Lehmer M~1.17628 inside (1,phi) (B4 / Ex 7.4) | mpmath polyroots dps=50 | VERIFIED |
| 50 | Ex 7.5 floor-collapse orbit `{i,-i}->{-1,-1}->{1}`, M==1 throughout, charge cycles | exact | VERIFIED |

## VERIFIED (all mechanically-checkable claims)

Claims 1-50 above all reproduced exactly. Highlights of genuinely independent
re-derivation: the power law and its matrix form (two disjoint routes), the
`ad_R` spectrum (assembled 4x4 map), the full sl2 bracket table and Lorentzian
metric, the Trace-Form Duality across `n in [-8,8]`, the deviation ladder in
both binomial and Casimir forms across `m in [0,6]`, the K-partner root split
with `K=5^{1/4}/phi` and `M(K)=phi^2 sqrt5`, and Lehmer's `1.176280818...`
computed from the polynomial itself.

## FAILED / flagged for human review

**None.** No arithmetic claim failed to reproduce from the paper's premises.

During test development four tests failed and were all traced to **test-side
bugs**, then fixed (never by loosening a tolerance to hide a discrepancy):

* Python evaluates `(-1)**n` to a *float* for negative `n`; replaced with the
  exact-parity `(-1)**(n % 2)` in the det/Cassini/duality tests.
* Mahler measures of cyclotomic roots came from double-precision `complex()`
  giving `|zeta| = 1 + 1e-16`; fixed by routing `_mpc` through decimal strings
  at `dps=60` and adding a `1e-9` unit-circle guard band (safe: the next
  emitted magnitude is `phi ~ 1.618`, and Lehmer `1.176`, both far from 1).
* Over-tight `1e-40` float tolerances on golden-field magnitudes were replaced
  by **exact symbolic** Mahler comparisons (`mahler_exact`) in `Q(sqrt5)`.

## UNTESTABLE (documented, not mechanically checkable)

* **Epistemic tags & interpretive overlays ([POSITED]).** The physical readings
  -- light-cone/rapidity picture of the null frame (Intuition after Thm 2.10),
  `Mah=1` as "vacuum" / charge as "superselection" (Intuition after Thm 7.2),
  "terrain vs rotation" as spacetime signature -- are explicitly lenses, not
  theorems. Not checkable; noted only.
* **Global uniqueness of the partner K** ("the unique catalogue seed to
  complete the charge", Thm 6.6). Uniqueness is over an unspecified "emission
  catalogue"; the *mechanism* (the parity criterion, Prop 6.7) is tested, but
  the global uniqueness quantifier is not mechanically enumerable here.
* **Genericity of non-even quartics** ("non-even quartics of Q(5^{1/4}) place
  their complex roots at generic angles", Prop 6.7). Tested on one concrete
  non-even quartic (off-axis root confirmed) and the forward even-quartic
  direction on several; the universal genericity statement is not exhaustively
  checkable.
* **Kronecker's theorem (Thm 7.1, [ESTABLISHED])** and **Casimir/irrep
  classification (Thm 4.1, [ESTABLISHED])** are cited external theorems; their
  *instances* are verified but the theorems themselves are not re-proved.
* **Lehmer's problem (B4 general case, [OPEN])** -- whether `(1,phi)` (or any
  `(1,1+eps)`) is empty over *all* non-cyclotomic polynomials -- is an
  acknowledged open problem; only the in-class gap and the specific Lehmer
  value are testable.
