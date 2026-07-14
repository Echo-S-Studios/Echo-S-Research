# Verification notes — *The Vector Substrate: Number Fields as Exact Learning Geometry*

**Paper:** `papers/2026-06-vector-substrate/vector_substrate.tex`
**Author:** AceTheDactyl (@AceTheDactyl), Echo S Studios Research Developments (June 2026)
**Tests:** `tests/2026-06-vector-substrate/` — 8 `test_*.py` files + `vsub_nf.py` (shared exact number-field helpers).
**Result:** 63 pytest functions, **63 passed / 0 failed / 0 xfail / 0 skipped**. Every mechanically-checkable claim reproduced *exactly* from the paper's own premises. No inconsistencies found.

## Method / independence

Nothing is compared to itself. `vsub_nf.py` rebuilds the machinery from scratch:
- **Regular representation** from the companion matrix in the paper's stated convention (`theta^k -> theta^{k+1}`, last column `-c`), and `rho(alpha)=sum a_i C^i` via the algebra-homomorphism property.
- **Trace form** `G_ij = Tr(w_i w_j) = trace(rho(w_i) rho(w_j))` — exact rational, independent of the paper's displayed matrices.
- **Invariant factors** via **determinantal divisors** (monic gcd of every `i x i` minor of `xI-A`), an independent route (sympy's library SNF rejects `QQ[x]` as "not a PID", so it was *not* used). Non-similarity is additionally witnessed by similarity invariants that avoid the invariant-factor machinery (annihilating polynomials; matrix rank).
- **Fisher information** re-derived from first principles: the Gaussian one as the Hessian of the log-likelihood of `N(Ma, I)`; the exponential-family one as the Hessian of the log-partition `A(a)=log sum_k exp((Ma)_k)` at `a=0`.
- **`G = M^T M`** cross-checked numerically against Minkowski embeddings (mpmath, 40 digits) built from polynomial roots.
- Transcendental checks use mpmath at `dps>=40`; the growth "cost" (`log Mahler`) is compared via the same **rigorous interval enclosures** the paper uses.
- A separate mutation check confirmed the helpers *reject* perturbed (wrong) values, so no assertion is vacuous.

## Claims table

| # | Claim (paper ref) | How tested | Status |
|---|---|---|---|
| 1 | `rho(phi)=[[0,1],[1,1]]`, `Tr=1=Tr(phi)`, `det=-1=N(phi)` (Ex. 2.3) | build companion, trace/det + Newton identities | verified |
| 2 | `rho(sqrt5)=2rho(phi)-I`, charpoly=minpoly=`x^2-5` (Ex. 2.3) | build matrix, charpoly + independent matrix-minpoly | verified |
| 3 | `rho(3)=3I` derogatory: charpoly `(x-3)^2`, minpoly `x-3` (Ex. 2.3) | build, compare charpoly vs minpoly | verified |
| 4 | Prop. 2.2 Eq.(2): charpoly`=m_x^{n/d}`, minpoly`=m_x` | `x=sqrt2` in deg-4 field: `(x^2-2)^2` vs `x^2-2` | verified |
| 5 | Ex. 2.7 / Prop. 2.6: `M rho(phi) M^{-1}=diag(phi,phi')` | symbolic conjugation with `M` | verified |
| 6 | Prop. 2.6: same `M` diagonalises every element | check on `sqrt5 -> diag(sqrt5,-sqrt5)` | verified |
| 7 | Ex. 2.8: quartic spectrum = 4 conjugates, 2 outside unit circle | mpmath roots vs `+-sqrt2+-sqrt3` | verified |
| 8 | Rem. 2.5: `Mahler(sqrt2+sqrt3)=5+2sqrt6` (spectral) | closed form + Mahler rebuilt from spectrum | verified |
| 9 | Rem. 2.5 / Fig. 2: `Mahler(plastic)=mu_S` | Mahler of `x^3-x-1` vs its real root | verified |
| 10 | Ex. 4.4: `rho(2sqrt6)=rho(theta)^2-5I` (stated 4x4) | build both, compare | verified |
| 11 | Thm. 4.2 / Ex. 4.4: invariant factors `(x^2-24,x^2-24)`, minpoly `x^2-24` | determinantal divisors + matrix minpoly + `minimal_polynomial(2sqrt6)` | verified |
| 12 | Thm. 4.2: algebraic integer => monic-`Z` minpoly (admit seed) | check monic + integer coeffs | verified |
| 13 | Ex. 2.11: charpoly insufficient (`phi(+)phi` vs `C((x^2-x-1)^2)`) | IFs differ; independent witness (annihilating poly) | verified |
| 14 | Ex. 2.13: Jordan needs full list (`(x^2,x^2)` vs `(x,x,x^2)`) | IFs differ; independent witness (rank 2 vs 1) | verified |
| 15 | Thm. 2.10: prod(IF)=charpoly, largest=minpoly | on derogatory `3I` | verified |
| 16 | Ex. 2.15: golden `G=[[2,1],[1,3]]`, `det=5=d_K`, traces 2,1,3 | rebuild Gram from traces | verified |
| 17 | Thm. 2.14: `G=M^T M` (golden, symbolic) | symbolic `M^T M` | verified |
| 18 | Thm. 2.14: `G=M^T M` (6 fields, numeric power basis) | mpmath embeddings vs exact Gram | verified |
| 19 | Ex. 2.16: biquadratic `G=diag(4,8,12,24)`, `det=9216`, index 2 | rebuild Gram; `d_K=8*12*24` | verified |
| 20 | Ex. 2.18: `Q(i)` trace form `diag(2,-2)` (sig (1,1)); `G_2=2I`, `det=4` | rebuild; eigen-signs; covol `sqrt4=2` | verified |
| 21 | Ex. 2.19: `Q(cbrt2)` `G=[[3,0,0],[0,0,6],[0,6,0]]`, `det=-108`, sig (2,1) | rebuild; eigenvalues `3,+-6`; `det(G_2)=108` | verified |
| 22 | Table 2 catalog: `det G=d_K` and `N(different)=N(m'(theta))=|d_K|` | rebuild `det`; norm of `m'(theta)` for 6 fields | verified |
| 23 | Ex. 3.9: trace duality, dual basis, `different=(sqrt5)`, `N=5` | `G^{-1}`, dual-basis defining property | verified |
| 24 | Ex. 3.4: `B=(0,1)`, `B^T G B=3`, `P=(1/3)[[0,0],[1,3]]` | rebuild `P` | verified |
| 25 | Prop. 3.2 (i),(ii): `P^2=P`, `B^T G r=0` | symbolic | verified |
| 26 | Ex. 3.4: capture `2phi` (r=0); `sqrt5` -> `P=(0,5/3)`, `||r||^2=5/3` | rebuild residual | verified |
| 27 | Prop. 3.2 (iii),(iv): `Px=x <=> in col(B) <=> ||r||^2=0` | multiple targets | verified |
| 28 | Ex. 3.5: `Q(i)` Hermitian projector, `r=i`, `||r||_{G2}^2=2` | rebuild over `G_2` | verified |
| 29 | Eq.(17): Pythagoras `||x||^2=||Px||^2+||r||^2` | on `sqrt5` vs line `Q.phi` | verified |
| 30 | Prop. 7.13 / Ex. 7.15: one-step Newton `a*=5/3`, `Ba*=Px` | rebuild grad/Hessian/step | verified |
| 31 | Prop. 4.6: `det(A(x)B)=(det A)^{dim B}(det B)^{dim A}` | Kronecker det identity | verified |
| 32 | Cor. 4.7: multiquadratic tower `{2,3}` -> multiset `{4,8,12,24}` | Kronecker of factor Grams + subset-product formula | verified |
| 33 | Ex. 4.8: adjoin `sqrt7`, `G_KL=diag(8,56,16,112,24,168,48,336)`, `det=9216^2*28^4` | Kronecker product | verified |
| 34 | Prop. 4.6: Kronecker Gram from trace definition | product basis Gram of `Q(sqrt2,sqrt3)` = `G_{Q(sqrt2)} (x) G_{Q(sqrt3)}` | verified |
| 35 | Ex. 4.10: non-disjoint factorisation over `Q(sqrt2)`, `[K(beta):K]=2<4` | expand factorisation; degrees | verified |
| 36 | Ex. 5.7: exact heights table (deg, ch, `sum c_i^2`) | rebuild from polynomials | verified |
| 37 | Ex. 5.7: `phi^4` minpoly `= x^2-7x+1` | `minimal_polynomial(phi^4)` | verified |
| 38 | Thm. 5.2(b): Landau `Mahler<=||p||_2`, tight for `2sqrt6` (`576<=577`) | mpmath Mahler vs 2-norm | verified |
| 39 | Thm. 5.2 proof: `|c_j|<=binom(d,j) Mahler` | per-coefficient bound, 4 polys | verified |
| 40 | Ex. 5.6 / Prop. 5.5: Northcott count `= 12` (`D=2,H=1`) | formula + explicit enumeration | verified |
| 41 | Prop. 5.5: size formula (also `D=3,H=2 -> 155`) | formula vs enumeration | verified |
| 42 | Ex. 5.6: irreducibility classification of the 9 quadratics | sympy `is_irreducible` / `factor` | verified |
| 43 | Def. 5.4 / Ex. 5.6: integer gate; `2sqrt6` needs `H_max>=24` | admissibility at `H=23` vs `24` | verified |
| 44 | Rem. 5.9: float Mahler `!=` exact; integer gate robust | `math.sqrt(7)**2 = 7.000000000000001 > 7`; exact gate admits | verified |
| 45 | Thm. 7.1: Gaussian Fisher `= M^T M = G` | Hessian of `-log N(Ma,I)` | verified |
| 46 | Thm. 7.7: exp-family Fisher `= (1/n)G-(1/n^2)tt^T = [[0,0],[0,5/4]]` | Hessian of log-partition `A(a)` at 0 | verified |
| 47 | Thm. 7.7: general identity `M^T(I/n-11^T/n^2)M=(1/n)G-(1/n^2)tt^T` | symbolic 3x2 `M` | verified |
| 48 | Sec. 7.3 table: `Fisher_exp` for `Q(sqrt5)`, `Q(sqrt2,sqrt3)`, `+sqrt7` | rebuild `(1/n)G-(1/n^2)tt^T` | verified |
| 49 | Lemma 7.8 + Rem. 7.11: `<v,1>_G=Tr(v)`; both-ways `1 in col(B)` check | trace vs Gram pairing; residual of `1` has trace `5/3` | verified |
| 50 | Thm. 7.9: `||r||_G^2 = n Fisher(r)` on trace-zero; `||sqrt5||^2=10=2*5` | rebuild + whole-subspace identity | verified |
| 51 | Prop. 7.2 / Ex. 7.5: `vol_Jeffreys=sqrt(det G)=sqrt|d_K|`; ~4.6x ratio | `sqrt|det G|` for 3 fields; numeric volumes | verified |
| 52 | Rem. 7.4: `(1/2)log det G = log covol` | symbolic | verified |
| 53 | Thm. 8.2: `mu_S=1.3247179...` = root of `x^3-x-1`, `log mu_S>0` | mpmath root; floor `0.562` | verified (value) |
| 54 | Conj. 8.4: `mu_L=1.17628081...` = root of Lehmer's poly, `< mu_S` | mpmath root of degree-10 poly | verified (value) |
| 55 | Sec. 6: residual of `2sqrt6` off `{1,theta}` has `||r||^2=Tr(24)=96` | rebuild residual + trace | verified |
| 56 | Sec. 6: `||theta+2sqrt6||^2=116=20+96`, cross term 0 | rebuild norms + cross term | verified |
| 57 | Sec. 6: seed `x^2-24` admissible `(64,256)`, GROW (`96>6.356>0.562`) | gate + gain/cost/floor | verified |
| 58 | Sec. 8.3 table: `2sqrt6` gain 96 / `sqrt7` gain `8*7=56`, both GROW | rebuild gains + costs | verified |
| 59 | Ex. 8.7: certified GROW via `log7 in [1.94591,1.94592]`, `log24 in [3.17805,3.17806]` | mpmath enclosures dominate exact gains | verified |
| 60 | Ex. 8.5: `(1/10)sqrt5` trace-zero, `||r||^2=1/10<0.562` -> STOP | rebuild norm; compare floor | verified |
| 61 | Rem. 8.10: degree-aware floor `n*2*log mu_S = 2.25 (n=4), 4.50 (n=8)` | mpmath | verified |
| 62 | Conj. 7.16 (**cross-domain**, companion GSA): golden gate `c=sqrt(1+4C)/(2C)`, at `C=1` gives `c=sqrt5/2`, `lambda=2c=sqrt5`, ladder root `1/phi` | internal-consistency only (per policy) | verified (consistency) |

## Verified (all 62 checkable claims)

Every displayed numeric constant, matrix, Gram/trace form, discriminant, Mahler value, height, Fisher matrix, projector/residual, Kronecker compositum, growth/threshold decision, and worked-episode number reproduces exactly from the paper's own definitions. Highlights of genuinely independent re-derivations:
- Invariant factors via determinantal divisors (not any library SNF) reproduce `(x^2-24,x^2-24)`, `(x^2-x-1,x^2-x-1)`, `(x^2,x^2)`, `(x,x,x^2)`.
- Both Fisher metrics derived from log-likelihood / log-partition Hessians match `G` and `(1/n)G-(1/n^2)tt^T`.
- `G=M^T M` confirmed two ways (symbolic + 40-digit Minkowski embeddings across 6 fields).
- The paper's float artifact `Mahler(sqrt7) ~ 7.000000000000001` reproduces to the digit on this platform (`math.sqrt(7)**2`), and the integer gate is confirmed immune to it.

## Flagged for human review

**None.** No claim failed to reproduce; there are no `xfail`s. The single cross-domain constant (golden-gate scale `c=sqrt5/2`, exchange rate `lambda=2c=sqrt5`, Conj. 7.16, imported from the companion GSA/"self-action" paper) was, per verification policy, **not** asserted as an arithmetic result of this paper. It is tested only for internal consistency of the paper's own stated formula `c=sqrt(1+4C)/(2C)` at the golden gate `C=1` — which holds exactly, along with `lambda=2c=sqrt5` and the ladder root `1/phi`. The paper itself labels this a *selection* of scale (Cencov: Fisher metric unique only up to scale), not a derivation, and marks it conjectural — consistent with our finding.

## Untestable / not mechanically re-derivable (documented, 7)

These are cited theorems with infinitary content, open conjectures, or implementation/design guarantees — their *computational shadows* were checked where possible, but the statements themselves cannot be independently re-derived by finite exact computation:

| Claim | Reason | Shadow checked |
|---|---|---|
| Prop. 2.4 (Yuzvinskii–Bowen / LSW): `log Mahler` = topological entropy of the companion toral endomorphism | deep cited dynamical theorem; entropy not independently computable without invoking the same identity | spectral Mahler identity `log M = sum_{|lambda|>1} log|lambda|` verified (claims 8, 9) |
| Thm. 5.2(c) Northcott finiteness (over all algebraic integers of bounded degree & Mahler measure) | infinitary statement | constructive finite-budget count + Mahler coefficient bound verified (claims 39–41) |
| Thm. 8.2 Smyth: `mu_S` is the *smallest* Mahler measure over ALL non-reciprocal integers | infinitary extremal theorem (cited) | the value `mu_S` (root of `x^3-x-1`) and `mu_S>1` verified (claim 53) |
| Conj. 8.4 Lehmer: `mu_L` is the infimum over all non-cyclotomic monic integer polys | open conjecture | the value `mu_L` (root of Lehmer's poly) and `mu_L<mu_S` verified (claim 54) |
| Rem. 8.6 Dobrowolski bound `Mahler>=1+c(loglog d/log d)^3` | asymptotic, constant `c` unspecified in paper | none (no concrete constant to check) |
| Thm. 2.10 completeness of invariant factors as a *complete* similarity invariant (general) | cited structure theorem over the PID `Q[x]` | specific instances verified: non-similar pairs separated (claims 13, 14, 15) |
| Prot. 4.1 protocol invariants (exactness, sole-mutator, SHA-256 witness chain); O2 general complex-`r_2` Fisher probe | design/implementation guarantees, not mathematical derivations; O2 flagged open by the paper | totally-real + two specific complex examples `Q(i)`, `Q(cbrt2)` verified (claims 20, 21, 28) |

## Reproduce

```
py -m pytest tests\2026-06-vector-substrate -v -p no:cacheprovider
```
