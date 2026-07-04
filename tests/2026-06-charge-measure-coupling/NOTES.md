# Verification notes: `2026-06-charge-measure-coupling`

**Paper:** *The Charge--Measure Coupling on a Spectral Semiring: Cyclicity,
Composition, and a Parity-Graded Mahler Floor* (v4, revised), AceTheDactyl /
Echo S Studios.
**Source:** `charge-measure-coupling-whitepaper-v4.tex`

## Method

Two characters were re-implemented from scratch in `_cmc_helpers.py` and used to
re-derive every numeric claim independently before comparing to the paper:

* **Character I (Mahler measure)** `M(p) = |lead| * prod_{|root|>1}|root|`, via
  `mpmath.polyroots` at 50 decimal digits.
* **Character II (charge group)** the least `n` with `alpha^n in R_{>0}` for every
  root (= least common denominator of `arg/(2pi)`); `None` = charge-inadmissible
  (irrational conjugate angle). Detection uses a `1e-12` tolerance against a
  `~1e-46` true-rational margin at 50 digits.

Symbolic identities use `sympy` (exact); finite-group claims are exact integer
computations in `Z/L`; two "computed" claims are replayed over a reduced but
honestly-documented finite window.

**No cross-domain / externally-seeded constants appear in this paper.** phi, the
plastic number mu_S, Lehmer's tau and phi^4 are all genuine roots of the integer
polynomials the paper names, so every constant is checked as an arithmetic
consequence of the paper's own premises (guardrail 4b not triggered; no xfail).

## Claims table

| # | Claim (paper ref) | How tested | Status |
|---|---|---|---|
| C1 | phi = (1+sqrt5)/2 is root of x^2-x-1 (Sec 1.4) | sympy exact | verified |
| C2 | phi' = -1/phi, phi*phi'=-1 (Sec 4.3) | sympy exact | verified |
| C3 | phi'<0, at argument pi (Lem 4.5) | sympy sign | verified |
| C4 | sqrt5 = phi+1/phi = phi-phi' (Sec 4.3) | sympy exact | verified |
| C5 | phi^2=(3+sqrt5)/2=phi+1~2.618 (Lem 2.6) | sympy + mpmath | verified |
| C6 | phi^4=(7+3sqrt5)/2=3phi+2~6.854 (Thm 6.5) | sympy + mpmath | verified |
| C7 | 2cos72=phi-1, 2cos144=-phi (Prop 6.1) | sympy exact | verified |
| C8 | pentagon cosines = roots of x^2+x-1, Galois-conj (Prop 6.1) | sympy exact | verified |
| C9 | mu_S=1.3247, real root of x^3-x-1 (Lem 2.5) | mpmath findroot | verified |
| C10 | mu_S is Pisot (conjugates inside unit disk) (Lem 2.5) | mpmath roots | verified |
| C11 | Lehmer tau=1.17628 in (1,phi) (Prop 8.2) | indep. Mahler | verified |
| C12 | beta_4 Mahler=1.72208 in (phi,2) (ledger F) | indep. Mahler | verified |
| C13 | phi=1.61803 (ledger H) | mpmath | verified |
| M1 | M(x^n-2)=2, n=3..7 (Thm 3.2) | indep. root product | verified |
| M2 | M(q_k)=phi, q_k=x^2k+x^k-1, k=2..5 (Thm 4.1) | indep. root product | verified |
| M3 | x^4+x^2-1=(x^2+phi)(x^2+phi') (Thm 4.1) | sympy expand | verified |
| M4 | M(beta_4)=1.72208... precise (ledger F) | indep. Mahler | verified |
| M5 | M(Lehmer)=1.17628... precise (Prop 8.2) | indep. Mahler | verified |
| M6 | M(L(x)(x-1))=tau, M(x-1)=1 (Prop 8.2) | indep. Mahler | verified |
| M7 | M((x-1)(x-2))=2 (ledger G) | indep. Mahler | verified |
| M8 | M(x^4-x^3+6x^2+4x+1)=phi^4 (ledger I) | indep. Mahler | verified |
| M9 | M(Phi_5*(x^2-3x+1))=phi^2 (ledger K) | indep. Mahler | verified |
| M10 | M(x^5-m)=m, m=2..6 (Thm 6.7) | indep. Mahler | verified |
| M11 | M((x^3-2)(x)(x^4-2))=128=2^7 (ledger D) | indep. product roots | verified |
| G1 | x^n-2: charge Z/n, all n charges (Thm 3.2) | indep. charge group | verified |
| G2 | x^4-2 cyclic Z/4, not Z/2xZ/2 (ledger B) | indep. charge group | verified |
| G3 | q_k: charge Z/2k, all 2k charges (Thm 4.1) | indep. charge group | verified |
| G4 | beta_4 charge-inadmissible (Lem 8.1) | indep. charge group=None | verified |
| G5 | Lehmer commutator inadmissible (ledger L) | indep. charge group=None | verified |
| G6 | pentagon quartic charge Z/5, irred, non-recip (ledger I) | charge + sympy irreducible | verified |
| G7 | Phi_5*(x^2-3x+1) charge Z/5, reciprocal (ledger K) | charge + palindrome | verified |
| G8 | x^5-2 charge Z/5, non-reciprocal (Thm 6.7) | charge + palindrome | verified |
| G9 | (x^3-2)(x)(x^4-2) charge Z/12=lcm(3,4) (Thm 3.5) | indep. charge group | verified |
| T1 | phi (x) phi charpoly=(x+1)^2(x^2-3x+1) (ledger M) | sympy from products | verified |
| T2 | M(phi (x) phi)=phi^2 (tropical), NOT phi^4 (Prop 3.3) | indep. Mahler | verified |
| T3 | off-circle tensor: factored law agrees, =2^7 (Thm 3.7) | indep. moduli>1 | verified |
| T4 | tensor charpoly = x^12-128 (ledger D) | numeric charpoly | verified |
| T5 | M(psi^k A)=M(A)^k (Sec 2.1 table) | indep. power roots | verified |
| T6 | x^6-2: psi^3->Z/2, psi^2->Z/3 (Thm 3.4) | indep. charge of powers | verified |
| T7 | charges add (sumset) under tensor (Sec 2.1/Thm 3.5) | indep. charges vs sumset | verified |
| S1 | every finite subgroup of Q/Z (i.e. of Z/L) is cyclic (Thm 3.1) | exact Z/L closure | verified |
| S2 | Z/p x Z/p is non-cyclic (unrealizable) (Thm 3.1 cor) | exact element orders | verified |
| S3 | CRT iso + Adams primary projectors psi^{n/p^e} (Thm 3.4) | exact Z/n maps | verified |
| S4 | (1/m)Z+(1/k)Z=(1/lcm)Z, order lcm (Thm 3.5) | exact Z/L closure | verified |
| S5 | safe composition <=> coprime; ker order=gcd (Thm 3.7) | exact | verified |
| F1 | real reciprocal unit: trace>=3 => r>=phi^2 (Lem 2.6) | sympy solve | verified |
| F2 | totally-positive gap: no M in (1,2), min>1 is 2 (Thm 4.6) | finite scan deg<=4,|c|<=6 | verified (window) |
| F3 | pentagon construction expands to minimiser (Thm 6.5) | sympy expand | verified |
| F4 | minimiser: +-72 pair mod phi^2, +-144 pair mod phi^-2 (Thm 6.5) | mpmath roots | verified |
| F5 | charge-Z/5 quartics: distinct M = {1, phi^4} (ledger J) | scan |c|<=6 + confirm | verified (window) |
| F6 | realification 2^{1/5}=1.1487 < mu_S=1.3247 (Thm 6.4) | mpmath | verified |
| F7 | emission gap M in {1}U[phi,inf) on all constructions (Prop 2.2) | indep. Mahler | verified (sanity) |
| X1 | beta_4 is Salem (1 out,1 in,2 on-circle) (Lem 8.1) | mpmath root moduli | verified |
| X2 | Lehmer Salem, deg 10 (1 out,1 in,8 on-circle) (Lem 8.1) | mpmath root moduli | verified |
| X3 | Lehmer polynomial has trace -1 (Prop 8.2) | sympy coeff | verified |
| X4 | L(x)(x-1)=x^11-x^9-x^8+x^3+x^2-1, trace 0 (Prop 8.2) | sympy expand | verified |
| X5 | companion matrix trace 0, integer, carries Lehmer's number (Prop 8.2) | numpy eig | verified |
| X6 | commutator: M=tau in (1,phi), charge bottom (Prop 8.2 i,iii) | indep. Mahler+charge | verified |

**58 distinct claims tested; all 58 verified. 0 failed / flagged / xfail.**

## VERIFIED (highlights)

* The whole consolidation ledger (App A, entries A--M) reproduces exactly in
  independent 50-digit arithmetic: charge groups, all-charge attainment, and
  Mahler measures 2, phi, phi^2, phi^4, 128, tau, 1.72208 all match.
* The v4 tensor-law *correction* is confirmed: `phi (x) phi` gives eigenvalues
  `{phi^2,-1,-1,phi'^2}`, charpoly `(x+1)^2(x^2-3x+1)`, and tropical measure
  `phi^2` -- NOT the retracted factored value `phi^4`.
* The parity mechanism checks out: the even construction `q_k` attains `phi` on a
  `Z/2k` lattice (contains pi); `x^n-2` gives the odd realised floor 2; the
  pi-ray obstruction `phi'` at argument pi is confirmed.
* The `Z/5` worked case checks out end to end: pentagon cosines, the Galois-coupled
  quartic construction expanding to `x^4-x^3+6x^2+4x+1` with measure `phi^4`, and
  the scan finding only `{1, phi^4}`.
* The commutator escape is confirmed concretely: `L(x)(x-1)` is a trace-zero
  integer polynomial whose companion matrix carries Lehmer's number (measure
  `1.17628` below the floor `phi`) and is charge-inadmissible.

## FAILED / FLAGGED for human review

None. No claim contradicted its own premises; no tolerance was loosened to hide a
discrepancy. `xfail_strict=false` in the shared `pytest.ini`, but no xfail markers
were needed.

## UNTESTABLE (documented, not mechanically decidable here)

These are cited external theorems or explicitly open/conjectural statements. The
paper itself tags them `[plausible]`/`[computed]` or attributes them to the
literature; they cannot be settled by a finite computation and are **not** treated
as errors.

1. **Kronecker's theorem** (Prop 2.1): `M=1 iff every root is a root of unity` --
   cited (Kronecker 1857). The `M=1` witnesses (Phi_5) are checked; the iff in full
   generality is not.
2. **Smyth's bound** (Lem 2.5): general non-reciprocal `M >= mu_S` -- cited
   (Smyth 1971). mu_S's value/Pisot property and specific witnesses are checked;
   the universal inequality is not re-proved.
3. **Universal emission floor phi** (Prop 2.2): tagged `[plausible]`; open. The
   *forced* content is only `M >= mu_S`; the value phi is conjectural. Verified as a
   sanity property on all constructions, not as a theorem.
4. **Even floor lower bound `mu(even) >= phi`** (Thm 4.1): `[plausible]`, inherits
   from (3); open.
5. **Odd realised floor `mu(odd)=2` for all odd n** (Thm 4.1): `[forced]` only for
   Z/3, `[computed]` beyond; the general equality is open. Finite windows verified.
6. **Salem gap `(phi,2)` empty for all odd n and all degrees** (Thm 8.3):
   `[forced]` for Z/3, `[plausible]` in general; open. Supported by the beta_4
   witness and the totally-positive / pentagon scans.
7. **Commutator realizability** (Prop 8.2): `trace-zero <=> commutator` and the
   integer refinement -- cited (Albert--Muckenhoupt 1957, Laffey--Reams 1994,
   Shoda 1936). The concrete witness (trace-zero integer companion carrying
   Lehmer) is fully checked; the general theorems are inputs.
8. **"mu_S is the smallest Pisot number"** (Lem 2.5): minimality is a cited fact;
   only "mu_S is Pisot" is verified.
9. **"beta_4 is the minimal degree-4 Salem number"** (ledger F): minimality not
   searched; "beta_4 is a Salem number with M=1.72208" is verified.
10. **Z/5 residual full-window emptiness** (Prop 6.8): `[computed]` over
    quartics |c|<=10, quintics |c|<=4, sextics |c|<=3. Only the quartic window
    (|c|<=6) is replayed here; the quintic/sextic windows are not.

## Reproduce

```
py -m pytest "C:\Users\acead\projects\Echo-S-Research\tests\2026-06-charge-measure-coupling" -v -p no:cacheprovider
```

Files: `test_constants.py`, `test_mahler_measure.py`, `test_charge_groups.py`,
`test_tensor_adams.py`, `test_group_structure.py`, `test_floor_bounds.py`,
`test_salem_commutator.py`, helper `_cmc_helpers.py`.
