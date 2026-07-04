# Verification notes — `2026-07-relational-charge`

**Paper:** *Relational Charge on the Spectral Semiring: Gauge Rigidity,
Coherence Types, and a Reference-Free Parity Floor* (AceTheDactyl / Echo S
Studios), `relational_charge_paper.tex`.

**Engines / method.** Every computational claim was independently re-derived in
exact arithmetic with `sympy` 1.14 (resultants, factorisation over Q, Sturm
root-counts, cyclotomic polynomials, Galois groups) and cross-checked, where
useful, with high-precision `mpmath` (dps=50). The paper's contact-signature
procedure trial-divides the ratio polynomial by every cyclotomic up to a
totient bound; **we instead factor the ratio polynomial over Q and identify the
cyclotomic factors** — a different, completeness-by-construction algorithm — and
additionally reproduce the small signatures with a purely numerical mpmath
engine (form every root ratio, detect roots of unity by their powers). Angles,
relational/absolute charge orders, and Mahler measures are recovered
independently from high-precision roots, never restated from the paper.

Result: **268 passed, 1 xfailed** (the single xfail is a documented ±1
bookkeeping-split issue in the census tally — see "Flagged for human review").

---

## Claims table

| Claim (paper location) | How tested | Status | Notes |
|---|---|---|---|
| Ratio poly `Rat_p` has degree `n^2`, root multiset = all ordered ratios (Def 6.3) | build via resultant; deg check; evaluate `Rat_p` at every root-ratio ≈0; leading coeff `(±p(0))^n` | verified | `test_ratio_object` |
| Two-route agreement: resultant = charpoly(`C_p⊗C_p^{-1}`) for `|p(0)|=1` (ledger G) | independent matrix construction, same signature | verified | β4,L,S6,S8 |
| Diagonal `Phi_1` multiplicity = `n` for squarefree `p` (Lem 6.6) | exact `(x-1)`-adic valuation | verified | |
| Totient bound `2φ(M)²≥M`, unique tight case `M=2` (Lem 4.5, ledger N) | independent sieve to 2·10⁵; cross-check vs `sympy.totient`; prime-power proof steps | verified | |
| Scan bound `2(deg Rat)²`; all quoted bounds (512,20000,2592,8192,41472,3200,131072) | arithmetic + `phi(M)≤d ⇒ M≤2d²` | verified | |
| Contact signatures A–H, P, U, W | factor-and-identify + numeric mpmath cross-check | verified | see list below |
| Gauge-blindness: `x^3±2` identical `{Φ1³,Φ3³}` (ledger A,B) | signatures equal | verified | |
| Shell signature coarser than type: `x^4+5x²±5` identical sig, differing Δ (ledger E,F) | sig equal; relational orders 2 vs 4 | verified | |
| β4, Lehmer relationally inert (Thm 7.9) | `{Φ1⁴}`, `{Φ1¹⁰}` | verified | |
| Twisted-shell non-inert witness `x^4+x²+2=q(x²)` (Ex 7.21, ledger U) | `p=q(x²)`, irreducible, sig `{Φ1⁴,Φ2⁴}`, `p(0)=2` | verified | |
| Cocycle / coboundary / rotation-invariance / trivialization (Prop 3.8 i–iv) | exact Q/Z identities on admissible objects | verified | `test_cocycle_gauge` |
| Class/relational groups cyclic `Z/m` (Prop 3.6) | difference set = `{j/m}`, order `m` | verified | |
| Rigidity `n∈{m,2m}`; anchors `x^m-2→n=m`, `x^3+2→n=2m` (Thm 4.1) | angles from roots; lcd vs difference-lcm | verified | |
| Group drop `x^4+5x²+5`: Eisenstein, purely-imaginary roots, abs Z/4 / rel Z/2, M=5 (Ex 4.6) | irreducibility, sign of `x²`-roots, angles, Mahler | verified | + `K` keeps Z/4 |
| Odd anchors full: `n` odd ⇒ Δ=Z/n (Lem 4.4) | absolute vs relational order | verified | |
| Adams `ψ^k` multiplies relative charge by `k`, `Δ(ψ^kO)=kΔ(O)` (Prop 4.8) | angle scaling in Q/Z; primary projectors | verified | |
| Tensor `lcm` law; CRT kernel cyclic order `gcd`, gen `(m/d,-k/d)`; safe⟺coprime (Prop 4.9) | subgroup arithmetic in Q/Z; surjection kernel | verified | |
| Golden internal relation `φ/φ'=-φ²`, `t_rel=½` (Lem 5.1, ledger J) | symbolic + angle | verified | |
| `q_k=x^{2k}+x^k-1`: real root each sign, angle set `(1/2k)Z/Z`, Δ=Z/2k, `M=φ` (Lem 5.2) | substitution `y²+y-1`; angles; Mahler | verified | k=1..6 |
| Sign twist `T(p)=(-1)^{deg}p(-x)` action; preserves M, Δ, monic; exchanges anchors at odd m (Lem 5.3, ledger L) | `T(x^3+2)=x^3-2`; Mahler/order invariance; coset arithmetic | verified | |
| Parity bit = evenness of `|Δ|` (Cor 5.7) | `½∈(1/m)Z/Z ⇔ m even` | verified | |
| Parity floor attainments: even by `q_{m/2}` (M=φ), odd at 2 by `x^m-2` (Thm 5.5) | relational order + Mahler | verified (attainments) | infimum **values** conditional on (EG), see below |
| Salem certificates: reciprocal + irreducible + trace-Sturm `(1,0,d-1)` for β4,S6,S8,L (Thm 7.10, ledger O) | own trace-poly construction + `count_roots` | verified | |
| `M(β4)=1.72208…`, `M(L)=1.17628…` (Thm 7.9) | Mahler via mpmath; root check | verified | to paper's stated digits |
| Plastic number certificate + inert; Smyth constant `θ0=1.32471…` (Ex 7.19, ledger W) | Sturm; modulus `θ0^{-1/2}`; sig `{Φ1³}` | verified | |
| Mixed ratio degree `(deg p)(deg q)`; no circle-locking, 6 pairs empty sig (Thm 8.2, Rem 8.4, ledger I,Q) | mixed resultant; empty signature; root multiset | verified | |
| Pinned modulus ⇒ inert (Thm 6.15) on Salem/Pisot; necessity via `x^4-2` torsion `{Φ1⁴,Φ2⁴,Φ4⁴}` | dominant-modulus gap + signature | verified (instances) | general proof structural, see below |
| `x^4-x+1`: irreducible, non-recip, 0 real roots; `charpoly(C⊗C)=S6²·(x^4+2x²-x+1)`; `Rat_p={Φ1⁴}`; nested `Rat_{Rat_p}` deg 256 `={Φ1²⁸}`; `Gal=S4`; `M=τ_{S6}` (Ex 6.20, ledger X) | Kronecker factor; nested resultant; `galois_group`; Mahler | verified | full-type showcase |
| `β4⊗β4` non-inert: `(x-1)`-mult 4, `deg gcd(F,F')=7`, 3 distinct positive reals; rational block of 6 (Ex 8.5, ledger S) | Kronecker charpoly, Sturm counts | verified | |
| Census: 729 family, 27 twist-fixed, 378 orbits, 37 Salem, 39 vanish@±1, tally sums 378; all 37 inert `{Φ1¹²}` (Thm 6.13, Rem 6.14, ledger T) | full enumeration + per-instance deg-144 scan | verified (robust invariants) | internal 257/45 split flagged below |
| Census exact `257/45` trace/reducible split (Rem 6.14) | classification | **xfail** | ±1 vs clean count-split 256/46; see below |
| Quadratic (EG) floors (App D): `n=-1` ⇒ M≥φ (=φ at |t|=1, by `x²-x-1`); `n=+1` ⇒ M≥φ²; no quadratic M in `(1,φ)`; `φ=(1+√5)/2`, `φ²=(3+√5)/2` (ledger R) | symbolic monotonicity + integer `(t,n)` sweep | verified | |

**Contact signatures verified (exact factor-route + numeric mpmath where feasible):**
A `x³-2 {Φ1³,Φ3³}`; B `x³+2 {Φ1³,Φ3³}`; C `x⁴-2 {Φ1⁴,Φ2⁴,Φ4⁴}`;
D `q2 {Φ1⁴,Φ2⁴}`; E `x⁴+5x²+5 {Φ1⁴,Φ2⁴}`; F `K {Φ1⁴,Φ2⁴}`; G `β4 {Φ1⁴}`;
H `Lehmer {Φ1¹⁰}`; P `S6 {Φ1⁶}, S8 {Φ1⁸}`; U `x²+x+2 {Φ1²}, x⁴+x²+2 {Φ1⁴,Φ2⁴}`;
W `plastic {Φ1³}`; X `x⁴-x+1 {Φ1⁴}` and nested `{Φ1²⁸}`; I/Q the six mixed
Salem pairs all empty.

---

## VERIFIED (summary)

All 84 distinct computational claims above reproduce **exactly** from the
paper's own premises: the ratio-object construction and its two independent
routes; the completeness/totient bound (unique tight case `M=2`); every
cyclotomic-contact signature in Appendix A that is a finite computation
(A–I, P, Q, U, W, X); rigidity and both anchors; the descent laws (Adams, `lcm`,
CRT kernel); the golden relation and the `q_k` / sign-twist parity machinery;
the Salem certificates and Mahler values (β4, Lehmer, S6, S8, plastic, x⁴-x+1);
the `x⁴-x+1` full-type showcase (Kronecker factorisation, nested `ν`-scan,
Galois `S4`, `M=τ_{S6}`); the `β4⊗β4` non-inert structure; the degree-12 census
combinatorics and the inertness of all 37 certified Salem instances; and the
quadratic emission-gap floors of Appendix D.

## FAILED / flagged for human review

**Census rejection tally, internal `257/45` split (Rem 6.14, `rem:censusscope`)
— xfail, strict.** The paper reports the 341 non-Salem twist-classes as
`39` (vanish at ±1) `/ 257` (fail trace-Sturm pattern `(1,0,5)`) `/ 45`
(reducible) `/ 0` (only-trace-irreducibility) `/ 37` Salem. An independent,
purely count-based Sturm-pattern test gives **256 / 46** instead of 257 / 45.
The two differ by exactly one polynomial:

  `x¹²−x¹¹−x¹⁰−x⁹−x⁷−x⁶−x⁵−x³−x²−x+1`

which is **reducible** — it has the cyclotomic factor `Φ₁₀ = x⁴−x³+x²−x+1`
(unit-circle roots) times a genuine Salem octic — yet whose trace polynomial has
root pattern `(1,0,5)` by count (the octic supplies the one root `>2`). Whether
this single edge case is bucketed under "trace-fail" or "reducible" is a
labelling convention; no clean mathematical rule singles it out (all 46
pattern-passing reducibles have reducible trace polynomials, so requiring an
irreducible trace poly moves *all* of them, giving 302 / 0, not 257 / 45).

**This is not a mathematical error.** Every robust invariant matches exactly:
`729` family / `27` twist-fixed / `378` orbits / `37` Salem / `39` vanish /
combined trace+reducible `= 302` (= 257+45 = 256+46) / sum `378`, and the
headline result — **all 37 certified Salem instances are relationally inert
`{Φ1¹²}`** — is fully reproduced. The exact `257` assertion is recorded as a
strict `xfail` so that any future reclassification is caught.

## UNTESTABLE (documented, not a finite computation)

- **(EG) emission gap `M∈{1}∪[φ,∞)` beyond degree 2** — an inherited,
  explicitly *conditional* hypothesis (`Prop 2.3` of the companion note `[CMC]`).
  Its quadratic case is elementary and **is** verified here (App D); the general
  case is not provable within this paper (requires `[CMC]`, DOI-archived). The
  paper tags this `[C]/[Pl]` and confines its consequence to two floor *values*.
- **Floor *values* as `φ` / `2` (Thm 5.5 lower bounds)** — the *attainments*
  (`M=φ` by `q_k`, `M=2` by `x^m-2`) are verified; identifying the infimum *as*
  `φ` rather than something in `[θ0,φ)` is conditional on (EG) and is a
  search/limit claim, not a finite check.
- **Modulus-pinning theorem (Thm 6.15) general proof** — the Galois-transport
  argument (transitive action moving a torsion ratio onto the pinned root) is
  structural. Its *consequences* are verified on 46 Salem/Pisot instances and 6
  mixed pairs; the universal statements Cor 6.16 (every Salem inert), Cor 6.17
  (no two Salem circle-lock), Cor 6.18 (Pisot) are general theorems corroborated,
  not exhausted, by computation.
- **Structure of partial coherence (Thm 4.3), Salem geometry (Cor 4.7)** —
  qualitative decomposition statements (rational block + mirrored offset
  classes); verified in spirit on the computed types, not as a single assertion.
- **Gauge/torsor interpretation and citation-provenance prose** — interpretive,
  not mathematical.
