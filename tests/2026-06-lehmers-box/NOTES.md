# Verification notes — *Lehmer's Box*

**Paper:** *Lehmer's Box: A Golden Floor and an Angle Lattice that Confine
Spectral Emission Away from Salem Numbers — Without Resolving Lehmer's Problem*
(AceTheDactyl / Echo S Studios). Source: `papers/2026-06-lehmers-box/lehmers_box.tex`.

Every checkable claim is **independently re-derived** (sympy exact / mpmath at
`dps>=40` / numpy linear algebra) and then compared to the paper's stated value —
no paper number is compared against itself. Shared re-derivation tools live in
`_helpers.py` (Mahler measure by the product formula, Salem classification from
root moduli, exact `Q(sqrt5)` sign). `conftest.py` only puts `_helpers` on the
path under pytest's importlib mode.

Final run: **45 passed, 1 xfailed, 0 failed, 0 errors.**

## Claims table

| # | Claim (paper location) | How tested | Status |
|---|---|---|---|
| 1 | `phi=(1+sqrt5)/2=1.6180339887…` (abstract) | mpmath, derive from sqrt5 | verified |
| 2 | `tau=(-1+sqrt5)/2=phi^-1` (Def 2.9) | tau·phi=1 exactly | verified |
| 3 | `mu_S=1.3247179572…` real root of `x^3-x-1`, smallest Pisot (§1.2, Lem 2.7) | Newton root + count real roots | verified |
| 4 | `log mu_S=0.281200`, `log phi=0.481212` nats (Cor 4.9) | mpmath log | verified |
| 5 | `Mah(L)=1.1762808182…` for Lehmer `L` (§1.2) | product-formula Mahler measure | verified |
| 6 | `L` is Salem, `Mah(L)<mu_S` (§1.2) | palindrome + root pattern (1,1,8) | verified |
| 7 | `beta_4=1.7220838057…`, dominant root of `x^4-x^3-x^2-x+1`, `>phi` (Cor 5.5) | mpmath roots | verified |
| 8 | `beta_4` is a degree-4 Salem, `Mah=beta_4` (Cor 5.5) | root pattern (1,1,2) | verified |
| 9 | Empty quadratic strip `(1,phi)`; 625 quadratics `\|b\|,\|c\|<=12` (Lem 3.1) | exhaustive scan, exact Mahler | verified |
| 10 | Minimiser `phi` attained only at disc 5 (`x^2±x-1`) (Lem 3.1) | argmin over scan | verified |
| 11 | `c=-1`: `Mah=(\|b\|+sqrt(b^2+4))/2`, min at `\|b\|=1` (Lem 3.1) | formula + monotonicity | verified |
| 12 | `c=1,\|b\|>=3`: dominant root `>= (3+sqrt5)/2 = phi^2` (Lem 3.1) | mpmath | verified |
| 13 | Operators preserve floor; semigroup `{phi,2,3,5,phi^4,beta^2}` (Lem 3.2, Rem 3.3) | product/square stay `>=phi^2` | verified |
| 14 | `gap = phi^4` (roots of `x^2-7x+1`) (Def 2.9) | `phi^4+phi^-4=7`, Mahler=`phi^4` | verified |
| 15 | Catalog arguments in `(pi/2)Z` (Lem 4.1) | each root real or pure imaginary | verified |
| 16 | `K`: `x^2=(-5±3sqrt5)/2 = 0.854 / -5.854` (Lem 4.1) | solve `y^2+5y-5=0` | verified |
| 17 | `(pi/2)Z` closed under `+` & doubling, not halving (Lem 4.2) | exact model `Z/4`; `pi/4` off lattice | verified |
| 18 | On-circle `+ (pi/2)Z ⇒ mu_4` fourth roots (Lem 4.3) | `e^{ik pi/2}∈{1,i,-1,-i}`, `z^4=1` | verified |
| 19 | Salem on-circle conjugates not roots of unity (Lem 4.4) | `L` irreducible, `≠Φ_11,Φ_22`; args off lattice | verified |
| 20 | Emission ⊆ Box: catalog in box (Thm 4.6) | floor + lattice membership | verified |
| 21 | Lehmer outside box on both walls (Thm 4.7) | `Mah<phi` and off-lattice on-circle | verified |
| 22 | A Salem above floor (`beta_4`) still lattice-excluded (Thm 4.7 / Fig 1) | `beta_4>phi` yet off-lattice pair | verified |
| 23 | Def 2.1: product form = Jensen integral form | both forms, no on-circle roots | verified |
| 24 | Lem 2.3(i): `Mah(p)>=1` monic | random monic integer polys | verified |
| 25 | Lem 2.3(ii): `Mah(pq)=Mah(p)Mah(q)` | convolve, compare | verified |
| 26 | Lem 2.3(iii): `Mah(p^[2])=Mah(p)^2` | build `p^[2]` from `p(x)p(-x)` | verified |
| 27 | Kronecker: `Mah=1` iff cyclotomic (Lem 2.4) | `Mah(Φ_n)=1`, `Mah(x^2-x-1)>1` | verified |
| 28 | Salem signature `(2,m-1)`, trace-form `(m+1,m-1)` (Prop 5.1) | embeddings + trace-form eigen-signs (`beta_4`,`L`) | verified |
| 29 | `[K:Q]=16` (Lem 5.2) | minpoly of `sqrt2+sqrt3+5^{1/4}` has degree 16 | verified |
| 30 | Real catalog differences lie in `K` (Lem 5.2) | `phi-psi=sqrt5=(5^{1/4})^2∈K` etc. | verified |
| 31 | 27-subfield signature census of `K` (Thm 5.3) | independent Galois comp. on `K(i)` | verified |
| 32 | Backing note "all **27 subgroups** enumerated" (Thm 5.3) | full subgroup count of `G` | **xfail — flagged** |
| 33 | `beta_4>phi` by exact sign of `m_{beta_4}(phi)=(1-sqrt5)/2<0` in `Q(sqrt5)` (Cor 5.5) | Galois-conjugate `(a,b)` + sign rule | verified |
| 34 | Only Salem-bearing subfields are the four `(2,1)` quartics (Cor 5.5) | filter census by signature shape | verified |
| 35 | Complex place off circle: `\|i beta\|=2.4195`, `\|5^{1/4}i\|=1.4953` (§5) | mpmath | verified |
| 36 | Trace-down identity `m_theta(x)=x^m T(x+1/x)` (Def 6.1) | symbolic reconstruction | verified |
| 37 | Salem ⇔ flip-straddle (one `T`-root `>2`, rest in `(-2,2)`) (Lem 6.2) | trace-down of `L`, `beta_4` | verified |
| 38 | Flip boundary `D=t^2-4` sign change at `t=±2` (Lem 6.2) | sign of `t^2-4` | verified |
| 39 | `rho(z)=z+1/z` sends `mu_4 → {2,0,-2}` (Prop 6.3) | exact; Salem straddle off this lattice | verified |
| 40 | Circulant eigenvalues `sum_k c_k ω^{jk}`; no Salem factor (Prop 7.1) | np vs formula; factor char-poly | verified |
| 41 | Shoda: commutator ⇒ trace 0; traceless Lehmer carrier exists (Lem 7.2) | `tr[X,Y]=0`; `L \| charpoly` of traceless block | verified (necessary dir.) |
| 42 | Self-action `spec(ad_R)={0,±sqrt5}` = difference set (Lem 7.4) | Kronecker `R⊗I−I⊗R^T` eigenvalues | verified |
| 43 | Guard ladder FORCED / FORCED_ABOVE_FLOOR / INVALID_CLOSURE (Prop 7.5) | `validate()` on catalog, `beta_4`, `L` | verified |
| 44 | Guard `beta<phi` by exact sign in `Q(sqrt5)`, no float (Prop 7.5) | same-sign / `a^2` vs `5b^2` rule | verified |

## Verified (43 claims)
Rows 1–31 and 33–44 above all reproduce exactly from the paper's own premises.
Highlights of the non-trivial independent re-derivations:

- **Empty strip (Lem 3.1):** the full 625-quadratic scan finds *no* Mahler
  measure in `(1,phi)`; the unique minimiser above 1 is exactly `phi` at
  discriminant 5 (`x^2±x-1`).
- **27-subfield census (Thm 5.3):** built `G=Gal(K(i)/Q)` from scratch as
  `C2×C2×D4` (order 32, associativity/inverse-checked), enumerated all subgroups,
  restricted to the 27 containing complex conjugation (= subfields of `K`), and
  computed each signature `(deg,r1,r2)` via the real-embedding count
  `r1=#{gH: g^{-1}cg∈H}`. The result matches the paper's table
  `{(1,1,0):1,(2,2,0):7,(4,2,1):4,(4,4,0):7,(8,4,2):6,(8,8,0):1,(16,8,4):1}`
  **exactly**, total 27, with `K` itself at `(8,4)` and four Salem `(2,1)` quartics.
- **Exact `Q(sqrt5)` guard (Cor 5.5 / Prop 7.5):** `m_{beta_4}(phi)=(1−sqrt5)/2`
  (negative ⇒ `beta_4>phi`, FORCED_ABOVE_FLOOR) and `L(phi)=(133+59 sqrt5)/2`
  (positive ⇒ `beta_L<phi`, INVALID_CLOSURE) — both decided with no floating point.

## Flagged for human review (1 claim)
- **Row 32 — backing note "all 27 subgroups enumerated" (Thm 5.3).**
  The **theorem is correct**: `K` has exactly 27 subfields with the stated
  signature census. But the parenthetical wording "all 27 **subgroups**
  enumerated" is imprecise: the Galois group `G≅C2×C2×D4` has **158 subgroups in
  total** (independently enumerated: 1 of order 1, 23 of order 2, 67 of order 4,
  51 of order 8, 15 of order 16, 1 of order 32). The number 27 is the count of
  subgroups **containing complex conjugation** `c` (equivalently `Gal(K(i)/K)`),
  which is precisely the set in Galois correspondence with the 27 subfields of
  `K`. Recommended reading: "the 27 subgroups fixing a subfield of `K`."
  Recorded as `xfail` in `test_signature_face.py::test_literal_total_subgroup_count_is_27`;
  the true total 158 is pinned by `test_total_subgroup_count_is_158`. This is a
  wording nuance in a small-font backing note, **not** a mathematical error in
  any derivation.

## Untestable / argument-only (documented, not mechanically checkable)
These carry no direct closed-form to re-derive; where possible a symbolic sanity
check or a concrete instance was added (noted), but the general statement is a
proof/logic claim:

1. **Smyth's bound (Lem 2.7)** `Mah(p)>=mu_S` for all non-reciprocal integer `p`
   — a cited deep theorem (Smyth 1971); only `mu_S` itself and Lehmer's `<mu_S`
   placement are checked, not the universal bound.
2. **Shoda converse (Lem 7.2)** "every traceless matrix *is* a commutator" — a
   cited theorem; only the necessary direction (`tr[X,Y]=0`) and eligibility of a
   Lehmer carrier are checked.
3. **Full induction Emission ⊆ Box over all of `S` (Thm 4.6/4.7)** — the operator
   base cases, closure lemmas and representative instances are verified; the
   word-length induction itself is an argument.
4. **Uniform self-action bound (Thm 5.7)** — assembles the difference-spectrum
   and census facts (both verified); the "uniform at every size" quantifier is a
   structural argument, not a finite computation.
5. **Logical position of the box theorem (Prop 8.1 (i)–(iv))** and **the
   sidestep (Cor 8.2)** — meta-logical statements ("(B) is consistent with both
   truth values of (P)", "(B) is unconditional"); not mechanically decidable.
6. **Open fronts (Rem 8.4)** — the free unbounded commutator and Lehmer's problem
   `(P)` are explicitly left open; nothing to verify.
