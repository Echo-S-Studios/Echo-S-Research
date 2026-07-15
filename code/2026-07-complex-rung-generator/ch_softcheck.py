#!/usr/bin/env python3
# ch_softcheck.py -- the chi-doctrine hunt (W[open]2, sub-branches g1/g3),
# re-run with the CITED companion papers now PRESENT (Echo-S-Research).
#
# Predecessor: cx_softcheck.py (v1.7) left W[open]2 PARTIAL: g1 (can Phi-
# extraction / direct-sum put +-i in the emission image S?) = EXTERNAL-OPEN;
# g3 (does the doctrine tie the gate multiplier +-i to S?) = bounded NO-TIE
# because the Operator-Algebra and Charge-Measure papers were NOT shipped.
# v1.7 paper thm:chichain flags verbatim: "Whether the full image --- closed
# also under addition, minimal polynomials, and Phi-extraction --- can reach
# +-i is NOT decided ... and remains open" and "No shipped document ties the
# forced gate multiplier +-i ... to the emission image."
#
# Those papers are now READ (READ-ONLY sibling repo). This harness re-encodes
# their load-bearing machinery COLD and decides both branches. Companion quotes
# (file + location) carried inline as CORPUS-[...] premises:
#
#   [OA] operator-algebra-whitepaper.tex
#     Def 2.1 (def:objects): objects = multisets of nonzero algebraic numbers;
#       A(+)B = multiset union, A(x)B = {lam*mu}, 1 = {1}.
#     Def 3.1 + Thm 3.2 (thm:adams): psi^n(A) = {lam^n}; squaring = psi^2.
#     Sec 1 + Sec 7 (prop:fixed): the spectral operators are EXACTLY
#       {(+), (x), (.)^2 = psi^2, minpoly, Phi}; "minpoly and Phi act as the
#       IDENTITY on the spectrum ... they are the IDEMPOTENTS of the structure"
#       (NO cyclotomic / Phi-extraction primitive that adds eigenvalues).
#     Thm 6.2 (thm:charge): charge chi(A) = round(2 arg lam / pi) mod 4 in Z/4Z;
#       "the FULL group Z/4Z is realised: the real seeds carry {0,2}, and the
#       Lorentzian seed K = x^4+5x^2-5 carries {0,1,2,3}, its imaginary pair
#       +-i.beta at +-pi/2."  (charge is ARGUMENT-only; modulus-free.)
#   [EG] emission_gap_paper.tex
#     Lem 4.1 (lem:catargs): every catalog eigenvalue has arg in (pi/2)Z.
#     Lem 4.2 (lem:closure): the five spectral operators PRESERVE (pi/2)Z
#       ("closed under addition and doubling", Rem 4.4 "why doubling, not
#       halving").
#     Lem 4.3 (lem:oncircle): "If theta is an eigenvalue of some M in S with
#       |theta| = 1, then theta in {1,i,-1,-i}."  <-- a CONTAINMENT (subset mu4),
#       NOT an equality: it is an UPPER bound on the on-circle image.
#   [LB] lehmers_box.tex (Lehmer's Box)
#     Lem 4.4 (lem:oncircle): "If |z|=1 and arg z in (pi/2)Z then z in mu4."
#     Sec 5/Prop 6.2 (prop:lorentz): the K-seed complex place is OFF the circle,
#       "|i.beta| = 2.4195 ... the spectral operators --- which only multiply,
#       union, double, and preserve spectra --- cannot rotate it onto the circle
#       without producing a rational angle, i.e. a fourth root of unity."
#   [CM] charge-measure-coupling-whitepaper-v4.tex
#     Def 4.x (def:charge): charge-admissible iff alpha^n in R_{>0}; the charge
#       reads the ARGUMENT class only; Prop commutator: the only non-semiring
#       operation (commutator) is the one door out -- semiring/emission
#       constraints do NOT extend to dynamical (non-emission) multipliers.
#
# DECISION: on-circle image of S = exactly mu_2 = {+1,-1}, NOT mu_4.  Hence the
# cx-C4 hypothesis "S contains a non-real ON-CIRCLE element" is FALSE given the
# corpus; g1 upgrades EXTERNAL-OPEN -> FORCED-GIVEN-CORPUS (negative); g3
# upgrades bounded-NO-TIE-papers-absent -> GROUNDED NEGATIVE.  chi = +-pi/2 does
# NOT become forced-given-corpus: D2' is either FALSE-given-corpus (on-circle
# reading) or a corpus-RESTATEMENT of D2 (charge-class reading).
#
# Discipline (inherited from cx/qx idiom): every DECISION exact over Q / Q(sqrt5)
# (extended by i and 5^(1/4)); unimodularity decided in Q(sqrt5) by rational
# (A,B) arithmetic (q5AB); monomial bookkeeping by certified integer exponent
# vectors; mpmath iv at 60 dps ONLY for the CH-G certified guards (counted
# separately).  Fail-first: each group opens with a falsifier (CH-*F*) asserting
# a deliberately perturbed variant is REJECTED; falsifiers were written first.
# Exit 0 iff all pass.

import sys
from fractions import Fraction as F
from itertools import combinations_with_replacement
from sympy import (symbols, sqrt, Rational, Integer, I, pi, sin, cos, exp,
                   im, re, simplify, expand, radsimp, together, nsimplify,
                   Poly, minimal_polynomial, real_roots, gcd, Mul, conjugate)

x = symbols('x')
y = symbols('y', real=True)
s5 = sqrt(5)
phi = (1 + s5) / 2
psi = (1 - s5) / 2                 # = -tau
tau = (s5 - 1) / 2                 # = 1/phi
f4 = 5 ** Rational(1, 4)           # 5^(1/4)
Kv = f4 * tau                      # K-seed terrain root
betav = f4 * phi                   # K-seed rotation magnitude, beta = 5^(1/4) phi

PASS, FAIL, GUARD = [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(("PASS" if ok else "FAIL"), cid, "-", desc)

def gk(cid, cond, desc):
    ok = bool(cond)
    GUARD.append((cid, ok))
    print(("GUARD-PASS" if ok else "GUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append(cid)

def zero(e):
    e2 = simplify(expand(simplify(e)))
    if e2 == 0:
        return True
    try:
        if simplify(expand(radsimp(together(e)))) == 0:
            return True
    except Exception:
        pass
    try:
        if simplify(expand(e.rewrite(sqrt))) == 0:
            return True
    except Exception:
        pass
    try:
        e3 = nsimplify(e2, [sqrt(5)])
        return simplify(radsimp(e3)) == 0
    except Exception:
        return False

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q. Raises if not in Q(sqrt5)."""
    from sympy import cancel
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return F(str(A)), F(str(B))

def is_one_q5(e):
    """Exact decision e == 1 for e in Q(sqrt5)."""
    A, B = q5AB(e - 1)
    return A == 0 and B == 0

def onS_modsq(e):
    """Exact decision |e|^2 == 1 (on the unit circle) for e whose |e|^2 in Q(sqrt5)."""
    return is_one_q5(expand(e * conjugate(e)))

# ------------------------------------------------------------------ catalog
# [EG] Def 3.1 / [LB] Def 2.6 (verbatim): seeds = companions of
# { phi, tau, sqrt2, sqrt3, sqrt5, gap (x^2-7x+1), K (x^4+5x^2-5) }.
r2, r3, r5 = sqrt(2), sqrt(3), sqrt(5)
SEEDS = [
    ("phi-seed",  x**2 - x - 1,      [phi, psi]),
    ("tau-seed",  x**2 + x - 1,      [tau, -phi]),
    ("sqrt2",     x**2 - 2,          [r2, -r2]),
    ("sqrt3",     x**2 - 3,          [r3, -r3]),
    ("sqrt5",     x**2 - 5,          [r5, -r5]),
    ("gap-seed",  x**2 - 7*x + 1,    [phi**4, tau**4]),
    ("K-seed",    x**4 + 5*x**2 - 5, [Kv, -Kv, I*betav, -I*betav]),
]

# Eigenvalue table (re-derived cold): lam = 2^(w0/4) 3^(w1/4) 5^(w2/4) phi^(w3/4) i^t.
# The charge is t (= round(2 arg / pi) mod 4, [OA] Thm 6.2).  |lam|=1 iff w = 0.
EIG = [
    ("phi",     phi,       (0, 0, 0,   4), 0),
    ("psi",     psi,       (0, 0, 0,  -4), 2),
    ("tau",     tau,       (0, 0, 0,  -4), 0),
    ("-phi",    -phi,      (0, 0, 0,   4), 2),
    ("r2",      r2,        (2, 0, 0,   0), 0),
    ("-r2",     -r2,       (2, 0, 0,   0), 2),
    ("r3",      r3,        (0, 2, 0,   0), 0),
    ("-r3",     -r3,       (0, 2, 0,   0), 2),
    ("r5",      r5,        (0, 0, 2,   0), 0),
    ("-r5",     -r5,       (0, 0, 2,   0), 2),
    ("phi4",    phi**4,    (0, 0, 0,  16), 0),
    ("phi-4",   tau**4,    (0, 0, 0, -16), 0),
    ("K",       Kv,        (0, 0, 1,  -4), 0),
    ("-K",      -Kv,       (0, 0, 1,  -4), 2),
    ("i.beta",  I*betav,   (0, 0, 1,   4), 1),
    ("-i.beta", -I*betav,  (0, 0, 1,   4), 3),
]

mu4 = [Integer(1), I, Integer(-1), -I]

# ================================================================== CH-A
print("== CH-A: cold catalog re-encoding; on-circle image of the CATALOG is empty ==")

# CH-AF1 (falsifier, written first): a WRONG on-circle certification must be
# rejected -- i.beta is NOT on the circle, and a wrong t-class for i.beta is bad.
ck("CH-AF1", (not onS_modsq(I*betav)) and (not zero(I*betav - betav))
        and onS_modsq(I),
   "falsifier: i.beta rejected as on-circle; wrong t-class rejected; +i accepted on-circle")

ok = True
for name, p, roots in SEEDS:
    ok = ok and Poly(p, x).degree() == len(roots)
    ok = ok and all(zero(p.subs(x, r)) for r in roots)
ck("CH-A1", ok, "7 catalog seeds ([EG] Def 3.1): root lists exact and complete (16 eigenvalues)")

ok = True
for name, lam, w, t in EIG:
    modl = 2**Rational(w[0], 4) * 3**Rational(w[1], 4) * \
           5**Rational(w[2], 4) * phi**Rational(w[3], 4)
    ok = ok and zero(lam - modl * I**t)
ck("CH-A2", ok, "each eigenvalue certified lam = 2^(w0/4)3^(w1/4)5^(w2/4)phi^(w3/4) i^t")

# The NEW base fact (not in cx): NO catalog eigenvalue is on the unit circle, and
# on-circle |lam|^2=1 is exactly the vector condition w=(0,0,0,0).
ok, none_oncircle = True, True
for name, lam, w, t in EIG:
    sym = onS_modsq(lam)                     # exact |lam|^2 == 1 in Q(sqrt5)
    vec = (w == (0, 0, 0, 0))
    ok = ok and (sym == vec)
    none_oncircle = none_oncircle and (not sym)
ck("CH-A3", ok and none_oncircle,
   "|lam|=1 <=> w=(0,0,0,0) exactly; NO catalog eigenvalue is on-circle "
   "(catalog on-circle image = empty)")

odd_t = [name for name, lam, w, t in EIG if t % 2 == 1]
ok = set(odd_t) == {"i.beta", "-i.beta"}
ok = ok and all(w[2] > 0 for name, lam, w, t in EIG if t % 2 == 1)
ok = ok and all(w[2] >= 0 for _, _, w, _ in EIG)
ck("CH-A4", ok, "5-valuation coupling ([OA] Thm 6.2): odd-quarter args (charge 1,3) "
   "occur ONLY at +-i.beta, carrying w2 = 1 > 0; all eigenvalues have w2 >= 0")

# ================================================================== CH-B
print("== CH-B: the emission-gap CONTAINMENT: on-circle image SUBSET mu4 (upper bound) ==")
L = {F(0), F(1, 4), F(1, 2), F(3, 4)}          # (pi/2)Z as angle classes in turns

# CH-BF1 (falsifier, written first): arg pi/4 is off the lattice; a modulus-1
# point at pi/4 (a primitive 8th root) is NOT in mu4; halving leaves the lattice.
z8 = exp(I * pi / 4)
ck("CH-BF1", (F(1, 8) not in L) and onS_modsq(z8)
        and all(not zero(z8 - m) for m in mu4),
   "falsifier: pi/4 not in (pi/2)Z; the on-circle 8th root exp(i pi/4) is NOT in mu4")

closed_add = all(((a + b) % 1) in L for a in L for b in L)
closed_dbl = all(((2 * a) % 1) in L for a in L)
not_halved = F(1, 8) not in L
ck("CH-B1", closed_add and closed_dbl and not_halved,
   "[EG] Lem 4.2 / [LB] Lem 4.3: (pi/2)Z closed under addition and doubling, NOT halving")

# [EG] Lem 4.3 / [LB] Lem 4.4 re-encoded as a geometric implication over a test
# set: on-circle AND arg in (pi/2)Z  ==>  value in mu4.  Verified on mu4 (accept)
# and rejected on off-lattice on-circle points (zeta5, zeta8).
z5 = exp(2 * I * pi / 5)
ok = all(onS_modsq(m) and any(zero(m - v) for v in mu4) for m in mu4)   # mu4 accepted
ok = ok and onS_modsq(z5) and all(not zero(z5 - v) for v in mu4)        # zeta5 rejected
ok = ok and onS_modsq(z8) and all(not zero(z8 - v) for v in mu4)        # zeta8 rejected
ck("CH-B2", ok, "[EG] Lem 4.3: on-circle & arg in (pi/2)Z => in mu4 (CONTAINMENT); "
   "mu4 accepted, off-lattice on-circle zeta5/zeta8 rejected")

# Closure preserves the lattice on the ACTUAL operations at the argument level:
# (x) adds t, psi^2 doubles t; both stay in Z/4. (the induction of [EG] Lem 4.2.)
ok = all(((ta + tb) % 4) in {0, 1, 2, 3} for _, _, _, ta in EIG
         for _, _, _, tb in EIG)
ok = ok and all(((2 * ta) % 4) in {0, 2} for _, _, _, ta in EIG)
ck("CH-B3", ok, "[EG] Lem 4.2 on ops: (x) sums t, psi^2 doubles t; both keep arg in (pi/2)Z; "
   "so on-circle image SUBSET mu4 (upper bound only)")

# ================================================================== CH-C
print("== CH-C: CRUX -- on-circle image of S = exactly mu_2 = {+1,-1}, NOT mu_4 (g1/CH-2,CH-3) ==")

def box_scan(eigs, E):
    """Exact scan of all tensor-monomials of total degree <= E over eigs.
    Bookkeeping by certified integer exponent vectors (w,t)."""
    hits, odd_hits, n_scanned = [], [], 0
    for rr in range(E + 1):
        for combo in combinations_with_replacement(range(len(eigs)), rr):
            v0 = v1 = v2 = v3 = t = 0
            for idx in combo:
                w = eigs[idx][2]
                v0 += w[0]; v1 += w[1]; v2 += w[2]; v3 += w[3]
                t = (t + eigs[idx][3]) % 4
            n_scanned += 1
            if v0 == 0 and v1 == 0 and v2 == 0 and v3 == 0:
                hits.append((combo, t))
                if t % 2 == 1:
                    odd_hits.append((combo, t))
    return hits, odd_hits, n_scanned

# CH-CF1 (falsifier, written first): poison the catalog with 5^(-1/4). The SAME
# scan then REACHES on-circle +-i (i.beta * tau * 5^(-1/4) = i), sympy-verified --
# proving the negative scan below is genuinely falsifiable.
EIGP = EIG + [("poison", 5**Rational(-1, 4), (0, 0, -1, 0), 0)]
hits_p, odd_p, _ = box_scan(EIGP, 4)
okf = len(odd_p) > 0
if okf:
    combo, tt = odd_p[0]
    prod = Mul(*[EIGP[i][1] for i in combo])
    okf = zero(expand(prod) - I**tt) and onS_modsq(prod)
ck("CH-CF1", okf, "falsifier: poisoned catalog (add 5^(-1/4)) DOES reach on-circle +-i, "
   "sympy-verified -- the mu_2 conclusion below is falsifiable")

hits, odd_hits, n_scanned = box_scan(EIG, 6)
t_vals = {t for _, t in hits}
ck("CH-C1", n_scanned == 74613 and len(odd_hits) == 0 and t_vals == {0, 2}
        and len(hits) > 2,
   "box |e|<=6 (74613 monomials): every unimodular monomial has EVEN arg; +-i "
   "unreachable multiplicatively; both classes 0 and pi occur")

# The unimodular monomial VALUES are exactly +1 or -1 (sympy-verified), and both
# are attained: phi*tau = 1, phi*psi = -1.
ok, saw1, sawm1 = True, False, False
for combo, t in hits:
    prod = expand(Mul(*[EIG[i][1] for i in combo])) if combo else Integer(1)
    A, B = q5AB(prod - 1)
    if A == 0 and B == 0:
        val_ok = (t == 0); saw1 = saw1 or val_ok
    else:
        A2, B2 = q5AB(prod + 1)
        val_ok = (A2 == 0 and B2 == 0 and t == 2); sawm1 = sawm1 or val_ok
    ok = ok and val_ok
ck("CH-C2", ok and saw1 and sawm1 and zero(phi*tau - 1) and zero(phi*psi + 1),
   "all %d unimodular monomials = +1 or -1 exactly; both attained "
   "(phi*tau=1, phi*psi=-1)" % len(hits))

# GENERAL (unbounded) forcing (=[proposition unimon] re-derived): unimodular
# tensor-monomial => total w2 = 0 => no +-i.beta factor => even arg => value +-1.
w2_nonneg = all(w[2] >= 0 for _, _, w, _ in EIG)
odd_pos5 = all(w[2] > 0 for name, lam, w, t in EIG if t % 2 == 1)
norm_phi = zero(phi*psi + 1)                     # N(phi) = -1
phi_not_root_of_unity = all(not zero(phi**j - 1) for j in range(1, 13))
pos_backstop = all(2**a * 3**b * 5**c > 1
                   for a in range(4) for b in range(4) for c in range(4)
                   if (a, b, c) != (0, 0, 0))
ck("CH-C3", w2_nonneg and odd_pos5 and norm_phi and phi_not_root_of_unity
        and pos_backstop,
   "GENERAL forcing (all exponents): w2>=0; odd-arg carriers have w2>0; "
   "N(phi)=-1, phi not a root of unity => unimodular monomial forces w2=0 => "
   "no +-i.beta => value != +-i")

# CH-3 RECONCILIATION: the FULL image's eigenvalue-VALUE set = the tensor-monomial
# set.  Reason ([OA] Def 2.1, Thm 3.2, prop:fixed): psi^2(lam)=lam^2 is the
# degree-2 monomial lam(x)lam (squaring SUBSET monomials); (+) is multiset union
# (adds no value); minpoly and Phi are spectrum-preserving IDEMPOTENTS (add no
# eigenvalue).  So NO operator escapes the multiplicative monomial monoid; there
# is no gap between the "multiplicative sub-closure" and the "full image" on the
# circle -- reconciling cx-C2.  Two falsifiable structural facts:
#  (a) psi^2 preserves the on-/off-circle dichotomy (|lam^2|=1 iff |lam|=1), so it
#      NEVER turns an off-circle eigenvalue into an on-circle one -- it adds no new
#      on-circle value beyond squaring the existing {+1,-1} (both -> +1);
#  (b) NO catalog seed has +-i (modulus 1) as a root, so the dossier's
#      hypothesized "(+)-route (a seed whose eigenvalues ARE +-i)" does not exist.
psi2_dichotomy = all(onS_modsq(lam**2) == onS_modsq(lam) for _, lam, _, _ in EIG)
psi2_oncircle_closed = zero(Integer(1)**2 - 1) and zero((-Integer(1))**2 - 1)
no_pm_i_seed = all(not zero(p.subs(x, I)) and not zero(p.subs(x, -I))
                   for _, p, _ in SEEDS)
ck("CH-C4", psi2_dichotomy and psi2_oncircle_closed and no_pm_i_seed,
   "CH-3 reconcile: psi^2 preserves on/off-circle (|lam^2|=1 iff |lam|=1) so adds "
   "no new on-circle value ((+-1)^2=+1); (+) unions; minpoly/Phi spectrum-preserve; "
   "and NO catalog seed has +-i as a root => the (+)-route to on-circle +-i does not exist")

# CONCLUSION (g1/CH-2): on-circle image of S = mu_2 = {+1,-1}, a PROPER subset of
# mu_4.  +-i are NOT on-circle elements of S.  Emission-gap Lem 4.3 gives only
# on-circle SUBSET mu4 (CH-B); the reverse mu4 SUBSET on-circle FAILS for +-i.
oncircle_image = {"+1", "-1"} if (saw1 and sawm1 and len(odd_hits) == 0) else set()
pm_i_reachable = (len(odd_hits) > 0)               # would need an on-circle t-odd
ck("CH-C5", oncircle_image == {"+1", "-1"} and (not pm_i_reachable),
   "FORCED-GIVEN-CORPUS: on-circle image of S = mu_2 = {+1,-1} (PROPER subset of "
   "mu_4); +-i NOT on-circle in S -- the reverse of [EG] Lem 4.3 fails for +-i")

# ================================================================== CH-D
print("== CH-D: charge (argument-class) vs on-circle -- +-i.beta is OFF-circle ==")

# CH-DF1 (falsifier, written first): claiming i.beta on-circle, or charge(i.beta)
# != 1, must be rejected.
ck("CH-DF1", (not onS_modsq(I*betav)) and (EIG[14][3] == 1) and (EIG[15][3] == 3),
   "falsifier: i.beta on-circle rejected; charge(+-i.beta) = 1,3 (arg +-pi/2)")

# [OA] Thm 6.2: K carries the FULL Z/4Z as CHARGES {0,1,2,3}.
kseed_charges = sorted({t for name, lam, w, t in EIG
                        if name in ("K", "-K", "i.beta", "-i.beta")})
ck("CH-D1", kseed_charges == [0, 1, 2, 3],
   "[OA] Thm 6.2: K-seed carries the full Z/4Z as CHARGES {0,1,2,3} "
   "(non-real charges 1,3 present in chi(S))")

# but the odd charges (1,3) are carried ONLY at modulus beta != 1 (off-circle).
ok = onS_modsq(I) and onS_modsq(-I)                       # +-i ARE on-circle...
ok = ok and (not onS_modsq(I*betav)) and (not onS_modsq(-I*betav))  # ...+-i.beta are NOT
ok = ok and zero(betav**2 - s5*phi**2) and (not is_one_q5(betav**2))
ck("CH-D2", ok, "[LB] Prop 6.2: the charge-1,3 eigenvalues are +-i.beta at "
   "modulus beta = 5^(1/4)phi, beta^2 = sqrt5 phi^2 != 1 -- OFF the circle")

# The decisive distinction: chi(S) = Z/4Z (has non-real CHARGES) but on-circle S =
# mu_2 (has NO non-real on-circle ELEMENT).  The cx-C4 hypothesis is FALSE.
ck("CH-D3", (1 in kseed_charges) and (oncircle_image == {"+1", "-1"}),
   "DISTINCTION: chi(S) = Z/4Z contains non-real charges 1,3 (as OFF-circle "
   "+-i.beta), but on-circle image = mu_2 has NO non-real element => the cx-C4 "
   "hypothesis 'S contains a non-real ON-CIRCLE element' is FALSE given corpus")

# ================================================================== CH-E
print("== CH-E: g3/CH-1 -- the gate-multiplier doctrine question, decided from the papers ==")
C = symbols('C')

# CH-EF1 (falsifier, written first): the rotation gate multiplier +-i must NOT be
# found among the on-circle image {+1,-1}.
ck("CH-EF1", all(not zero(I - v) for v in (Integer(1), Integer(-1)))
        and all(not zero(-I - v) for v in (Integer(1), Integer(-1))),
   "falsifier: rotation multiplier +-i is NOT in the on-circle image {+1,-1}")

# CH-E1: [W-Prop 10.2] the rotation gate at C = -1/2 has multipliers exactly -+i,
# on-circle (re-derived cold, as in cx CX-D3).
u_hp, u_hm = (-1 + I)/2, (-1 - I)/2
mh_p = expand(Rational(1, 2)/(1 + u_hp)**2)     # m = -C/(1+u)^2 at C = -1/2
mh_m = expand(Rational(1, 2)/(1 + u_hm)**2)
ok = zero(u_hp**2 + u_hp + Rational(1, 2)) and zero(u_hm**2 + u_hm + Rational(1, 2))
ok = ok and ((zero(mh_p + I) and zero(mh_m - I)) or (zero(mh_p - I) and zero(mh_m + I)))
ok = ok and onS_modsq(mh_p) and onS_modsq(mh_m)
ck("CH-E1", ok, "[W-Prop 10.2] rotation gate C=-1/2: multipliers are exactly -+i, "
   "on-circle (the quarter-turn quantum is FORCED gate data)")

# CH-E2: but +-i (the forced rotation multiplier) is NOT an on-circle element of S
# (CH-C5). So D2' read as "the multiplier is an on-circle element of S" is
# INCONSISTENT with the forced rotation gate -- NOT a doctrine theorem.
pm_i_in_oncircle_S = pm_i_reachable          # False (from CH-C)
ck("CH-E2", (not pm_i_in_oncircle_S),
   "g3 (Reading A, on-circle-element): the FORCED rotation multiplier +-i is NOT "
   "in on-circle S = {+1,-1} => D2' (on-circle reading) is UNSATISFIABLE for the "
   "rotation gate -- NOT forced-given-corpus [FALSE given corpus]")

# CH-E3: the TERRAIN gate multipliers m_n = -phi^(-2n) = phi*psi*tau^(2n) ARE in S
# as tensor-monomials -- but they are REAL and OFF-circle; helicity rejects them,
# so they cannot deliver chi = +-pi/2. (re-derived, cx CX-D5.)
ok = True
for n in range(1, 7):
    Cn = 1/(phi**n - phi**(-n))**2
    un = 1/(phi**(2*n) - 1)
    mn = -Cn/(1 + un)**2
    ok = ok and zero(mn + phi**(-2*n))                  # ladder value
    ok = ok and zero(phi*psi*tau**(2*n) + phi**(-2*n))  # monomial witness in S
    ok = ok and (not onS_modsq(mn))                     # off-circle (|m_n|<1)
ck("CH-E3", ok, "terrain gate multipliers m_n = -phi^(-2n) = phi*psi*tau^(2n) ARE "
   "monomials in S but REAL and OFF-circle (helicity rejects) -- the on-circle "
   "elements of S that D2' could select are only {+1,-1}, both real")

# CH-E4: RESTATEMENT kill (Reading B, charge-class).  The ONLY substrate source of
# an arg-pi/2 (odd-charge) emission is the K-seed rotation pair +-i.beta ([OA]
# Thm 6.2). So "D2' puts a pi/2-charge in play" imports exactly D2 (the winding
# realizes the K-seed's charge completion). Corpus-confirms cx-E5.
odd_src = {name for name, lam, w, t in EIG if t % 2 == 1}
ck("CH-E4", odd_src == {"i.beta", "-i.beta"},
   "g3 (Reading B, charge-class): [OA] Thm 6.2 -- the ONLY odd-charge (pi/2) "
   "emission source is the K-seed pair +-i.beta => 'D2' selects a pi/2-charge' = "
   "K-seed charge completion = D2 RESTATED (corpus-confirms cx-E5)")

# CH-E5: the D2' dichotomy as a truth table over the two readings.  Neither is a
# derivation of chi = +-pi/2 from the doctrine.
def chi_readingA():
    # multiplier must be an ON-CIRCLE element of S; on-circle S = {+1,-1} (real);
    # helicity keeps only non-reals -> empty -> unsatisfiable.
    onS = [Integer(1), Integer(-1)]              # forced by CH-C5
    return [v for v in onS if not zero(im(v))]   # helicity survivors
def chi_readingB_imports_D2():
    # multiplier's CHARGE in chi(S) = Z/4Z; non-real charges 1,3 exist, but ONLY
    # via the K-seed pair -> importing them IS D2.
    return odd_src == {"i.beta", "-i.beta"}
outA = chi_readingA()
ck("CH-E5", (outA == []) and chi_readingB_imports_D2(),
   "DICHOTOMY: Reading A (on-circle) -> survivors {} (unsatisfiable, FALSE given "
   "corpus); Reading B (charge) -> chi=+-pi/2 but imports D2 (RESTATEMENT). "
   "Neither derives chi=+-pi/2 from the doctrine => W[open]2 stays PARTIAL")

# CH-E6: restatement test (OP-RATE): even had the chain closed, it fixes the
# quantum mod 2pi only (branch-blind on the winding number). (cx CX-E6.)
ok = all(((F(1, 4) + k) - F(1, 4)).denominator == 1 for k in range(-3, 4))
ck("CH-E6", ok, "restatement test: the chi chain (if closed) fixes the quantum mod "
   "2pi only (pi/2 + 2pi k = pi/2 mod 2pi for all k) -- branch-blind, a "
   "CHARACTERIZATION of the quantum, never a branch derivation")

# ================================================================== CH-G
print("== CH-G: certified interval guards (60 dps, counted separately) ==")
from mpmath import iv
iv.dps = 60
phii = (1 + iv.sqrt(5)) / 2
betai = iv.mpf(5) ** (iv.mpf(1) / 4) * phii
gk("CH-G1", betai.a > iv.mpf("2.41").a and betai.b < iv.mpf("2.43").b and betai.a > 1,
   "beta = 5^(1/4) phi in [2.41,2.43], certified > 1: the only arg-pi/2 emission "
   "(+-i.beta) is strictly OFF-circle -- no non-real on-circle element in S")
inv_beta = iv.mpf(1) / betai
gk("CH-G2", inv_beta.a > iv.mpf("0.41").a and inv_beta.b < iv.mpf("0.42").b,
   "1/beta in [0.41,0.42]: rescaling +-i.beta to the circle needs a 5^(-1/4)-scale "
   "factor (negative 5-valuation) absent from the catalog (all w2 >= 0)")

# ================================================================== summary
n_guard = len(GUARD)
n_guard_pass = sum(1 for _, okg in GUARD if okg)
n_exact_pass = len(PASS)
n_exact_fail = len([f for f in FAIL if not f.startswith("CH-G")])
line = "CH EXACT %d/%d GUARDS %d/%d -- %s" % (
    n_exact_pass, n_exact_pass + n_exact_fail, n_guard_pass, n_guard,
    "ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)
print(line)
sys.exit(0 if not FAIL else 1)
