#!/usr/bin/env python3
# cx_softcheck.py -- chi-selection dossier harness (HANDOFF-v1.0 SS4, W[open]2).
# Written cold this session from whitepaper SS5.1-5.2 (catalog, Lemmas 5.4-5.6,
# Remark 5.3), Prop 8.2, Prop 10.2, Thm 11.1, Cor 11.2, Decl 12.2 (D2), and the
# ECHO-S-RESEARCH_UNIFIED.md catalog/operator definition (the only shipped
# corpus-synthesis source; Operator-Algebra and Charge-Measure papers NOT shipped).
#
# Discipline: every DECISION exact over Q / Q(sqrt5) (extended by i and 5^(1/4)
# where needed); unimodularity decided in Q(sqrt5) by rational (A,B) coefficient
# arithmetic (q5AB), monomial bookkeeping by integer exponent vectors that are
# themselves certified symbolically per eigenvalue (CX-A2).  mpmath iv at 60 dps
# is used ONLY for the certified interval guards (CX-G), counted separately.
# Fail-first: each group opens with falsifier guards (CX-*F*) asserting that
# deliberately perturbed variants are REJECTED; the falsifiers were written first.
#
# Groups:
#   CX-A  catalog replication + axis/angle-confinement (Lemma 5.4 content)   (6)
#   CX-B  angle-lattice arithmetic: pentagon exclusion by image (C1 / g2)    (6)
#   CX-C  unimodularity: tensor-monomials never reach +-i; on-circle
#         value set of the shipped closure of S is exactly {+1,-1} (C2,C3)   (8)
#   CX-D  gate layer: obstruction links (Prop 8.2 / Thm 11.1 / Prop 10.2)
#         + terrain/rotation multiplier dichotomy (C4 links, C6)             (7)
#   CX-E  chi-selection chain truth table over mu_4 (C4)                     (6)
#   CX-G  certified 60-dps interval guards (counted separately)              (2)
# Exit 0 iff all pass.

import sys
from fractions import Fraction as F
from itertools import combinations_with_replacement
from sympy import (symbols, sqrt, Rational, Integer, I, pi, sin, cos, exp,
                   im, re, simplify, expand, radsimp, together, nsimplify,
                   Poly, minimal_polynomial, real_roots, gcd, fibonacci,
                   Mul, conjugate)

x = symbols('x')
y = symbols('y', real=True)
s5 = sqrt(5)
phi = (1 + s5) / 2
psi = (1 - s5) / 2                 # = -tau, the golden conjugate
tau = (s5 - 1) / 2                 # = 1/phi
f4 = 5 ** Rational(1, 4)           # 5^(1/4)
Kv = f4 * tau                      # K-seed terrain root, K = 5^(1/4)/phi
betav = f4 * phi                   # K-seed rotation magnitude, beta = 5^(1/4) phi
r2, r3, r5 = sqrt(2), sqrt(3), sqrt(5)

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

# ------------------------------------------------------------------ catalog
# Whitepaper SS5.1 (verified verbatim): seeds = companions of
# { phi, tau, sqrt2, sqrt3, sqrt5, gap (x^2-7x+1), K (x^4+5x^2-5) }.
SEEDS = [
    ("phi-seed",  x**2 - x - 1,      [phi, psi]),
    ("tau-seed",  x**2 + x - 1,      [tau, -phi]),
    ("sqrt2",     x**2 - 2,          [r2, -r2]),
    ("sqrt3",     x**2 - 3,          [r3, -r3]),
    ("sqrt5",     x**2 - 5,          [r5, -r5]),
    ("gap-seed",  x**2 - 7*x + 1,    [phi**4, tau**4]),
    ("K-seed",    x**4 + 5*x**2 - 5, [Kv, -Kv, I*betav, -I*betav]),
]

# Eigenvalue table: (name, symbol, w, t) with |lam|^4 = 2^w0 3^w1 5^w2 phi^w3
# and lam = |lam| * i^t (t in Z/4).  Every entry is CERTIFIED in CX-A2.
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

def on_axes(p):
    """Exact test: every root of the rational-coefficient squarefree poly p is
    real or purely imaginary (i.e. all eigenvalue arguments in (pi/2)Z)."""
    deg = Poly(p, x).degree()
    n_real = len(real_roots(Poly(p, x)))
    e = expand(p.subs(x, I*y))
    ree, ime = expand(re(e)), expand(im(e))
    if ime == 0:
        g = ree
    elif ree == 0:
        g = ime
    else:
        g = gcd(Poly(ree, y), Poly(ime, y)).as_expr()
    n_imag = 0 if g.is_number else len([r for r in real_roots(Poly(g, y))
                                        if r != 0])
    return n_real + n_imag == deg

# ================================================================== CX-A
print("== CX-A: catalog replication + angle confinement (Lemma 5.4) ==")

# CX-AF1 (falsifier, written first): the plastic seed x^3 - x - 1 must FAIL the
# axis test (its complex pair is neither real nor purely imaginary: Re p(iy) = -1
# identically), and a wrong argument-class certification must be REJECTED.
plastic = x**3 - x - 1
wrong_cert = zero(I*betav - betav)   # claims t=0 for i.beta -- must be False
ck("CX-AF1", (not on_axes(plastic)) and (not wrong_cert),
   "falsifiers: plastic seed rejected by axis test; wrong t-class for i.beta rejected")

ok = True
for name, p, roots in SEEDS:
    ok = ok and Poly(p, x).degree() == len(roots)
    ok = ok and all(zero(p.subs(x, r)) for r in roots)
ck("CX-A1", ok, "7 catalog seeds: root lists exact and complete (16 eigenvalues)")

ok = True
for name, lam, w, t in EIG:
    modl = 2**Rational(w[0], 4) * 3**Rational(w[1], 4) * \
           5**Rational(w[2], 4) * phi**Rational(w[3], 4)
    ok = ok and zero(lam - modl * I**t)
ck("CX-A2", ok, "each eigenvalue certified: lam = 2^(w0/4) 3^(w1/4) 5^(w2/4) phi^(w3/4) i^t")

ck("CX-A3", all(on_axes(p) for _, p, _ in SEEDS),
   "angle confinement replicated: every catalog seed has all roots on the axes")

odd_t = [name for name, lam, w, t in EIG if t % 2 == 1]
ok = set(odd_t) == {"i.beta", "-i.beta"}
ok = ok and all(w[2] > 0 for name, lam, w, t in EIG if t % 2 == 1)
ok = ok and all(w[0] >= 0 and w[1] >= 0 and w[2] >= 0 for _, _, w, _ in EIG)
ck("CX-A4", ok, "5-valuation coupling: odd-quarter args occur ONLY at +-i.beta, "
   "which carry positive 5-exponent; all non-unit exponents nonnegative")

ck("CX-A5", zero(Kv*phi - f4) and zero(betav*tau - f4) and zero(Kv*betav - s5)
        and zero(betav - Kv*phi**2),
   "moduli generators: K phi = beta tau = 5^(1/4); K beta = sqrt5; beta = K phi^2")

# ================================================================== CX-B
print("== CX-B: angle-lattice arithmetic -- pentagon exclusion by image (C1/g2) ==")
L = {F(0), F(1, 4), F(1, 2), F(3, 4)}          # image angle classes, in turns

# CX-BF1 (falsifier, written first): the membership test must ACCEPT lattice
# points, REJECT pentagon points, and flip on the perturbed lattice (1/5)Z.
L5 = {F(k, 5) for k in range(5)}
ck("CX-BF1", (F(1, 4) in L) and (F(1, 2) in L) and (F(2, 5) not in L)
        and (F(2, 5) in L5) and (F(1, 4) not in L5),
   "falsifier: lattice membership test distinguishes (1/4)Z from (1/5)Z")

ck("CX-B1", all((F(k, 5) % 1) not in L and (4 * F(k, 5)).denominator != 1
                for k in (1, 2, 3, 4)),
   "pentagon classes k/5 (k=1..4) not in (1/4)Z mod 1")

closed_add = all(((a + b) % 1) in L for a in L for b in L)
closed_dbl = all(((2 * a) % 1) in L for a in L)
not_halved = F(1, 8) not in L                   # halving 1/4 leaves the lattice
ck("CX-B2", closed_add and closed_dbl and not_halved,
   "Lemma 5.5 content: lattice closed under addition and doubling, NOT halving")

z5 = exp(2 * I * pi / 5)
Phi5 = x**4 + x**3 + x**2 + x + 1
mu4 = [Integer(1), I, Integer(-1), -I]
ok = Poly(minimal_polynomial(z5, x), x) == Poly(Phi5, x)
ok = ok and zero(z5 * conjugate(z5) - 1)
ok = ok and all(not zero(Phi5.subs(x, lam)) for lam in mu4)
ck("CX-B3", ok, "zeta_5 exact: minpoly = Phi_5, unimodular, Phi_5 has no root in mu_4")

ck("CX-B4", set(F(k, 5) % 1 for k in range(5)) & L == {F(0)},
   "pentagon angle classes meet the image lattice only at 0: zeta_5 not in S")

m0 = (tau + I * sqrt((5 + s5) / 2)) / 2         # root of m^2 - tau m + 1
ok = zero(2 * cos(2 * pi / 5) - tau) and zero(2 * cos(4 * pi / 5) + phi)
ok = ok and zero(expand((x**2 - tau*x + 1) * (x**2 + phi*x + 1)) - Phi5)
ok = ok and Poly(minimal_polynomial(m0, x), x) == Poly(Phi5, x)
ck("CX-B5", ok, "pentagon gates: traces 2cos(2pi/5) = tau, 2cos(4pi/5) = -phi; "
   "multipliers are primitive 5th roots (pair product = Phi_5)")

# ================================================================== CX-C
print("== CX-C: unimodularity -- tensor-monomials never reach +-i (C2,C3/g1) ==")

def box_scan(eigs, E):
    """Exact scan of all tensor-monomials of total degree <= E over eigs.
    Returns (all unimodular hits, unimodular hits with odd quarter-arg, count)."""
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

# CX-CF1 (falsifier, written first): poison the catalog with 5^(-1/4).  The SAME
# scan must then REACH +-i (e.g. i.beta * tau * 5^(-1/4) = i), sympy-verified --
# proving the scan detects reachability when it exists.
EIGP = EIG + [("poison", 5**Rational(-1, 4), (0, 0, -1, 0), 0)]
hits_p, odd_p, _ = box_scan(EIGP, 4)
okf = len(odd_p) > 0
if okf:
    combo, tt = odd_p[0]
    prod = Mul(*[EIGP[i][1] for i in combo])
    okf = zero(expand(prod) - I**tt)
ck("CX-CF1", okf, "falsifier: poisoned catalog (add 5^(-1/4)) DOES reach +-i, "
   "sympy-verified -- the negative scan below is falsifiable")

hits, odd_hits, n_scanned = box_scan(EIG, 6)
t_vals = {t for _, t in hits}
ck("CX-C1", n_scanned == 74613 and len(odd_hits) == 0 and t_vals == {0, 2}
        and len(hits) > 2,
   "box |e|<=6 (74613 monomials): every unimodular monomial has even quarter-arg;"
   " +-i unreachable; both classes 0 and pi occur")

ok, saw1, sawm1 = True, False, False
for combo, t in hits:
    prod = expand(Mul(*[EIG[i][1] for i in combo])) if combo else Integer(1)
    A, B = q5AB(prod - 1)
    if A == 0 and B == 0:
        val_ok = (t == 0)
        saw1 = saw1 or val_ok
    else:
        A2, B2 = q5AB(prod + 1)
        val_ok = (A2 == 0 and B2 == 0 and t == 2)
        sawm1 = sawm1 or val_ok
    ok = ok and val_ok
ck("CX-C2", ok and saw1 and sawm1,
   "all %d unimodular hits sympy-verified exactly equal +1 or -1" % len(hits))

ck("CX-C3", zero(phi*tau - 1) and zero(phi*psi + 1),
   "on-circle witnesses IN S: 1 in spec(phi (x) tau), -1 in spec(phi (x) phi)")

# Independent route: decide |monomial|^2 == 1 exactly in Q(sqrt5) via q5AB for
# ALL monomials of degree <= 3 and cross-check against the vector method.
ok = True
for rr in range(4):
    for combo in combinations_with_replacement(range(len(EIG)), rr):
        prod = Mul(*[EIG[i][1] for i in combo]) if combo else Integer(1)
        m2 = expand(prod * conjugate(prod))
        uni_sym = is_one_q5(m2)
        v = [sum(EIG[i][2][j] for i in combo) for j in range(4)]
        uni_vec = (v == [0, 0, 0, 0])
        ok = ok and (uni_sym == uni_vec)
ck("CX-C4", ok, "cross-route |e|<=3 (969 monomials): Q(sqrt5)-exact |lam|^2 = 1 "
   "decision agrees with the certified exponent-vector method everywhere")

# GENERAL (unbounded) forcing of +-i unreachability via 5-valuation.
# Premises, each checked exactly:
#  (P1) every catalog eigenvalue has 5-exponent w2 >= 0 (no eigenvalue supplies
#       a negative power of 5);
#  (P2) the ONLY eigenvalues with odd quarter-argument (t odd) are +-i.beta,
#       and they carry w2 = 1 > 0.
# For a tensor-monomial: total w2 = sum of factor w2 (all >= 0); total t = sum of
# factor t (mod 4).  Unimodular => modulus^4 = 2^w0 3^w1 5^w2 phi^w3 = 1 => (norm
# to Q, N(phi) = -1) the rational part 2^w0 3^w1 5^w2 = 1 with w0,w1,w2 >= 0, so
# w2 = 0 => (nonneg summands) NO +-i.beta factor => total t even => value != +-i.
# The harness certifies the premises and the norm identity; the conclusion is
# then a finite piece of logic (nonnegativity) valid for ALL exponents, not just
# the box.  CX-C1 corroborates on |e| <= 6.
w2_nonneg = all(w[2] >= 0 for _, _, w, _ in EIG)
odd_carriers = [name for name, lam, w, t in EIG if t % 2 == 1]
odd_pos5 = all(w[2] > 0 for name, lam, w, t in EIG if t % 2 == 1)
norm_phi = zero(phi*psi + 1)                     # N_{Q(phi)/Q}(phi) = phi*psi = -1
phi_not_root_of_unity = all(not zero(phi**j - 1) for j in range(1, 13))
# positivity backstop: any positive nonneg 2,3,5-exponent makes the rational > 1
pos_backstop = all(2**a * 3**b * 5**c > 1
                   for a in range(4) for b in range(4) for c in range(4)
                   if (a, b, c) != (0, 0, 0))
ck("CX-C5", w2_nonneg and odd_carriers == ["i.beta", "-i.beta"] and odd_pos5
        and norm_phi and phi_not_root_of_unity and pos_backstop,
   "GENERAL forcing: 5-valuation >= 0 on all eigenvalues; odd-arg carriers = "
   "{+-i.beta} with w2 > 0; N(phi) = -1, phi not a root of unity => unimodular "
   "tensor-monomial forces w2 = 0 => no +-i.beta factor => value != +-i (all exps)")

# Closure semantics (Lemma 5.5 verbatim: (x) multiplies, psi^2 squares, (+)
# unions, minpoly/Phi-extraction select subsets): conjugates of catalog
# eigenvalues are eigenvalues of the SAME seed, so conjugate-closure adds no
# values; squares of monomials are monomials by vector doubling (structural).
ok = True
for name, p, roots in SEEDS:
    for r in roots:
        ok = ok and any(zero(conjugate(r) - r_) for r_ in roots)
ck("CX-C6", ok, "each seed's root multiset is conjugation-closed: minpoly/Phi "
   "extraction introduces no values outside the monomial set")

# Honest scope on g1.  The MULTIPLICATIVE (tensor-monomial) sub-closure of S has
# on-circle value set exactly {+1,-1} [FORCED, CX-C5].  Reaching +-i requires a
# NON-multiplicative operator (Lemma 5.5: minpoly / Phi-extraction / (+)).  The
# corpus-exhibited candidate host is Remark 5.3's content polynomial
# x^4-1 = Phi_1 Phi_2 Phi_4 [forced SC-I13], whose factor Phi_4 = x^2+1 is rooted
# at +-i (the "rotation fiber"); its real face is exactly {+-1} (the "terrain").
# Whether x^4-1 (or a +-i-bearing object) lies in the SHIPPED operator-closure of
# S is BOUNDED: the Operator-Algebra / Charge-Measure papers that fix the closure
# are absent from the file set, so g1 cannot be decided here -- report it as such.
Phi1, Phi2, Phi4 = x - 1, x + 1, x**2 + 1
mult_onS = saw1 and sawm1 and len(odd_hits) == 0
r53_factor = zero(expand(Phi1*Phi2*Phi4) - (x**4 - 1))
r53_fiber = zero(Phi4.subs(x, I)) and zero(Phi4.subs(x, -I)) \
    and zero(I*conjugate(I) - 1)                 # +-i unimodular, roots of Phi_4
r53_terrain = sorted(int(r) for r in real_roots(Poly(x**4 - 1, x))) == [-1, 1]
ck("CX-C7", mult_onS and r53_factor and r53_fiber and r53_terrain,
   "g1 (bounded): tensor sub-closure on-circle = {+1,-1} [forced]; +-i needs "
   "Phi-extraction -- host x^4-1 = Phi_1 Phi_2 Phi_4 (Rem 5.3), Phi_4 = x^2+1 "
   "rooted at +-i; closure-membership BOUNDED (operator-algebra paper absent)")

# ================================================================== CX-D
print("== CX-D: gate layer -- obstruction links + multiplier dichotomy ==")
C = symbols('C')
s = symbols('s', positive=True)
th = symbols('th', real=True)
Dd = sqrt(1 + 4*C)
up, um = (-1 + Dd)/2, (-1 - Dd)/2
mp_, mm_ = -C/(1 + up)**2, -C/(1 + um)**2

# CX-DF1 (falsifier, written first): perturbed flip identities must be REJECTED.
ck("CX-DF1", (not zero(mp_*mm_ + 1)) and (not zero(mp_ + mm_ + 1/C - 2)),
   "falsifiers: m+m- = -1 and m + 1/m = -1/C + 2 both rejected symbolically")

ck("CX-D1", zero(mp_*mm_ - 1) and zero(mp_ + mm_ + 1/C + 2),
   "Prop 8.2 links: m+m- = 1 and m + 1/m = -1/C - 2, generic symbolic C")

m_s = s * exp(I * th)
ok = zero(im(m_s + 1/m_s) - (s - 1/s) * sin(th))
ok = ok and [r for r in real_roots(Poly(x**2 - 1, x)) if r > 0] == [1]
ck("CX-D2", ok,
   "Thm 11.1 algebra: Im(m + 1/m) = (s - 1/s) sin th; s - 1/s = 0 has unique "
   "positive root s = 1 (real C => multiplier real or unimodular)")

u_hp, u_hm = (-1 + I)/2, (-1 - I)/2
mh_p = expand(Rational(1, 2)/(1 + u_hp)**2)     # m = -C/(1+u)^2 at C = -1/2
mh_m = expand(Rational(1, 2)/(1 + u_hm)**2)
ok = zero(u_hp**2 + u_hp + Rational(1, 2)) and zero(u_hm**2 + u_hm + Rational(1, 2))
ok = ok and ((zero(mh_p + I) and zero(mh_m - I)) or
             (zero(mh_p - I) and zero(mh_m + I)))
ok = ok and not (zero(mh_p - 1) or zero(mh_p + 1))
ck("CX-D3", ok, "Prop 10.2: at C = -1/2 the multipliers are exactly -+i "
   "(and not +-1): the quarter-turn quantum is gate data")

ck("CX-D4", zero(-1/(-tau**2) - 2 - tau) and zero(-1/(-phi**2) - 2 + phi),
   "pentagon gates C = -phi^-2, -phi^2: traces tau and -phi "
   "(with CX-B5: multipliers are primitive 5th roots of unity)")

ok = True
for n in range(1, 7):
    Cn = 1/(phi**n - phi**(-n))**2
    un = 1/(phi**(2*n) - 1)
    mn = -Cn/(1 + un)**2
    ok = ok and zero(mn + phi**(-2*n))                  # Thm 9.1 ladder value
    ok = ok and zero(phi*psi*tau**(2*n) + phi**(-2*n))  # monomial witness in S
ok = ok and zero(Rational(1, 5)/(1 + 1/(phi**4 - 1))**2 - tau**4)  # m_2 = -gap
ck("CX-D5", ok, "dichotomy, terrain half: every ladder multiplier m_n = -phi^(-2n)"
   " = phi*psi*tau^(2n) is a tensor-monomial value (emission-realized); n = 1..6")

ck("CX-D6", len(odd_hits) == 0 and all((F(k, 5) % 1) not in L for k in (1, 2, 3, 4)),
   "dichotomy, rotation half: +-i (CX-C1) and the pentagon multipliers (CX-B1/B4)"
   " are NOT values of S -- the rotation face is emission-transcendent")

# ================================================================== CX-E
print("== CX-E: the chi-selection chain (C4) ==")

def order_mu4(lam):
    for j in (1, 2, 3, 4):
        if zero(lam**j - 1):
            return j
    return 0

# CX-EF1 (falsifier, written first): the PERTURBED chain (excluding non-reals
# instead of reals) leaves {+1,-1}, whose elements have order <= 2: they cannot
# generate Z/4Z -- the perturbed selection is rejected by the charge anatomy.
pert = [lam for lam in mu4 if zero(im(lam))]
ck("CX-EF1", all(order_mu4(lam) < 4 for lam in pert) and len(pert) == 2,
   "falsifier: perturbed chain (keep reals) yields {+1,-1}, no Z/4Z generator")

ok = all(zero(lam * conjugate(lam) - 1) for lam in mu4)
ok = ok and all(zero(I**k - mu4[k]) for k in range(4))
ok = ok and sorted(F(k, 4) % 1 for k in range(4)) == sorted(L)
ck("CX-E1", ok, "D2''-class: on-circle + arg in (pi/2)Z <=> mu_4 = {1, i, -1, -i}"
   " (one on-circle representative per lattice class)")

surv = [lam for lam in mu4 if not zero(im(lam))]
ok = len(surv) == 2 and any(zero(lam - I) for lam in surv) \
     and any(zero(lam + I) for lam in surv)
ck("CX-E2", ok, "helicity link (Thm 11.1/Cor 11.2 anatomy): real multipliers "
   "carry no rotation; survivors = {+i, -i}, i.e. chi = +-pi/2")

ok = ({k for k in range(4) if not zero(im(mu4[k]))} ==
      {k for k in range(4) if order_mu4(mu4[k]) == 4})
ck("CX-E3", ok, "de-axiomatization identity: within mu_4, non-real = order-4 = "
   "generator of Z/4Z -- D2's generator clause is a theorem given helicity")

# The chain, made honest as a truth table over the g1 disjunction (is +-i in S?).
# on_circle_S(present) := the on-circle elements of S that D2' may select.
#  - present = True  (Phi_4-extraction lands in S): set = {+1,-1,+i,-i}
#  - present = False (multiplicative closure only):  set = {+1,-1}
# In both branches, D2' + helicity (Thm 11.1: real multipliers carry no rotation)
# removes {+1,-1}; the outcome is the SURVIVOR set.
def chi_outcome(pm_i_present):
    onS = [Integer(1), Integer(-1)] + ([I, -I] if pm_i_present else [])
    return [lam for lam in onS if not zero(im(lam))]   # helicity keeps non-reals
out_yes = chi_outcome(True)
out_no = chi_outcome(False)
ok = (len(out_yes) == 2 and any(zero(l - I) for l in out_yes)
      and any(zero(l + I) for l in out_yes)          # +-i => chi = +-pi/2 forced
      and out_no == [])                              # no +-i => D2' unsatisfiable
ck("CX-E4", ok, "chi truth table: IF +-i in S then chi = +-pi/2 [forced given D2'];"
   " IF NOT then D2' has no helicity-carrying multiplier (unsatisfiable) -> fall "
   "back to D2''. Outcome hinges entirely on the bounded g1 question (CX-C7)")

# KILL: the ONLY substrate carrier of an odd-quarter (+-pi/2) emission argument is
# the K-seed's rotation pair +-i.beta (CX-A4). So "+-i in S" -- the hypothesis
# that makes the chain non-vacuous -- is sourced from exactly the K-seed's Z/4Z
# charge, i.e. Declaration D2 ("the winding realizes the seed's charge
# completion") restated at the eigenvalue level. Motivating D2' this way IS D2.
odd_src = {name for name, lam, w, t in EIG if t % 2 == 1}
kseed_rot = {"i.beta", "-i.beta"}
ck("CX-E5", odd_src == kseed_rot,
   "RESTATEMENT kill (axiom-replacement claim): the substrate's ONLY odd-quarter "
   "emission argument comes from the K-seed's rotation pair +-i.beta; so 'D2' "
   "puts +-i in S' is the K-seed charge-completion = D2 restated (chain still ships)")

ok = all(((F(1, 4) + k) - F(1, 4)).denominator == 1 for k in range(-3, 4))
ck("CX-E6", ok, "restatement test: the chain fixes the quantum mod 2pi only "
   "((pi/2 + 2pi k) = pi/2 mod 2pi for all k): it selects the QUANTUM "
   "(OP-RATE sub-target i), branch-blind on the winding number")

# ================================================================== CX-G
print("== CX-G: certified interval guards (60 dps, counted separately) ==")
from mpmath import iv
iv.dps = 60
phii = (1 + iv.sqrt(5)) / 2
betai = iv.mpf(5) ** (iv.mpf(1) / 4) * phii
gk("CX-G1", betai.a > iv.mpf("2.41").a and betai.b < iv.mpf("2.43").b
        and betai.a > 1,
   "beta = 5^(1/4) phi in [2.41, 2.43], certified > 1: the quarter-turn "
   "eigenvalue is strictly off-circle")
sep = iv.pi / 2 - 2 * iv.pi / 5
gk("CX-G2", sep.a > iv.mpf("0.31").a,
   "certified pi/2 - 2pi/5 = pi/10 > 0.31: pentagon quantum separated from "
   "the lattice point (corroborates the exact Fraction decision)")

# ================================================================== summary
n_guard = len(GUARD)
n_guard_pass = sum(1 for _, okg in GUARD if okg)
n_exact_pass = len(PASS)
n_exact_fail = len([f for f in FAIL if not f.startswith("CX-G")])
line = "CX EXACT %d/%d GUARDS %d/%d -- %s" % (
    n_exact_pass, n_exact_pass + n_exact_fail, n_guard_pass, n_guard,
    "ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)
print(line)
sys.exit(0 if not FAIL else 1)
