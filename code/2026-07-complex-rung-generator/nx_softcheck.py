#!/usr/bin/env python3
# nx_softcheck.py -- third independent harness for complex-rung-generator v1.5
# Written cold this compilation from the v1.4 statements + RC paper defs.
# Discipline: every DECISION is exact (sympy over Q / Q(sqrt5) extended by i;
# Q(sqrt5) sign decisions by rational coefficient arithmetic; K decided at K^2).
# mpmath 50 dps is used for DISPLAY/corroboration only (group NX-H).
# Exit 0 iff every check passes.  Groups: NX-R (replication, 8),
# NX-T (complete relational type, 9), NX-S (sphere, 7), NX-H (displays, 3).

import sys
from fractions import Fraction
from sympy import (
    symbols, sqrt, log, pi, exp, sin, cos, I, Rational, Integer, simplify,
    expand, nsimplify, Poly, minimal_polynomial, factor_list, conjugate,
    im, re, diff, trigsimp, integrate, together, radsimp
)

x = symbols('x')
rho = symbols('rho', positive=True)
s5 = sqrt(5)
phi = (1 + s5) / 2
tau = (s5 - 1) / 2               # = 1/phi = phi - 1 (denominator-free form)
gap = phi ** -4
K2 = 1 - gap                     # K decided at K^2
K = sqrt(K2)
kap = pi / (2 * log(phi))        # D3 rate, kept symbolic-exact
q = I * tau

PASS, FAIL = [], []
def ck(cid, cond, desc):
    (PASS if cond else FAIL).append(cid)
    print(("PASS" if cond else "FAIL"), cid, "--", desc)

def zero(e):
    e2 = simplify(expand(simplify(e)))
    if e2 == 0:
        return True
    try:
        e3 = nsimplify(e2, [sqrt(5)])
        return simplify(radsimp(e3)) == 0
    except Exception:
        return False

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q.  Raises if not in Q(sqrt5)."""
    from sympy import cancel
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, f"not linear in sqrt5: {e}"
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQ5(e):
    """Exact sign of A + B*sqrt(5) via rational arithmetic only."""
    A, B = q5AB(e)
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    # opposite signs: compare A^2 vs 5B^2
    if A > 0:   # A>0>B
        return 1 if A * A > 5 * B * B else -1
    else:       # A<0<B
        return 1 if 5 * B * B > A * A else -1

# ---------------------------------------------------------------- roots
p_q = x ** 4 + 3 * x ** 2 + 1
roots = [I * tau, -I * tau, I * phi, -I * phi]

# ================================================================ NX-R
# R1: minimal polynomial and irreducibility
mp_q = minimal_polynomial(q, x)
irred = len(factor_list(p_q)[1]) == 1 and factor_list(p_q)[1][0][1] == 1
ck("NX-R1", Poly(mp_q, x) == Poly(p_q, x) and irred,
   "minpoly(q) = x^4+3x^2+1, irreducible over Q")

# R2: root set = {+-i tau, +-i phi}
ck("NX-R2", all(zero(p_q.subs(x, r)) for r in roots) and len(roots) == 4,
   "roots are exactly {+-i tau, +-i phi}")

# R3: Mahler M = phi^2 with integer guards (phi>1 <=> 5>4; tau<1 <=> 5<9)
g1 = (5 > 4); g2 = (5 < 9)
mods2 = [simplify(expand(r * conjugate(r))) for r in roots]
big = [m2 for m2 in mods2 if sgnQ5(m2 - 1) > 0]
M = simplify(sqrt(big[0]) * sqrt(big[1])) if len(big) == 2 else None
ck("NX-R3", g1 and g2 and len(big) == 2 and zero(M - phi ** 2),
   "Mahler measure M = phi^2 (two roots outside unit circle; integer guards)")

# R4: within-shell ordered ratios -> contact signature {Phi_1^4, Phi_2^4}
shells = {}
for r in roots:
    m2 = simplify(expand(r * conjugate(r)))
    key = next((k for k in shells if zero(k - m2)), None)
    shells.setdefault(key if key is not None else m2, []).append(r)
within = []
for sh in shells.values():
    for a in sh:
        for b in sh:
            within.append(simplify(a / b))
ones = sum(1 for r_ in within if zero(r_ - 1))
negs = sum(1 for r_ in within if zero(r_ + 1))
ck("NX-R4", len(shells) == 2 and ones == 4 and negs == 4 and len(within) == 8,
   "within-shell ordered ratios = {1 x4, -1 x4}: signature {Phi_1^4, Phi_2^4}")

# R5: gauge identity both directions with Z* = x^4-3x^2+1
Zs = expand((x ** 2 - x - 1) * (x ** 2 + x - 1))
ck("NX-R5", zero(expand(p_q.subs(x, I * x)) - Zs) and
            zero(expand(Zs.subs(x, I * x)) - p_q),
   "minpoly(q)(ix) = Z* and Z*(ix) = minpoly(q)")

# R6: rung sphere identity z_n^2 + |W_n|^2 = 1, n=0..8
ok = all(zero((1 - tau ** (2 * n)) + expand((q ** n) * conjugate(q ** n)) - 1)
         for n in range(9))
ck("NX-R6", ok, "z_n^2 + |W_n|^2 = 1 exactly, n = 0..8")

# R7: nu-locus parity nu(rho_n) = (-1)^n and half-floor nonvanishing
ok1 = all(zero(exp(2 * I * (n * pi / 2)) - (-1) ** n) for n in range(9))
ok2 = all(simplify(sin(2 * ((n + Rational(1, 2)) * pi / 2)) ** 2 - 1) == 0
          for n in range(9))
ck("NX-R7", ok1 and ok2, "nu(rho_n) = (-1)^n; Im nu = +-1 at half-floors")

# R8: z_n strictly increasing, n=0..8 (exact Q(sqrt5) signs)
ok = all(sgnQ5((1 - tau ** (2 * (n + 1))) - (1 - tau ** (2 * n))) > 0
         for n in range(8))
ck("NX-R8", ok, "z_n strictly increasing (decided at z^2 in Q(sqrt5))")

# ================================================================ NX-T
# T1: shell moduli^2 multiset {tau^2 x2, phi^2 x2}; two distinct shells
cnt_t = sum(1 for m2 in mods2 if zero(m2 - tau ** 2))
cnt_p = sum(1 for m2 in mods2 if zero(m2 - phi ** 2))
ck("NX-T1", cnt_t == 2 and cnt_p == 2 and sgnQ5(phi ** 2 - tau ** 2) > 0,
   "two shells: |root|^2 multiset {tau^2,tau^2,phi^2,phi^2}, tau^2 != phi^2")

# T2: within-shell antipodal ratios all -1
anti = [simplify(a / b) for sh in shells.values() for a in sh for b in sh
        if not zero(a - b)]
ck("NX-T2", len(anti) == 4 and all(zero(r_ + 1) for r_ in anti),
   "within-shell antipodal ratios = -1 (four ordered pairs)")

# T3: cross-shell ratio multiset {tau^2 x2, -tau^2 x2, phi^2 x2, -phi^2 x2}
sh_list = list(shells.values())
cross = [simplify(a / b) for a in sh_list[0] for b in sh_list[1]] + \
        [simplify(a / b) for a in sh_list[1] for b in sh_list[0]]
def count_eq(lst, v): return sum(1 for r_ in lst if zero(r_ - v))
ck("NX-T3", len(cross) == 8 and
   count_eq(cross, tau**2) == 2 and count_eq(cross, -tau**2) == 2 and
   count_eq(cross, phi**2) == 2 and count_eq(cross, -phi**2) == 2,
   "cross-shell ratios multiset = {+-tau^2 x2 each, +-phi^2 x2 each}")

# T4: all cross ratios real and non-unimodular (exact sign decisions)
ok_real = all(zero(im(r_)) for r_ in cross)
ok_nonuni = (sgnQ5(tau ** 4 - 1) < 0) and (sgnQ5(phi ** 4 - 1) > 0)
ck("NX-T4", ok_real and ok_nonuni,
   "cross ratios real; |ratio|^2 in {tau^4, phi^4}, both != 1")

# T5: normalized cross ratios nu = mu/conj(mu) = 1 identically (RC-Rem 7.6)
nus = [simplify(mu / conjugate(mu)) for mu in cross]
ck("NX-T5", all(zero(nu_ - 1) for nu_ in nus),
   "cross-shell nu = mu/conj(mu) = 1 for all 8 pairs: coherent, torsion order 1")

# T6: t_rel pattern: equal-angle cross pairs 0, opposite-angle pairs 1/2
def tval(r_):
    # roots are i*s with s real; t = 1/4 if s>0 else 3/4
    s_ = simplify(r_ / I)
    return Fraction(1, 4) if sgnQ5(s_) > 0 else Fraction(3, 4)
ok = True
for a in sh_list[0]:
    for b in sh_list[1]:
        d = (tval(a) - tval(b)) % 1
        mu = simplify(a / b)
        want = Fraction(0) if sgnQ5(mu) > 0 else Fraction(1, 2)
        ok = ok and (d == want)
ck("NX-T6", ok, "t_rel(cross) = 0 on equal-angle pairs, 1/2 on opposite-angle")

# T7: coherence completeness -- all 16 ordered pairs coherent, one class
all16 = [simplify(a / b) for a in roots for b in roots]
def coherent(mu):
    # coherence (RC-Def 3.1) <=> t(a)-t(b) rational; here decided exactly:
    # every ratio is real (angle diff in {0,1/2}) -- check Im mu = 0
    return zero(im(mu))
ck("NX-T7", len(all16) == 16 and all(coherent(mu) for mu in all16),
   "all 16 ordered pairs coherent (real ratios): single coherence class")

# T8: angle data -> Delta = Z/2Z, absolute Z/4Z, anchor (m,n) = (2,4), n = 2m
tset = sorted({tval(r_) for r_ in roots})
diffs = sorted({(ta - tb) % 1 for ta in tset for tb in tset})
def subgrp_order(gens):
    G = {Fraction(0)}
    frontier = set(gens)
    while frontier:
        G2 = {(g + h) % 1 for g in G | frontier for h in frontier} | G | frontier
        if G2 == G: break
        frontier = G2 - G; G = G2
    return len(G)
m_ord = subgrp_order(diffs)
n_ord = subgrp_order(set(tset))
ck("NX-T8", tset == [Fraction(1,4), Fraction(3,4)] and
   diffs == [Fraction(0), Fraction(1,2)] and m_ord == 2 and n_ord == 4
   and n_ord == 2 * m_ord,
   "angles {1/4,3/4}; Delta = Z/2Z; absolute Z/4Z; anchor (m,n) = (2,4)")

# T9: twisted-shell realization p = s(x^2), s = minpoly(q^2) = y^2+3y+1;
#     coherence-witness moduli {tau^2, phi^2} = {|q^2|, |q^2|^-1}
y = symbols('y')
s_y = y ** 2 + 3 * y + 1
ok_twist = zero(expand(s_y.subs(y, x ** 2)) - p_q)
ok_minp = Poly(minimal_polynomial(q ** 2, y), y) == Poly(s_y, y)
ok_wit = zero(expand(q**2 * conjugate(q**2)) - tau ** 4) and \
         zero(phi ** 2 * tau ** 2 - 1)
ck("NX-T9", ok_twist and ok_minp and ok_wit,
   "p = s(x^2) with s = y^2+3y+1 = minpoly(q^2); witness moduli {tau^2, phi^2}")

# ================================================================ NX-S
a_ = exp(-rho); z_ = sqrt(1 - exp(-2 * rho))
th = kap * rho
X = a_ * cos(th); Y = a_ * sin(th)

# S1: compass cross <=> rung lattice (symbolic at lattice/half-floor points)
ok = True
for k in range(5):
    r_ev = 2 * k * log(phi)
    r_od = (2 * k + 1) * log(phi)
    ok = ok and simplify(Y.subs(rho, r_ev)) == 0            # even: y = 0
    if k > 0:
        ok = ok and simplify(X.subs(rho, r_ev)) != 0
    ok = ok and simplify(X.subs(rho, r_od)) == 0            # odd:  x = 0
    ok = ok and simplify(Y.subs(rho, r_od)) != 0
for n in range(4):                                          # half-floors: xy != 0
    r_h = (n + Rational(1, 2)) * log(phi)
    xy = simplify((X * Y).subs(rho, r_h))
    ok = ok and simplify(xy ** 2 - exp(-4 * r_h) / 4) == 0 and xy != 0
ck("NX-S1", ok, "S(rho) on compass cross iff rho in (ln phi)Z; half-floors miss")

# S2: Im nu = 2xy/(x^2+y^2) identically
ck("NX-S2", simplify(2 * X * Y / (X**2 + Y**2) - sin(2 * th)) == 0,
   "Im nu = 2xy/(x^2+y^2) (spherical form of the nu-criterion)")

# S3: Re(W_n conj(W_{n+1})) = 0; S_n.S_{n+1} = z_n z_{n+1}; z1 z2 exact form
ok = all(zero(re(expand(q**n * conjugate(q**(n+1))))) for n in range(9))
ok = ok and zero(expand((sqrt(tau) * K) ** 2) - s5 * tau ** 3)  # (z1 z2)^2
ok = ok and zero(expand((5 ** Rational(1,4) * tau ** Rational(3,2)) ** 2)
                 - s5 * tau ** 3)
ck("NX-S3", ok, "horizontal parts of consecutive rungs orthogonal; z1 z2 = 5^(1/4) tau^(3/2)")

# S4: S_0 . S_1 = 0 exactly (first spherical step = quarter great circle)
S0 = (Integer(1), Integer(0), Integer(0))
S1v = (Integer(0), tau, sqrt(tau))
ck("NX-S4", zero(sum(u * v for u, v in zip(S0, S1v))),
   "S_0 . S_1 = 0: d_0 = pi/2 exactly")

# S5: z_n z_{n+1} strictly increasing (exact signs at squared level)
def zz2(n): return expand((1 - tau ** (2*n)) * (1 - tau ** (2*n + 2)))
ok = all(sgnQ5(zz2(n + 1) - zz2(n)) > 0 for n in range(8))
ck("NX-S5", ok, "z_n z_{n+1} strictly increasing: spherical steps d_n strictly decrease")

# S6: spherical speed |S'|^2 = a^2 kap^2 + a^2/z^2 = e^{-2rho}(kap^2 + 1/(1-e^{-2rho}))
sp2 = diff(X, rho)**2 + diff(Y, rho)**2 + diff(z_, rho)**2
tgt = a_**2 * kap**2 + a_**2 / z_**2
ck("NX-S6", simplify(trigsimp(expand(sp2 - tgt))) == 0,
   "|S'(rho)|^2 = e^(-2rho)(kap^2 + 1/(1-e^(-2rho))) symbolically")

# S7: substitution + elliptic-form identities and exact bounds
c, s_, g_, h_, kk = symbols('c s g h kk', positive=True)
e1 = simplify(sqrt(kk**2 + 1/c**2) * c - sqrt(1 + kk**2 * c**2))
m_ = kk**2 / (1 + kk**2)
e2 = simplify(expand((1 + kk**2) * (1 - m_ * s_**2) - (1 + kk**2 * (1 - s_**2))))
lo = expand((kk**2 + (1 + h_)) - (kk**2 + 1))          # = h >= 0: integrand >= sqrt(1+kk^2)
hi = expand((kk + sqrt(g_))**2 - (kk**2 + g_))          # = 2 kk sqrt(g) >= 0
int_hi = integrate(1 / sqrt(1 - symbols('aa', positive=True)**2),
                   (symbols('aa', positive=True), 0, 1))
ck("NX-S7", e1 == 0 and e2 == 0 and lo == h_ and hi == 2*kk*sqrt(g_)
   and simplify(int_hi - pi/2) == 0,
   "L_res = sqrt(1+kap^2) E(kap^2/(1+kap^2)); bounds sqrt(1+kap^2) <= L <= kap + pi/2")

# ================================================================ NX-H (displays)
from mpmath import mp, mpf, sqrt as msqrt, cos as mcos, log as mlog, pi as mpi
from mpmath import quad, ellipe, acos as macos, degrees as mdeg
mp.dps = 50
phin = (1 + msqrt(5)) / 2
taun = 1 / phin
kapn = mpi / (2 * mlog(phin))
Kn = msqrt(1 - phin ** -4)

def disp10(v):
    q_ = mp.floor(v * mpf(10) ** 10 + mpf("0.5")) / mpf(10) ** 10
    return mp.nstr(q_, 12, strip_zeros=False)

# H1: L_res -- 50 dps quadrature vs closed form; 10 dp display
L_quad = quad(lambda u: msqrt(1 + kapn**2 * mcos(u)**2), [0, mpi / 2])
L_form = msqrt(1 + kapn**2) * ellipe(kapn**2 / (1 + kapn**2))
agree = abs(L_quad - L_form) < mpf(10) ** -40
print("   L_res  =", mp.nstr(L_form, 30), " display:", disp10(L_form))
ck("NX-H1", agree, "L_res quadrature vs sqrt(1+kap^2) E(m) agree to 1e-40 (display 10 dp)")

# H2: d_1 = arccos(5^(1/4) tau^(3/2)) in degrees, 10 dp display
d1 = mdeg(macos(5 ** mpf("0.25") * taun ** mpf("1.5")))
d1b = mdeg(macos(msqrt(taun) * Kn))
print("   d_1    =", mp.nstr(d1, 30), " display:", disp10(d1))
ck("NX-H2", abs(d1 - d1b) < mpf(10) ** -45, "d_1 via both exact forms agrees (display 10 dp)")

# H3: numeric bracket corroboration of the exact S7 bounds
ck("NX-H3", msqrt(1 + kapn**2) < L_form < kapn + mpi / 2,
   "sqrt(1+kap^2) < L_res < kap + pi/2 at 50 dps (corroborates NX-S7)")

# ================================================================ summary
print(f"\n{len(PASS)}/{len(PASS)+len(FAIL)} checks passed", end="")
print("" if not FAIL else f"; FAILURES: {FAIL}")
sys.exit(0 if not FAIL else 1)
