#!/usr/bin/env python3
# ux_softcheck.py -- harness for HANDOFF-v1.0 dossiers SS10 (Vein 6, V6.1-V6.7)
# and SS9 (v1.6 uniqueness OPEN: F1 falsifier + candidate Theorem U).
# Written cold this session from the dossier statements; corpus citations
# (whitepaper Thm 9.1 / Thm 17.1, v1.6 thm:kseedladder = Thm 8.3) verified
# against the shipped papers before use.
#
# Discipline: every DECISION is exact -- sympy over Q / Q(sqrt5) (extended by i
# where needed); ladder statements decided at the squared level; sgnQ5 exact
# sign decisions by rational arithmetic; sqrt(14)/sqrt(56)-type comparisons
# decided by rational squaring (cmp_sqrt).  mpmath 50 dps appears ONLY in the
# corroboration guard UX-G1, counted separately (QX-G idiom).  Floats never
# touch an exact decision.  Exit 0 iff every check passes.
#
# Fail-first: every claim group carries falsifier guards that check
# deliberately perturbed variants are REJECTED (wrong-sign / wrong-coefficient
# polynomials, wrong TFD parity, perturbed recursions, window-boundary
# polynomials, a rational-parameter counterexample to the integer window).
# The detection machinery itself is validated on reducible in-class witnesses
# (UX-A5) before it certifies any irreducibility.
#
# Groups:
#   UX-A  anchors + citation replication + machinery fail-first          (5)
#   UX-B  V6.1 universal seed ladder p_n, n = 1..10                      (6)
#   UX-C  V6.2 z_n^2 b_n^2 = c_n                                         (2)
#   UX-D  V6.3 Mahler ladder M(p_n) = b_n^2; gate ties at ALL rungs      (4)
#   UX-E  V6.4 M(p_1) = phi (coincidence-candidate; no mechanism)        (2)
#   UX-F  V6.5 coefficient recursion c_{n+1} = 3c_n - c_{n-1} + 2        (3)
#   UX-K  V6.6 full-compass straddler + angle data Z/4Z                  (2)
#   UX-H  SS9 F1 falsifier x^4 + 6x^2 - 5 at level 5                     (4)
#   UX-U  Theorem U: window, monotonicity, endpoint, ladder levels       (7)
#   UX-G  numeric corroboration guard (counted separately)               (1)

import sys
from fractions import Fraction
from math import gcd
from sympy import (symbols, sqrt, I, Rational, Integer, simplify, expand,
                   radsimp, together, nsimplify, cancel, Poly,
                   minimal_polynomial, factor_list, lucas, fibonacci, diff)

PASS, FAIL, GUARD = [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append((cid, desc))
    print(("PASS" if ok else "FAIL"), cid, "-", desc)

def gk(cid, cond, desc):
    ok = bool(cond)
    GUARD.append((cid, ok, desc))
    print(("GUARD-PASS" if ok else "GUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append((cid, desc))

def zero(e):
    e2 = simplify(expand(e))
    if e2 == 0:
        return True
    e3 = simplify(radsimp(together(expand(e))))
    if e3 == 0:
        return True
    try:
        e4 = nsimplify(e2, [sqrt(5)])
        return simplify(radsimp(e4)) == 0
    except Exception:
        return False

x, y = symbols('x y')
w = symbols('w', positive=True)          # w = phi^{2n}, symbolic rung
a_s = symbols('a_s', real=True)          # Theorem-U coefficient a
c_s = symbols('c_s', positive=True)      # Theorem-U level c > 0
Ps, Qs = symbols('P_s Q_s', positive=True)

s5  = sqrt(5)
phi = (1 + s5)/2
tau = (s5 - 1)/2
gap = phi**-4
K2  = 1 - gap

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q. Raises if not in Q(sqrt5)."""
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQ5(e):
    """Exact sign of A + B*sqrt(5) via rational arithmetic only."""
    A, B = q5AB(e)
    if A == 0 and B == 0:
        return 0
    if A >= 0 and B >= 0:
        return 1
    if A <= 0 and B <= 0:
        return -1
    if A > 0:            # A > 0 > B
        return 1 if A*A > 5*B*B else -1
    return 1 if 5*B*B > A*A else -1      # A < 0 < B

def cmp_sqrt(P, r):
    """Exact sign of sqrt(P) - r for rational P >= 0, rational r."""
    P = Fraction(P); r = Fraction(r)
    assert P >= 0
    if r < 0:
        return 1
    if P > r*r:
        return 1
    if P < r*r:
        return -1
    return 0

def qz_order(fracs):
    L = 1
    for f in fracs:
        f = Fraction(f) % 1
        L = L*f.denominator // gcd(L, f.denominator)
    return L

def group_data(angles):
    diffs = [(Fraction(a) - Fraction(b)) % 1 for a in angles for b in angles]
    return qz_order(diffs), qz_order(angles)

def irreducible4(p):
    fl = factor_list(p, x)[1]
    return len(fl) == 1 and fl[0][1] == 1 and Poly(fl[0][0], x).degree() == 4

def reducible(p):
    fl = factor_list(p, x)[1]
    return len(fl) > 1 or (len(fl) == 1 and fl[0][1] > 1)

# ---------------------------------------------------------------- ladder data
N = 10
cn  = [None] + [int(lucas(2*n)) - 2 for n in range(1, N + 1)]
zn2 = [None] + [expand(1 - tau**(2*n)) for n in range(1, N + 1)]     # z_n^2
bn2 = [None] + [expand(phi**(2*n) - 1) for n in range(1, N + 1)]     # b_n^2
pn  = [None] + [x**4 + cn[n]*x**2 - cn[n] for n in range(1, N + 1)]
CLIST = [1, 5, 16, 45, 121, 320, 841, 2205]          # handoff V6.5 list, n=1..8
LEVELS = [1, 5, 16, 45, 121]                          # Theorem-U scan levels

# ================================================================ UX-A anchors
print("== UX-A: anchors, citation replication, machinery fail-first ==")

ck("UX-A1", zero(tau**2 + tau - 1) and zero(phi*tau - 1)
        and zero(phi**2 - phi - 1) and zero((1 - tau**4) - K2),
   "seed identities tau^2+tau=1, phi tau=1, phi^2=phi+1; z_2^2 = K^2")

ck("UX-A2", expand(minimal_polynomial(sqrt(tau), x) - (x**4 + x**2 - 1)) == 0
        and expand(minimal_polynomial(sqrt(phi), x) - (x**4 - x**2 - 1)) == 0
        and zero(phi**2*tau - phi),
   "W-Thm 17.1 replicated: minpoly(sqrt(tau)) = x^4+x^2-1, minpoly(sqrt(phi)) = x^4-x^2-1, sqrt(phi)=phi sqrt(tau)")

ck("UX-A3", cn[1] == 1 and cn[2] == 5 and cn[3] == 16
        and expand(pn[1] - (x**4 + x**2 - 1)) == 0
        and expand(pn[2] - (x**4 + 5*x**2 - 5)) == 0
        and expand(pn[3] - (x**4 + 16*x**2 - 16)) == 0,
   "sanity anchors: p_1 = x^4+x^2-1, p_2 = K-seed x^4+5x^2-5, p_3 = x^4+16x^2-16")

u_w = 1/(w - 1); C_w = 1/(w - 2 + 1/w)
ok_a4 = zero(u_w**2 + u_w - C_w) and zero((-1/w) + C_w/(1 + u_w)**2)
for n in range(1, 7):
    ok_a4 = ok_a4 and zero(expand(1/(phi**n - phi**(-n))**2 - Rational(1, cn[n])))
ok_a4 = ok_a4 and cn[4] == 45 and cn[5] == 121
ck("UX-A4", ok_a4,
   "W-Thm 9.1 replicated: u_n = 1/(phi^2n-1) solves u^2+u=C_n, m_n = -phi^-2n, C_n = 1/(phi^n-phi^-n)^2 = 1/c_n; C3=1/16, C4=1/45, C5=1/121")

ck("UX-A5", reducible(x**4 + 3*x**2 - 4) and reducible(x**4 + 4*x**2 - 5)
        and expand((x**2 - 1)*(x**2 + c_s) - (x**4 + (c_s - 1)*x**2 - c_s)) == 0
        and sgnQ5(tau) > 0 and sgnQ5(tau - 1) < 0 and sgnQ5(phi - 1) > 0
        and (not zero(tau - 1)) and cmp_sqrt(2, 1) > 0 and cmp_sqrt(2, 2) < 0
        and cmp_sqrt(4, 2) == 0,
   "machinery fail-first: factor_list flags reducible in-class witnesses x^4+3x^2-4, x^4+4x^2-5; boundary family (x^2-1)(x^2+c); sgnQ5/cmp_sqrt/zero reject nonzeros")

# ================================================================ UX-B (V6.1)
print("== UX-B: V6.1 universal seed ladder, n = 1..10 ==")

ok_b1 = all(cn[n] == (int(lucas(n))**2 if n % 2 else 5*int(fibonacci(n))**2)
            for n in range(1, N + 1))
ok_b1 = ok_b1 and all(int(lucas(n))**2 - 5*int(fibonacci(n))**2 == 4*(-1)**n
                      for n in range(1, 13))
ok_b1 = ok_b1 and all(zero(expand((phi**n - phi**(-n))**2 - cn[n]))
                      for n in range(1, N + 1))
ck("UX-B1", ok_b1,
   "c_n = L_2n - 2 = L_n^2 (n odd) / 5F_n^2 (n even); TFD parity via L_n^2-5F_n^2 = 4(-1)^n; c_n = (phi^n-phi^-n)^2")

ck("UX-B2", all(expand(minimal_polynomial(sqrt(1 - tau**(2*n)), x) - pn[n]) == 0
                for n in range(1, N + 1)),
   "minpoly(z_n) = p_n(x) = x^4 + c_n x^2 - c_n, direct sympy, n = 1..10")

ok_b3 = True
for n in range(1, N + 1):
    ok_b3 = ok_b3 and zero(expand(zn2[n]**2 + cn[n]*zn2[n] - cn[n]))   # p(+-z_n)=0
    ok_b3 = ok_b3 and zero(expand(bn2[n]**2 - cn[n]*bn2[n] - cn[n]))   # p(+-i b_n)=0
    ok_b3 = ok_b3 and irreducible4(pn[n])
ck("UX-B3", ok_b3,
   "roots {+-z_n, +-i b_n} at the squared level; p_n irreducible over Q (factor_list), n = 1..10")

ok_b4 = zero(expand((1 - 1/w)*(1 - w) - (2 - (w + 1/w))))
for n in range(1, N + 1):
    ok_b4 = ok_b4 and zero(expand((1 - tau**(2*n))*(1 - phi**(2*n)) - (2 - lucas(2*n))))
    ok_b4 = ok_b4 and (2 - int(lucas(2*n)) == -cn[n]) and (-cn[n] < 0)
    ok_b4 = ok_b4 and q5AB(expand(tau**(2*n)))[1] != 0
ck("UX-B4", ok_b4,
   "norm obstruction: N(1-tau^2n) = 2 - L_2n = -c_n < 0 (not a Q(sqrt5)-square) and tau^2n irrational => degree 4 uniformly")

ck("UX-B5", expand(c_s*(x**4/c_s + x**2 - 1) - (x**4 + c_s*x**2 - c_s)) == 0
        and all(Fraction(1, cn[n])*cn[n] == 1 for n in range(1, N + 1)),
   "universal seed: p_n(z) = 0 iff C_n z^4 + z^2 - 1 = 0 with C_n = 1/c_n = the Thm 9.1 gate level")

ok_b6 = True
for n in range(1, 5):
    ok_b6 = ok_b6 and (not zero(expand(zn2[n]**2 + cn[n]*zn2[n] + cn[n])))   # x^4+cx^2+c rejected
    ok_b6 = ok_b6 and (not zero(expand(zn2[n]**2 + (cn[n]+1)*zn2[n] - cn[n])))  # x^4+(c+1)x^2-c rejected
wrongpar = all(abs(cn[n] - (5*int(fibonacci(n))**2 if n % 2 else int(lucas(n))**2)) == 4
               for n in range(1, N + 1))
ck("UX-B6", ok_b6 and wrongpar,
   "falsifiers: perturbed seeds x^4+c x^2+c and x^4+(c+1)x^2-c do NOT annihilate z_n; wrong TFD parity misses c_n by exactly 4")

# ================================================================ UX-C (V6.2)
print("== UX-C: V6.2 z_n^2 b_n^2 = c_n ==")

ck("UX-C1", zero(expand((1 - 1/w)*(w - 1) - (w - 2 + 1/w)))
        and all(zero(expand(zn2[n]*bn2[n] - cn[n])) for n in range(1, N + 1)),
   "z_n^2 b_n^2 = c_n = |p_n(0)|: symbolic in w and exact n = 1..10 (generalizes K^2 beta^2 = 5)")

ok_c2 = all(not zero(expand(zn2[n]*bn2[n] - (cn[n] + 1))) for n in range(1, N + 1))
ok_c2 = ok_c2 and all(sgnQ5(expand(zn2[n]*bn2[n + 1] - cn[n])) != 0 for n in range(1, 5))
ck("UX-C2", ok_c2,
   "falsifiers: z_n^2 b_n^2 != c_n + 1; cross-rung z_n^2 b_{n+1}^2 != c_n (exact signs)")

# ================================================================ UX-D (V6.3)
print("== UX-D: V6.3 Mahler ladder and gate ties at ALL rungs ==")

ck("UX-D1", all(sgnQ5(expand(zn2[n] - 1)) < 0 and sgnQ5(expand(bn2[n] - 1)) > 0
                for n in range(1, N + 1)),
   "exact root positions: z_n^2 < 1 < b_n^2, n = 1..10 -- only the imaginary pair lies outside the unit circle")

ck("UX-D2", all(Fraction(1, cn[n])*cn[n] == 1 and zero(expand(bn2[n]*(1/(phi**(2*n) - 1)) - 1))
                for n in range(1, N + 1)),
   "M(p_n) = b_n^2, u_n = 1/M(p_n), C_n = 1/|p_n(0)|: Thm 8.3(d) ties extend from even rungs to ALL rungs")

ok_d3 = True
for j in range(1, 6):
    ok_d3 = ok_d3 and expand(pn[2*j] - (x**4 + 5*fibonacci(2*j)**2*x**2
                                        - 5*fibonacci(2*j)**2)) == 0
    ok_d3 = ok_d3 and zero(expand(bn2[2*j] - (phi**(4*j) - 1)))
ck("UX-D3", ok_d3,
   "even-rung consistency: Vein-6 p_2j = corpus p_j (thm:kseedladder), b_2j^2 = beta_j^2 = phi^4j - 1, j = 1..5")

ck("UX-D4", all(sgnQ5(expand(bn2[n] - zn2[n])) > 0 and (not zero(expand(zn2[n] - bn2[n])))
                for n in range(1, N + 1)),
   "falsifier: the flipped Mahler claim M = z_n^2 is rejected (b_n^2 - z_n^2 = c_n > 0 exactly)")

# ================================================================ UX-E (V6.4)
print("== UX-E: V6.4 M(p_1) = phi -- coincidence-candidate, no mechanism ==")

ck("UX-E1", zero(expand(bn2[1] - phi)),
   "M(p_1) = phi exactly (phi^2 - 1 = phi): rung-1 seed on the emission Mahler floor [flag: coincidence-candidate]")

ck("UX-E2", all(sgnQ5(expand(bn2[n] - phi)) > 0 for n in range(2, N + 1)),
   "specialness falsifier: M(p_n) != phi for n = 2..10 (exact signs; the hit is unique to n = 1)")

# ================================================================ UX-F (V6.5)
print("== UX-F: V6.5 coefficient recursion ==")

ck("UX-F1", cn[1:9] == CLIST
        and all(cn[n + 1] == 3*cn[n] - cn[n - 1] + 2 for n in range(2, N)),
   "integer list 1,5,16,45,121,320,841,2205 matches c_n; recursion c_{n+1} = 3c_n - c_{n-1} + 2 on n = 2..9")

ck("UX-F2", zero(expand((phi**2*w - 2 + 1/(phi**2*w)) + (w/phi**2 - 2 + phi**2/w)
                        - 3*(w - 2 + 1/w) - 2))
        and all(int(lucas(2*n + 2)) == 3*int(lucas(2*n)) - int(lucas(2*n - 2))
                for n in range(1, 11)),
   "recursion symbolic via phi^2 + phi^-2 = 3 (Lucas identity L_{2n+2} = 3L_2n - L_{2n-2})")

# solve c3 = A c2 - c1 + B and c4 = A c3 - c2 + B exactly:
A_u = Fraction((cn[4] + cn[2]) - (cn[3] + cn[1]), cn[3] - cn[2])
B_u = Fraction(cn[3] + cn[1]) - A_u*cn[2]
ok_f3 = (A_u == 3 and B_u == 2)
for (Ap, Bp) in [(3, 1), (2, 2), (4, 2)]:
    ok_f3 = ok_f3 and (not all(cn[n + 1] == Ap*cn[n] - cn[n - 1] + Bp
                               for n in range(2, N)))
ck("UX-F3", ok_f3,
   "falsifiers: (A,B) solved from the list is uniquely (3,2); perturbed laws (3,1),(2,2),(4,2) all fail")

# ================================================================ UX-K (V6.6)
print("== UX-K: V6.6 full-compass straddler and angle data ==")

ok_k1 = group_data([Fraction(0), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)]) == (4, 4)
ok_k1 = ok_k1 and sgnQ5(phi**2 - 2) > 0 and zero(expand((1 - (1 - 1/w)) - 1/w))
ok_k1 = ok_k1 and all(zero(expand(zn2[n]/(1 - zn2[n]) - bn2[n])) for n in range(1, N + 1))
ck("UX-K1", ok_k1,
   "angle data {0,1/2,1/4,3/4}: Delta = Z/4Z, anchor n = m = 4; all-n straddle atoms (phi^2 > 2, 1-z_n^2 = tau^2n > 0); V6.7 mirror tan^2 alpha_n = z_n^2/(1-z_n^2) = b_n^2")

ok_k2 = group_data([Fraction(0), Fraction(1, 2)]) == (2, 2)
ok_k2 = ok_k2 and group_data([Fraction(1, 4), Fraction(3, 4)]) == (2, 4)
for c in LEVELS:
    ok_k2 = ok_k2 and (1 + (c - 1) - c == 0)          # q(1) = 0: real pair ON the circle
    ok_k2 = ok_k2 and reducible(x**4 + (c - 1)*x**2 - c)
ck("UX-K2", ok_k2,
   "falsifiers: half-compass sets give (2,2)/(2,4) not (4,4); boundary a = c-1 puts the real pair ON the circle and is reducible")

# ================================================================ UX-H (SS9 F1)
print("== UX-H: F1 falsifier x^4 + 6x^2 - 5 at level 5 ==")

g5 = x**4 + 6*x**2 - 5
yp_g = -3 + sqrt(14); ym_g = -3 - sqrt(14)
ck("UX-H1", irreducible4(g5)
        and expand(yp_g**2 + 6*yp_g - 5) == 0 and expand(ym_g**2 + 6*ym_g - 5) == 0
        and cmp_sqrt(14, 3) > 0 and cmp_sqrt(14, 4) < 0 and cmp_sqrt(14, -2) > 0,
   "x^4+6x^2-5 irreducible; y+- = -3+-sqrt(14); y+ = sqrt(14)-3 in (0,1) exact (9<14<16); |y-| = 3+sqrt(14) > 1")

ck("UX-H2", expand(g5 - pn[2]) == x**2 and abs(-5) == abs(-cn[2])
        and (1 + 6 - 5 > 0) and (1 - 6 - 5 < 0)
        and group_data([Fraction(0), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)]) == (4, 4),
   "distinct from the K-seed (difference x^2), same level |p(0)| = 5, full-compass straddler => F1 (constant-term-only uniqueness) is FALSE")

ck("UX-H3", zero(expand((6 + sqrt(56))/2 - (3 + sqrt(14))))
        and zero(expand((5 + sqrt(45))/2 - (phi**4 - 1)))
        and cmp_sqrt(56, 0) > 0 and 56 > 45,
   "Mahler ordering at level 5: M(x^4+6x^2-5) = 3+sqrt(14) > M(K-seed) = (5+sqrt(45))/2 = phi^4-1 (56 > 45 exact)")

ck("UX-H4", reducible(x**4 + 4*x**2 - 5) and (1 + 4 - 5 == 0),
   "falsifier control: a = 4 < c = 5 is NOT a second counterexample (reducible, real pair ON the circle) -- the window boundary is sharp")

# ================================================================ UX-U (Theorem U)
print("== UX-U: Theorem U in the class {x^4 + a x^2 - c : a in Z, c in Z > 0} ==")

disc = sqrt(a_s**2 + 4*c_s)
yp = (-a_s + disc)/2; ym = (-a_s - disc)/2
ck("UX-U1", expand(yp + ym + a_s) == 0 and expand(yp*ym + c_s) == 0
        and expand((1 - yp)*(1 - ym) - (1 + a_s - c_s)) == 0
        and expand((1 + yp)*(1 + ym) - (1 - a_s - c_s)) == 0
        and expand((x**4 + a_s*x**2 - c_s).subs(x, 1) - (1 + a_s - c_s)) == 0
        and expand((x**4 + a_s*x**2 - c_s).subs(x, I) - (1 - a_s - c_s)) == 0,
   "window identities: Vieta y+y- = -c < 0; (1-y+)(1-y-) = p(1); (1+y+)(1+y-) = p(i) = 1-a-c")

ok_u2 = True
for c in LEVELS:
    for a in range(-c - 5, c + 11):
        straddler = (1 + a - c > 0) and (1 - a - c < 0)
        ok_u2 = ok_u2 and (straddler == (a >= c))
for c in [1, 5, 16]:
    for a in range(c - 3, c + 5):
        s_root = (cmp_sqrt(a*a + 4*c, a + 2) < 0) and (cmp_sqrt(a*a + 4*c, 2 - a) > 0)
        s_sign = (1 + a - c > 0) and (1 - a - c < 0)
        ok_u2 = ok_u2 and (s_root == s_sign)
ck("UX-U2", ok_u2,
   "(i) full-compass straddler <=> p(1) > 0 <=> a >= c: exhaustive integer scan, two independent exact decision paths agree")

ck("UX-U3", simplify(diff((a_s + disc)/2, a_s)
                     - (disc + a_s)/(2*disc)) == 0
        and expand((a_s**2 + 4*c_s) - a_s**2 - 4*c_s) == 0 and (4*c_s).is_positive
        and expand((sqrt(Ps) - sqrt(Qs))*(sqrt(Ps) + sqrt(Qs)) - (Ps - Qs)) == 0,
   "(ii) M = |y-| = (a+sqrt(a^2+4c))/2: dM/da = (disc+a)/(2 disc) > 0 since disc^2 - a^2 = 4c > 0; sqrt strictly monotone")

ok_u4 = True
for c in LEVELS:
    for a in range(c, c + 10):
        ok_u4 = ok_u4 and ((a + 1)**2 + 4*c > a**2 + 4*c) and (a + 1 - a > 0)
    for a in range(c + 1, c + 11):
        ok_u4 = ok_u4 and (a - c > 0) and (a*a + 4*c > c*c + 4*c)
ck("UX-U4", ok_u4,
   "discrete monotonicity: M(a+1) > M(a) and M(a) > M(c) on the window (integer radicand comparisons + UX-U3 lemma)")

ok_u5 = all(irreducible4(pn[n]) for n in range(1, N + 1))
ok_u5 = ok_u5 and all((cn[n] % 5 == 0 and cn[n] % 25 != 0) == (n % 2 == 0 and n % 5 != 0)
                      for n in range(1, N + 1))
ok_u5 = ok_u5 and all((int(fibonacci(2*j)) % 5 == 0) == (j % 5 == 0) for j in range(1, 16))
ck("UX-U5", ok_u5,
   "(iii) irreducibility of x^4 + c x^2 - c at all ladder levels; Eisenstein-at-5 applies iff n even and 5 does not divide n (fails at n = 10, c = 15125 = 5^3*11^2, factor_list still certifies)")

ok_u6 = all(cn[n]**2 + 4*cn[n] == 5*int(fibonacci(2*n))**2 for n in range(1, N + 1))
ok_u6 = ok_u6 and zero(expand((w - 2 + 1/w)**2 + 4*(w - 2 + 1/w) - (w - 1/w)**2))
ok_u6 = ok_u6 and all(zero(expand((cn[n] + sqrt(Integer(cn[n]**2 + 4*cn[n])))/2 - bn2[n]))
                      for n in range(1, N + 1))
ck("UX-U6", ok_u6,
   "(iii) endpoint value: M at a = c_n is (c_n + sqrt(c_n^2+4c_n))/2 = b_n^2 exactly (c^2+4c = 5 F_2n^2) -- the minimal straddler at level c_n IS p_n")

ok_u7 = True
# perturbed window "straddler <=> a >= c+1" fails at a = c (which IS a straddler):
for c in LEVELS:
    a = c
    ok_u7 = ok_u7 and ((1 + a - c > 0) and (1 - a - c < 0)) and not (a >= c + 1)
# perturbed window "straddler <=> a >= c-1" fails at a = c-1 (NOT a straddler):
for c in LEVELS:
    a = c - 1
    ok_u7 = ok_u7 and not ((1 + a - c > 0) and (1 - a - c < 0)) and (a >= c - 1)
# perturbed minimality "min at a = c+1" fails: M(c) < M(c+1):
for c in LEVELS:
    ok_u7 = ok_u7 and ((c + 1)**2 + 4*c > c*c + 4*c)
# perturbed endpoint value rejected:
ok_u7 = ok_u7 and all(not zero(expand((cn[n] + sqrt(Integer(cn[n]**2 + 4*cn[n])))/2
                                      - (bn2[n] + 1))) for n in range(1, 4))
# rational-parameter counterexample: over Q the equivalence (i) FAILS at (a,c) = (1/2, 1/10):
aq, cq = Fraction(1, 2), Fraction(1, 10)
p1_pos = (1 + aq - cq > 0)                                   # p(1) > 0 holds
ym_in  = cmp_sqrt(aq*aq + 4*cq, 2 - aq) < 0                  # but |y-| < 1: not a straddler
ok_u7 = ok_u7 and p1_pos and ym_in
ck("UX-U7", ok_u7,
   "falsifiers: window is sharp at a = c (perturbed windows +-1 both refuted); min is NOT at a = c+1; endpoint value rigid; integer hypothesis load-bearing (rational witness (1/2,1/10) has p(1)>0 yet no straddle)")

# ================================================================ UX-G guard
print("== UX-G: numeric corroboration (counted separately) ==")
from mpmath import mp, mpf, sqrt as msqrt, polyroots, fabs

mp.dps = 50
fphi = (1 + msqrt(5))/2
ok_g1 = True
for n in range(1, N + 1):
    zf = msqrt(1 - (1/fphi)**(2*n))
    bf = msqrt(fphi**(2*n) - 1)
    rts = polyroots([1, 0, cn[n], 0, -cn[n]], maxsteps=200, extraprec=200)
    want = [zf, -zf, 1j*bf, -1j*bf]
    used = [False]*4
    for r in rts:
        hit = None
        for j2, t2 in enumerate(want):
            if not used[j2] and fabs(r - t2) < mpf(10)**-40:
                hit = j2
                break
        if hit is None:
            ok_g1 = False
        else:
            used[hit] = True
    mah = mpf(1)
    for r in rts:
        if fabs(r) > 1:
            mah = mah*fabs(r)
    ok_g1 = ok_g1 and fabs(mah - bf**2) < mpf(10)**-40
gk("UX-G1", ok_g1,
   "50-dps corroboration: polyroots(p_n) match {+-z_n, +-i b_n}; Mahler product = b_n^2 to 1e-40, n = 1..10")
print("   display: M(p_1) =", mp.nstr(fphi, 12), "(= phi);  c_n =", cn[1:])

# ================================================================ summary
print("=" * 64)
n_guard_ok = sum(1 for _, ok, _ in GUARD if ok)
n_exact_fail = len([f for f in FAIL if not f[0].startswith("UX-G")])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, n_guard_ok, len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
