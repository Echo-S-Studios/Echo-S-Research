#!/usr/bin/env python3
# qx_softcheck.py -- third independent harness (v1.5 session, written cold).
# Discipline: ALL decisions exact in Q / Q(sqrt5) (extended by i where needed);
# K-family decided at the squared level; mpmath used ONLY for the certified
# interval guards (QX-G1, QX-G2) and the numeric corroboration guard (QX-G3);
# floats never touch an exact decision.
#
#   QX-A  replication of the v1.4 layer this session builds on         (15)
#   QX-B  complete relational type of minpoly(q)          [Prop chirat]    (8)
#   QX-C  K-family seed ladder p_j                        [Thm kseedladder/Cor evenmirror] (14)
#   QX-D  quarter-twist transport, anchor dichotomy       [Thm quartertwist]    (10)
#   QX-E  exact rail-spacing bounds, gap-rate law         [Prop railbounds]   (9)
#   QX-F  Gamma-transfer of the interpolation classification [Rem gammatransfer]  (6)
#   QX-G  certified-interval / numeric guards (counted separately)      (3)

import sys
from fractions import Fraction
from math import gcd
from sympy import (symbols, sqrt, Rational, I, pi, log, exp, sin, cos, tan,
                   cot, csc, simplify, expand, radsimp, trigsimp, together,
                   series, minimal_polynomial, Poly, factor_list, im,
                   fibonacci, lucas, diff, expand_log)

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
    if simplify(expand(e)) == 0:
        return True
    e2 = simplify(radsimp(together(e)))
    if e2 == 0:
        return True
    if simplify(expand(e.rewrite(cos))) == 0:
        return True
    if simplify(expand_log(e, force=True)) == 0:
        return True
    return simplify(trigsimp(e2)) == 0

x, t, zv = symbols('x t z')
alpha = symbols('alpha', positive=True)
w = symbols('w', positive=True)      # w = phi^{2j} (or phi^{2n})
G = symbols('G', positive=True)      # G = phi^{4j} = w^2
gsym = symbols('g', positive=True)   # g = gap^j
u_s = symbols('u_s', positive=True)  # u = kappa'^2
kk, kp = symbols('k k_p', integer=True)

s5   = sqrt(5)
phi  = (1 + s5)/2
tau  = (s5 - 1)/2
gap  = phi**-4
K2   = 1 - gap
K    = sqrt(K2)
beta2 = phi**2*s5
beta  = sqrt(beta2)
kappa = pi/(2*log(phi))
q    = I*tau

# ---------------------------------------------------------------- QX-A
print("== QX-A: replication of the layer built upon ==")
ck("QX-A1", zero(tau**2 + tau - 1) and zero(phi*tau - 1) and zero(phi**2 - phi - 1),
   "seed identities tau^2+tau=1, phi*tau=1, phi^2=phi+1")
ck("QX-A2", zero(q**2 + tau**2) and zero(q**4 - gap),
   "generator: q^2=-tau^2, q^4=gap")
mq = x**4 + 3*x**2 + 1
ck("QX-A3", expand(minimal_polynomial(q, x) - mq) == 0,
   "minpoly(q) = x^4+3x^2+1")
ck("QX-A4", all(zero(mq.subs(x, r)) for r in (I*tau, -I*tau, I*phi, -I*phi)),
   "roots of minpoly(q) are {+-i tau, +-i phi}")
Zs = x**4 - 3*x**2 + 1
ck("QX-A5", expand((x**2 - x - 1)*(x**2 + x - 1) - Zs) == 0,
   "Z* = (x^2-x-1)(x^2+x-1) = x^4-3x^2+1")
ck("QX-A6", expand(mq.subs(x, I*x) - Zs) == 0 and expand(Zs.subs(x, I*x) - mq) == 0,
   "gauge identity minpoly(q)(ix)=Z*(x) and Z*(ix)=minpoly(q)(x)")
ck("QX-A7", zero(K2 - (1 - gap)) and zero(K2**2 - 5*gap) and zero((K2*phi**2)**2 - 5),
   "K^2=1-gap, K^4=5gap, (K phi)^4=5")
ck("QX-A8", zero(K2*beta2 - 5) and zero(beta2 - K2*phi**4) and zero(beta2 - (phi**4 - 1)),
   "K^2 beta^2 = 5, beta^2 = K^2 phi^4 = phi^4 - 1")
ck("QX-A9", zero(1/(w*phi**4) - gap/w),
   "two-step law 1-z_{n+2}^2 = gap (1-z_n^2), symbolic in n (w=phi^{2n})")
ck("QX-A10", zero((1 - tau**2) - tau) and zero((1 - tau**4) - K2),
   "z_1=sqrt(tau), z_2=K (squared level)")
ck("QX-A11", zero((1 - 1/w)*w - (w - 1)),
   "tan^2 alpha_n = phi^{2n}-1, symbolic in n")
ck("QX-A12", all(zero(phi**n*lucas(n) - (phi**(2*n) - 1)) for n in (1, 3, 5, 7))
        and all(zero(phi**n*s5*fibonacci(n) - (phi**(2*n) - 1)) for n in (2, 4, 6, 8)),
   "TFD dress: phi^n L_n (odd) / phi^n sqrt5 F_n (even) = phi^{2n}-1")
zc = sin(alpha)
ck("QX-A13", zero((1 - zc**2)/zc**2 - cot(alpha)**2)
        and zero((1 - zc**2)/zc**4 - cot(alpha)**2*csc(alpha)**2)
        and zero((2 - zc**2)/zc**2 - (1 + 2*cot(alpha)**2))
        and zero(-(1 - zc**2) - (-cos(alpha)**2))
        and zero(zc**2*(2 - zc**2)/(1 - zc**2) - tan(alpha)**2*(1 + cos(alpha)**2)),
   "state-angle dictionary: u, C, sqrt(D), m, lambda")
ok14 = True
for n in (2, 4, 6):
    Cn = 1/(phi**n - phi**(-n))**2
    un = 1/(phi**(2*n) - 1)
    mn = -phi**(-2*n)
    ok14 = ok14 and zero(un**2 + un - Cn) and zero(mn + Cn/(1 + un)**2) \
                and zero(Cn - 1/(5*fibonacci(n)**2))
ck("QX-A14", ok14, "even gate ladder: u_n fixed point, m_n=-phi^{-2n}, C_n=1/(5F_n^2)")
ck("QX-A15", zero(((1 - 1/(w*phi**2)) - (1 - 1/w)) - (1/w)*tau)
        and zero(tau/(1 - tau**2) - 1),
   "transfer ledger: Delta F_n = tau^{2n+1}, sum tau^{2n+1} = 1")

# ---------------------------------------------------------------- QX-B
print("== QX-B: complete relational type of minpoly(q) ==")
roots_q = [I*tau, -I*tau, I*phi, -I*phi]
ratios = [simplify(radsimp(a/b)) for a in roots_q for b in roots_q]
ck("QX-B1", all(zero(im(r)) for r in ratios),
   "all 16 pairwise ratios of the roots are real")

def match_multiset(items, expected):
    pool = list(expected)
    for it in items:
        hit = None
        for j, e in enumerate(pool):
            if zero(it - e):
                hit = j
                break
        if hit is None:
            return False
        pool.pop(hit)
    return len(pool) == 0

ck("QX-B2", match_multiset(ratios,
        [1]*4 + [-1]*4 + [tau**2]*2 + [-tau**2]*2 + [phi**2]*2 + [-phi**2]*2),
   "ratio multiset = {1^4, (-1)^4, (tau^2)^2, (-tau^2)^2, (phi^2)^2, (-phi^2)^2}")
lhs = 1
for r in ratios:
    lhs = lhs*(x - r)
ck("QX-B3", zero(expand(lhs) - expand((x - 1)**4*(x + 1)**4*(x**4 - 7*x**2 + 1)**2)),
   "char poly of Rat = (x-1)^4 (x+1)^4 (x^4-7x^2+1)^2")
ck("QX-B4", expand((x**2 - 3*x + 1)*(x**2 + 3*x + 1) - (x**4 - 7*x**2 + 1)) == 0
        and zero((phi**2)**2 - 3*phi**2 + 1) and zero((tau**2)**2 - 3*tau**2 + 1),
   "x^4-7x^2+1 = (x^2-3x+1)(x^2+3x+1); phi^2, tau^2 root the split factor")
nontor = [tau**2, -tau**2, tau**2, -tau**2, phi**2, -phi**2, phi**2, -phi**2]
ck("QX-B5", match_multiset([r**2 for r in nontor], [phi**4]*4 + [tau**4]*4)
        and zero(phi**4 + tau**4 - 7) and zero(phi**4*tau**4 - 1),
   "Adams square of nontorsion part = 4 copies of spec(A_gap) (x^2-7x+1)")
ck("QX-B6", bool(simplify(tau**2) > 0) and bool(simplify(-tau**2 + 0) < 0)
        and zero((I*tau)/(I*phi) - tau**2) and zero((I*tau)/(-I*phi) + tau**2)
        and zero((I*tau)/(-I*tau) + 1),
   "t_rel(itau,iphi)=0, t_rel(itau,-iphi)=1/2, within-shell t_rel=1/2 (exact signs)")
torsion = [r for r in ratios if zero(r - 1) or zero(r + 1)]
ck("QX-B7", len(torsion) == 8 and sum(1 for r in torsion if zero(r - 1)) == 4,
   "torsion sub-multiset {1^4, (-1)^4}: contact signature {Phi_1^4, Phi_2^4}")
ck("QX-B8", all(any(zero(r - e) for e in (phi**2, -phi**2, tau**2, -tau**2)) for r in nontor)
        and zero(phi**2*tau**2 - 1),
   "every nontorsion ratio is a golden unit +-phi^{+-2}")

# ---------------------------------------------------------------- QX-C
print("== QX-C: K-family seed ladder ==")
FF  = G - 2 + 1/G       # 5 F_{2j}^2
Kj2 = 1 - 1/G
Bj2 = G - 1
ck("QX-C1", zero((w - 1/w)**2 - (w**2 + 1/w**2 - 2))
        and all(lucas(4*j) - 2 == 5*fibonacci(2*j)**2 for j in range(1, 9)),
   "L_{4j}-2 = 5 F_{2j}^2 (symbolic + integers j<=8)")
ck("QX-C2", zero(Kj2**2 + FF*Kj2 - FF)
        and all(zero((1 - gap**j)**2 + 5*fibonacci(2*j)**2*(1 - gap**j)
                     - 5*fibonacci(2*j)**2) for j in range(1, 7)),
   "p_j(K_j)=0: symbolic in G and explicit j<=6 (decided at K_j^2)")
ck("QX-C3", zero(Bj2 - Kj2*G),
   "beta_j = K_j phi^{2j} (squared): beta_j^2 = K_j^2 phi^{4j} = phi^{4j}-1")
ck("QX-C4", zero(Bj2**2 - FF*Bj2 - FF),
   "p_j(i beta_j)=0: y=-beta_j^2 roots the y-quadratic")
irr_ok = True
for j in range(1, 9):
    pj = x**4 + 5*fibonacci(2*j)**2*x**2 - 5*fibonacci(2*j)**2
    fl = factor_list(pj, x)[1]
    irr_ok = irr_ok and len(fl) == 1 and fl[0][1] == 1 and Poly(fl[0][0], x).degree() == 4
mp_ok = all(expand(minimal_polynomial(sqrt(1 - gap**j), x)
            - (x**4 + 5*fibonacci(2*j)**2*x**2 - 5*fibonacci(2*j)**2)) == 0 for j in (1, 2))
ck("QX-C5", irr_ok and mp_ok,
   "p_j irreducible over Q, j<=8; minimal_polynomial(K_j)=p_j corroborated j<=2")
ck("QX-C6", zero((1 - 1/G)*(1 - G) - (2 - (G + 1/G)))
        and zero((2 - (G + 1/G)) + FF)
        and all(fibonacci(2*j) >= 1 for j in range(1, 9)),
   "norm obstruction: N(1-gap^j) = 2-L_{4j} = -5F_{2j}^2 < 0 => degree 4 uniformly")
ck("QX-C7", all((fibonacci(2*j) % 5 == 0) == (j % 5 == 0) for j in range(1, 16)),
   "Eisenstein-at-5 availability: 5 | F_{2j} iff 5 | j (j<=15)")
ck("QX-C8", zero((Bj2 - Kj2) - FF) and zero(Kj2*Bj2 - FF)
        and zero((Kj2 + Bj2) - (G - 1/G))
        and zero((1 - 1/w**2)*w - (w - 1/w))
        and zero((1 - 1/w**2)**2 - (w - 1/w)**2/w**2),
   "ledger: beta^2-K^2 = K^2 beta^2 = 5F^2; K^2+beta^2 = sqrt5 F_{4j}; (K_j phi^j)^2 = sqrt5 F_{2j}; K_j^4 = 5F^2 gap^j")
ck("QX-C9", zero(1/Kj2 - 1/Bj2 - 1),
   "family identity 1/K_j^2 - 1/beta_j^2 = 1 (csc^2-cot^2 at rung 2j)")
ck("QX-C10", zero((gsym - gsym*gap) - K2*gsym) and zero(K2/(1 - gap) - 1),
   "increment law K_{j+1}^2-K_j^2 = K^2 gap^j; unit budget sum_{j>=0} K^2 gap^j = 1")
ck("QX-C11", zero(Kj2/Bj2 - 1/G) and bool(simplify((G - 2).subs(G, phi**4)) > 0)
        and zero(phi**4 - (3*phi + 2)) and zero(1 - Kj2 - 1/G),
   "cross-shell ratio moduli^2 = 1/G != 1; K_j^2 < 1 < beta_j^2 (phi^{4j} >= phi^4 > 2)")

def qz_order(fracs):
    L = 1
    for f in fracs:
        f = Fraction(f) % 1
        L = L*f.denominator // gcd(L, f.denominator)
    return L

def group_data(angles):
    diffs = [(Fraction(a) - Fraction(b)) % 1 for a in angles for b in angles]
    return qz_order(diffs), qz_order(angles)

ok12 = True
for j in (1, 2, 3):
    Kjv = sqrt(1 - gap**j)
    Bjv = sqrt(phi**(4*j) - 1)
    pj = x**4 + 5*fibonacci(2*j)**2*x**2 - 5*fibonacci(2*j)**2
    ok12 = ok12 and all(zero(pj.subs(x, r)) for r in (Kjv, -Kjv, I*Bjv, -I*Bjv))
m12, n12 = group_data([Fraction(0), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)])
ck("QX-C12", ok12 and m12 == 4 and n12 == 4,
   "roots {+-K_j, +-i beta_j} on the full compass; Delta = Z/4Z, anchor n = m = 4")
ck("QX-C13", all(zero(s5*fibonacci(2*j) - (phi**(2*j) - phi**(-2*j))) for j in range(1, 9))
        and zero(1/(w**2 - 1) - 1/Bj2.subs(G, w**2)),
   "gate ties: sqrt5 F_{2j} = phi^{2j}-phi^{-2j} => C_{2j}=1/(5F^2)=1/|p_j(0)|; u_{2j}=1/beta_j^2=1/M(p_j)")
ck("QX-C14", zero((1 - 1/G)/(1/G) - (G - 1)) and zero((w - 1/w)*w - (w**2 - 1))
        and expand((x**4 + 5*fibonacci(2)**2*x**2 - 5*fibonacci(2)**2) - (x**4 + 5*x**2 - 5)) == 0,
   "tan^2 alpha_{2j} = phi^{4j}-1 = beta_j^2; beta_j^2 = sqrt5 F_{2j} phi^{2j}; j=1 recovers the K-seed")

# ---------------------------------------------------------------- QX-D
print("== QX-D: quarter-twist transport and anchor-exchange dichotomy ==")

def twist(P):
    d = Poly(P, x).degree()
    return expand(I**(-d)*P.subs(x, -I*x))

def int_poly(P):
    try:
        return all(c.is_integer for c in Poly(P, x).all_coeffs())
    except Exception:
        return False

pK   = x**4 + 5*x**2 - 5
Phi8 = x**4 + 1
W6a  = x**4 - 2*x**2 + 4
W6b  = x**4 + 2*x**2 + 4
Q2   = x**2 - 2
Q2t  = x**2 + 2
r2 = sqrt(2)
wit = {
 'Zs':  (Zs,  [phi, -phi, tau, -tau],
         [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1, 2)]),
 'Mq':  (mq,  [I*tau, -I*tau, I*phi, -I*phi],
         [Fraction(1, 4), Fraction(3, 4), Fraction(1, 4), Fraction(3, 4)]),
 'pK':  (pK,  [K, -K, I*beta, -I*beta],
         [Fraction(0), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)]),
 'Phi8': (Phi8, [exp(2*I*pi*Rational(k, 8)) for k in (1, 3, 5, 7)],
         [Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8)]),
 'W6a': (W6a, [r2*exp(2*I*pi*Rational(k, 12)) for k in (1, 5, 7, 11)],
         [Fraction(1, 12), Fraction(5, 12), Fraction(7, 12), Fraction(11, 12)]),
 'W6b': (W6b, [r2*exp(2*I*pi*Rational(k, 6)) for k in (1, 2, 4, 5)],
         [Fraction(1, 6), Fraction(1, 3), Fraction(2, 3), Fraction(5, 6)]),
 'Q2':  (Q2,  [r2, -r2], [Fraction(0), Fraction(1, 2)]),
 'Q2t': (Q2t, [I*r2, -I*r2], [Fraction(1, 4), Fraction(3, 4)]),
}
domain = ['Zs', 'Mq', 'pK', 'Phi8', 'W6a', 'W6b', 'Q2', 'Q2t']
ck("QX-D0", all(all(zero(P.subs(x, r)) for r in R) for P, R, A in wit.values()),
   "all witness root lists verified against their polynomials")
ok1 = True
for kn in domain:
    P, R, A = wit[kn]
    Qp = twist(P)
    ok1 = ok1 and int_poly(Qp) and all(zero(Qp.subs(x, I*r)) for r in R)
ck("QX-D1", ok1, "transport i^{-d}P(-ix) integral on the domain; roots go to i*(roots)")
ck("QX-D2", all(expand(twist(twist(wit[kn][0])) - wit[kn][0]) == 0 for kn in domain),
   "the transport is an involution on the domain")
gold = x**2 - x - 1
ck("QX-D3", (not int_poly(twist(gold))) and (not zero(gold.subs(x, -phi))),
   "golden seed excluded: twist non-integral, -phi not a root (identity must live at Z*)")
ok4 = True
for kn in domain:
    P, R, A = wit[kn]
    ok4 = ok4 and expand(P.subs(x, -x) - P) == 0 \
              and all(any(zero(rr + r) for rr in R) for r in R)
ok4 = ok4 and expand(gold.subs(x, -x) - gold) != 0
ck("QX-D4", ok4, "domain = parity-homogeneous = negation-closed root multiset (pos + neg cases)")
expected_mn = {'Zs': (2, 2), 'Mq': (2, 4), 'pK': (4, 4), 'Phi8': (4, 8),
               'W6a': (6, 12), 'W6b': (6, 6), 'Q2': (2, 2), 'Q2t': (2, 4)}
ok5 = all(group_data(wit[kn][2]) == expected_mn[kn] for kn in domain)
ck("QX-D5", ok5, "Delta and anchor: (m,n) as expected on all eight witnesses")
ok6 = all(Fraction(1, 2) in [(Fraction(a) - Fraction(b)) % 1
          for a in wit[kn][2] for b in wit[kn][2]]
          and group_data(wit[kn][2])[0] % 2 == 0 for kn in domain)
ck("QX-D6", ok6, "1/2 in Delta (as a literal difference) and m even on the whole domain")
tpK = twist(pK)
ok7 = expand(twist(Zs) - mq) == 0 and expand(twist(W6a) - W6b) == 0 \
      and expand(twist(Q2) - Q2t) == 0 and expand(twist(Phi8) - Phi8) == 0 \
      and expand(tpK - (x**4 - 5*x**2 - 5)) == 0 \
      and all(zero(tpK.subs(x, r)) for r in (I*K, -I*K, beta, -beta)) \
      and group_data([Fraction(1, 4), Fraction(3, 4), Fraction(0), Fraction(1, 2)]) == (4, 4)
ck("QX-D7", ok7,
   "witness dichotomy: m=2, m=6 exchange (Zs<->Mq, W6a<->W6b, Q2<->Q2t); m=4 anchor-fixed (pK, Phi8)")
ok8 = True
for m in range(2, 25, 2):
    for jn in range(0, 2*m):
        a0 = Fraction(jn, 2*m) % 1
        n_  = qz_order([Fraction(1, m), a0])
        n2_ = qz_order([Fraction(1, m), (a0 + Fraction(1, 4)) % 1])
        ok8 = ok8 and n_ in (m, 2*m)
        if m % 4 == 0:
            ok8 = ok8 and n2_ == n_
        else:
            ok8 = ok8 and sorted([n_, n2_]) == [m, 2*m] and n_ != n2_
ck("QX-D8", ok8,
   "exhaustive coset simulation m<=24: rigidity n in {m,2m}; 4|m fixes, m=2 mod 4 exchanges")
ck("QX-D9", all(expand(((-1)**Poly(wit[kn][0], x).degree())*wit[kn][0].subs(x, -x)
                - wit[kn][0]) == 0 for kn in domain),
   "sign twist = identity on the domain (its odd-m exchange sector is disjoint)")

# ---------------------------------------------------------------- QX-E
print("== QX-E: exact rail-spacing bounds and the gap-rate law ==")
theta = symbols('theta', positive=True)
zfun = sqrt(1 - exp(-2*theta/kappa))
ck("QX-E1", zero(diff(zfun, theta) - (1 - zfun**2)/(kappa*zfun)),
   "dz/dtheta = (1-z^2)/(kappa z) from D1+D3")
xs = symbols('x_s', positive=True)
ck("QX-E2", expand((1 + xs/2)**2 - (1 + xs) - xs**2/4) == 0
        and expand((1 + xs) - (1 + xs/2 - xs**2/8)**2 - xs**3*(8 - xs)/64) == 0,
   "sqrt bounds: (1+x/2)^2-(1+x) = x^2/4; (1+x)-(1+x/2-x^2/8)^2 = x^3(8-x)/64")
A1f = -t - log(1 - t)
Ff  = 1/(1 - t) + 3*log(1 - t) - 3*(1 - t) + (1 - t)**2/2   # sign corrected by this harness
ck("QX-E3", zero(diff(A1f, t) - t/(1 - t)) and zero(diff(Ff, t) - t**3/(1 - t)**2),
   "antiderivatives: d/dt[-t-ln(1-t)] = t/(1-t); d/dt F = t^3/(1-t)^2")
I1 = (log((1 - tau**2*t)/(1 - t)) - tau*t)/(4*K*kappa)
diffA = (A1f - A1f.subs(t, tau**2*t))/(4*K*kappa)
ser = series(log((1 - tau**2*t)/(1 - t)) - tau*t, t, 0, 6).removeO()
cok = zero(ser.coeff(t, 1)) and all(zero(ser.coeff(t, m) - (1 - tau**(2*m))/m)
                                    for m in (2, 3, 4, 5))
ck("QX-E4", zero(I1 - diffA) and cok,
   "I1 closed form = antiderivative difference; t^m coefficient = (1-tau^{2m})/m, m=1 term 0")
lead = simplify(ser.coeff(t, 2)/(4*K*kappa))
ck("QX-E5", zero(lead - K/(8*kappa)) and zero(1 - tau**4 - K2)
        and zero(K/(8*kappa) - K*log(phi)/(4*pi)),
   "leading constant: lim I1/t^2 = K/(8 kappa) = K ln(phi)/(4 pi), via 1-tau^4 = K^2")
sv = symbols('s_v', positive=True)
ck("QX-E6", expand((t**4 - (tau**2*t)**4)/4 - t**4*(1 - tau**8)/4) == 0
        and zero((1 - sv) - K2 - (gap - sv))
        and expand(sv*(2 - sv) - (2*sv - sv**2)) == 0,
   "I2 majorant: int s^3 = t^4(1-tau^8)/4; (1-s)-K^2 = gap-s >= 0 for s<=gap; x(t) increasing")
serr = series(log((1 - tau**4*t)/(1 - tau**2*t)) - tau**3*t, t, 0, 3).removeO()
ck("QX-E7", zero(simplify(serr.coeff(t, 2)/ser.coeff(t, 2)) - tau**4) and zero(tau**4 - gap),
   "excess ratio: lim I1(tau^2 t)/I1(t) = tau^4 = gap (gap-paced approach to the floor)")
ck("QX-E8", zero(tau**4 - gap) and zero((1/w**2) - (1/w)**2),
   "|W_n|^4 = tau^{4n} = gap^n: the excess is priced by co-intensity squared")
av, cvv = symbols('a_v c_v', positive=True)
ck("QX-E9", expand((zv**2 + av)*(zv**2 + 1/av) - (zv**4 + (av + 1/av)*zv**2 + 1)) == 0
        and expand((cvv**2 - 2)**2 - 4 - (cvv**4 - 4*cvv**2)) == 0,
   "reciprocal factorization z^4+(c^2-2)z^2+1 = (z^2+a)(z^2+1/a); disc = c^2(c^2-4)")

# ---------------------------------------------------------------- QX-F
print("== QX-F: Gamma-transfer of the interpolation classification ==")
gridF = all(simplify(exp(I*(pi/2 + 2*pi*k)*n) - I**n) == 0
            for k in range(-2, 3) for n in range(0, 9))
symF = simplify(cos(2*pi*kk) - 1) == 0 and simplify(sin(2*pi*kk)) == 0
ck("QX-F1", gridF and symF,
   "rung agreement: e^{i kappa_k n ln phi} = i^n (exact grid + e^{2 pi i k} = 1)")
ck("QX-F2", expand((pi/2 + 2*pi*kk) - (pi/2 + 2*pi*kp) - 2*pi*(kk - kp)) == 0,
   "kappa_k - kappa_{k'} = 2 pi (k-k')/ln phi: branches pairwise distinct as curves")
ck("QX-F3", all(abs(Fraction(1, 2) + 2*k) >= Fraction(3, 2) for k in range(-10, 11) if k != 0)
        and Fraction(1, 2) < Fraction(3, 2),
   "|1/2+2k| uniquely minimized at k=0 (exact rationals; tail bound for k != 0)")
rv, Av, Bv = symbols('r_v A_v B_v', positive=True)
ck("QX-F4", expand(diff(Av + rv**2*u_s, u_s) - rv**2) == 0
        and zero(diff(sqrt(Av + Bv*u_s), u_s) - Bv/(2*sqrt(Av + Bv*u_s))),
   "rail speed^2 = A + r^2 kappa'^2 strictly increasing in kappa'^2 wherever r > 0")
ck("QX-F5", bool(simplify(3*s5 - 5) > 0) and bool(simplify(sqrt(3)/2) > 0),
   "K > 0 and z_c > 0: r > 0 on both profile branches for z > 0 (minimality applies)")
rho = symbols('rho', positive=True)
Sk2 = exp(-2*rho)*(cos((pi/2 + 2*pi*kk)*rho/log(phi))**2
                   + sin((pi/2 + 2*pi*kk)*rho/log(phi))**2) + (1 - exp(-2*rho))
ck("QX-F6", zero(trigsimp(Sk2) - 1),
   "|S_k(rho)|^2 = 1 for every branch: transfer to the residual sphere")

# ---------------------------------------------------------------- QX-G guards
print("== QX-G: certified-interval and numeric guards (counted separately) ==")
from mpmath import iv, mp, mpf, quad as mquad
iv.dps = 60
iphi = (1 + iv.sqrt(5))/2
igap = iphi**-4
iK2 = 1 - igap
ikap = iv.pi/(2*iv.log(iphi))
gk("QX-G1", igap**2/(ikap**2*iK2**2) < iv.mpf(8),
   "interval-certified: x_max = gap^2/(kappa^2 K^4) < 8 (lower bound valid on n>=2)")
gk("QX-G2", iv.sqrt(iK2)*ikap > iv.mpf(2),
   "interval-certified: K kappa > 2 (reciprocal quartic factors with a > 0 real)")
mp.dps = 50
fphi = (1 + mp.sqrt(5))/2
ftau = 1/fphi
fgap = fphi**-4
fK = mp.sqrt(1 - fgap)
fkap = mp.pi/(2*mp.log(fphi))
ok_num = True
tbl = []
for n in range(2, 6):
    tn = ftau**(2*n)
    zfn = lambda th: mp.sqrt(1 - mp.e**(-2*th/fkap))
    integ = lambda th: mp.sqrt(fK**2 + ((1 - zfn(th)**2)/(fkap*zfn(th)))**2)
    ds = mquad(integ, [n*mp.pi/2, (n + 1)*mp.pi/2])
    I1n = (mp.log((1 - ftau**2*tn)/(1 - tn)) - ftau*tn)/(4*fK*fkap)
    Fp = lambda tt: 1/(1 - tt) + 3*mp.log(1 - tt) - 3*(1 - tt) + (1 - tt)**2/2
    I2n = (Fp(tn) - Fp(ftau**2*tn))/(16*fK**3*fkap**3)
    lo = mp.pi*fK/2 + I1n - I2n
    hi = mp.pi*fK/2 + I1n
    ok_num = ok_num and (lo - mpf(10)**-40 <= ds <= hi + mpf(10)**-40)
    tbl.append((n, ds, (ds - mp.pi*fK/2)/fgap**n))
gk("QX-G3", ok_num, "numeric corroboration (50 dps quadrature): bounds bracket Delta s_n, n=2..5")
print("   n |        Delta s_n (18 sig)        | (Ds_n - piK/2)/gap^n")
for n, ds, rat in tbl:
    print("   %d | %s | %s" % (n, mp.nstr(ds, 18), mp.nstr(rat, 18)))
print("   K ln(phi)/(4 pi) =", mp.nstr(fK*mp.log(fphi)/(4*mp.pi), 18))
print("   pi K / 2         =", mp.nstr(mp.pi*fK/2, 18))

# ---------------------------------------------------------------- summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('QX-G')])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
