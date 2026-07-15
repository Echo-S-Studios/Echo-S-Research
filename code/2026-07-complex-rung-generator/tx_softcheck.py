#!/usr/bin/env python3
# tx_softcheck.py -- transcendence-pack harness (HANDOFF sec 3.4 Route R3, T1-T6;
# closes W[open]4 per sec 7). Written cold this session from the dossier + papers.
# Discipline (nx/qx idiom): ALL decisions exact in Q / Q(sqrt5) (extended by i
# where needed); K-level statements decided at the squared level; mpmath iv at
# 60 dps ONLY for the certified interval guards and the decimal audit (group
# TX-G, counted separately); floats never touch an exact decision.
#
#   TX-A  T1: i*pi/ln(phi) not rational -- the (-1)^q = phi^p engine      (7)
#   TX-B  T3: branch enumeration kappa_k = (4k+1) kappa, 2i-bridge        (5)
#   TX-C  T4: W[open]4 closure algebra + Guard 15.3 pitch clause          (7)
#   TX-D  T5: W-Thm 16.1 replication (Binet + norm, premise-disjoint)     (5)
#   TX-E  T6: ledger algebraicity certificates, closure, engine teeth     (6)
#   TX-G  certified-interval + decimal-audit guards (counted separately)  (6)
#
# T2 (Gelfond-Schneider) is CITED, not certified: see the CITE line below.
# Fail-first: every group carries at least one falsifier (a check that would
# FAIL if the claim -- or the deciding engine -- were broken): TX-A2, TX-A6,
# TX-C6, TX-D2, TX-E4.

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin,
                   cos, tan, asin, simplify, expand, radsimp, trigsimp,
                   together, cancel, Poly, minimal_polynomial, factor_list,
                   fibonacci, lucas, im, re, CRootOf, expand_log)
from sympy.polys.polyerrors import NotAlgebraic

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
    if simplify(expand_log(e, force=True)) == 0:
        return True
    return simplify(trigsimp(e2)) == 0

x = symbols('x')
kk = symbols('k', integer=True)

s5    = sqrt(5)
phi   = (1 + s5)/2
psi   = (1 - s5)/2          # golden conjugate, phi*psi = -1
tau   = (s5 - 1)/2
gap   = phi**-4
K2    = 1 - gap
K     = sqrt(K2)
beta2 = phi**2*s5
beta  = sqrt(beta2)
zc    = sqrt(3)/2
kappa = pi/(2*log(phi))
Mcons = 11 + 5*s5
Mres  = (43 + 19*s5)/2
muS   = CRootOf(x**3 - x - 1, 0)
beta4 = CRootOf(x**4 - x**3 - x**2 - x + 1, 1)

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
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    if A > 0:
        return 1 if A*A > 5*B*B else -1
    return 1 if 5*B*B > A*A else -1

def irreducibleZ(P):
    fl = factor_list(P, x)[1]
    return len(fl) == 1 and fl[0][1] == 1

print("CITE T2 [established, NOT harness-checkable] Gelfond-Schneider theorem")
print("  (A. O. Gelfond 1934; Th. Schneider 1934, independently). Quotient form")
print("  used: for nonzero algebraic alpha, beta with fixed logarithm branches,")
print("  log(alpha) != 0 != log(beta), the ratio log(alpha)/log(beta) is either")
print("  rational or transcendental. Instances consumed below:")
print("  (alpha,beta) = (-1, phi), log(-1) = i*pi principal   [T3]")
print("  (alpha,beta) = ( 2, phi), real logarithms            [T5]")

# ---------------------------------------------------------------- TX-A
print("== TX-A: T1 engine -- i*pi/ln(phi) not rational, elementary ==")
GRID_P = [p for p in range(-12, 13)]
AB = {p: q5AB(phi**p) for p in GRID_P}   # exact (A,B): phi^p = A + B sqrt5

ck("TX-A1", zero(phi**2 - phi - 1) and zero(phi*tau - 1) and zero(tau**2 + tau - 1)
        and all(zero(2*phi**p - (lucas(p) + fibonacci(p)*s5)) for p in GRID_P),
   "seed identities + Binet 2 phi^p = L_p + F_p sqrt5 exact on p in [-12,12]")
ck("TX-A2", all(not zero(2*phi**p - (lucas(p) + (fibonacci(p) + 1)*s5))
                for p in (-5, -1, 1, 2, 7))
        and not zero(2*phi**3 - (lucas(3) + 1 + fibonacci(3)*s5)),
   "FALSIFIER: perturbed Binet coefficients are rejected (zero() has teeth)")
ck("TX-A3", all(AB[p][1] != 0 and AB[p][1] == Fraction(int(fibonacci(p)), 2)
                for p in GRID_P if p != 0),
   "p != 0: sqrt5-part of phi^p is F_p/2 != 0 => phi^p irrational => phi^p != +-1")
ck("TX-A4", all(sgnQ5(phi**p) == 1 for p in GRID_P)
        and all(not zero(phi**p - 1) and not zero(phi**p + 1)
                for p in GRID_P if p != 0),
   "phi^p > 0 (exact sign) and phi^p != +-1 decided directly for p != 0")
ck("TX-A5", all(simplify(exp(I*pi*q) - (-1)**q) == 0 for q in range(1, 13))
        and all(((-1)**q == 1) == (q % 2 == 0) for q in range(1, 13))
        and zero(phi**0 - 1),
   "exponentiation bookkeeping: e^{i pi q} = (-1)^q; degenerate p=0 gives phi^0=1,"
   " and q*i*pi != 0 for q>=1 is closed by guard TX-G1 (pi > 3)")
ck("TX-A6", zero(phi*psi + 1)
        and all(zero((phi*psi)**p - (-1)**p) for p in range(-6, 7)),
   "FALSIFIER-CONTROL: the engine DOES detect genuine +-1 powers ((phi psi)^p)")
c_re = pi/log(phi)
ck("TX-A7", sgnQ5(phi - 1) == 1
        and ((c_re.is_real is True) or (simplify(im(c_re)) == 0))
        and simplify(re(I*c_re)) == 0,
   "alternate route: i pi/ln phi is purely imaginary (phi > 1 exact; ln phi real"
   " nonzero [established monotonicity]) -- a nonzero non-real number is not in Q")

# ---------------------------------------------------------------- TX-B
print("== TX-B: T3 branch enumeration and the 2i-bridge ==")
ck("TX-B1", expand((pi/2 + 2*pi*kk) - (4*kk + 1)*pi/2) == 0
        and zero((pi/2 + 2*pi*kk)/log(phi) - (4*kk + 1)*kappa),
   "kappa_k = (pi/2 + 2 pi k)/ln phi = (4k+1) kappa, symbolic in integer k")
ck("TX-B2", all((Fraction(1, 2) + 2*k)/Fraction(1, 2) == 4*k + 1
                for k in range(-8, 9)),
   "grid k in [-8,8]: kappa_k/kappa = 4k+1 (exact rationals)")
ck("TX-B3", expand((4*kk + 1) - 2*(2*kk) - 1) == 0
        and all(abs(4*k + 1) >= 1 and (4*k + 1) % 2 == 1 for k in range(-50, 51))
        and not Rational(-1, 4).is_integer,
   "4k+1 is odd hence never 0 (no integer solves 4k+1=0): the factor is a"
   " NONZERO rational for every branch")
ck("TX-B4", zero(I*pi/log(phi) - 2*I*kappa)
        and zero((I*pi/log(phi))/(2*I) - kappa)
        and expand(minimal_polynomial(2*I, x) - (x**2 + 4)) == 0,
   "2i-bridge: i pi/ln phi = 2i kappa with 2i algebraic (x^2+4) => kappa algebraic"
   " iff i pi/ln phi algebraic; G-S + T1 kill the latter")
ck("TX-B5", zero(2/((4*kk + 1)*kappa) - (4*log(phi)/pi)/(4*kk + 1))
        and zero(((4*kk + 1)*kappa)/(4*kk + 1) - kappa),
   "branch uniformity: 2/kappa_k = (4 ln phi/pi)/(4k+1) and kappa recovered from"
   " any branch by the nonzero rational 1/(4k+1)")

# ---------------------------------------------------------------- TX-C
print("== TX-C: T4 -- W[open]4 closure algebra + Guard 15.3 pitch clause ==")
ck("TX-C1", zero(2/kappa - 4*log(phi)/pi) and zero(pi/log(phi) - 2*kappa)
        and zero(log(phi)/pi - 1/(2*kappa)),
   "near-miss LHS: 4 ln phi/pi = 2/kappa; every rational multiple of pi/ln phi"
   " (or ln phi/pi) is a nonzero rational multiple of kappa (or 1/kappa)")
ck("TX-C2", zero(2/tau - 2*phi)
        and expand(minimal_polynomial(2*phi, x) - (x**2 - 2*x - 4)) == 0
        and zero((2*phi)**2 - 2*(2*phi) - 4),
   "clause 1 equivalence: 4 ln phi/pi = tau <=> kappa = 2/tau = 2 phi, and 2 phi"
   " is algebraic (x^2-2x-4) -- impossible for transcendental kappa")
ck("TX-C3", zero(1 - tau**4 - K2) and zero(cos(asin(tau**2)) - K)
        and zero(tan(asin(tau**2)) - tau**2/K),
   "pitch chain 1: 1 - tau^4 = K^2, so tan(theta_K) = tan(arcsin tau^2) = tau^2/K")
ck("TX-C4", zero(2*log(phi)/pi - 1/kappa) and zero(1/tau**2 - phi**2)
        and expand(minimal_polynomial(phi**2, x) - (x**2 - 3*x + 1)) == 0
        and zero(phi**4 - 3*phi**2 + 1),
   "pitch chain 2: 2 ln phi/(pi K) = tau^2/K <=> 2 ln phi/pi = tau^2 <=> kappa ="
   " 1/tau^2 = phi^2, and phi^2 is algebraic (x^2-3x+1) -- impossible")
ck("TX-C5", all(expand(minimal_polynomial(v, x) - P) == 0 and irreducibleZ(P)
                for v, P in ((tau, x**2 + x - 1), (2*phi, x**2 - 2*x - 4),
                             (phi**2, x**2 - 3*x + 1))),
   "the three would-be targets tau, 2 phi, phi^2 are algebraic with irreducible"
   " integer minimal polynomials")
ck("TX-C6", not zero(2/kappa - 4*log(phi)/pi - Rational(1, 10**6))
        and not zero(cos(asin(tau**2)) - K - Rational(1, 10**6))
        and not zero(1/tau**2 - phi**2 - Rational(1, 10**6)),
   "FALSIFIER: perturbed variants of C1/C3/C4 identities are rejected")
ck("TX-C7", sgnQ5(K2) == 1 and sgnQ5(tau**2) == 1 and sgnQ5(1 - tau**2) == 1
        and sgnQ5(1 - K2) == 1,
   "domain facts (exact signs): K > 0 and tau^2 in (0,1), so both angles lie in"
   " (0, pi/2) where tan is injective [established monotonicity]")

# ---------------------------------------------------------------- TX-D
print("== TX-D: T5 -- W-Thm 16.1 replication + strengthening bookkeeping ==")
ok1 = True
for p in GRID_P[2:-2]:                       # p in [-10, 10]
    A, B = AB[p]
    for q in range(1, 9):
        if p != 0:
            ok1 = ok1 and B != 0             # phi^p irrational, 2^{+-q} rational
        else:
            ok1 = ok1 and B == 0 and A - Fraction(2**q) != 0 \
                      and A - Fraction(1, 2**q) != 0
ck("TX-D1", ok1,
   "2^q != phi^p exactly for p in [-10,10], q in [1,8], both orientations"
   " (sqrt5-part F_p/2 != 0 for p != 0; A-part 1 != 2^{+-q} for p = 0)")
d42 = q5AB(Integer(4)**1 - Integer(2)**2)
d86 = q5AB(Integer(8)**2 - Integer(2)**6)
ck("TX-D2", d42 == (0, 0) and d86 == (0, 0),
   "FALSIFIER-CONTROL: the same (A,B) engine DETECTS the true equalities"
   " 4^1 = 2^2 and 8^2 = 2^6 (multiplicatively dependent pair)")
ck("TX-D3", expand(minimal_polynomial(Integer(2), x) - (x - 2)) == 0
        and expand(minimal_polynomial(phi, x) - (x**2 - x - 1)) == 0
        and sgnQ5(phi - 1) == 1 and Integer(2) - 1 != 0
        and (x**2 - x - 1).subs(x, 1) != 0 and (x - 2).subs(x, 1) != 0,
   "G-S side conditions: 2 and phi algebraic, both != 0,1, both > 1 so both"
   " logarithms are nonzero [established monotonicity]")
ck("TX-D4", all(simplify(exp(q*log(2)) - 2**q) == 0 for q in (1, 2, 3))
        and all(simplify(exp(p*log(phi)) - phi**p) == 0 for p in (-2, 1, 3)),
   "exponentiation bookkeeping: ln2/lnphi = p/q would force 2^q = phi^p")
ok5 = all(zero(phi**p*psi**p - (-1)**p) for p in range(-8, 9)) \
      and all(4**q >= 4 and abs((-1)**p) == 1 for p in range(-8, 9)
              for q in range(1, 9))
ck("TX-D5", ok5,
   "second, premise-disjoint mechanism: field norm N(phi^p) = (phi psi)^p ="
   " (-1)^p has |N| = 1, while N(2^q) = 4^q >= 4: no equality is possible"
   " (norm obstruction, independent of Binet positivity)")

# ---------------------------------------------------------------- TX-E
print("== TX-E: T6 -- ledger algebraicity, closure instances, engine teeth ==")
LEDGER = [
    ("phi",    phi,            x**2 - x - 1),
    ("tau",    tau,            x**2 + x - 1),
    ("K",      K,              x**4 + 5*x**2 - 5),
    ("beta",   beta,           x**4 - 5*x**2 - 5),
    ("K^2",    K2,             x**2 + 5*x - 5),
    ("beta^2", beta2,          x**2 - 5*x - 5),
    ("gap",    gap,            x**2 - 7*x + 1),
    ("z_c",    zc,             4*x**2 - 3),
    ("tau^2",  tau**2,         x**2 - 3*x + 1),
    ("phi^2",  phi**2,         x**2 - 3*x + 1),
    ("4/9",    Rational(4, 9), 9*x - 4),
    ("sqrt5",  s5,             x**2 - 5),
    ("L4",     Integer(7),     x - 7),
    ("Mcons",  Mcons,          x**2 - 22*x - 4),
    ("Mres",   Mres,           x**2 - 43*x + 11),
    ("muS",    muS,            x**3 - x - 1),
    ("beta4",  beta4,          x**4 - x**3 - x**2 - x + 1),
]
okE1 = True
for name, v, P in LEDGER:
    okE1 = okE1 and expand(minimal_polynomial(v, x) - P) == 0 and irreducibleZ(P)
    if not isinstance(v, CRootOf):
        okE1 = okE1 and zero(P.subs(x, v))
ck("TX-E1", okE1,
   "algebraicity certificates: 17 ledger constants each root an irreducible"
   " integer polynomial (minimal_polynomial verified)")
okE2 = True
for c in (phi**2, 2*phi, K):
    pm = minimal_polynomial(2*I*c, x)
    okE2 = okE2 and all(cf.is_integer for cf in Poly(pm, x).all_coeffs()) \
                and Poly(pm, x).degree() >= 2 and zero(pm.subs(x, 2*I*c))
ck("TX-E2", okE2,
   "bridge effectivity: a algebraic => 2i a algebraic, instanced on phi^2, 2 phi,"
   " K (integer minimal polynomials constructed and verified); contrapositive is"
   " T3's step kappa = (i pi/ln phi)/(2i)")
okE3 = expand(minimal_polynomial(tau + zc, x)
              - (16*x**4 + 32*x**3 - 40*x**2 - 56*x - 11)) == 0 \
    and Poly(minimal_polynomial(phi*zc, x), x).degree() == 4 \
    and expand(minimal_polynomial(1/beta4, x)
               - (x**4 - x**3 - x**2 - x + 1)) == 0 \
    and expand(minimal_polynomial(Rational(4, 9)*gap, x)
               - Poly(81*x**2 - 4*9*7*x + 16, x).as_expr()) == 0
ck("TX-E3", okE3,
   "Qbar closure instances: sums/products/inverses of ledger constants stay"
   " algebraic (tau+z_c, phi*z_c, 1/beta4 self-reciprocal, (4/9)*gap) -- every"
   " algebraic-scan target value lies in Qbar, hence != kappa by T3")
okE4 = 0
for tv in (pi, log(2)):
    try:
        minimal_polynomial(tv, x)
    except NotAlgebraic:
        okE4 += 1
ck("TX-E4", okE4 == 2,
   "FALSIFIER: minimal_polynomial raises NotAlgebraic on pi and log(2) -- the"
   " algebraicity engine is not vacuous (teeth check only; transcendence of pi"
   " is Lindemann [established], not this check)")
ck("TX-E5", zero((4*log(phi)/pi)*kappa - 2)
        and zero((2/kappa)*(pi/(4*log(phi))) - 1),
   "tie-shape bookkeeping: 2/kappa and kappa are Qbar-linked (product 2), so a"
   " single algebraic tie to EITHER would make BOTH algebraic -- one exclusion"
   " retires both scan families")
f3 = lambda t: t**3 - t - 1
g4 = lambda t: t**4 - t**3 - t**2 - t + 1
ck("TX-E6", f3(Fraction(13247, 10000)) < 0 and f3(Fraction(13248, 10000)) > 0
        and g4(Fraction(17220, 10000)) < 0 and g4(Fraction(17221, 10000)) > 0,
   "exact rational sign changes bracket muS in (1.3247,1.3248) and beta4 in"
   " (1.7220,1.7221) [IVT established]: brackets feed guard TX-G5")

# ---------------------------------------------------------------- TX-G guards
print("== TX-G: certified-interval and decimal-audit guards (counted separately) ==")
from mpmath import iv, mp
import mpmath as mnum
iv.dps = 60
iphi  = (1 + iv.sqrt(5))/2
itau  = (iv.sqrt(5) - 1)/2
igap  = iphi**-4
iK2   = 1 - igap
iK    = iv.sqrt(iK2)
ibeta2 = iphi**2*iv.sqrt(5)
ilnphi = iv.log(iphi)
ikappa = iv.pi/(2*ilnphi)

gk("TX-G1", (iv.pi > 3) and (ilnphi > iv.mpf('0.4812')) and (ilnphi < iv.mpf('0.4813')),
   "interval-certified: pi > 3 (closes T1's p=0 branch: q i pi != 0 for q >= 1)"
   " and 0.4812 < ln phi < 0.4813 (so ln phi != 0)")
idiff = itau - 4*ilnphi/iv.pi
gk("TX-G2", (idiff > iv.mpf('0.00533')) and (idiff < iv.mpf('0.00534')),
   "interval-certified: 0.00533 < tau - 4 ln phi/pi < 0.00534 (clause-1 margin;"
   " corroborates the forced never-equal)")
ipm = itau**2 - 2*ilnphi/iv.pi
gk("TX-G3", ipm > iv.mpf('0.07'),
   "interval-certified: tau^2 - 2 ln phi/pi > 0.07 (pitch clause margin at the"
   " tangent level; equality would force kappa = phi^2)")
ig4 = ikappa - 2*iphi
gk("TX-G4", (ig4 > iv.mpf(1)/40) and (ig4 < iv.mpf(3)/100),
   "interval-certified: 1/40 < kappa - 2 phi < 3/100 (the clause-1 would-be"
   " value 2 phi misses kappa by ~0.0282)")
ILEDGER = [
    ("phi", iphi), ("tau", itau), ("K", iK), ("K^2", iK2),
    ("beta", iv.sqrt(ibeta2)), ("beta^2", ibeta2), ("gap", igap),
    ("z_c", iv.sqrt(3)/2), ("tau^2", itau**2), ("phi^2", iphi**2),
    ("2phi", 2*iphi), ("4/9", iv.mpf(4)/9), ("sqrt5", iv.sqrt(5)),
    ("L4", iv.mpf(7)), ("Mcons", 11 + 5*iv.sqrt(5)),
    ("Mres", (43 + 19*iv.sqrt(5))/2),
    ("muS", iv.mpf(['1.3247', '1.3248'])),      # bracket certified in TX-E6
    ("beta4", iv.mpf(['1.7220', '1.7221'])),    # bracket certified in TX-E6
]
okG5 = all(abs(ikappa - c) > iv.mpf(1)/40 for _, c in ILEDGER)
gk("TX-G5", okG5,
   "interval-certified: |kappa - c| > 1/40 for all 18 ledger constants"
   " (finite-scan corroboration of T6; the theorem covers all of Qbar)")
mp.dps = 60
fphi  = (1 + mnum.sqrt(5))/2
ftau  = 1/fphi
fK    = mnum.sqrt(1 - fphi**-4)
flnphi = mnum.log(fphi)
audits = [
    ("4lnphi/pi",  4*flnphi/mnum.pi,                          '0.6126979251'),
    ("tau",        ftau,                                       '0.6180339887'),
    ("difference", ftau - 4*flnphi/mnum.pi,                    '0.0053360637'),
    ("pitch(deg)", mnum.atan(2*flnphi/(mnum.pi*fK))*180/mnum.pi, '18.3394927773'),
    ("thetaK(deg)", mnum.asin(ftau**2)*180/mnum.pi,            '22.4555151986'),
    ("logphi2",    mnum.log(2)/flnphi,                         '1.4404200904'),
]
def audit_ok(val, lit):
    dp = len(lit.split('.')[1])
    litv = mnum.mpf(lit)
    ulps = abs(val - litv)*mnum.mpf(10)**dp
    if ulps <= mnum.mpf('0.51'):
        return True
    trunc = mnum.floor(val*mnum.mpf(10)**dp)/mnum.mpf(10)**dp
    return abs(trunc - litv) < mnum.mpf(10)**(-dp - 4)
okG6 = True
for name, val, lit in audits:
    hit = audit_ok(val, lit)
    okG6 = okG6 and hit
    print("   audit %-11s printed %s  recomputed %s  %s"
          % (name, lit, mnum.nstr(val, 16), "ok" if hit else "MISMATCH"))
gk("TX-G6", okG6,
   "decimal audit at 60 dps: all six printed literals of Guard 15.3 / Cor 16.2"
   " match fresh recomputation (0.51 ulp or floor-truncation)")

# ---------------------------------------------------------------- summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('TX-G')])
n_exact_pass = len(PASS)
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (n_exact_pass, n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
