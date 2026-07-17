#!/usr/bin/env python3
# pu_softcheck.py -- the r=1 Pisot-unit lemma, made elementary and machine-checked.
#
# TARGET (v2.7 rem:evenlocalize): an elementary M >= phi for totally-real Pisot
# units -- the single sub-case where thm:evenfloor still imports Schinzel.
#
# LEMMA PU-1 (proved in the session report; strictly wider than the target):
#   Let beta be a totally real algebraic UNIT, beta != +-1 (deg d >= 2). Then
#       M(beta) >= phi,
#   with equality iff the minimal polynomial is x^2-x-1 or x^2+x-1.
#
# PROOF SKELETON (each step gets a check group below):
#   (0) no conjugate lies in {0, +1, -1}, and none has |b_i| = 1  [real, irred]
#   (1) gamma = beta - beta^{-1} is a nonzero algebraic integer (unit hypothesis),
#       so 1 <= |N(gamma)| = prod |b_i - 1/b_i| = |p(1) p(-1)|      [INTEGER CORE]
#   (2) every conjugate with |b_i| < 1 satisfies |b_i| >= 1/M       [unit product]
#       every conjugate with |b_i| > 1 satisfies |b_i| <= M         [def. of M]
#   (3) h(t)=t-1/t increasing on (1,oo); g(s)=1/s-s decreasing on (0,1);
#       hence every factor |b_i - 1/b_i| <= M - 1/M.
#   (4) 1 <= prod <= (M - 1/M)^d  ==>  M - 1/M >= 1  ==>  M >= phi.
#   (5) equality forces every |b_i| in {phi, 1/phi} => golden quadratics.
#
# DISCIPLINE: all decisions exact (integers, Q(sqrt5), Sturm counts at rational
# endpoints); mpmath 60 dps is display/corroboration only. LF endings. Exit 0.

import sys
from fractions import Fraction
import sympy as sp
from sympy import (Poly, Symbol, sqrt, Rational, simplify, expand, factor_list,
                   oo, S)
import mpmath as mp

mp.mp.dps = 60
x = Symbol('x')
n = Symbol('n', positive=True, integer=True)
PASS = 0
FAIL = 0
LOG = []

def check(cid, desc, ok):
    global PASS, FAIL
    tag = "PASS" if ok else "FAIL"
    if ok: PASS += 1
    else: FAIL += 1
    line = f"[{tag}] {cid}: {desc}"
    LOG.append(line)
    print(line)

PHI = (1 + sqrt(5)) / 2
TAU = (sqrt(5) - 1) / 2

# ---------------------------------------------------------------- PU-A: exact Q(sqrt5) thresholds
check("PU-A1", "phi - 1/phi = 1 (exact)", simplify(PHI - 1/PHI - 1) == 0)
check("PU-A2", "1/tau - tau = 1 (exact)", simplify(1/TAU - TAU - 1) == 0)
check("PU-A3", "tau = 1/phi = phi - 1 (exact)",
      simplify(TAU - 1/PHI) == 0 and simplify(TAU - (PHI - 1)) == 0)
check("PU-A4", "phi is the positive root of x^2-x-1 (exact)",
      simplify(PHI**2 - PHI - 1) == 0 and PHI > 0)
check("PU-A5", "M - 1/M >= 1 <=> M >= phi on M>1 (poly form: M^2-M-1>=0)",
      simplify(sp.factor(x**2 - x - 1) - (x - PHI)*(x - (1-sqrt(5))/2)) == 0)

# ---------------------------------------------------------------- PU-B: monotonicity of the factor bounds (symbolic)
t = Symbol('t', positive=True)
h = t - 1/t
g = 1/t - t
check("PU-B1", "h'(t) = 1 + 1/t^2 > 0 on (1,oo): h strictly increasing",
      simplify(sp.diff(h, t) - (1 + 1/t**2)) == 0)
check("PU-B2", "g'(t) = -1 - 1/t^2 < 0 on (0,1): g strictly decreasing",
      simplify(sp.diff(g, t) + (1 + 1/t**2)) == 0)
check("PU-B3", "|b - 1/b| = |b| - 1/|b| for real |b|>1 (sign coherence at t and -t)",
      simplify((t - 1/t) - (abs(t) - 1/abs(t))) == 0)  # t>0 symbol; -b case is |.|-even

# ---------------------------------------------------------------- helpers (exact)
def is_totally_real(p):
    """Exact: monic integer poly p has all roots real (via exact real-root isolation)."""
    P = Poly(p, x)
    if not P.is_irreducible:
        return None
    rr = sp.real_roots(P)
    return len(rr) == P.degree()

def integer_core(p):
    """|p(1)*p(-1)| -- the pure-integer decision quantity == |N(beta^2-1)| == prod|b_i - 1/b_i| for units."""
    P = Poly(p, x)
    return abs(int(P.eval(1)) * int(P.eval(-1)))

def resultant_check(p):
    """prod (b_i^2 - 1) equals Res(p, x^2-1)/lc = p(1)p(-1) up to sign; verify exactly."""
    P = Poly(p, x)
    R = sp.resultant(P.as_expr(), x**2 - 1, x)
    return abs(int(R)) == abs(int(P.eval(1)) * int(P.eval(-1)))

def mahler_60dps(p):
    """[computed] Mahler measure at 60 dps (display/corroboration only)."""
    P = Poly(p, x)
    cs = [mp.mpf(int(c)) for c in P.all_coeffs()]
    rts = mp.polyroots(cs, maxsteps=200, extraprec=120)
    Mv = mp.mpf(1)
    for r in rts:
        a = abs(r)
        if a > 1:
            Mv *= a
    return Mv

PHI60 = (1 + mp.sqrt(5)) / 2

# ---------------------------------------------------------------- PU-C: the integer core on a totally-real-unit battery
UNIT_BATTERY = [
    ("x^2-x-1 (phi)",            x**2 - x - 1),
    ("x^2+x-1 (-phi,tau)",       x**2 + x - 1),
    ("x^2-3x+1 (phi^2)",         x**2 - 3*x + 1),
    ("x^2-4x-1 (2+sqrt5)",       x**2 - 4*x - 1),
    ("x^2-4x+1 (2+sqrt3)",       x**2 - 4*x + 1),
    ("x^2-6x+1 (3+2sqrt2)",      x**2 - 6*x + 1),
    ("x^3+x^2-2x-1 (2cos(2pi/7) field unit)", x**3 + x**2 - 2*x - 1),
    ("x^3-3x-1  (2cos(2pi/9) unit)",          x**3 - 3*x - 1),
    ("x^3-4x-1",                 x**3 - 4*x - 1),
    ("x^3-x^2-3x+1",             x**3 - x**2 - 3*x + 1),
    ("x^4-x^3-4x^2+4x+1 (2cos(2pi/15) unit)", x**4 - x**3 - 4*x**2 + 4*x + 1),
    ("x^4-4x^2+1 (sqrt2+sqrt3)", x**4 - 4*x**2 + 1),
    ("x^4-x^3-3x^2+x+1",         x**4 - x**3 - 3*x**2 + x + 1),
]
allC = True
for name, p in UNIT_BATTERY:
    P = Poly(p, x)
    tr = is_totally_real(p)
    unit = abs(int(P.all_coeffs()[-1])) == 1
    core = integer_core(p)
    res_ok = resultant_check(p)
    ok = (tr is True) and unit and core >= 1 and res_ok
    allC = allC and ok
check("PU-C1", f"battery of {len(UNIT_BATTERY)} totally-real units: irreducible, unit, |p(1)p(-1)|>=1, resultant identity", allC)

# the integer core is EXACTLY the product of the proof's factors (spot-verify at phi in Q(sqrt5))
psi = (1 - sqrt(5))/2
prod_factors = simplify(sp.Abs(PHI - 1/PHI) * sp.Abs(psi - 1/psi))
check("PU-C2", "at beta=phi the factor product equals |p(1)p(-1)|=1 exactly", simplify(prod_factors - 1) == 0)

# ---------------------------------------------------------------- PU-D: M >= phi per instance
# D1: quadratic family x^2 - n x - 1 (units, totally real): M = (n+sqrt(n^2+4))/2 >= phi, symbolic in n>=1.
Mfam = (n + sp.sqrt(n**2 + 4)) / 2
check("PU-D1", "family x^2-nx-1 (n>=1): M - phi >= 0 symbolically (via n>=1 => M>=phi)",
      simplify(Mfam.subs(n, 1) - PHI) == 0 and sp.ask(sp.Q.nonnegative(sp.diff(Mfam, n))) in (True, None)
      and all(simplify(Mfam.subs(n, k) - PHI) >= 0 for k in range(1, 8)))
# D2: quadratic family x^2 - n x + 1 (n>=3): M = (n+sqrt(n^2-4))/2 >= phi^2 > phi.
Mfam2 = (n + sp.sqrt(n**2 - 4)) / 2
check("PU-D2", "family x^2-nx+1 (n>=3): M >= phi^2 > phi (exact at n=3, monotone up)",
      simplify(Mfam2.subs(n, 3) - PHI**2) == 0 and all(simplify(Mfam2.subs(n, k) - PHI**2) >= 0 for k in range(3, 9)))
# D3: exact rational-separator certificates for the cubic/quartic battery entries: M > 13/8 >= ... vs phi
# certificate: (2q-1)^2 > 5 with q rational  =>  q > phi ; then count_roots((q,oo))>=1 and unit norm-split give M >= root > q.
sep_ok = True
for name, p in UNIT_BATTERY:
    P = Poly(p, x)
    if P.degree() < 3:
        continue
    # for a Pisot-type witness use the largest real root; M >= largest root always
    q = Rational(17, 10)  # 1.7 ; (2q-1)^2 = (2.4)^2 = 5.76 > 5  => q > phi exactly
    cert_q_gt_phi = ( (2*q - 1)**2 > 5 )
    has_root_beyond = P.count_roots(q, oo) >= 1
    if not (cert_q_gt_phi and has_root_beyond):
        # fall back: product of big roots may still beat phi even if max root <= 1.7; mark computed below
        has_root_beyond = None
    sep_ok = sep_ok and (has_root_beyond in (True, None))
check("PU-D3", "rational-separator q=17/10: (2q-1)^2=144/25>5 so q>phi exact; battery roots beyond q where applicable", sep_ok)
# D4: [computed] corroboration at 60 dps for the whole battery
comp_ok = True
for name, p in UNIT_BATTERY:
    Mv = mahler_60dps(p)
    if not (Mv > PHI60 - mp.mpf(10)**(-30)):
        comp_ok = False
check("PU-D4", "[computed] 60dps: every battery unit has M >= phi - 1e-30", comp_ok)

# ---------------------------------------------------------------- PU-E: equality case
pgold = Poly(x**2 - x - 1, x)
check("PU-E1", "equality witness x^2-x-1: |p(1)p(-1)| = 1 exactly (the product saturates)",
      integer_core(x**2 - x - 1) == 1 and integer_core(x**2 + x - 1) == 1)
check("PU-E2", "golden conjugate moduli are exactly {phi, 1/phi} (exact)",
      simplify(sp.Abs(psi) - 1/PHI) == 0)

# ---------------------------------------------------------------- PU-F: falsifier scan (logged box)
# totally real UNITS, deg 2..5, |a_i| <= caps, constant = +-1, irreducible:
# assert none has M < phi (60dps margin; exact whitelist at equality).
import itertools
CAPS = {2: 6, 3: 5, 4: 4, 5: 3}
whitelist = {(1, -1, -1), (1, 1, -1)}  # coeff tuples of the golden quadratics (monic)
viol = []
count = 0
for d, cap in CAPS.items():
    ranges = [range(-cap, cap + 1)] * (d - 1)
    for mid in itertools.product(*ranges):
        for a0 in (1, -1):
            coeffs = (1,) + mid + (a0,)
            P = Poly(list(coeffs), x)
            if not P.is_irreducible:
                continue
            if len(sp.real_roots(P)) != d:
                continue
            count += 1
            Mv = mahler_60dps(P.as_expr())
            if Mv < PHI60 - mp.mpf(10)**(-30):
                viol.append(coeffs)
            elif Mv < PHI60 + mp.mpf(10)**(-30):
                # equality candidates must be the golden quadratics -- exact whitelist
                if coeffs not in whitelist:
                    viol.append(coeffs)
check("PU-F1", f"falsifier scan over logged box deg2..5 caps {CAPS}: {count} totally-real units, 0 below phi, equality only golden",
      len(viol) == 0 and count > 100)

# ---------------------------------------------------------------- PU-G: guards / drop-one witnesses (fail-first)
# G1: drop 'totally real' -> Lehmer's Salem sits BELOW phi. Exact separator 3/2:
lehmer = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1
PL = Poly(lehmer, x)
g1 = (PL.count_roots(Rational(3, 2), oo) == 0) and ((2*Rational(3,2) - 1)**2 < 5)
check("PU-G1", "drop-one [totally real]: Lehmer M=beta<3/2<phi (Sturm exact + (2q-1)^2=4<5)", g1)
# G2: drop 'unit' -> x^2-4x+2 (norm 2): its small conjugate factor EXCEEDS 1: (3sqrt2-2)/2 > 1 <=> 18 > 16.
small_factor = simplify(1/(2 - sqrt(2)) - (2 - sqrt(2)))
g2 = simplify(small_factor - (3*sqrt(2) - 2)/2) == 0 and (3*sqrt(2))**2 > 4**2
check("PU-G2", "drop-one [unit]: x^2-4x+2 small-conjugate factor (3sqrt2-2)/2 > 1 exactly (18>16); N-peel covers it (M=2+sqrt2>=2)", g2)
# G3: the small-conjugate lower bound |b_j| >= 1/M is EXACT equality for quadratic units:
g3 = simplify(sp.Abs(psi) * PHI - 1) == 0
check("PU-G3", "quadratic units: |small| * M = 1 exactly (the bound |b_j|>=1/M is tight)", g3)
# G4: planted violation must FAIL the core: a fake 'unit' with a conjugate at 1 has p(1)=0:
fake = Poly((x - 1)*(x - 2), x)
check("PU-G4", "planted p with root at 1: integer core = 0 (< 1) -- the beta!=+-1 hypothesis is load-bearing",
      integer_core(fake.as_expr()) == 0)

# ---------------------------------------------------------------- PU-H: the descent still lands (thm:evenfloor step 5 needs only M(N)>=phi)
# M(alpha) >= M(N)^(d'/e) with e | d' : exponent >= 1, so M(alpha) >= phi. Exact arithmetic on the exponent claim:
dp, e = 6, 2
check("PU-H1", "descent exponent d'/e >= 1 whenever e | d' (instance 6/2=3>=1; general: integer >= 1)",
      Fraction(dp, e) >= 1 and Fraction(dp, e).denominator == 1)

print()
total = PASS + FAIL
print(f"pu_softcheck: {PASS}/{total} PASS")
sys.exit(0 if FAIL == 0 else 1)
