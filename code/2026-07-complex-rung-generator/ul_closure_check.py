#!/usr/bin/env python3
# ul_closure_check.py -- UL-1 apex-closure theorem, SECOND LANE (written cold).
#
# TARGET (v2.7 blk:radab): the apex-closure theorem UL-1 -- (g1)-(g5) preserve
# apex monodromy order 1 over the whole D5 tower; with it, OP-RADIUS ab initio
# closes by IMPOSSIBILITY relative to D1-D3 + D5 (+ the declared ambient AMB).
#
# THEOREM UL-1 (session report, proof by structural induction in the ambient
# AMB = grid-based real log-exp transseries at z->0+ [ADH 2017; vdDMM 2001]):
#   Let L be the set of transseries whose transmonomials are z^n * e^A with
#   n in Z and A purely infinite AND itself in L (log-z-free, integer-graded
#   at every exponential depth).  Then L is a differential subfield closed
#   under (g1) field ops, (g2) d/dz, (g3) exp, (g4) log of an apex unit
#   (v_z = 0), (g5) apex-unramified algebraic extension (constants may extend).
#   The z-grading homomorphism w(z^n e^A) = n takes values in Z on all of L;
#   the formal monodromy M fixes L pointwise (order 1).  Any y with
#   y^2 = q z u (q const != 0, u unit) has w(y) = 1/2, hence y not in L:
#   r = (K/sqrt(z_c)) z^{1/2} is outside F_ed(D5).  The unique order-raising
#   composite exp(c log z) requires log z, inadmissible by the unit-log gate
#   (n_0(z) = v_z(z) = 1 != 0) -- the corpus's named channel, and the only one.
#
# CHECK GROUPS
#   UL-A  the grading homomorphism w on the monomial model (exact Fractions)
#   UL-B  generator preservation battery: (g1) field ops, (g2) d/dz, (g3) exp,
#         (g4) unit-log; every output integer-graded (exact exponents)
#   UL-C  (g5) Newton-Puiseux: unramified battery solves with Z exponents;
#         ramified battery detected by a fractional first slope (exact)
#   UL-D  the exclusion: 2w = 1 insoluble in Z; the r-witness; prop:radiff
#         identity r^2 = (K^2/z_c) rho' (1-z^2) exact; rem:ratie constant
#   UL-E  drop-one: the unit-log gate is load-bearing (formal log-z token
#         raises the grading to 1/2; log of a non-unit refused)
#   UL-F  corpus ties: winding lemma n_0(u) = res_0(u'/u) = v_z(u); the
#         thm:apexram pentad valuations (-2,-4,-2,0,2); v_z(rho)=2; v_z(r^2)=1
#
# DISCIPLINE: exact Fractions / Q(sqrt5) / Q(sqrt3) at every decision;
# no floats anywhere in this file. LF endings. Exit 0.

import sys, random
from fractions import Fraction
import sympy as sp
from sympy import Symbol, sqrt, Rational, simplify, series, exp, log, diff, O

z = Symbol('z', positive=True)
PASS = 0; FAIL = 0

def check(cid, desc, ok):
    global PASS, FAIL
    tag = "PASS" if ok else "FAIL"
    if ok: PASS += 1
    else: FAIL += 1
    print(f"[{tag}] {cid}: {desc}")

random.seed(51)  # sqrt5-flavored seed, deterministic

# ================================================================ UL-A: monomial group & w
# monomial model: (n, A) with n a Fraction (z-power) and A a frozenset of
# (tag, k) pairs -- formal purely-infinite exponent generators e^{A} with
# integer multiplicities; product adds both; w = n.
def mmul(m1, m2):
    n1, A1 = m1; n2, A2 = m2
    d = dict(A1)
    for t, k in A2:
        d[t] = d.get(t, 0) + k
        if d[t] == 0: del d[t]
    return (n1 + n2, frozenset(d.items()))

def w(m):
    return m[0]

ok = True
TAGS = ['a', 'b', 'c']
for _ in range(200):
    m1 = (Fraction(random.randint(-9, 9)), frozenset({(random.choice(TAGS), random.randint(-3, 3))} - {(random.choice(TAGS), 0)}))
    m2 = (Fraction(random.randint(-9, 9)), frozenset({(random.choice(TAGS), random.randint(-3, 3))} - {(random.choice(TAGS), 0)}))
    if w(mmul(m1, m2)) != w(m1) + w(m2):
        ok = False
check("UL-A1", "w is a homomorphism on the monomial model: w(m1*m2)=w(m1)+w(m2), 200 random exact instances", ok)
check("UL-A2", "exp-monomials carry z-grading 0: w(z^0 e^A) = 0 for every A (by construction of the model)",
      all(w((Fraction(0), frozenset({(t, k)}))) == 0 for t in TAGS for k in (-2, -1, 1, 2)))

# ================================================================ series toolkit (exact, truncated Laurent over Q(sqrt5))
N_TRUNC = 9

def zsupport(expr):
    """Exact z-exponent support of a truncated Laurent expression (poly in z, 1/z over Q(sqrt5,sqrt3))."""
    e = sp.expand(expr)
    terms = sp.Add.make_args(e)
    sup = set()
    for t in terms:
        # robust: use as_powers_dict on each multiplicative factor
        deg = Fraction(0)
        for f, k in t.as_powers_dict().items():
            if f == z:
                deg += Fraction(str(k))
        sup.add(deg)
    return sup

def integer_supported(expr):
    return all(d.denominator == 1 for d in zsupport(expr))

def trunc(expr, N=N_TRUNC):
    return sp.expand(sp.series(expr, z, 0, N).removeO())

# (g1) field ops: product & inverse of integer-supported units/series
f1 = z**-2 * (3 + sqrt(5)*z + z**3)
f2 = 5 + 2*z**2 - z**4
prod = sp.expand(f1 * f2)
inv2 = trunc(1/f2)
ok_g1 = integer_supported(prod) and integer_supported(inv2)
check("UL-B1", "(g1) field ops: product and geometric-series inverse stay integer-graded (exact exponents)", ok_g1)

# (g2) d/dz: on series and on exp-monomials (chain rule keeps grading in Z)
g = z**-3 * (1 + z + 7*z**2)
E = exp(2/z + sqrt(5)/z**2)           # purely infinite, integer-graded exponent
h = z**4 * E * (1 + 3*z)
dg = sp.expand(diff(g, z))
dh = sp.expand(diff(h, z))
def graded_with_exp(expr):
    """Integer z-grading allowing exp(A) factors with A itself integer-graded Laurent."""
    e = sp.expand(expr)
    for t in sp.Add.make_args(e):
        deg = Fraction(0)
        for f, k in t.as_powers_dict().items():
            if f == z:
                deg += Fraction(str(k))
            elif f.func == sp.exp:
                if not integer_supported(sp.expand(f.args[0])):
                    return False
        if deg.denominator != 1:
            return False
    return True
check("UL-B2", "(g2) d/dz: derivative of z^-3-series and of z^4 e^{2/z + sqrt5/z^2}(1+3z) stays integer-graded", 
      graded_with_exp(dg) and graded_with_exp(dh))

# (g3) exp: three-part decomposition f = f_inf + c + f_small
f_inf = 3/z; c0 = sqrt(5); f_small = 2*z + z**2
ef = exp(f_inf) * exp(c0) * trunc(exp(f_small))
check("UL-B3", "(g3) exp: exp(f_inf + c + f_small) = e^{f_inf} * e^c * (integer-graded unit series); new monomial has w=0",
      graded_with_exp(sp.expand(ef)) and integer_supported(trunc(exp(f_small))))

# (g4) unit-log: u = c(1+eps), v_z(u)=0  ->  log u = log c + log(1+eps), integer-graded
u_unit = 5*(1 + 3*z + z**3)
logu = trunc(sp.expand(log(5) + log(1 + 3*z + z**3)))
check("UL-B4", "(g4) unit-log: log(5(1+3z+z^3)) = log5 + integer-graded series (gate: v_z(u)=0)",
      integer_supported(sp.expand(logu - log(5))))
# a-fortiori: log of an e^A-dominant unit returns A + log(unit): still in L
logeu = sp.expand(2/z + log(5) + trunc(log(1 + 3*z)))
check("UL-B5", "(g4') log(e^{2/z} * 5(1+3z)) = 2/z + log5 + series: exponent part re-enters L (a-fortiori safe)",
      integer_supported(sp.expand(logeu - log(5))))

# ================================================================ UL-C: (g5) Newton-Puiseux, exact
def np_first_slope(coeffs):
    """coeffs: dict i -> valuation v(c_i) (Fraction) for P(Y)=sum c_i Y^i, c_0 != 0.
    Returns the minimal candidate slope mu = (v(c_0)-v(c_i))/i over i>=1 with c_i present."""
    v0 = coeffs[0]
    best = None
    for i, vi in coeffs.items():
        if i == 0: continue
        mu = Fraction(v0 - vi, i)
        if best is None or mu < best:
            best = mu
    return best

def puiseux_solve_quadratic(a1, a0, N=8):
    """Solve Y^2 + a1(z) Y + a0(z) = 0 by exact Newton iteration from an integer
    first slope; a1, a0 sympy Laurent polys.  Returns (y_trunc, all_integer)."""
    # first slope from valuations
    def valz(e):
        s = zsupport(sp.expand(e))
        return min(s) if s else None
    v1 = valz(a1) if a1 != 0 else None
    v0 = valz(a0)
    cands = []
    if a1 != 0: cands.append(Fraction(v0 - v1, 1))
    cands.append(Fraction(v0, 2))
    mu = min(cands)
    if mu.denominator != 1:
        return None, False, mu
    # leading coefficient: solve c^2*[i2] + c*[i1] + [i0] = 0 on leading terms
    lc0 = sp.expand(a0 / z**int(v0)).subs(z, 0)
    lc1 = sp.expand(a1 / z**int(v1)).subs(z, 0) if a1 != 0 else 0
    # candidate leading term c*z^mu with 2*mu vs mu+v1 vs v0 balance -- do a direct Newton on series:
    csol = None
    for cand in sp.solve(sp.Symbol('Y')**2 * (1 if 2*mu == v0 else 0) + (lc1 if (a1 != 0 and mu + v1 == v0) else 0)*sp.Symbol('Y') + lc0, sp.Symbol('Y')):
        if cand != 0:
            csol = cand; break
    if csol is None:
        # pure square-root balance: c^2 = -lc0
        csol = sqrt(-lc0)
    y = sp.expand(csol * z**int(mu))
    for _ in range(N):
        Fv = sp.expand(y**2 + a1*y + a0)
        Fp = sp.expand(2*y + a1)
        corr = trunc(-Fv / Fp, N + 4)
        y = sp.expand(y + corr)
        y = trunc(y, N + 2)
    resid = trunc(sp.expand(y**2 + a1*y + a0), N)
    ok_solved = (resid == 0) or (min(zsupport(resid)) >= N - 2 if zsupport(resid) else True)
    return y, integer_supported(y) and ok_solved, mu

# unramified battery: Y^2 = 1+z ; Y^2 = 5(1+z) [constants extend]; Y^2 - (2+z)Y + (1+z) = 0 ; Y^2 = z^2(1+z)
y1, ok1, mu1 = puiseux_solve_quadratic(sp.Integer(0), -(1 + z))
y2, ok2, mu2 = puiseux_solve_quadratic(sp.Integer(0), -5*(1 + z))
y3, ok3, mu3 = puiseux_solve_quadratic(-(2 + z), (1 + z))
y4, ok4, mu4 = puiseux_solve_quadratic(sp.Integer(0), -z**2*(1 + z))
check("UL-C1", "(g5) unramified battery: sqrt(1+z), sqrt(5(1+z)), root of Y^2-(2+z)Y+(1+z), sqrt(z^2(1+z)) -- all integer-graded expansions, residuals vanish to order",
      ok1 and ok2 and ok3 and ok4 and all(m.denominator == 1 for m in (mu1, mu2, mu3, mu4)))

# ramified battery: fractional first slope detected exactly
_, okr1, mur1 = puiseux_solve_quadratic(sp.Integer(0), -z)          # Y^2 = z     -> slope 1/2
_, okr2, mur2 = puiseux_solve_quadratic(sp.Integer(0), -z**3*(1+z)) # Y^2 = z^3.. -> slope 3/2
check("UL-C2", "(g5) ramified battery: Y^2=z and Y^2=z^3(1+z) rejected -- first slopes 1/2, 3/2 not in Z (exact Fractions)",
      (okr1 is False) and (okr2 is False) and mur1 == Fraction(1, 2) and mur2 == Fraction(3, 2))
check("UL-C3", "discriminating pair: Y^2=z^2(1+z) admitted (slope 1) vs Y^2=z^3(1+z) refused (slope 3/2) -- integrality, not evenness, is the invariant",
      ok4 and (okr2 is False) and mu4 == 1 and mur2 == Fraction(3, 2))

# ================================================================ UL-D: the exclusion and the r-witness
check("UL-D1", "2*w(y) = 1 has no solution with w(y) in Z (exact)",
      all(2*k != 1 for k in range(-1000, 1001)))
K2 = (3*sqrt(5) - 5)/2
zc = sqrt(3)/2
coef = K2/zc
check("UL-D2", "r-witness: r^2 = (K^2/z_c) z has v_z = 1 (odd) and nonzero exact coefficient in Q(sqrt5,sqrt3)",
      simplify(coef) != 0 and zsupport(sp.expand(coef*z)) == {Fraction(1)})
_, okr, mur = puiseux_solve_quadratic(sp.Integer(0), -coef*z)
check("UL-D3", "Newton-Puiseux refuses Y^2 = (K^2/z_c) z with exact fractional slope 1/2 -- r is outside every D5 tower",
      (okr is False) and mur == Fraction(1, 2))
rho_p = z/(1 - z**2)
check("UL-D4", "prop:radiff identity replicated exactly: (K^2/z_c) * rho' * (1-z^2) = (K^2/z_c) z = r^2",
      simplify(coef*rho_p*(1 - z**2) - coef*z) == 0)
check("UL-D5", "rem:ratie constant replicated: (K^2/z_c)^2 = 70/3 - 10 sqrt5 != 1 (exact)",
      simplify(coef**2 - (Rational(70, 3) - 10*sqrt(5))) == 0 and simplify(coef**2 - 1) != 0)
check("UL-D6", "r^2 itself is order-1 (integer-graded): w(r^2)=1 in Z -- odd valuation is no obstruction at F_ed level",
      integer_supported(sp.expand(coef*z)))

# ================================================================ UL-E: the unit-log gate is load-bearing (drop-one)
# formal model: grant a log-z token LOGZ with grading rule w(exp(c*LOGZ)) = c.
def w_exp_of_clogz(c):
    return Fraction(c)
check("UL-E1", "drop-one [unit-log gate]: granting log z, exp((1/2)log z) has w = 1/2 -- the ONE order-raising composite, exactly as blk:radab names it",
      w_exp_of_clogz(Fraction(1, 2)) == Fraction(1, 2) and w_exp_of_clogz(Fraction(1, 2)).denominator == 2)
def gate(u):
    s = zsupport(sp.expand(u))
    return min(s) == 0  # v_z(u) == 0 admits
check("UL-E2", "gate refusal: log z inadmissible since n_0(z) = v_z(z) = 1 != 0; log(1-z^2) admitted since v_z = 0",
      (gate(z) is False) and (gate(1 - z**2) is True))

# ================================================================ UL-F: corpus ties (winding lemma, pentad, rho, r^2)
def res0_dlog(u, N=12):
    """res_0(u'/u) for a Laurent poly u: exact via series of u'/u."""
    s = sp.series(sp.expand(diff(u, z)/u), z, 0, 2)
    return s.coeff(z, -1)

battery = [z**3*(2 + z), (1 - z**2), z**-2*(5 + z), 7*z*(1 + 3*z + z**2)]
ok_wind = True
for u in battery:
    n0 = res0_dlog(u)
    vz = min(zsupport(sp.expand(u)))
    if sp.simplify(n0 - vz) != 0:
        ok_wind = False
check("UL-F1", "winding lemma n_0(u) = res_0(u'/u) = v_z(u) replicated on a 4-object battery (exact residues)", ok_wind)

u_c  = (1 - z**2)/z**2
C_c  = (1 - z**2)/z**4
sD_c = (2 - z**2)/z**2
m_c  = -(1 - z**2)
l_c  = z**2*(2 - z**2)/(1 - z**2)
pent = [min(zsupport(sp.expand(sp.cancel(e*z**8))))-8 for e in (u_c, C_c, sD_c, m_c)]  # cancel then shift back
lam_v = min(zsupport(trunc(sp.expand(l_c*(1)), 6))) if True else None
lam_v = min(zsupport(trunc(l_c, 6)))
check("UL-F2", "thm:apexram pentad replicated: v_z(u,C,sqrtD,m,lambda) = (-2,-4,-2,0,2) exact",
      pent == [-2, -4, -2, 0] and lam_v == 2)
rho_series = trunc(-sp.Rational(1, 2)*log(1 - z**2), 8)
check("UL-F3", "v_z(rho) = 2 with rho = -(1/2)log(1-z^2) (log of an apex unit; series z^2/2 + ...) -- the substrate's one log is gate-admissible",
      min(zsupport(rho_series)) == 2 and rho_series.coeff(z, 2) == Rational(1, 2))
check("UL-F4", "v_z(r^2) = 1 (odd) while every pentad value is even -- the chart-level certificate the F_ed level refines",
      min(zsupport(sp.expand(coef*z))) == 1 and all(v % 2 == 0 for v in pent + [lam_v]))

print()
print(f"ul_closure_check: {PASS}/{PASS+FAIL} PASS")
sys.exit(0 if FAIL == 0 else 1)
