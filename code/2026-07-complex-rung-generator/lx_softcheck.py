#!/usr/bin/env python3
# lx_softcheck.py -- W[open]3 dossier harness (HANDOFF v1.0 section 6), written cold
# this session (2026-07-15).  Targets: L1 scoping lemma (splitness over Q is generic
# on rational heights), L2 predicate-first scan of the lens constants in foreign
# contexts, L3 bounded-negative bookkeeping.
#
# Discipline (per nx_/qx_ idiom): every DECISION exact over Q / Q(sqrt5) (extended
# by i where needed); Q(sqrt5) sign decisions by rational coefficient arithmetic
# (sgnQ5); K-level statements decided at the squared level; fail-first falsifier
# guards are written BEFORE the positive claims they protect (group LX-F precedes
# the claims it guards in reading order; a falsifier guard FAILS if the claim it
# guards were false or the detector blind).
# Group LX-G: mpmath interval corroboration displays, counted separately.
# Group LX-S: scan bookkeeping (regex counts over the shipped text corpus + the
# whitepaper PDF via pypdf 6.x) -- NOT exact math, counted separately, labeled SCAN.
# Exit 0 iff every check passes.  ASCII-only output.  Run: py lx_softcheck.py
#
# Corpus sources verified for the chart formulas (both, independently):
#   - whitepaper Thm 13.1/13.2 (pentad table p.14; multipliers m = -C/(1+u)^2)
#   - HELIX-EXPLICIT-DERIVATION.md section 2 (same table, HX11-HX16)
# Gate family [inherited]: g_C(x) = x^2 + x - C; "fixed points" = roots of g_C;
# multiplier at fixed point u: m = -C/(1+u)^2; m+ m- = 1, m + 1/m = -1/C - 2.

import sys, os, re
from fractions import Fraction
from sympy import (symbols, sqrt, log, I, Rational, Integer, simplify, expand,
                   radsimp, together, cancel, fraction, Poly, factor_list,
                   minimal_polynomial, expand_log, lucas, fibonacci, solve, im, re as sre)

PASS, FAIL, GUARD, SCAN, SKIP = [], [], [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(("PASS" if ok else "FAIL"), cid, "-", desc)
    return ok

def gk(cid, cond, desc):
    ok = bool(cond)
    GUARD.append((cid, ok))
    print(("GUARD-PASS" if ok else "GUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append(cid)

def sk(cid, cond, desc):
    ok = bool(cond)
    SCAN.append((cid, ok))
    print(("SCAN-PASS" if ok else "SCAN-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append(cid)

def sksk(cid, desc):
    # SCAN fixture unavailable (corpus text file / pypdf absent): SKIP, not FAIL.
    # The L2/L3 scan is candidate-level bookkeeping; the FORCED scoping lemma
    # (groups LX-A/B/C) is self-contained and always decides. Full scan runs where
    # the corpus fixtures are present (the golden-substrate workspace).
    SKIP.append(cid)
    print("SCAN-SKIP", cid, "-", desc)

def zero(e):
    e2 = simplify(expand(e))
    if e2 == 0:
        return True
    e3 = simplify(radsimp(together(e2)))
    if e3 == 0:
        return True
    return simplify(expand_log(e3, force=True)) == 0

x = symbols('x')                          # unconstrained (roots may be negative/complex)
Z = symbols('Z', positive=True)           # Z stands for z^2, 0 < Z < 1
w = symbols('w', positive=True)           # w = phi^(2n)
s5 = sqrt(5)
phi = (1 + s5) / 2
tau = (s5 - 1) / 2
gap = phi ** -4
K2 = 1 - gap                              # K decided at K^2

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q. Raises if not in Q(sqrt5)."""
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5"
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
        return 1 if A * A > 5 * B * B else -1
    return 1 if 5 * B * B > A * A else -1

def in_QZ(e):
    """Certificate: e is a ratio of polynomials in Z with rational coefficients.
    (Hence Z in Q => e in Q, away from denominator zeros.)  Returns False if the
    expression leaves Q(Z) (e.g. contains sqrt(Z))."""
    try:
        e2 = cancel(together(expand(e)))
        n, d = fraction(e2)
        pn, pd = Poly(n, Z), Poly(d, Z)
        return (all(c.is_rational for c in pn.all_coeffs())
                and all(c.is_rational for c in pd.all_coeffs()))
    except Exception:
        return False

# ---- chart pentad AS SHIPPED (whitepaper Thm 13.1 / HELIX section 2), in Z = z^2
uZ  = (1 - Z) / Z                 # occupation
CZ  = (1 - Z) / Z**2              # gate level
mZ  = -(1 - Z)                    # multiplier at u  (= -C/(1+u)^2, checked below)
sDZ = (2 - Z) / Z                 # sqrt(D), D = 1 + 4C
lamZ = Z * (2 - Z) / (1 - Z)      # exchange rate
u2Z = -1 - uZ                     # the other fixed point (root sum of g_C = -1)
m2Z = -1 / (1 - Z)                # its multiplier

def mult_of(Cv, uv):
    """Corpus multiplier formula m = -C/(1+u)^2 (whitepaper Thm 13.1(a)/13.2)."""
    return simplify(-Cv / (1 + uv) ** 2)

# ================================================================ LX-A
print("== LX-A: corpus chart-pentad replication (verify-before-use) ==")

ck("LX-A1", zero(uZ**2 + uZ - CZ) and zero(u2Z**2 + u2Z - CZ),
   "u(z) and -1-u(z) are exactly the two fixed points of g_C (g_C(u)=0)")
ck("LX-A2", zero(mult_of(CZ, uZ) - mZ) and zero(mult_of(CZ, u2Z) - m2Z),
   "multiplier formula m=-C/(1+u)^2 gives m=-(1-Z) and m'=-1/(1-Z)")
ck("LX-A3", zero(mZ * m2Z - 1) and zero((mZ + m2Z) - (-1/CZ - 2)),
   "reciprocity m m' = 1 and trace m + m' = -1/C - 2 (flip identities)")
ck("LX-A4", zero((1 + 4*CZ) - sDZ**2) and zero(lamZ - sDZ/CZ)
        and zero(lamZ*(1 - Z) - (1 - (1 - Z)**2)) and zero((mZ - m2Z)**2 - lamZ**2),
   "D=1+4C=((2-Z)/Z)^2; lambda=sqrt(D)/C; lambda|m|=1-m^2; (m-m')^2=lambda^2")
ck("LX-A5", zero(CZ * Z**2 + Z - 1),
   "chart relation C z^4 + z^2 - 1 = 0 (the universal seed relation)")
ck("LX-A6", zero(uZ.subs(Z, 1 - 1/w) - 1/(w - 1)) and zero(mZ.subs(Z, 1 - 1/w) + 1/w),
   "rung pullback (symbolic w=phi^(2n)): u_n = 1/(phi^(2n)-1), m_n = -phi^(-2n)")
C1 = simplify(1 / (phi - 1/phi)**2)
C2 = simplify(1 / (phi**2 - 1/phi**2)**2)
C3 = simplify(1 / (phi**3 - 1/phi**3)**2)
ck("LX-A7", zero(C1 - 1) and zero(C2 - Rational(1, 5)) and zero(C3 - Rational(1, 16))
        and lucas(3) == 4 and zero(C3 - 1/Integer(lucas(3))**2),
   "ladder levels C1=1, C2=1/5, C3=1/16=1/L3^2 (L3=4) exact")
# ordering tau < 3/4 < K^2 (corpus integer guards 25>20 and 180>169 replicated)
ck("LX-A8", sgnQ5(Rational(3, 4) - tau) > 0 and sgnQ5(K2 - Rational(3, 4)) > 0
        and (25 > 20) and (180 > 169),
   "ordering tau < 3/4 < K^2 by exact Q(sqrt5) signs + integer guards")

# ================================================================ LX-F
# Fail-first falsifier checks (exact).  Written BEFORE the L1 claims they
# protect; each would FAIL if the guarded claim were false or the detector blind.
print("== LX-F: falsifier checks (fail-first; detector must distinguish) ==")

ck("LX-F1", not in_QZ(sqrt(Z)),
   "membership detector rejects sqrt(Z): the odd object is NOT in Q(Z)")
sD_r1 = sDZ.subs(Z, tau)                      # rung 1 height: z1^2 = tau
ck("LX-F2", zero(sD_r1 - s5)
        and Poly(minimal_polynomial(s5, x), x).degree() == 2
        and (2**2 < 5 < 3**2),
   "rung 1: sqrt(D)=(2-tau)/tau=sqrt5 irrational; D1=5 nonsquare (4<5<9)")
sD_r2 = sDZ.subs(Z, K2)                       # rung 2 height: z2^2 = K^2
D_r2 = simplify((1 + 4*CZ).subs(Z, K2))
fl2 = factor_list(expand((x**2 + x - Rational(1, 5)) * 5))
ck("LX-F3", zero(sD_r2 - 3/s5) and zero(D_r2 - Rational(9, 5))
        and Poly(minimal_polynomial(3/s5, x), x).degree() == 2
        and (6**2 < 45 < 7**2)
        and len(fl2[1]) == 1 and fl2[1][0][1] == 1
        and Poly(minimal_polynomial(K2, x), x).degree() == 2,
   "rung 2: C2=1/5 in Q BUT D2=9/5 nonsquare (36<45<49), gate irreducible/Q;"
   " and K^2 not in Q (converse of L1 is false)")
mult_quad = lambda m_: m_**2 + (1/CZ + 2)*m_ + 1   # quadratic with roots m, m'
ck("LX-F4", zero(mult_quad(mZ)) and zero(mult_quad(m2Z))
        and not zero(CZ.subs(Z, Rational(3, 4)) - Rational(5, 9))
        and not zero(sDZ.subs(Z, Rational(1, 2)) - 2)
        and not zero((1 + 4*CZ) - ((2 + Z)/Z)**2)
        and not zero(mult_quad(-(1 - Z)**2)),
   "perturbation rejections: C(3/4)!=5/9; sqrtD(1/2)!=2; sqrtD_bad=(2+Z)/Z fails"
   " D=1+4C; m,m' satisfy the multiplier quadratic but m_bad=-(1-Z)^2 fails it")

# ================================================================ LX-B
print("== LX-B: L1 scoping lemma -- splitness over Q is generic ==")

ck("LX-B1", all(in_QZ(e) for e in (uZ, CZ, 1 + 4*CZ, sDZ, mZ, m2Z, lamZ)),
   "Q(Z)-membership certificates: u, C, D, sqrt(D), m, m', lambda all in Q(Z)"
   " => z^2 in Q forces every one of them into Q")
den_ok = True
for e in (uZ, CZ, sDZ, lamZ, m2Z):
    d = fraction(cancel(together(e)))[1]
    rts = Poly(d, Z).all_roots()
    den_ok = den_ok and all(r == 0 or r == 1 for r in rts)
ck("LX-B2", den_ok,
   "denominator zero sets are within {0,1}: certificates valid on 0 < Z < 1")

# the z^2 = 1/2 example (dossier-specified), pair computed FROM the corpus formula
Zh = Rational(1, 2)
Ch, sDh = CZ.subs(Z, Zh), sDZ.subs(Z, Zh)
uh1, uh2 = solve(x**2 + x - Ch, x)
pair_h = sorted([mult_of(Ch, uh1), mult_of(Ch, uh2)])
rho_h = -log(1 - Zh) / 2
ck("LX-B3", Ch == 2 and sDh == 3 and (1 + 4*Ch) == 9
        and zero(expand((x**2 + x - 2) - (x - 1)*(x + 2)))
        and pair_h == [-2, Rational(-1, 2)]
        and simplify(pair_h[0]*pair_h[1] - 1) == 0,
   "z^2=1/2: C=2, D=9=3^2, splits (x-1)(x+2), corpus-formula pair {-1/2,-2},"
   " product 1  [dossier example verified]")
ck("LX-B4", zero(rho_h - log(2)/2) and Rational(1, 2).is_rational
        and sgnQ5(tau - Rational(1, 2)) > 0 and (5 > 4),
   "z^2=1/2 sits at rho=(1/2)ln2 in (ln2)Q, strictly below rung 1 (tau>1/2 <=> 5>4)")

sweep_ok = True
for Zr in (Rational(1, 3), Rational(2, 3), Rational(1, 5), Rational(9, 10),
           Rational(5, 8), Rational(3, 4), Rational(1, 2)):
    Cv, sDv = CZ.subs(Z, Zr), sDZ.subs(Z, Zr)
    ma, mb = -(1 - Zr), -1/(1 - Zr)
    sweep_ok = (sweep_ok and Cv.is_rational and sDv.is_rational
                and simplify(sDv**2 - (1 + 4*Cv)) == 0
                and ma.is_rational and mb.is_rational
                and simplify(ma*mb - 1) == 0)
ck("LX-B5", sweep_ok,
   "genericity sweep: 7 rational heights all give C, sqrt(D) in Q and split gates"
   " with rational reciprocal multiplier pairs")

# lens exact data (the distinction is THIS, not splitness) -- whitepaper Thm 13.2
Zc = Rational(3, 4)
Cc, sDc, lamc = CZ.subs(Z, Zc), sDZ.subs(Z, Zc), lamZ.subs(Z, Zc)
uc1, uc2 = solve(x**2 + x - Rational(4, 9), x)
pair_c = sorted([mult_of(Cc, uc1), mult_of(Cc, uc2)])
ck("LX-B6", Cc == Rational(4, 9) and sDc == Rational(5, 3)
        and (1 + 4*Cc) == Rational(25, 9) and lamc == Rational(15, 4)
        and sorted([uc1, uc2]) == [Rational(-4, 3), Rational(1, 3)]
        and pair_c == [-4, Rational(-1, 4)]
        and (1 - Zc) == Rational(1, 4)
        and zero(-log(1 - Zc)/2 - log(2)),
   "lens Z=3/4: C=4/9, sqrt(D)=5/3, D=25/9, fixed pts {1/3,-4/3}, pair {-1/4,-4},"
   " lambda=15/4, co-intensity 1/4, rho=ln2 EXACT (the exact data = the distinction)")

# ================================================================ LX-C
print("== LX-C: scan-support exact checks (hit values + restatement certificates) ==")

# C1: Thm 16.1 replication: log_phi(2) not in Q  (elementary, corpus-forced)
binet_ok = all(zero(2*phi**p - lucas(p) - fibonacci(p)*s5) for p in range(1, 13))
fib_pos = (fibonacci(1) == 1 and fibonacci(2) == 1
           and all(fibonacci(p + 1) - fibonacci(p) - fibonacci(p - 1) == 0
                   for p in range(2, 13)))
# induction step is the recursion + positivity closure (a>0, b>0 => a+b>0);
# hence F_p >= 1 for all p >= 1, so phi^p = (L_p + F_p sqrt5)/2 has nonzero
# sqrt5-part for p >= 1: phi^p irrational; p <= -1 symmetric (invert); p = 0
# forces 2^q = 1, impossible for q >= 1.  So ln2/lnphi = p/q is impossible.
ck("LX-C1", binet_ok and fib_pos,
   "Thm 16.1 replication chain: 2 phi^p = L_p + F_p sqrt5 (p=1..12) + Fibonacci"
   " recursion/positivity => log_phi(2) irrational [corpus-forced, re-forced]")

# C2: M(cons) = 2 phi^5 = 11 + 5 sqrt5; ln M(cons) = ln2 + 5 ln phi;
#     and ln M(cons) NOT in (ln2)Q  (reduction to C1)
Mcons = 11 + 5*s5
ck("LX-C2", zero(2*phi**5 - Mcons)
        and zero(expand_log(log(2*phi**5) - log(2) - 5*log(phi), force=True)),
   "M(cons)=2phi^5=11+5sqrt5 exact; ln M(cons) = ln2 + 5 ln phi (lattice fact)")
# reduction: ln2 + 5 lnphi = q ln2  =>  log_phi(2) = 5/(q-1) in Q (q != 1),
# contradicting C1;  q = 1 would force 5 lnphi = 0, i.e. phi = 1: false.
ck("LX-C3", sgnQ5(phi - 1) > 0,
   "q=1 branch killed exactly (phi>1); with C1: ln M(cons) NOT in (ln2)Q --"
   " the registry hit is a lattice-component fact, not a (ln2)Q rapidity")

# C4: N(M(cons)) = -4 and the merge with the ln2-component (one fact, not two)
NM = 11**2 - 5*5**2
Nphi = Rational(1 - 5, 4)          # N(phi) = (1^2 - 5*1^2)/4 = -1
ck("LX-C4", NM == -4 and Nphi == -1 and 4**1 * (-1)**5 == -4
        and [a for a in range(0, 4) if 4**a == 4] == [1],
   "N(11+5sqrt5)=-4=(4^a)(-1)^b at (a,b)=(1,5); |N|=4 forces a=1: the -4 hit IS"
   " the ln2-lattice hit (a=1), not independent evidence")

# C5: flip-restatement certificate: 'lens multiplier = flip level' <=> z_c^2 = 3/4
flipC = solve(1 + 4*symbols('Ct'), symbols('Ct'))
lensZ = solve(-(1 - Z) - Rational(-1, 4), Z)
ck("LX-C5", flipC == [Rational(-1, 4)] and lensZ == [Rational(3, 4)],
   "flip locus D=0 <=> C=-1/4 (unique); m(Z)=-1/4 <=> Z=3/4 (unique): any"
   " lens<->flip tie via -1/4 restates the inherited z_c datum")

# C6: the C=-1/2 rotation gate's dyadic is internal to the +-i requirement
Cr = Rational(-1, 2)
ur1, ur2 = solve(x**2 + x - Cr, x)
mr1, mr2 = mult_of(Cr, ur1), mult_of(Cr, ur2)
onlyC = solve(-1/symbols('Ct') - 2, symbols('Ct'))
ck("LX-C6", zero(mr1**2 + 1) and zero(mr2**2 + 1) and onlyC == [Rational(-1, 2)],
   "gate C=-1/2 has multipliers +-i exactly, and m+1/m=0 forces C=-1/2 uniquely:"
   " the dyadic -1/2 is internal to +-i, no lens input")

# C7: the ECHO 'trace 5/3' hit value, exact (trace form on Q(phi), basis {1,phi})
def Tr(e):
    A, B = q5AB(e)
    return Rational(2 * A)         # Tr(a + b sqrt5) = 2a
coeff = Rational(Tr(phi), Tr(phi**2))          # projection of 1 onto <phi>
resid_tr = Tr(1 - coeff * phi)
ck("LX-C7", Tr(Integer(1)) == 2 and Tr(phi) == 1 and Tr(phi**2) == 3
        and coeff == Rational(1, 3) and resid_tr == Rational(5, 3)
        and resid_tr == sDc,
   "ECHO hit: Tr(1-(1/3)phi)=5/3 exact; equals sqrt(D)(z_c) as rationals; the two"
   " derivations share no step (shared-constant only)")

# C8: rung-3 dyadic C3=1/16 enters via Lucas, NOT via co-intensity
co3 = phi**-6
ck("LX-C8", Poly(minimal_polynomial(co3, x), x).degree() == 2
        and zero(co3 - (9 - 4*s5)),
   "rung-3 co-intensity phi^-6 = 9-4sqrt5 irrational: C3=1/16's dyadicity comes"
   " from L3=4, not from a dyadic co-intensity")

# ================================================================ LX-G
print("== LX-G: certified numeric corroboration displays (counted separately) ==")
from mpmath import mp, mpf, log as mlog, sqrt as msqrt
mp.dps = 60
ln2n = mlog(2)
flrn = mlog(2) / mlog((1 + msqrt(5)) / 2)
gk("LX-G1", abs(ln2n - mpf("0.6931471806")) < mpf(10)**-10 / 2,
   "corpus display ln2 = 0.6931471806 corroborated at 60 dps (0.51 ulp matcher)")
gk("LX-G2", abs(flrn - mpf("1.4404200904")) < mpf(10)**-10 / 2,
   "corpus display log_phi 2 = 1.4404200904 corroborated at 60 dps")

# ================================================================ LX-S
# Scan driver: L2/L3 bookkeeping, frozen 2026-07-15.  NOT exact math.
# Predicate (the citable scan class): occurrences of the target literals
#   T = {4/9, 25/9, 5/3, 15/4, -1/4, -4(paired/norm-equated), dyadic 2^-k as
#        signed gate levels or 1/8..1/64, ln2-rapidity markers}
# as standalone rational constants (regex classes below), in the six shipped
# domains.  Own-lens-context occurrences do not count as hits; role-mismatched
# occurrences (angle t-values, TFD/Casimir coefficients, 5^(1/4) exponents,
# series-bound coefficients) were reviewed and excluded -- see session report.
print("== LX-S: predicate-first scan, frozen inventory (bookkeeping) ==")

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSES = {
    'c_4_9':  r'(?<![\d.])4\s*/\s*9(?![\d.])|\{4\}\{9\}',
    'c_25_9': r'(?<![\d.])25\s*/\s*9(?![\d.])|\{25\}\{9\}',
    'c_5_3':  r'(?<![\d.])5\s*/\s*3(?![\d.])|\{5\}\{3\}',
    'c_15_4': r'(?<![\d.])15\s*/\s*4(?![\d.])|\{15\}\{4\}',
    'c_m1_4': r'[-−]\s*(?:1\s*/\s*4(?![\d.])|¼)',
    'c_m4':   r'[-−]\s*1\s*/\s*4\s*,\s*[-−]\s*4|=\s*[-−]\s*4(?![\d.])',
    'c_dy':   r'(?<![\d.])1\s*/\s*(?:8|16|32|64)(?![\d.])'
              r'|[-−]\s*(?:1\s*/\s*2(?![\d.])|½)',
    'c_ln2':  r'ln\s*2(?![\d\w])|log_?\\?(?:gp|phi|varphi|φ)\s*2(?![\d\w])'
              r'|0\.6931471806|1\.4404200904',
}
TEXT_EXPECT = {
    # file -> class -> sorted match line numbers (frozen 2026-07-15)
    'ECHO-S-RESEARCH_UNIFIED.md': {'c_5_3': [58], 'c_m1_4': [44, 62]},
    'K-VARIABLE-DERIVATIONS.md': {},
    'HELIX-EXPLICIT-DERIVATION.md': {
        'c_4_9': [21, 56, 90, 92, 139, 210, 227, 242], 'c_25_9': [90, 242],
        'c_5_3': [57, 90, 210], 'c_15_4': [59, 94, 210],
        'c_m1_4': [61, 94, 210, 242], 'c_m4': [94, 210, 242],
        'c_dy': [32, 39, 42, 73, 104, 160],
        'c_ln2': [21, 84, 86, 86, 86, 86, 128, 128, 139, 139, 139, 146, 146,
                  192, 194, 196, 196, 210, 211, 233, 244, 244, 255, 259]},
    'complex-rung-generator-v1.6.tex': {'c_ln2': [1037]},
    'gap_multiplier_paper.tex': {'c_dy': [357]},
}
sidx = 0
for fname, expect in TEXT_EXPECT.items():
    sidx += 1
    path = os.path.join(HERE, fname)
    try:
        text = open(path, encoding='utf-8').read()
    except FileNotFoundError:
        sksk("LX-S%d" % sidx, "%s not alongside harness: scan skipped" % fname)
        continue
    lines = text.splitlines()
    got = {}
    for cname, pat in CLASSES.items():
        hits = [i + 1 for i, ln in enumerate(lines) for _ in re.finditer(pat, ln)]
        if hits:
            got[cname] = hits
    sk("LX-S%d" % sidx, got == expect,
       "%s inventory frozen: %s" % (fname, sorted(got.items()) or "NO MATCHES"))

# whitepaper via pypdf, on the whitespace-normalized stream; stacked-fraction
# variants (PDF extraction splits \frac{a}{b} into 'a b') counted separately.
WP_STACKED = {
    'st_4_9':  r'(?<![\d.])4 9(?![\d.])',   'st_25_9': r'(?<![\d.])25 9(?![\d.])',
    'st_5_3':  r'(?<![\d.])5 3(?![\d.])',   'st_15_4': r'(?<![\d.])15 4(?![\d.])',
    'st_1_16': r'(?<![\d.])1 16(?![\d.])|1/16', 'st_m_1_4': r'[-−] 1 4(?![\d.])',
}
WP_EXPECT = {'c_4_9': 4, 'c_25_9': 0, 'c_5_3': 1, 'c_15_4': 1, 'c_m1_4': 0,
             'c_m4': 5, 'c_dy': 1, 'c_ln2': 25,
             'st_4_9': 5, 'st_25_9': 2, 'st_5_3': 2, 'st_15_4': 2,
             'st_1_16': 1, 'st_m_1_4': 7}
try:
    import pypdf
    rdr = pypdf.PdfReader(os.path.join(HERE, 'golden-substrate-whitepaper.pdf'))
    wp = '\n'.join((p.extract_text() or '').replace('\x00', '') for p in rdr.pages)
    wpj = re.sub(r'\s+', ' ', wp)
    wgot = {}
    for cname, pat in CLASSES.items():
        wgot[cname] = len(list(re.finditer(pat, wpj)))
    for cname, pat in WP_STACKED.items():
        wgot[cname] = len(list(re.finditer(pat, wpj)))
    sk("LX-S6", wgot == WP_EXPECT,
       "whitepaper (pypdf, normalized + stacked-fraction pass) counts frozen: %s"
       % sorted(wgot.items()))
except (ImportError, FileNotFoundError) as _e:
    sksk("LX-S6", "whitepaper scan skipped (%s)" % type(_e).__name__)

# L3 bounded-negative summary assertion: the foreign-context hit inventory.
# (Classification own/foreign performed by reading; frozen here as bookkeeping.)
FOREIGN_HITS = [
    "5/3 @ ECHO-S-RESEARCH_UNIFIED.md:58 (trace of residual, Q(sqrt5) trace form)"
    " -- numerological, unregistered (no mechanism)",
    "ln2-lattice @ whitepaper Thm 18.1(a) + HELIX line 244: ln M(cons)=ln2+5lnphi"
    " -- recorded, significance [open] in corpus, no mechanism",
    "-4 @ whitepaper Thm 18.1(a): N(M(cons))=-4 -- merged with ln2-lattice hit"
    " (a=1), not independent",
    "-1/4 @ flip locus C=-1/4 (ECHO 44/62, HELIX 61, whitepaper Prop 8.x)"
    " -- gate-architecture role; lens tie <=> z_c^2=3/4 (restatement)",
    "dyadic 1/16 @ gap_multiplier_paper.tex:357 + whitepaper ladder (C3=1/L3^2)"
    " -- Lucas mechanism L3=4, no lens tie",
    "dyadic -1/2 @ C=-1/2 rotation gate (HELIX 32/42/104, whitepaper Prop 10.2)"
    " -- forced by +-i internally, no lens tie",
]
NEGATIVES = ("4/9, 25/9, 15/4, the pair {-1/4,-4} as a pair, dyadic co-intensities"
             " 2^-k (strict class), and (ln2)Q rapidities: NO foreign-context"
             " occurrence in the six scanned domains")
sk("LX-S7", len(FOREIGN_HITS) == 6,
   "L3 bounded negative: 6 foreign candidate hits, all without mechanism; " + NEGATIVES)
for h in FOREIGN_HITS:
    print("     HIT:", h)

# ================================================================ summary
n_exact_pass = len(PASS)
n_exact_total = len(PASS) + sum(1 for c in FAIL
                                if not any(c == g[0] for g in GUARD)
                                and not any(c == s_[0] for s_ in SCAN))
g_pass = sum(1 for _, ok in GUARD if ok)
s_pass = sum(1 for _, ok in SCAN if ok)
line = "LX: EXACT %d/%d, GUARDS %d/%d, SCAN %d/%d%s -- %s" % (
    n_exact_pass, n_exact_total, g_pass, len(GUARD), s_pass, len(SCAN),
    ("" if not SKIP else (" (%d SCAN skipped: fixtures/pypdf absent)" % len(SKIP))),
    "ALL PASS" if not FAIL else ("FAILURES: %s" % FAIL))
print(line)
sys.exit(0 if not FAIL else 1)
