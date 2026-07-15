#!/usr/bin/env python3
# rx_softcheck.py -- OP-RADIUS dossier harness (HANDOFF-v1.0 section 5), written cold
# this session.  Routes: RAD-0 (parity no-go), RAD-2 (area-transfer reformulation),
# RAD-4 (rotating-equipotential test).
#
# Discipline (qx/nx idiom): every DECISION exact in Q / Q(sqrt5) (extended by sqrt3
# and by symbols where needed); K-level statements decided at the squared level;
# mpmath iv interval arithmetic at 60 dps used ONLY for the certified guards
# (RX-G group, counted separately).  Interval DISJOINTNESS is a rigorous proof of
# non-equality; interval overlap proves nothing (the falsifier guard RX-G1 pins
# this honestly).  Floats never touch an exact decision.
#
# Corpus anchors (verified against the shipped files this session):
#   - W-Thm 13.1 (whitepaper p.13-14): u=(1-z^2)/z^2, m=-(1-z^2), C=(1-z^2)/z^4,
#     sqrtD=(2-z^2)/z^2, lambda=z^2(2-z^2)/(1-z^2).   [NOTE: the handoff dossier's
#     "u = z^2/(1-z^2)" is inverted; corpus u is the reciprocal.  Verified here.]
#   - D1 (W-Decl 12.1): 1-z^2 = e^{-2 rho}, rho(z) = -(1/2) ln(1-z^2).
#   - D3 (W-Decl 12.3): theta = kappa rho, kappa = pi/(2 ln phi).
#   - inherited radius (W-(7)/(9), [open] 1; register defect D-1 corrected wording):
#     r(z) = K sqrt(z/z_c) on [0, z_c], r = K on [z_c, 1], z_c = sqrt(3)/2.
#   - v1.6 thm:stateangle: same pentad as trig dictionary (replicated RX-A2).
#   - W-Thm 15.1 / Prop 15.2a: z=r point algebra, kink r'(z_c-) = K/(2 z_c).
#
# Groups:
#   RX-A  corpus replication / definition pinning                          (4)
#   RX-P  RAD-0 parity no-go (fail-first falsifiers precede acceptances)   (9)
#   RX-T  RAD-2 area-transfer law (fail-first falsifiers first)            (4)
#   RX-E  RAD-4 rotating equipotential + clock ledger                      (7)
#   RX-G  certified-interval guards, counted separately                    (4)
#
# Exit 0 iff every exact check passes and every guard certifies.

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   tan, cot, csc, simplify, expand, radsimp, trigsimp, together,
                   cancel, fraction, Poly, diff, solve, integrate, expand_log)

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

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q.  Raises if not in Q(sqrt5)."""
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, sqrt(5))
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(sqrt(5)) if p.degree() == 1 else Integer(0)
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

z = symbols('z')
h = symbols('h')
alpha = symbols('alpha', positive=True)
w = symbols('w')
kk = symbols('k', integer=True)
omega, gg = symbols('omega g', positive=True)
c1, c2, c3, c4, c5, c6 = symbols('c1:7')

s5   = sqrt(5)
phi  = (1 + s5)/2
tau  = (s5 - 1)/2
gap  = phi**-4
K2   = 1 - gap
K    = sqrt(K2)
zc   = sqrt(3)/2
kappa = pi/(2*log(phi))

# chart pentad as height functions (W-Thm 13.1, corpus-verified)
u_  = (1 - z**2)/z**2
C_  = (1 - z**2)/z**4
sD  = (2 - z**2)/z**2
m_  = -(1 - z**2)
lam = z**2*(2 - z**2)/(1 - z**2)
rho   = -log(1 - z**2)/2          # D1
theta = kappa*rho                 # D2+D3
r2  = K2*z/zc                     # inherited horn branch, squared: r^2 = K^2 z/z_c
PENTAD = [("u", u_), ("C", C_), ("sqrtD", sD), ("m", m_), ("lambda", lam)]

# ------------------------------------------------------- parity detectors
def rat_numden(e):
    return fraction(cancel(together(expand(e))))

def even_cert(e):
    """Syntactic certificate e in Q(...)(z^2): canonical num and den have only
    even z-degrees.  Sound AND complete for rational functions: sigma z->-z maps
    the gcd-free pair to itself up to a shared sign; a shared odd sign would make
    z divide both -- impossible gcd-free."""
    num, den = rat_numden(e)
    pn, pd = Poly(num, z), Poly(den, z)
    return (all(mm[0] % 2 == 0 for mm in pn.monoms())
            and all(mm[0] % 2 == 0 for mm in pd.monoms()))

def sub_even(e):
    return simplify(cancel(together(e.subs(z, -z) - e))) == 0

def sub_odd(e):
    return simplify(cancel(together(e.subs(z, -z) + e))) == 0

# ================================================================ RX-A
print("== RX-A: corpus replication / definition pinning ==")
ck("RX-A1", zero(phi*tau - 1) and zero(K2 - (1 - gap)) and zero(K2**2 - 5*gap)
        and zero(K2 - (3*s5 - 5)/2) and zero(zc**2 - Rational(3, 4))
        and zero((1 - zc**2) - Rational(1, 4)),
   "seed identities: phi*tau=1, K^2=1-gap=(3sqrt5-5)/2, K^4=5gap, z_c^2=3/4")

zA = sin(alpha)
ck("RX-A2", zero(u_.subs(z, zA) - cot(alpha)**2)
        and zero(C_.subs(z, zA) - cot(alpha)**2*csc(alpha)**2)
        and zero(sD.subs(z, zA) - (1 + 2*cot(alpha)**2))
        and zero(m_.subs(z, zA) + cos(alpha)**2)
        and zero(lam.subs(z, zA) - tan(alpha)**2*(1 + cos(alpha)**2))
        and zero(exp(-2*rho.subs(z, zA)) - cos(alpha)**2),
   "pentad = state-angle dictionary (W-Thm 13.1 / v1.6 thm:stateangle); e^(-2 rho)=cos^2 (rho=-ln cos, decided at the squared level, alpha in (0,pi/2))")

z1 = sqrt(tau)
tab1 = (zero(u_.subs(z, z1) - tau) and zero(C_.subs(z, z1) - 1)
        and zero(sD.subs(z, z1) - s5) and zero(m_.subs(z, z1) + tau**2)
        and zero(lam.subs(z, z1) - s5))
tabc = (zero(u_.subs(z, zc) - Rational(1, 3)) and zero(C_.subs(z, zc) - Rational(4, 9))
        and zero(sD.subs(z, zc) - Rational(5, 3)) and zero(m_.subs(z, zc) + Rational(1, 4))
        and zero(lam.subs(z, zc) - Rational(15, 4)))
tabK = (zero(u_.subs(z, K) - 1/(phi**4 - 1)) and zero(C_.subs(z, K) - Rational(1, 5))
        and zero(sD.subs(z, K) - 3/s5) and zero(m_.subs(z, K) + gap)
        and zero(lam.subs(z, K) - 3*s5))
ck("RX-A3", tab1 and tabc and tabK,
   "W-Thm 13.1 evaluation table replicated at z1=sqrt(tau), z_c, z2=K")

ck("RX-A4", zero((K/(2*zc))**2 - (3*s5 - 5)/6) and sgnQ5(6*s5 - 13) > 0
        and zero(rho.subs(z, zc) - log(2))
        and expand(K2*z/zc - z**2 + z*(z - K2/zc)) == 0,
   "profile pins: kink slope^2 K^2/3=(3sqrt5-5)/6; K^2>z_c^2 (180>169); rho(z_c)=ln2; z=r factorization")

# ================================================================ RX-P  (RAD-0)
print("== RX-P: RAD-0 parity no-go ==")
# P1 -- FAIL-FIRST falsifiers: the detectors must reject planted odd/mixed
# functions and accept planted even ones (incl. hidden cancellations).
odd_plant   = z**3/(1 - z**2)
mixed_plant = z + z**2
hidden_odd  = (z**5 + z**3)/z**2          # = z^3 + z after cancellation
hidden_even = (z**7 + z**5)/z**3          # = z^4 + z^2 after cancellation
ck("RX-P1", (not even_cert(odd_plant)) and (not sub_even(odd_plant)) and sub_odd(odd_plant)
        and (not even_cert(mixed_plant)) and (not sub_even(mixed_plant)) and (not sub_odd(mixed_plant))
        and (not even_cert(hidden_odd)) and (not sub_even(hidden_odd))
        and even_cert(hidden_even) and sub_even(hidden_even),
   "falsifier: parity detectors reject planted odd/mixed (incl. cancelled), accept even")

ck("RX-P2", all(even_cert(p) and sub_even(p) for _, p in PENTAD)
        and zero(u_ - ((1 - w)/w).subs(w, z**2))
        and zero(C_ - ((1 - w)/w**2).subs(w, z**2))
        and zero(sD - ((2 - w)/w).subs(w, z**2))
        and zero(m_ - (-(1 - w)).subs(w, z**2))
        and zero(lam - (w*(2 - w)/(1 - w)).subs(w, z**2)),
   "every pentad member is a rational function of z^2 (certificate + explicit R(z^2) form)")

ck("RX-P2b", cancel(1/(1 + u_) - z**2) == 0 and cancel(u_/(1 + u_) - (1 - z**2)) == 0,
   "sharpness: z^2 = 1/(1+u) -- the pentad's rational closure EQUALS Q(sqrt5)(z^2)")

ck("RX-P3", simplify(rho - (-log(1 - w)/2).subs(w, z**2)) == 0 and sub_even(rho)
        and sub_even(theta)
        and sub_even(((pi/2 + 2*pi*kk)/log(phi))*rho)
        and zero((pi/2 + 2*pi*kk)/log(phi) - (4*kk + 1)*kappa),
   "rho = g(z^2), theta = kappa_k g(z^2) even for EVERY branch k (parity layer is branch-blind)")

ck("RX-P4", sub_odd(r2) and cancel((r2 + r2.subs(z, -z))/2) == 0
        and (not even_cert(r2)) and simplify(r2) != 0,
   "r^2 = K^2 z/z_c is odd (r^2(-z) = -r^2(z)), nonzero: even part vanishes identically")

wit = r2.subs(z, Rational(1, 2)) - r2.subs(z, Rational(-1, 2))
ck("RX-P5", zero(wit - (3*s5 - 5)/sqrt(3)) and sgnQ5(3*s5 - 5) > 0,
   "exact witness z0=1/2: r^2(z0)-r^2(-z0) = (3sqrt5-5)/sqrt3 > 0 (sign in Q(sqrt5))")

E1 = (c1*u_**3 - c2*C_*sD + c3*m_*lam**2)/(c4 + c5*u_*C_ + c6*lam)
E2 = (u_ - C_)/(m_ + lam*sD) + (C_**2*m_**3)/(u_ + 1)
E3 = (c1*theta**2 + c2*rho**3*u_ - c3*lam/(1 + rho**2))/(c4 + rho*theta)
ck("RX-P6", even_cert(E1) and sub_even(E1) and even_cert(E2) and sub_even(E2)
        and sub_even(E3),
   "closure: adversarial rational combos of the pentad (+rho,theta, generic coeffs) stay even")

r4 = cancel(r2**2)
ck("RX-P7", even_cert(r4) and sub_even(r4)
        and cancel((r2*u_ + (r2*u_).subs(z, -z))/2) == 0
        and cancel(r2*zc/K2 - z) == 0,
   "Z/2 sharpening: r^4 IS chart-even; even*r^2 stays odd; F(r^2) = Q(sqrt5)(z) (z = r^2 z_c/K^2)")

cands_even = [z**2, K2*z**2/zc**2, lam, u_*C_, sD*m_]
ck("RX-P8", all(simplify(cancel(together(r2 - cand))) != 0 for cand in cands_even),
   "direct corroboration: natural even chart candidates all differ from r^2 exactly")

ck("RX-P9", (not even_cert(z)) and sub_odd(z)
        and Rational(1, 2) - Rational(-1, 2) != 0,
   "the odd generator hunt is well-posed: z itself is odd -- z, sqrt(z) are bridge candidates")

# ================================================================ RX-T  (RAD-2)
print("== RX-T: RAD-2 area-transfer law ==")
ck("RX-T1", simplify(cancel(pi*r2 - (K2/(2*zc))*2*pi*z**2)) != 0
        and simplify(cancel(pi*r2 - (K2/(3*zc))*2*pi*z)) != 0
        and simplify(cancel(pi*r2 - (K2/(2*zc))*4*pi*z)) != 0,
   "falsifier: perturbed zone laws (z^2 zone, wrong constant, doubled area) all rejected")

ck("RX-T2", cancel((2*pi)**2*(1 - h**2)*(1 + (diff(sqrt(1 - h**2), h))**2) - (2*pi)**2) == 0
        and zero(integrate(2*pi, (h, 0, z)) - 2*pi*z),
   "Archimedes zone area forced: (2 pi sqrt(1-h^2))^2 (1+x'(h)^2) = (2 pi)^2; A_zone(z) = 2 pi z")

ck("RX-T3", zero(pi*r2 - (K2/(2*zc))*(2*pi*z))
        and zero(diff(pi*r2, z) - pi*K2/zc)
        and zero(pi*K2/zc - pi*(3*s5 - 5)/sqrt(3))
        and cancel(solve(pi*w - (K2/(2*zc))*2*pi*z, w)[0] - K2*z/zc) == 0,
   "pi r(z)^2 = (K^2/(2 z_c)) A_zone(z); deposition rate d(pi r^2)/dz = pi K^2/z_c constant; both directions")

ck("RX-T4", zero(pi*r2.subs(z, zc) - pi*K2),
   "lens continuity: transferred area at z_c equals the cap disc pi K^2 (C0 kink survives)")

# ================================================================ RX-E  (RAD-4)
print("== RX-E: RAD-4 rotating-equipotential test ==")
g_bad = solve(omega**2/(2*gg) - (zc/K2 + 1), gg)[0]
ck("RX-E1", simplify(cancel(g_bad - omega**2*K2/(2*zc))) != 0,
   "falsifier: perturbed horn coefficient z_c/K^2+1 yields g_eff != omega^2 K^2/(2 z_c)")

g_sol = solve(omega**2/(2*gg) - zc/K2, gg)[0]
ck("RX-E2", simplify(cancel(g_sol - omega**2*K2/(2*zc))) == 0
        and simplify(g_sol*z - omega**2*(K2*z/zc)/2) == 0,
   "horn z=(z_c/K^2)r^2 matches rigid-rotation surface z=(omega^2/2g)r^2 iff g_eff=omega^2 K^2/(2 z_c)")

g1 = kappa**2*K2/(2*zc)
ck("RX-E3", zero(g1 - pi**2*K2/(8*zc*log(phi)**2))
        and zero(((4*kk + 1)*kappa)**2*K2/(2*zc) - (4*kk + 1)**2*g1),
   "clock omega=kappa (per unit rho, D3): g_1 = pi^2 K^2/(8 z_c ln^2 phi); branch k scales by (4k+1)^2")

# T1 (handoff 3.4): i pi / ln phi is irrational.  Proof shape p/q => (-1)^q = phi^p.
det_hits = [(p, q) for q in range(1, 5) for p in range(-8, 9)
            if Integer(2)**p == Integer(4)**q]                 # log4/log2 = 2 must be FOUND
scan_ok = True
for p in range(-8, 9):
    if p == 0:
        continue
    val = expand(phi**p) if p > 0 else expand(tau**(-p))
    scan_ok = scan_ok and sgnQ5(val - 1) != 0 and sgnQ5(val + 1) != 0
wpos = symbols('w_pos', positive=True)
step = (wpos*(phi - 1)).is_positive                            # induction step w>0 => w(phi-1)>0
ck("RX-E4", len(det_hits) > 0 and scan_ok and sgnQ5(phi - 1) > 0 and step is True,
   "T1 forced: detector finds log4/log2 rational (falsifier); (-1)^q = phi^p impossible (phi^p != +-1, p != 0; phi>1 monotone; p=0 needs pi=0, killed by RX-G0)")

ck("RX-E5", sgnQ5(3*s5 - 5) > 0 and zero(g1/kappa**2 - K2/(2*zc))
        and zero(K2/(2*zc) - (3*s5 - 5)/(2*sqrt(3))),
   "g_1 = kappa^2 * (K^2/(2 z_c)) with K^2/(2 z_c) = (3sqrt5-5)/(2sqrt3) algebraic nonzero: g_1 algebraic iff kappa algebraic (Gelfond-Schneider kills it)")

dsdth = sqrt(K2 + ((1 - z**2)/(kappa*z))**2)
ck("RX-E6", zero(dsdth.subs(z, 1) - K) and zero((1/K)**2*K2/(2*zc) - 1/(2*zc)),
   "arc-length clock: ds/dtheta -> K at the top, omega_s = 1/K; g_2 = 1/(2 z_c) -- K^2 cancels: RESTATEMENT of the horn coefficient, not a tie")

dthdz = kappa*z/(1 - z**2)
ck("RX-E7", simplify(cancel(diff(dthdz, z))) != 0
        and zero(dthdz.subs(z, zc) - 2*sqrt(3)*kappa)
        and zero((2*sqrt(3)*kappa)**2*K2/(2*zc) - 12*g1),
   "per-height clock: d theta/dz non-constant (no rigid clock exists); at z_c it is 2sqrt3 kappa, g = 12 g_1 -- still kappa^2 * algebraic")

# ================================================================ RX-G guards
print("== RX-G: certified-interval guards (counted separately) ==")
from mpmath import iv, mp

iv.dps = 60
iphi  = (1 + iv.sqrt(5))/2
ilnph = iv.log(iphi)
ikap  = iv.pi/(2*ilnph)
igap  = iphi**-4
iK2   = 1 - igap
iK    = iv.sqrt(iK2)
izc   = iv.sqrt(3)/2
ig1   = ikap**2*iK2/(2*izc)

def sep(a, b):
    """Certified disjointness of two 60-dps intervals: a rigorous proof of !=."""
    return bool(a.b < b.a) or bool(b.b < a.a)

gk("RX-G0", bool(iv.pi.a > 0),
   "interval-certified: pi > 0 (closes the p=0 branch of T1 in RX-E4)")

gk("RX-G1", (not sep(ig1, ig1)) and (not sep(iv.pi*iK/2, iv.pi*iv.sqrt(1 - iphi**-4)/2)),
   "falsifier: the separator does NOT separate equal quantities (two constructions of pi K/2)")

ledger = [
    ("1",                 iv.mpf(1)),
    ("K",                 iK),
    ("K^2",               iK2),
    ("phi",               iphi),
    ("tau",               1/iphi),
    ("gap",               igap),
    ("5",                 iv.mpf(5)),
    ("sqrt5",             iv.sqrt(5)),
    ("z_c",               izc),
    ("sqrt3",             iv.sqrt(3)),
    ("1/(2 z_c)",         1/(2*izc)),
    ("phi^2",             iphi**2),
    ("2 phi^2",           2*iphi**2),
    ("21/4",              iv.mpf(21)/4),
    ("ln phi",            ilnph),
    ("1/ln phi",          1/ilnph),
    ("pi",                iv.pi),
    ("pi^2",              iv.pi**2),
    ("2 pi",              2*iv.pi),
    ("pi phi",            iv.pi*iphi),
    ("kappa",             ikap),
    ("kappa^2",           ikap**2),
    ("kappa^2/2",         ikap**2/2),
    ("pi kappa/2",        iv.pi*ikap/2),
    ("kappa K",           ikap*iK),
    ("pi K/2",            iv.pi*iK/2),
    ("K ln phi/(4 pi)",   iK*ilnph/(4*iv.pi)),
    ("kappa + pi/2",      ikap + iv.pi/2),
]
all_sep = all(sep(ig1, val) for _, val in ledger)
gk("RX-G2", all_sep,
   "NO-TIE certified: g_1 = kappa^2 K^2/(2 z_c) differs from EVERY listed ledger constant (disjoint 60-dps intervals)")

# L_res exclusion via the corpus-forced bracket (NX-S7): L_res <= kappa + pi/2 < g_1
gk("RX-G3", bool(ig1.a > (ikap + iv.pi/2).b),
   "L_res exclusion: g_1 > kappa + pi/2 >= L_res (corpus bound NX-S7), hence g_1 != L_res")

# mismatch table (display only)
mp.dps = 35
fphi = (1 + mp.sqrt(5))/2
flnp = mp.log(fphi)
fkap = mp.pi/(2*flnp)
fK2  = 1 - fphi**-4
fzc  = mp.sqrt(3)/2
fg1  = fkap**2*fK2/(2*fzc)
print("   RAD-4 mismatch table (30 sig figs, display only):")
print("   g_1 = kappa^2 K^2 / (2 z_c) =", mp.nstr(fg1, 30))
fled = [("kappa", fkap), ("kappa^2", fkap**2), ("kappa^2/2", fkap**2/2),
        ("pi kappa/2", mp.pi*fkap/2), ("pi", mp.pi), ("pi^2", mp.pi**2),
        ("2 phi^2", 2*fphi**2), ("21/4", mp.mpf(21)/4),
        ("pi K/2", mp.pi*mp.sqrt(fK2)/2), ("kappa + pi/2", fkap + mp.pi/2),
        ("ln phi", flnp), ("1/(2 z_c)", 1/(2*fzc))]
for name, val in fled:
    print("     g_1 - (%s) = %s" % (name, mp.nstr(fg1 - val, 20)))
print("   g_2 (arc-length clock 1/K) = 1/(2 z_c) =", mp.nstr(1/(2*fzc), 30), " [RESTATEMENT]")

# ================================================================ summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('RX-G')])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
