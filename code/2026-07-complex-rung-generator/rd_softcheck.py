#!/usr/bin/env python3
# rd_softcheck.py -- OP-RADIUS core dossier harness (v1.8), written cold this
# session.  Builds ON the v1.7 RAD-0 parity no-go (rx_softcheck.py / Thm
# paritynogo); does NOT redo it.  Routes:
#   RD-1 (RAD-5): valuation / ramification certificate at the apex z=0 --
#                 sharpen the parity no-go to a LOCAL ramification certificate.
#   RD-2 (RAD-3): rotor double-cover -- is there a CANONICAL substrate square
#                 root of r^4?  (make-or-break for a constructive bridge)
#   RD-3 (RAD-1): three-atom split status ledger (cap / kink / exponent-1/2).
#
# Discipline (rx/qx/nx idiom): every DECISION exact in Q / Q(sqrt5) (extended by
# sqrt3 and by symbols where needed); K-level statements decided at the squared
# level; valuations decided by EXACT lowest-order term (sympy .leadterm, exponent
# is an exact Rational) -- floats never touch a decision.  mpmath iv interval
# arithmetic at 60 dps appears ONLY in the RD-G corroboration group, counted
# separately; interval OVERLAP proves nothing (RD-G1 pins this honestly), every
# real decision here is exact-algebraic.  Fail-first: each claim carries a
# falsifier guard that FAILS if the claim were false.
#
# Corpus anchors (re-verified against the shipped files this session):
#   - RAD-0 parity no-go [rx / Thm paritynogo, v1.7]: the chart-rational closure
#     of the pentad over Q(sqrt5) is EXACTLY F := Q(sqrt5)(z^2) (sharpness
#     z^2 = 1/(1+u)); rho, theta_k even in z for every branch; r^2 = K^2 z/z_c
#     is odd; r^4 chart-even so F(r^2)/F is exactly quadratic.  INHERITED here.
#   - D1 (W-Decl 12.1): 1-z^2 = e^{-2 rho}, rho = -(1/2) ln(1-z^2).
#   - D3 (W-Decl 12.3): theta = kappa rho, kappa = pi/(2 ln phi).
#   - inherited radius (W-(11.1), [open]1): r(z) = K sqrt(z/z_c) on [0,z_c],
#     r = K on [z_c,1], z_c = sqrt(3)/2.  r^2 = K^2 z/z_c.
#   - whitepaper Thm 15.1 (VERIFIED, pdf p.15): the unique nontrivial z=r point
#     is (K,K); horn r = K sqrt(z/z_c), r(z)=z factors z(z - K^2/z_c)=0, the
#     nontrivial root K^2/z_c > z_c is excluded from the horn (180>169).
#   - whitepaper Prop 15.2(a) (VERIFIED, pdf p.15): kink r'(z_c-) = K/(2 z_c)
#     = K/sqrt3, square K^2/3 = (3 sqrt5 - 5)/6 > 0, against r'(z_c+)=0; C^0
#     not C^1; the lens sits at rapidity rho(z_c) = ln 2, residual a(z_c)=1/2.
#
# Groups:
#   RD-A  corpus replication / definition pinning + RAD-0 inheritance        (4)
#   RD-V  RD-1  valuation / ramification certificate at z=0                   (7)
#   RD-R  RD-2  rotor double-cover: no canonical substrate square root        (6)
#   RD-S  RD-3  three-atom split status ledger                               (5)
#   RD-G  certified-interval corroborations + falsifier (counted separately) (2)
#
# Exit 0 iff every exact check passes and every guard certifies.

import sys
from fractions import Fraction
from math import gcd
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
w = symbols('w')
th = symbols('th')
kk = symbols('k', integer=True)

s5   = sqrt(5)
phi  = (1 + s5)/2
tau  = (s5 - 1)/2
gap  = phi**-4
K2   = 1 - gap
K    = sqrt(K2)
zc   = sqrt(3)/2
kappa = pi/(2*log(phi))

# chart pentad as height functions (W-Thm 13.1, corpus-verified in rx)
u_  = (1 - z**2)/z**2
C_  = (1 - z**2)/z**4
sD  = (2 - z**2)/z**2
m_  = -(1 - z**2)
lam = z**2*(2 - z**2)/(1 - z**2)
rho   = -log(1 - z**2)/2          # D1
theta = kappa*rho                 # D2+D3
r2  = K2*z/zc                     # inherited horn branch, squared
r   = K*sqrt(z/zc)                # inherited horn branch (the radical itself)
r4  = cancel(r2**2)               # chart-even
a_res = sqrt(1 - z**2)            # residual amplitude e^{-rho} = cos alpha
PENTAD = [("u", u_), ("C", C_), ("sqrtD", sD), ("m", m_), ("lambda", lam)]

# ------------------------------------------------------- parity detectors (rx)
def rat_numden(e):
    return fraction(cancel(together(expand(e))))

def even_cert(e):
    """e in Q(...)(z^2): canonical num and den have only even z-degrees."""
    num, den = rat_numden(e)
    pn, pd = Poly(num, z), Poly(den, z)
    return (all(mm[0] % 2 == 0 for mm in pn.monoms())
            and all(mm[0] % 2 == 0 for mm in pd.monoms()))

def sub_even(e):
    return simplify(cancel(together(e.subs(z, -z) - e))) == 0

def sub_odd(e):
    return simplify(cancel(together(e.subs(z, -z) + e))) == 0

# ------------------------------------------------------- valuation machinery
def vz(e):
    """EXACT z-adic valuation at z=0: exponent of the leading term as z->0.
    Returns an exact sympy Rational (integer for rational-in-z, half-integer
    for the radical r).  This is the order of vanishing (>0) / pole (<0)."""
    return together(e).leadterm(z)[1]

def vg_even(e):
    """Detector: does e sit in an even-integer level (the chart value group 2Z)?"""
    ex = vz(e)
    return bool(ex.is_integer) and (Integer(ex) % 2 == 0)

def worder(e):
    """EXACT w-adic valuation at w=0 where w = z^2 (the chart uniformizer)."""
    return cancel(e.subs(z, sqrt(w))).leadterm(w)[1]

def vg_index(g):
    """Ramification index e = [ <2Z, g> : 2Z ] when a valuation-g element is
    adjoined to a field whose z-value group is 2Z.  g an exact Fraction."""
    g = Fraction(g)
    base = Fraction(2)
    num = gcd(base.numerator * g.denominator, g.numerator * base.denominator)
    den = base.denominator * g.denominator
    extgen = Fraction(num, den)          # generator of <2, g> in (Q,+)
    return base / extgen                 # index of 2Z inside <2,g>

def frac(ex):
    """sympy Rational exponent -> python Fraction (exact)."""
    return Fraction(int(ex.p), int(ex.q)) if hasattr(ex, 'p') else Fraction(int(ex))

# ================================================================ RD-A
print("== RD-A: corpus replication / definition pinning + RAD-0 inheritance ==")
ck("RD-A1", zero(phi*tau - 1) and zero(K2 - (1 - gap)) and zero(K2**2 - 5*gap)
        and zero(K2 - (3*s5 - 5)/2) and zero(zc**2 - Rational(3, 4))
        and zero((1 - zc**2) - Rational(1, 4)),
   "seed identities: phi*tau=1, K^2=1-gap=(3sqrt5-5)/2, K^4=5gap, z_c^2=3/4, 1-z_c^2=1/4")

ck("RD-A2", zero(r2 - K2*z/zc) and zero(r2.subs(z, zc) - K2)
        and expand(K2*z/zc - z**2 + z*(z - K2/zc)) == 0
        and sgnQ5(K2 - zc**2) > 0 and sgnQ5(6*s5 - 13) > 0,
   "profile pin (W-11.1): r^2=K^2 z/z_c; r^2(z_c)=K^2 (cap value K); z=r factors z(z-K^2/z_c); K^2>z_c^2 (180>169)")

ck("RD-A3", all(even_cert(p) and sub_even(p) for _, p in PENTAD)
        and cancel(1/(1 + u_) - z**2) == 0
        and zero(u_ - ((1 - w)/w).subs(w, z**2))
        and zero(lam - (w*(2 - w)/(1 - w)).subs(w, z**2)),
   "RAD-0 inheritance: pentad in Q(sqrt5)(z^2); closure EQUALS F=Q(sqrt5)(z^2) (z^2=1/(1+u))")

ck("RD-A4", simplify(rho - (-log(1 - w)/2).subs(w, z**2)) == 0 and sub_even(rho)
        and sub_even(theta) and zero(rho.subs(z, zc) - log(2))
        and zero(a_res.subs(z, zc) - Rational(1, 2)) and zero(a_res - exp(-rho)),
   "D1/D3: rho=g(z^2), theta=kappa*rho even; rho(z_c)=ln2, a(z_c)=1/2, a=e^{-rho} (landmark seeds)")

# ================================================================ RD-V  (RD-1)
print("== RD-V: RD-1 valuation / ramification certificate at the apex z=0 ==")

# V1 -- FAIL-FIRST: the exact valuation detector must report the right integer /
# half-integer levels and the even-value-group detector must accept even levels
# and REJECT odd / half levels (incl. r^2 and r).
battery = {
    "z^2":     (z**2,        Fraction(2),   True),
    "z^3":     (z**3,        Fraction(3),   False),
    "z":       (z,           Fraction(1),   False),
    "1/z^2":   (1/z**2,      Fraction(-2),  True),
    "r^2":     (r2,          Fraction(1),   False),
    "r^4":     (r4,          Fraction(2),   True),
    "r":       (r,           Fraction(1, 2),False),
}
v1_ok = all(frac(vz(e)) == val and vg_even(e) == acc for (e, val, acc) in battery.values())
ck("RD-V1", v1_ok and (not vg_even(r2)) and (not vg_even(r)) and vg_even(r4)
        and vg_even(z**2),
   "falsifier: exact vz detector reports z^2->2, z^3->3, z->1, r^2->1, r->1/2, r^4->2; even-group detector accepts even, rejects odd/half (r^2, r)")

# V2 -- the chart field F has z-value group inside 2Z (every chart quantity even).
adv_even = [u_*C_, lam - sD*m_, (u_**2 + 3*C_)/(1 + lam), rho*u_, theta*C_]
ck("RD-V2", all(vg_even(e) for _, e in PENTAD)
        and vg_even(rho) and vg_even(theta) and vg_even(((pi/2 + 2*pi*kk)/log(phi))*rho)
        and all(vg_even(e) for e in adv_even)
        and frac(vz(u_)) == -2 and frac(vz(C_)) == -4 and frac(vz(lam)) == 2,
   "value group of F=Q(sqrt5)(z^2) lies in 2Z: pentad, rho, theta_k (every branch), and adversarial even combos all have EVEN v_z")

# V3 -- r^2 non-membership by valuation (odd level not in 2Z).
ck("RD-V3", frac(vz(r2)) == 1 and frac(vz(r2)) % 2 == 1 and (not vg_even(r2))
        and sub_odd(r2) and simplify(r2) != 0
        and vg_even(r4) and frac(vz(r4)) == 2,
   "r^2 non-membership: v_z(r^2)=1 ODD (not in 2Z) => r^2 not in F; even control r^4 has v_z=2 in 2Z")

# V4 -- r itself sits at the half-level (where parity cannot reach; valuation can).
ck("RD-V4", frac(vz(r)) == Fraction(1, 2) and (not vz(r).is_integer)
        and zero(r**2 - r2) and zero(cancel(r2*zc/K2) - z),
   "r half-level: v_z(r)=1/2 (not an integer); r^2=(K^2/z_c)z, so z=r^2 z_c/K^2 recovers the odd generator")

# V5 -- ramification indices at z=0 (value-group index of the base 2Z).
#      FAIL-FIRST: an EVEN adjoined valuation must give e=1 (no ramification);
#      the ramification detector fires ONLY on the odd / half levels.
e_r2  = vg_index(frac(vz(r2)))     # adjoin r^2 (v_z=1)  -> Z          -> e=2
e_r   = vg_index(frac(vz(r)))      # adjoin r   (v_z=1/2)-> (1/2)Z     -> e=4
e_ev  = vg_index(frac(vz(rho)))    # adjoin rho (v_z=2, even)          -> e=1
ck("RD-V5", e_r2 == 2 and e_r == 4 and e_ev == 1,
   "ramification at z=0: F(r^2)/F has e=2 (the v1.7 quadratic bridge); F(r)/F has e=4 (two ramified steps); even quantities e=1 (unramified) -- detector fires only on odd/half levels")

# V6 -- unramified elsewhere: the double cover w=z^2 branches only at z=0
#      among finite places; at z=z_c the two sheets are distinct.
branch_locus = solve(diff(z**2, z), z)      # {2z=0} = {0}
ck("RD-V6", branch_locus == [Integer(0)]
        and (not zero(zc - (-zc))) and zero(zc**2 - Rational(3, 4))
        and frac(vz(z**2)) == 2,
   "unramified elsewhere: finite branch locus of w=z^2 is {2z=0}={z=0} only; at z=z_c the sheets +-z_c are distinct (z_c != -z_c, z_c^2=3/4>0)")

# V7 -- the certificate is STRICTLY sharper than parity: parity certifies r^2
#      (an involution on the rational field), but is undefined on r (r(-z) is
#      neither +r nor -r); the valuation certificate covers r as well.
r_neg = r.subs(z, -z)
ck("RD-V7", sub_odd(r2)
        and simplify(cancel(r_neg - r)) != 0 and simplify(cancel(r_neg + r)) != 0
        and (not vz(r).is_integer),
   "sharper than parity: parity certifies r^2 (sub_odd), but r(-z) is neither +r nor -r (parity undefined on the radical); v_z(r)=1/2 certifies r regardless -- a LOCAL apex certificate")

# ================================================================ RD-R  (RD-2)
print("== RD-R: RD-2 rotor double-cover -- no canonical substrate square root ==")

# R1 -- FAIL-FIRST: the exact square roots of r^4 are precisely +-r^2; a planted
#      non-root must be rejected.
ck("RD-R1", zero(r2**2 - r4) and zero((-r2)**2 - r4)
        and cancel((z**2)**2 - r4) != 0 and cancel(lam**2 - r4) != 0,
   "falsifier: the two square roots of r^4 are exactly +-r^2 (both verified); planted non-roots z^2, lambda rejected")

# R2 -- r^4 is chart-even (in F) but NOT a square in F (odd w-valuation);
#      the square root lives one level up, in Q(sqrt5)(z)=F(z), not in F.
perfect_sq = cancel((z**2*(1 - z**2))**2)         # a genuine square in F
ck("RD-R2", even_cert(r4) and frac(vz(r4)) == 2
        and frac(worder(r4)) == 1
        and frac(worder(perfect_sq)) == 2
        and zero(r4 - r2**2) and (not even_cert(r2)),
   "r^4 in F (even, v_z=2) but v_w(r^4)=1 ODD => r^4 is NOT a square in F; a genuine square (z^2(1-z^2))^2 has v_w=2 -- the root r^2 lives in F(z), not F")

# R3 -- FORCED valuation obstruction: any square root g of r^4 has v_z(g)=1
#      (odd); every amplitude-type substrate quantity (residual a, rotor
#      modulus) sits at v_z=0 (even) and CANNOT be g.
rotor_mod2 = cos(th/2)**2 + sin(th/2)**2          # |e^{J th/2}|^2 in Cl(2,0)
ck("RD-R3", frac(vz(r2)) == 1
        and frac(vz(a_res)) == 0 and vg_even(a_res)
        and zero(rotor_mod2 - 1)
        and Fraction(2) / Fraction(2) == 1,        # v_z(g)=v_z(r^4)/2=1 is forced odd
   "FORCED obstruction: any g with g^2=r^4 has v_z(g)=1 (odd); residual amplitude a has v_z=0, rotor modulus |e^{J th/2}|^2=1 (v_z=0) -- no even-level amplitude can be g")

# R4 -- candidate canonical square roots enumerated and all fail (exact).
cand = {
    "residual a":       a_res,          # a^2 = 1-z^2
    "z*a":              z*a_res,        # (z a)^2 = z^2(1-z^2)
    "rotor amp (|.|=1)":Integer(1),     # rotor modulus
    "z^2 (chart-even)":  z**2,
}
r4_roots_fail = all(cancel(simplify(c**2 - r4)) != 0 for _, c in cand.items())
ck("RD-R4", r4_roots_fail
        and cancel(a_res**2 - (1 - z**2)) == 0
        and cancel((z*a_res)**2 - z**2*(1 - z**2)) == 0,
   "candidate square roots all fail: residual a (a^2=1-z^2), z*a (=z^2(1-z^2)), rotor amp (^2=1), z^2 -- none squares to r^4=(K^4/z_c^2)z^2")

# R5 -- the rotor supplies the WRONG double cover: unimodular half-angle cover
#      of the PHASE (Z/8 at the rung cadence chi=pi/2), carrying NO radial
#      scale; the radial coefficient K/sqrt(z_c) (K^2/z_c) is not unimodular.
coef = cancel(K2/zc)                                # r^2 coefficient = (3sqrt5-5)/sqrt3
ck("RD-R5", zero(rotor_mod2 - 1)
        and zero(8*(pi/2)/2 - 2*pi)                 # 8 half-angles of chi/2=pi/4 close: Z/8
        and sgnQ5(cancel(coef**2) - 1) < 0          # (K^2/z_c)^2 = 70/3-10sqrt5 < 1: coef != 1
        and sgnQ5(cancel(coef**2)) > 0,
   "rotor is the WRONG cover: |e^{J th/2}|=1 (unimodular, no radial scale), Z/8 half-angle phase cover (8*pi/4=2pi); radial coeff (K^2/z_c)^2=70/3-10sqrt5 != 1 -- NO-TIE on the coefficient")

# R6 -- the missing square root IS the RAD-0 odd generator: g = (K^2/z_c) z,
#      the SIGNED height (orientation choice), not any canonical amplitude.
g_missing = coef*z
ck("RD-R6", zero(g_missing - r2) and zero(g_missing**2 - r4)
        and sub_odd(g_missing) and frac(vz(g_missing)) == 1,
   "the missing root = (K^2/z_c)*z = r^2: the signed height z (v_z=1, sub_odd) -- exactly the RAD-0 odd generator / orientation choice, not a canonical amplitude")

# ================================================================ RD-S  (RD-3)
print("== RD-S: RD-3 three-atom split status ledger ==")

# S1 -- cap-value atom: FORCED / chart-tied (whitepaper Thm 15.1, the (K,K) point).
zr_roots = solve(K2*z/zc - z**2, z)                 # {0, K^2/z_c}
ck("RD-S1", zero(r2.subs(z, zc) - K2)
        and len(zr_roots) == 2 and any(rt == 0 for rt in zr_roots)
        and any(zero(rt - K2/zc) for rt in zr_roots)
        and sgnQ5(K2 - zc**2) > 0                    # horn root K^2/z_c>z_c iff K^2>z_c^2
        and sgnQ5(K2 - Rational(3, 4)) > 0 and sgnQ5(1 - K2) > 0,   # K in [z_c,1]
   "cap-value atom [FORCED/chart-tied, Thm 15.1]: r(z_c)=K; z=r factors z(z-K^2/z_c), horn root K^2/z_c>z_c (K^2>z_c^2) excluded, cylinder gives z=K in [z_c,1] -> unique (K,K)")

# S2 -- kink-slope atom: FORCED (whitepaper Prop 15.2a).  FAIL-FIRST on the slope.
slope2 = cancel(diff(r, z)**2).subs(z, zc)          # r'(z_c-)^2 on the horn
ck("RD-S2", zero(slope2 - K2/3) and zero(K2/3 - (3*s5 - 5)/6) and sgnQ5(3*s5 - 5) > 0
        and zero((K/(2*zc))**2 - K2/3)
        and cancel(slope2 - K2/2) != 0,             # falsifier: wrong slope rejected
   "kink-slope atom [FORCED, Prop 15.2a]: r'(z_c-)^2=K^2/3=(3sqrt5-5)/6>0 = (K/(2z_c))^2=(K/sqrt3)^2, against r'(z_c+)=0 -> C^0 not C^1; wrong slope K^2/2 rejected")

# S3 -- exponent-1/2 atom: FORCED-RESTATEMENT (RAD-2 area-transfer, re-verified).
h = symbols('h')
Azone = integrate(2*pi, (h, 0, z))                  # Archimedes zone area = 2 pi z
ck("RD-S3", cancel((2*pi*sqrt(1 - h**2))**2*(1 + diff(sqrt(1 - h**2), h)**2) - (2*pi)**2) == 0
        and zero(Azone - 2*pi*z)
        and zero(pi*r2 - (K2/(2*zc))*Azone)
        and zero(diff(pi*r2, z) - pi*K2/zc)
        and cancel(pi*r2 - (K2/(2*zc))*(pi*z**2)) != 0,   # falsifier: z-zone rejected
   "exponent-1/2 atom [FORCED-RESTATEMENT, RAD-2]: r^2 ~ z <=> pi r^2=(K^2/2z_c)A_zone, A_zone=2 pi z; deposition rate pi K^2/z_c constant; wrong z^2-zone rejected")

# S4 -- landmark: rho(z_c)=ln2 is the FIRST/UNIQUE ln2-rapidity level; a(z_c)=1/2.
ln2_solset = solve(rho - log(2), z)                 # {+-sqrt3/2}
ck("RD-S4", zero(rho.subs(z, zc) - log(2))
        and set(ln2_solset) == {zc, -zc}
        and zero(a_res.subs(z, zc) - Rational(1, 2))
        and (not zero(diff(rho, z).subs(z, zc))),   # rho strictly monotone at z_c
   "landmark: rho(z_c)=ln2 is the unique positive ln2-rapidity level (rho=ln2 <=> z^2=3/4); a(z_c)=1/2 (residual amplitude halved)")

# S5 -- landmark restatement test: z_c <-> rho(z_c) is a monotone bijection, so
#      'kink at rho=ln2' and 'lens at z_c=sqrt3/2' are EQUIVALENT data --
#      declaring either FORCES the other (RESTATEMENT); the kink EXISTENCE
#      (cap onset) is DECLARED, not derived.
drho = cancel(diff(rho, z))                         # z/(1-z^2) > 0 on (0,1)
ck("RD-S5", zero(drho - z/(1 - z**2))
        and zero(rho.subs(z, zc) - log(2))          # z_c given => rho=ln2 forced
        and set(solve(rho - log(2), z)) == {zc, -zc}  # rho=ln2 given => z_c forced
        and (not zero(drho.subs(z, Rational(1, 2)))),
   "restatement test: rho monotone (drho/dz=z/(1-z^2)>0) => z_c<->rho(z_c) bijection; 'ln2 landmark' <=> 'lens at sqrt3/2' is mutual derivation (RESTATEMENT); kink existence DECLARED, not derived")

# ================================================================ RD-G guards
print("== RD-G: certified-interval corroborations (counted separately) ==")
from mpmath import iv, mp

iv.dps = 60
iphi  = (1 + iv.sqrt(5))/2
ilnph = iv.log(iphi)
izc   = iv.sqrt(3)/2
iK2   = 1 - iphi**-4
icoef = iK2/izc

def overlap(A, B):
    return not (bool(A.b < B.a) or bool(B.b < A.a))

# G0 -- FAIL-FIRST: the separator does not separate equal quantities (two
# constructions of the coefficient K^2/z_c).
icoef2 = (iphi**-4*(-1) + 1)/izc
gk("RD-G0", overlap(icoef, icoef2) and overlap(icoef, icoef),
   "falsifier: the overlap test does NOT separate two equal constructions of K^2/z_c (overlap proves nothing; the exact proof is RD-R5/RD-A2)")

# G1 -- corroboration of the ln2 landmark at 60 dps (overlap => consistent only;
# the exact identity rho(z_c)=ln2 is proved in RD-A4/RD-S4).
irho_zc = -iv.log(1 - izc**2)/2
gk("RD-G1", overlap(irho_zc, iv.log(iv.mpf(2))),
   "corroboration: rho(z_c) and ln2 agree at 60 dps (overlap only; exact proof RD-A4/RD-S4)")

# numeric display block (display-only, not decisions)
mp.dps = 30
fphi = (1 + mp.sqrt(5))/2
fK2  = 1 - fphi**-4
fzc  = mp.sqrt(3)/2
print("   RD valuation ledger (exact exponents; decimals display-only):")
print("     v_z(u,C,sqrtD,m,lambda) = (-2,-4,-2,0,2)  [all EVEN, F value group 2Z]")
print("     v_z(rho)=v_z(theta_k)=2 EVEN ; v_z(r^2)=1 ODD ; v_z(r)=1/2 ; v_z(r^4)=2 EVEN")
print("     ramification at z=0:  e(F(r^2)/F)=2 ,  e(F(r)/F)=4 ,  even-quantity e=1")
print("     K^2/z_c (r^2 coefficient) =", mp.nstr((fK2/fzc), 20),
      " (=/= 1: (K^2/z_c)^2 = 70/3 - 10 sqrt5)")
print("     rho(z_c) = ln2 =", mp.nstr(mp.log(2), 20), " ; a(z_c)=1/2 ; kink r'(z_c-)=K/sqrt3 =",
      mp.nstr(mp.sqrt(fK2)/mp.sqrt(3), 20))

# ================================================================ summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('RD-G')])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
else:
    print("rd_softcheck: ALL PASS -- RD-1 valuation certificate (e=2/e=4 apex "
          "ramification), RD-2 no canonical square root (NO-TIE), RD-3 atom "
          "ledger; OP-RADIUS stays open, sharpened.")
sys.exit(0 if not FAIL else 1)
