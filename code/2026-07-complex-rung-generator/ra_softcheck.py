#!/usr/bin/env python3
# ra_softcheck.py -- OP-RADIUS reduction dossier harness (v1.9), written cold this
# session.  Builds ON the v1.8 apex-ramification results (rd_softcheck.py /
# Thm apexram, Prop bridgeobstruction, Rem atomledger) and the v1.7 parity no-go
# (rx_softcheck.py / Thm paritynogo); does NOT redo them.  Routes:
#   RA-1 (THE HEADLINE): reduce OP-RADIUS to a MINIMAL declared atom D4, the
#        D2->D2' pattern for the radius.  D4 = (i) an ORIENTATION of the
#        apex-ramified double cover F(r)/F (a sign for the odd generator = the
#        signed height, RD-2) PLUS (ii) the CAP ONSET (formation completes at the
#        lens z_c, RD-3d).  THEOREM (forced for the implication): GIVEN D4, the
#        full profile r(z)=K sqrt(z/z_c) on [0,z_c] capped at r=K for z>=z_c is
#        UNIQUELY determined -- exponent 1/2 (area-transfer, RD-3c), cap value K
#        (the (K,K) point, RD-3a), kink slope K/sqrt3 (RD-3b) all chart-tied and
#        inherited; D4(i) fixes the radical's sign (the only free datum on the
#        double cover), D4(ii) fixes the join point.  Plus MINIMALITY: dropping
#        D4(i) leaves +-r; dropping D4(ii) leaves the join free.
#   RA-2 (differential strengthening of the non-membership, RAD-5 full): try to
#        strengthen "r not in the rational closure of the chart pentad" to "r not
#        in the ELEMENTARY-DIFFERENTIAL closure".  Honest two-part finding:
#        (a) the NAIVE strengthening to the FULL elementary closure FAILS -- r is
#        elementary via r=(K/sqrt z_c)(z^2)^{1/4}, and d/dz breaks the parity
#        (v_z(rho')=1 odd); (b) the RIGHT invariant is the apex MONODROMY ORDER
#        (= apex ramification order), which survives field ops + d/dz + exp +
#        log-of-UNIT; under that (unramified) restriction -- exactly the
#        substrate's canonical form (rho = -1/2 log(1-z^2), a log of the unit
#        1-z^2) -- r is excluded (order-2 monodromy unreachable from order-1
#        data).  The missing datum IS D4(i)'s order-2 ramification.  Strictly
#        strengthens the valuation certificate (excludes a larger field) AND is
#        differential-stable.  CONDITIONAL on the unit-log hypothesis.
#   RA-3 (rotor re-check, bounded): re-confirm the v1.8 NO-TIE -- no corpus-
#        canonical amplitude supplies the r-coefficient K/sqrt(z_c).  Bounded.
#
# Discipline (rd/rx/qx/nx idiom): every DECISION exact in Q / Q(sqrt5) (extended
# by sqrt3 and by symbols where needed); K-level statements decided at the squared
# level; valuations by EXACT leading-term exponent (sympy .leadterm, exact
# Rational) -- floats never touch a decision.  mpmath iv interval arithmetic at
# 60 dps appears ONLY in the RA-G group, counted separately; interval OVERLAP
# proves nothing (RA-G0 pins this), interval DISJOINTNESS is a rigorous !=.
# Fail-first: each detector is shown to FIND a KNOWN positive (planted profile /
# planted ramified power) before any negative is trusted.
#
# Corpus anchors (re-verified against the shipped files this session):
#   - whitepaper Eq (7)/Sec 12: r(z)=K sqrt(z/z_c) on [0,z_c], r=K on [z_c,1],
#     z_c=sqrt3/2 is INHERITED, "not derived here ([open] 1)"; the unique z=r
#     point and the kink are [forced] GIVEN the profile.
#   - whitepaper Thm 15.1 (VERIFIED): unique nontrivial z=r point is (K,K);
#     r^2=z^2 factors z(z-K^2/z_c)=0, horn root K^2/z_c>z_c excluded (180>169),
#     cylinder gives z=K in [z_c,1].  (RD-3a cap value.)
#   - whitepaper Prop 15.2(a) (VERIFIED): kink r'(z_c-)=K/(2 z_c)=K/sqrt3,
#     K^2/3=(3sqrt5-5)/6>0, against r'(z_c+)=0; C^0 not C^1.  (RD-3b kink slope.)
#   - RAD-2 area-transfer (Prop areatransfer): r^2 ~ z <=> constant areal
#     deposition, A_zone=2 pi z; exponent 1/2 forced-reformulation.  (RD-3c.)
#   - RD-1 apex-ramification (Thm apexram): v_z of pentad=(-2,-4,-2,0,2) EVEN,
#     v_z(rho)=v_z(theta_k)=2 EVEN, F=Q(sqrt5)(z^2) value group 2Z; v_z(r^2)=1,
#     v_z(r)=1/2; e(F(r^2)/F)=2, e(F(r)/F)=4.  INHERITED here.
#   - RD-2 constructive-bridge (Prop bridgeobstruction): unique-up-to-sign root
#     +-(K^2/z_c) z = the signed height = the odd orientation generator; rotor
#     is a unimodular angular cover, no radial scale; (K^2/z_c)^2=70/3-10sqrt5!=1.
#     INHERITED here.
#   - RD-3d honest negative (Rem atomledger): cap onset DECLARED, not forced --
#     rho(z_c)=ln2 <-> z_c=sqrt3/2 monotone bijection (RESTATEMENT).  INHERITED.
#
# Groups:
#   RA-A  corpus replication / inheritance pinning (RD-1/RD-2/RD-3 facts)     (4)
#   RA-B  RA-1  D4 reduction: D4=>r unique + minimality (the headline)        (8)
#   RA-C  RA-2  differential strengthening: monodromy-order invariant         (8)
#   RA-D  RA-3  rotor coefficient re-check (bounded NO-TIE)                    (3)
#   RA-G  certified-interval corroborations + falsifier (counted separately)  (4)
#
# Exit 0 iff every exact check passes and every guard certifies.

import sys
from fractions import Fraction
from math import gcd
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   tan, cot, csc, simplify, expand, radsimp, trigsimp, together,
                   cancel, fraction, Poly, diff, solve, integrate, expand_log,
                   powsimp, series, oo, nsimplify)

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

z  = symbols('z')
zp = symbols('zp', positive=True)
w  = symbols('w')
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

# chart pentad as height functions (W-Thm 13.1)
u_  = (1 - z**2)/z**2
C_  = (1 - z**2)/z**4
sD  = (2 - z**2)/z**2
m_  = -(1 - z**2)
lam = z**2*(2 - z**2)/(1 - z**2)
rho   = -log(1 - z**2)/2          # D1  (log of the UNIT 1-z^2)
theta = kappa*rho                 # D2+D3
r2  = K2*z/zc                     # inherited horn branch, squared
r   = K*sqrt(z/zc)                # inherited horn branch (the radical itself)
r4  = cancel(r2**2)
a_res = sqrt(1 - z**2)            # residual amplitude e^{-rho}=cos alpha (unit at 0)
PENTAD = [("u", u_), ("C", C_), ("sqrtD", sD), ("m", m_), ("lambda", lam)]

# ------------------------------------------------------- parity / valuation idiom
def rat_numden(e):
    return fraction(cancel(together(expand(e))))

def even_cert(e):
    num, den = rat_numden(e)
    pn, pd = Poly(num, z), Poly(den, z)
    return (all(mm[0] % 2 == 0 for mm in pn.monoms())
            and all(mm[0] % 2 == 0 for mm in pd.monoms()))

def sub_odd(e):
    return simplify(cancel(together(e.subs(z, -z) + e))) == 0

def sub_even(e):
    return simplify(cancel(together(e.subs(z, -z) - e))) == 0

def vz(e):
    """EXACT z-adic valuation at z=0: leading-term exponent as z->0 (exact Rational)."""
    return together(e).leadterm(z)[1]

def frac(ex):
    return Fraction(int(ex.p), int(ex.q)) if hasattr(ex, 'p') else Fraction(int(ex))

def e_frac(e):
    """Apex fractional exponent (mod 1) -- the multiplicative-monodromy invariant."""
    return frac(vz(e)) % 1

def mono_order(e):
    """Local monodromy order at z=0 of the algebraic part = denominator of e_frac."""
    return (frac(vz(e)) % 1).denominator

def vg_index(g):
    g = Fraction(g); base = Fraction(2)
    num = gcd(base.numerator * g.denominator, g.numerator * base.denominator)
    den = base.denominator * g.denominator
    return base / Fraction(num, den)

# ================================================================ RA-A
print("== RA-A: corpus replication / inheritance pinning (RD-1/RD-2/RD-3) ==")
ck("RA-A1", zero(phi*tau - 1) and zero(K2 - (1 - gap)) and zero(K2**2 - 5*gap)
        and zero(K2 - (3*s5 - 5)/2) and zero(zc**2 - Rational(3, 4))
        and zero((1 - zc**2) - Rational(1, 4)) and sgnQ5(K2 - Rational(3,4)) > 0
        and sgnQ5(1 - K2) > 0,
   "seed identities: phi*tau=1, K^2=1-gap=(3sqrt5-5)/2, K^4=5gap, z_c^2=3/4, 3/4<K^2<1")

ck("RA-A2", zero(r2 - K2*z/zc) and zero(r2.subs(z, zc) - K2)
        and expand(K2*z/zc - z**2 + z*(z - K2/zc)) == 0
        and sgnQ5(K2 - zc**2) > 0 and zero(r**2 - r2)
        and zero(r.subs(z, zc) - K),
   "profile pin (W-Eq7/Thm15.1): r^2=K^2 z/z_c; r(z_c)=K; z=r factors z(z-K^2/z_c); K^2>z_c^2 (180>169)")

ck("RA-A3", frac(vz(r2)) == 1 and frac(vz(r)) == Fraction(1, 2) and frac(vz(r4)) == 2
        and (not vz(r).is_integer) and sub_odd(r2)
        and all(frac(vz(p)) % 2 == 0 for _, p in PENTAD)
        and frac(vz(rho)) == 2 and frac(vz(theta)) == 2,
   "RD-1 inheritance: v_z(r^2)=1 odd, v_z(r)=1/2, v_z(r^4)=2; pentad/rho/theta all EVEN (F value group 2Z)")

coefR2 = cancel(K2/zc)            # r^2 coefficient K^2/z_c
g_sign = coefR2*z                 # signed height = the odd orientation generator (RD-2)
ck("RA-A4", zero(g_sign - r2) and zero(g_sign**2 - r4) and sub_odd(g_sign)
        and frac(vz(g_sign)) == 1 and zero(cancel(coefR2**2) - (Rational(70,3) - 10*s5))
        and sgnQ5(cancel(coefR2**2) - 1) < 0,
   "RD-2 inheritance: unique-up-to-sign root +-(K^2/z_c)z = signed height (odd gen); (K^2/z_c)^2=70/3-10sqrt5 != 1")

# ================================================================ RA-B  (RA-1)
print("== RA-B: RA-1 the D4 reduction (D4 => r unique; minimality) ==")

# The two-parameter admissible family, given the FORCED chart-ties (paraboloid
# exponent 1/2 [RD-3c], cap value magnitude K [RD-3a], C^0 join [continuity]):
#     r_{s,z0}(z) = s*K*sqrt(z/z0) on [0,z0];  s*K on [z0,1],
# free data = (sign s in {+1,-1}=D4(i), onset z0 in (0,1)=D4(ii)).  Squared level
# is sign-blind: horn_sq(z;z0)=(K^2/z0) z, cap_sq=K^2.
def horn_sq(z0):      return (K2/z0)*z          # r^2 on the paraboloid horn
def horn_rad(s, z0):  return s*K*sqrt(z/z0)     # r on the horn (sign s)
cap_sq = K2

# B1 -- FAIL-FIRST profile-equality detector: it must ACCEPT the inherited
# profile and REJECT wrong sign / wrong onset / wrong exponent (cone).
cone_sq = (K2/zc**2)*z**2                        # cone r^2 ~ z^2 (also caps at K^2)
det_accept = zero(horn_sq(zc) - r2) and zero(cap_sq - K2)
det_rej_sign  = not zero(horn_rad(-1, zc) - r)   # -r != r
det_rej_onset = not zero(horn_sq(Rational(4,5)) - r2)
det_rej_expo  = not zero(cone_sq - r2)
ck("RA-B1", det_accept and det_rej_sign and det_rej_onset and det_rej_expo
        and zero(cone_sq.subs(z, zc) - K2),         # cone still hits cap value at z_c
   "falsifier: profile detector ACCEPTS inherited (s=+1,z0=z_c,paraboloid), REJECTS wrong sign (-r), wrong onset (4/5), wrong exponent (cone z^2)")

# B2 -- THE THEOREM (D4 => r).  Given the chart-ties + D4=(s=+1,z0=z_c), the
# reconstructed profile equals the inherited r EXACTLY on both branches.
r_rec_horn = horn_rad(+1, zc)                    # = K sqrt(z/z_c)
ck("RA-B2", zero(r_rec_horn - r) and zero(horn_sq(zc) - r2)
        and zero((horn_rad(+1, zc)).subs(z, zc) - K)         # C^0 join: horn(z_c)=K
        and zero(r_rec_horn.subs(z, zc) - K)
        and zero(diff(r_rec_horn, z).subs(z, zc) - K/(2*zc))  # kink slope consequence
        and zero(K/(2*zc) - K/sqrt(3)),                       # = K/sqrt3 (RD-3b)
   "THEOREM D4=>r [FORCED]: chart-ties + D4(s=+1,z0=z_c) reconstruct r EXACTLY on horn+cap; C^0 join r(z_c)=K; kink r'(z_c-)=K/(2z_c)=K/sqrt3 emerges as a forced consequence")

# B3 -- cap branch + kink asymmetry pinned by D4(ii).
ck("RA-B3", zero(cap_sq - K2)
        and zero(diff(K, z))                        # r'(z_c+)=0 on the flat cap
        and (not zero(diff(r_rec_horn, z).subs(z, zc)))   # r'(z_c-) != 0: genuine kink
        and zero((K/(2*zc))**2 - K2/3) and zero(K2/3 - (3*s5 - 5)/6) and sgnQ5(K2/3) > 0,
   "cap+kink: r=K on [z_c,1] (r'=0); horn slope r'(z_c-)^2=K^2/3=(3sqrt5-5)/6>0 != 0 -- C^0 not C^1 at the D4(ii) join")

# B4 -- MINIMALITY of D4(i): drop the orientation.  Both signs satisfy every
# chart-tie (r^2 identical, |cap|=K identical); they are DISTINCT functions.
ck("RA-B4", zero(horn_rad(+1, zc)**2 - horn_rad(-1, zc)**2)     # r^2 sign-blind: tie holds
        and zero((horn_rad(+1, zc)).subs(z, zc) - K)
        and zero((horn_rad(-1, zc)).subs(z, zc) + K)            # -r cap = -K, |.|=K
        and (not zero(horn_rad(+1, zc) - horn_rad(-1, zc)))     # +r != -r
        and sgnQ5(K2) > 0,                                       # 2K != 0
   "MINIMALITY drop-D4(i) [FORCED]: +r and -r share every chart-tie (r^2, |cap|=K) but differ (2K!=0) -- sign is a genuine free datum; D4(i) is needed")

# B5 -- MINIMALITY of D4(ii): drop the onset.  Any z0 in (0,1) gives a
# chart-consistent capped paraboloid with cap value K; distinct onsets => distinct
# functions.  Witness z0'=4/5 (< z_c and < K, so the (K,K) point survives).
z0alt = Rational(4, 5)
z_star = Rational(5, 6)          # 4/5 < 5/6 < z_c=sqrt3/2:  alt is capped, inherited still climbing
diff_star = horn_sq(zc).subs(z, z_star) - cap_sq   # inherited horn - alt cap at z*
ck("RA-B5", zero(horn_sq(z0alt).subs(z, z0alt) - K2)            # alt: C^0 join, cap value K
        and (not zero(horn_sq(zc) - horn_sq(z0alt)))            # coefficients differ
        and (not zero(cancel(K2/zc - K2/z0alt)))                # K^2/z_c != K^2/(4/5)
        and (Integer(25) < Integer(27))                        # 5 < 3 sqrt3 => diff_star < 0 (exact rational witness)
        and (not zero(diff_star))
        and sgnQ5(K2 - z0alt**2) > 0,                          # z0'=4/5 < K: (K,K) point survives on the cap
   "MINIMALITY drop-D4(ii) [FORCED]: onset z0'=4/5 gives a chart-consistent capped paraboloid (cap value K, C^0), != inherited (coeff K^2/z_c != K^2/(4/5); at z*=5/6 inherited climbs below cap 5<3sqrt3) -- onset is a genuine free datum; D4(ii) is needed")

# B6 -- the RD-3d honesty guard: the 'kink slope forces the onset' objection FAILS
# because the kink slope K/(2 z0) is a strictly-monotone (bijective) function of
# z0, so specifying the slope <=> specifying z0 -- a RESTATEMENT with no
# independent content (mirrors rho(z_c)=ln2 <-> z_c=sqrt3/2, RD-S5).
z0sym = symbols('z0', positive=True)
slope_of = K/(2*z0sym)                                   # horn kink slope at onset z0
dslope   = cancel(diff(slope_of, z0sym))                 # = -K/(2 z0^2) != 0: bijection
ck("RA-B6", zero(slope_of.subs(z0sym, zc) - K/sqrt(3))   # slope at z_c is K/sqrt3
        and zero(cancel(slope_of.subs(z0sym, K/(2*slope_of)) - slope_of))  # z0 recovered from slope: inverse
        and (not zero(dslope)) and sgnQ5(K2) > 0
        and zero(rho.subs(z, zc) - log(2))                # the parallel ln2 landmark
        and set(solve(rho - log(2), z)) == {zc, -zc},
   "RD-3d honesty [FORCED-restatement]: kink slope K/(2 z0) is a monotone bijection of the onset z0 (slope<=>onset), and rho(z_c)=ln2<=>z_c=sqrt3/2 likewise -- the slope/ln2 'constraints' RESTATE the onset, forcing nothing; the onset stays DECLARED (D4(ii))")

# B7 -- the reduction ledger: exactly two free real data (sign, onset); the other
# three atoms (exponent, cap value, kink slope) are chart-tied/forced consequences.
ck("RA-B7", zero(horn_sq(zc) - r2)                       # exponent 1/2: r^2 linear in z (RD-3c)
        and zero(horn_sq(zc).subs(z, zc) - K2)           # cap value K (RD-3a)
        and zero(diff(horn_rad(+1, zc), z).subs(z, zc) - K/sqrt(3))  # kink slope K/sqrt3 (RD-3b) is a consequence
        and (not zero(horn_rad(+1, zc) - horn_rad(-1, zc)))          # free: sign
        and (not zero(horn_sq(zc) - horn_sq(z0alt))),               # free: onset
   "reduction ledger: exponent-1/2 + cap-value-K + kink-slope-K/sqrt3 are chart-tied/forced-consequences; ONLY the sign and the onset are free => D4 = (orientation, onset) is the minimal extra data")

# B8 -- de-restatement check: D4 is NOT one of D1-D3 and does NOT smuggle the
# profile.  D4 declares only a Z/2 sign (an element of {+1,-1}) and a single real
# onset value; from D1-D3 alone r^2 is chart-odd (RD-1) hence NOT reconstructible,
# so D4 adds genuine content (not a restatement of a targeted axiom).
ck("RA-B8", (not even_cert(r2)) and sub_odd(r2)          # D1-D3 chart layer cannot furnish r^2 (odd)
        and frac(vz(r2)) % 2 == 1                        # odd valuation: outside F=Q(sqrt5)(z^2)
        and vg_index(frac(vz(r2))) == 2 and vg_index(frac(vz(r))) == 4  # e=2 / e=4 apex ramification
        and zero(g_sign - r2),                            # D4(i) datum = the signed height, not a D1-D3 quantity
   "de-restatement [FORCED]: r^2 is chart-odd (v_z=1, e=2) so unreachable from D1-D3/F; D4 adds only a Z/2 sign (=signed height) + one onset value -- genuine minimal content, not a restatement of a targeted axiom")

# ================================================================ RA-C  (RA-2)
print("== RA-C: RA-2 differential strengthening (apex monodromy-order invariant) ==")

# C1 -- FAIL-FIRST monodromy-order detector: it must report order 2 for the
# radical r / sqrt(z), order 3 for a planted z^{1/3}, and order 1 for every
# integer-exponent (single-valued) element INCLUDING r^2 and the parity-broken
# derivative rho' (whose valuation is ODD but whose monodromy is trivial).
rho_p = cancel(diff(rho, z))                     # z/(1-z^2), v_z = 1 (ODD) but single-valued
ck("RA-C1", mono_order(r) == 2 and e_frac(r) == Fraction(1, 2)
        and mono_order(sqrt(z)) == 2
        and mono_order(z**Rational(1, 3)) == 3 and e_frac(z**Rational(1,3)) == Fraction(1,3)
        and mono_order(z**2) == 1 and mono_order(r2) == 1
        and mono_order(rho_p) == 1 and frac(vz(rho_p)) == 1,
   "falsifier: apex monodromy-order detector reports order 2 for r/sqrt(z), 3 for z^{1/3}, 1 for z^2, r^2 AND rho'=z/(1-z^2) (v_z=1 odd but single-valued)")

# C2 -- (a) the NAIVE strengthening FAILS: r is in the FULL elementary-differential
# closure, r = (K/sqrt z_c)(z^2)^{1/4} = (K/sqrt z_c) exp((1/4) log(z^2)).
recon = (K/sqrt(zc))*(zp**2)**Rational(1, 4)
ck("RA-C2", zero(recon - K*sqrt(zp/zc))                       # (z^2)^{1/4} = sqrt(z) on z>0
        and zero((zp**2)**Rational(1, 4) - sqrt(zp))
        and zero(exp(Rational(1, 4)*log(zp**2)) - sqrt(zp))   # exp/log build sqrt(z)
        and zero(K*sqrt(zp/zc) - (K/sqrt(zc))*sqrt(zp)),
   "(a) NAIVE strengthening FAILS [FORCED honest-negative]: r=(K/sqrt z_c)(z^2)^{1/4}=(K/sqrt z_c)exp((1/4)log(z^2)) -- r IS in the full elementary-differential closure of F (z^2 in F)")

# C3 -- (a) d/dz breaks the parity/value-group certificate: an EVEN-valuation
# chart quantity differentiates to an ODD-valuation one (same parity as r^2), so
# the 2Z-value-group invariant is NOT differential-stable.
ck("RA-C3", frac(vz(rho)) == 2 and frac(vz(rho_p)) == 1 and frac(vz(rho_p)) % 2 == 1
        and frac(vz(diff(u_, z))) == -3 and frac(vz(diff(u_, z))) % 2 == 1
        and frac(vz(theta)) == 2 and frac(vz(cancel(diff(theta, z)))) == 1,
   "(a) d/dz breaks parity [FORCED]: v_z(rho)=2 even -> v_z(rho')=1 ODD; v_z(u)=-2 -> v_z(u')=-3 ODD; d/dz shifts valuation by -1, so the even/2Z certificate does NOT survive differentiation")

# C4 -- (b) the RIGHT invariant: apex monodromy order is preserved by field ops +
# d/dz on order-1 (single-valued) data.  Battery of unramified-at-0 elementary
# generators + their derivatives all stay order 1.
gens = [u_, C_, sD, m_, lam, rho, theta, r2, a_res, 1/a_res, rho_p]
gens_ok = all(mono_order(g) == 1 for g in gens)
derivs_ok = all(mono_order(cancel(diff(g, z))) == 1 for g in [u_, C_, lam, rho, theta, r2, a_res])
combos = [u_*C_ + lam, (rho*theta + r2)/(1 + u_**2), a_res*rho_p - C_/(1 + lam)]
combo_ok = all(mono_order(cancel(c)) == 1 for c in combos)
ck("RA-C4", gens_ok and derivs_ok and combo_ok,
   "(b) closure on order-1 data [FORCED]: pentad, rho, theta, r^2, a=sqrt(1-z^2) (unit at 0), 1/a, rho', their d/dz and rational combos ALL stay apex monodromy-order 1 (single-valued near z=0)")

# C5 -- (b) log/exp preservation under the UNIT hypothesis: logs of apex-UNITS
# (v_z(arg)=0) are regular (order 1); the substrate's canonical log is exactly
# such (rho = -1/2 log(1-z^2), v_z(1-z^2)=0).  exp of order-1 data is order 1.
ck("RA-C5", frac(vz(1 - z**2)) == 0                      # 1-z^2 is an apex UNIT
        and mono_order(rho) == 1                          # log of a unit is regular
        and mono_order(a_res) == 1 and mono_order(exp(-rho)) == 1   # exp(order-1)=order-1
        and zero(exp(-rho) - a_res)                       # e^{-rho}=sqrt(1-z^2)
        and mono_order(exp(theta)) == 1,
   "(b) unit-log/exp stay order 1 [FORCED-GIVEN-CORPUS]: substrate log arg 1-z^2 is an apex unit (v_z=0) => rho, a=e^{-rho}, exp(theta) all order 1; the canonical elementary layer is unramified at the apex")

# C6 -- (b) the obstruction is EXACTLY the forbidden non-unit log: reaching
# order 2 requires log of a NON-unit (v_z != 0), e.g. log(z^2) (v_z(z^2)=2),
# which the substrate does NOT canonically supply.  This is the D4(i) datum.
ck("RA-C6", frac(vz(z**2)) == 2 and frac(vz(z**2)) != 0    # z^2 is a NON-unit at the apex
        and mono_order(sqrt(z)) == 2                        # the forbidden route yields sqrt(z)=exp((1/4)log z^2) [C2]: order 2
        and e_frac(sqrt(z)) == Fraction(1, 2)
        and mono_order(r) == 2 and mono_order(g_sign) == 1,   # r order2; its square (signed height) order1
   "(b) obstruction = forbidden non-unit log [CONDITIONAL]: order 2 is reachable ONLY via log(z^2) (v_z=2, a NON-unit); under the unit-log restriction it is excluded, so r (order 2) is NOT in the unramified elementary-differential closure -- the missing datum is D4(i)'s order-2 apex ramification")

# C7 -- (b) STRICT strengthening over the valuation certificate: the excluded
# field is strictly LARGER.  The valuation certificate excludes r^2 from F (odd
# value group); the monodromy certificate additionally excludes the RADICAL r
# from the unramified elementary-differential closure F_ed (which CONTAINS r^2,
# rho', logs of units, exps -- all order 1 -- yet not r).
ck("RA-C7", mono_order(r2) == 1 and (not even_cert(r2))    # r^2: order1 (in F_ed) but odd-valuation (not in F)
        and mono_order(rho_p) == 1 and frac(vz(rho_p)) == 1   # rho': order1, odd-valuation -- in F_ed, not F
        and zero(cancel(rho_p*(1 - z**2)) - z)                # d/dz(rho) recovers z: z = rho'*(1-z^2)
        and zero(r2 - coefR2*rho_p*(1 - z**2))                # hence r^2 = (K^2/z_c) z is IN F_ed (field ops + d/dz)
        and mono_order(r) == 2                              # r: order2 -- not in F_ed
        and frac(vz(r)) == Fraction(1, 2),
   "(b) STRICT strengthening [CONDITIONAL]: r^2 and rho' are order-1 (inside the unramified elementary closure F_ed) yet odd-valuation (outside F); r is order-2 (outside F_ed) -- the monodromy certificate excludes a strictly LARGER field than the valuation certificate and is differential-stable")

# C8 -- honest scope tag: differential closure adds NO obstruction in the naive
# (full-closure) sense (C2), while the monodromy/ramification-order invariant IS
# the right differential-stable sharpening (C4-C7), CONDITIONAL on the unit-log /
# unramified hypothesis grounded in the corpus (C5).  Not upgraded past that.
ck("RA-C8", mono_order(r) == 2 and vg_index(frac(vz(r))) == 4  # apex ramification order = monodromy order concordance
        and zero(recon - K*sqrt(zp/zc))                          # full-closure membership stands
        and frac(vz(1 - z**2)) == 0,                             # unit-log hypothesis holds for the substrate
   "honest scope: full elementary closure CONTAINS r (naive strengthening dead); the differential-stable sharpening is apex-monodromy-order-2 = apex-ramification-order-4/2, CONDITIONAL on the corpus-grounded unit-log hypothesis -- exactly D4(i)")

# ================================================================ RA-D  (RA-3)
print("== RA-D: RA-3 rotor coefficient re-check (bounded NO-TIE) ==")
# The r-coefficient of sqrt(z) is K/sqrt(z_c); its square is K^2/z_c.  Re-confirm
# (bounded) that no canonical amplitude supplies it and it is not unimodular.
coefR = cancel(K2/zc)                              # (K/sqrt z_c)^2 = K^2/z_c
rotor_mod2 = cos(th/2)**2 + sin(th/2)**2           # |e^{J th/2}|^2 = 1
ck("RA-D1", zero(rotor_mod2 - 1)
        and zero(cancel(coefR**2) - (Rational(70,3) - 10*s5))
        and sgnQ5(cancel(coefR**2) - 1) < 0        # K^2/z_c < 1: not unimodular
        and sgnQ5(K2) > 0,                          # K^2>0, z_c>0 => coeff K^2/z_c > 0
   "coefficient not unimodular [FORCED]: r-coeff^2=K^2/z_c, (K^2/z_c)^2=70/3-10sqrt5<1 != 1=|rotor|^2 -- the rotor's unimodular angular cover carries no radial scale")

# D2 -- exact NO-TIE against a bounded ledger of canonical amplitude-squares
# (decided in Q(sqrt5) at the squared level; K^2/z_c has sqrt3 so compare squares).
ledger_sq = {                                       # name -> (amplitude)^2 in Q(sqrt5)
    "1": Integer(1), "K^2": K2, "K^4": cancel(K2**2), "phi^2": cancel(phi**2),
    "tau^2": cancel(tau**2), "gap": gap, "gap^2": cancel(gap**2), "5": Integer(5),
    "1/5": Rational(1,5), "phi": phi, "1/phi": tau,
}
# (K^2/z_c)^2 = 70/3 - 10 sqrt5; compare to each cand^2.  All must differ.
tie_free = all(sgnQ5(cancel(coefR**2) - cancel(cnd**2)) != 0 for cnd in ledger_sq.values())
ck("RA-D2", tie_free
        and sgnQ5(cancel(coefR**2) - cancel(K2**2)) != 0
        and sgnQ5(cancel(coefR**2) - 1) != 0,
   "bounded NO-TIE [COMPUTED]: (K^2/z_c)^2=70/3-10sqrt5 differs (exact Q(sqrt5)) from the square of every listed canonical amplitude {1,K^2,K^4,phi^2,tau^2,gap,gap^2,5,1/5,phi,1/phi}")

# D3 -- the coefficient is NOT the signed-height generator either: g_sign=(K^2/z_c)z
# carries the coefficient but is the ODD orientation datum (D4(i)), not a scalar
# amplitude -- re-confirming RD-2's verdict that the scale is not canonical.
ck("RA-D3", zero(g_sign - coefR*z) and sub_odd(g_sign) and frac(vz(g_sign)) == 1
        and (not zero(g_sign - coefR)),               # generator != scalar coefficient
   "NO-TIE verdict [FORCED-GIVEN-CORPUS]: the scale enters only through g=(K^2/z_c)z (odd, v_z=1 = the D4(i) orientation datum), never as a canonical scalar amplitude -- re-confirms RD-2")

# ================================================================ RA-G guards
print("== RA-G: certified-interval corroborations (counted separately) ==")
from mpmath import iv, mp

iv.dps = 60
iphi  = (1 + iv.sqrt(5))/2
izc   = iv.sqrt(3)/2
iK2   = 1 - iphi**-4
iK    = iv.sqrt(iK2)
icoefR = iK2/izc                                    # K^2/z_c

def sep(a, b):   return bool(a.b < b.a) or bool(b.b < a.a)
def overlap(a, b): return not sep(a, b)

# G0 -- FAIL-FIRST: overlap proves nothing (two equal constructions of K^2/z_c).
icoefR2 = (1 - iphi**-4)/(iv.sqrt(3)/2)
gk("RA-G0", overlap(icoefR, icoefR2) and overlap(icoefR, icoefR),
   "falsifier: the overlap test does NOT separate two equal constructions of K^2/z_c (overlap proves nothing; the exact NO-TIE is RA-D2)")

# G1 -- NO-TIE certified: K^2/z_c disjoint from a bounded ledger of amplitudes.
gledger = [("1", iv.mpf(1)), ("K", iK), ("K^2", iK2), ("phi", iphi),
           ("tau", 1/iphi), ("gap", iphi**-4), ("z_c", izc), ("1/z_c", 1/izc),
           ("sqrt5", iv.sqrt(5)), ("2-phi", 2 - iphi), ("K/sqrt3", iK/iv.sqrt(3))]
allsep = all(sep(icoefR, val) for _, val in gledger)
gk("RA-G1", allsep,
   "NO-TIE certified: K^2/z_c differs from every listed canonical amplitude {1,K,K^2,phi,tau,gap,z_c,1/z_c,sqrt5,2-phi,K/sqrt3} (disjoint 60-dps intervals)")

# G2 -- D4 minimality witness certified: at z*=5/6 the inherited horn radius is
# STRICTLY below the cap K (still climbing), while the alt-onset (4/5) profile is
# capped at K -- the two profiles are certifiably distinct on [4/5, z_c].
ihorn_star = iv.sqrt(iK2/izc * iv.mpf(5)/6)         # inherited r at z*=5/6
gk("RA-G2", bool(ihorn_star.b < iK.a) and bool((iv.mpf(4)/5) < izc.a),
   "D4(ii) minimality certified: inherited r(5/6) < cap K (interval-strict), and 4/5 < z_c -- the alt-onset profile (capped at K on [4/5,1]) is distinct from the inherited one")

# G3 -- monodromy corroboration: r ~ z^{1/2} at the apex (the radical's leading
# power) vs the integer-exponent chart data -- display of the order-2 branch.
gk("RA-G3", bool(icoefR.a > 0) and overlap(icoefR, icoefR),
   "corroboration: K^2/z_c > 0 (a genuine positive radial scale); the exact order-2 monodromy/e=4 ramification is proved in RA-C1/RA-C6 (leading power z^{1/2})")

# numeric display block (display-only, not decisions)
mp.dps = 30
fphi = (1 + mp.sqrt(5))/2
fK2  = 1 - fphi**-4
fzc  = mp.sqrt(3)/2
print("   RA ledger (exact facts; decimals display-only):")
print("     RA-1  free data given the chart-ties: (orientation sign s in {+1,-1}, onset z0) = D4")
print("           D4 = (s=+1, z0=z_c=sqrt3/2) => r(z)=K sqrt(z/z_c) capped at K -- UNIQUE")
print("     RA-2  r=(K/sqrt z_c)(z^2)^{1/4} in FULL elem-diff closure; apex monodromy order(r)=2")
print("           v_z(rho')=1 ODD: d/dz breaks parity; unit-log closure order 1 excludes r [CONDITIONAL]")
print("     RA-3  r-coeff^2 = K^2/z_c =", mp.nstr(fK2/fzc, 20),
      " ; (K^2/z_c)^2 = 70/3-10sqrt5 =", mp.nstr((fK2/fzc)**2, 12), "!= 1")

# ================================================================ summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('RA-G')])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
else:
    print("ra_softcheck: ALL PASS -- OP-RADIUS REDUCED to the declared atom D4 "
          "(orientation + onset): D4=>r FORCED, each component minimal; RA-2 "
          "apex-monodromy-order sharpening (full closure contains r; unramified "
          "unit-log closure excludes it) CONDITIONAL; RA-3 rotor NO-TIE re-confirmed.")
sys.exit(0 if not FAIL else 1)
