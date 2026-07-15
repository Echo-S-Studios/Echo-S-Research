#!/usr/bin/env python3
# ex_softcheck.py -- rail closed-form harness (HANDOFF-v1.0 section 8 dossier;
# written cold for the E1/E2 deliverables, 2026-07-15).
#
# E1: explicit incomplete-elliptic closed form for the above-lens rail length.
#     Corpus anchors verified against complex-rung-generator-v1.6.tex:
#       prop:trichotomy  Delta s_n = int sqrt(K^2 + (dz/dtheta)^2) dtheta,
#                        theta in [n pi/2, (n+1) pi/2]
#       QX-E1            dz/dtheta = (1-z^2)/(kappa z)   [from D1+D3]
#       QX-E9            z^4 + (K^2 kappa^2 - 2) z^2 + 1 = (z^2+a)(z^2+1/a),
#                        a + 1/a = K^2 kappa^2 - 2
#     Closed form certified here (FORCED, symbolic differentiation certificate):
#       G(z) = -z sqrt((z^2+a)/(z^2+1/a)) + sqrt(a) E(theta|m)
#              - a^(-3/2) F(theta|m) + (a+1) a^(-3/2) Pi(1+1/a; theta|m)
#       theta = arctan(z sqrt(a)),  m = 1 - 1/a^2,
#       Delta s_n = G(z_{n+1}) - G(z_n),  z_n = sqrt(1 - tau^{2n}).
#
# E2: non-elementarity certificate ingredients: K kappa > 2 (interval) => a
#     real, a > 1 => quartic kernel has 4 distinct roots (genus-1); the F/E/Pi
#     coefficients are exactly nonzero; Pi characteristic non-degenerate;
#     residues at z = +-1 are -+ K kappa / 2 (nonzero).  Established anchor
#     (cited, not re-proven): Liouville's theorem on integration in finite
#     terms -- Rosenlicht, "Integration in finite terms", Amer. Math. Monthly
#     79 (1972) 963-972; algorithmic form: Trager, "Integration of algebraic
#     functions", PhD thesis, MIT, 1984.  Nonvanishing elliptic terms over a
#     genus-1 kernel with non-degenerate parameters => no elementary
#     antiderivative.  Tag: FORCED-GIVEN-ESTABLISHED.
#
# Discipline (qx_softcheck idiom): every DECISION exact in sympy over
# Q / Q(sqrt5) / Q(a, sqrt(a)); algebraic-number nonzero decisions via
# minimal_polynomial (exact); mpmath iv at 60 dps ONLY for transcendental
# strict inequalities and mp at 50 dps ONLY for numeric corroboration --
# guards counted separately (EX-G).  Floats never touch an exact decision.
# Fail-first: EX-A3, EX-B7, EX-G5 are falsifier guards written before the
# passing certificates; each FAILS if the corresponding claim were false.
# Exit 0 iff every check passes.
#
#   EX-A  integrand derivation in the z-chart (replicates QX-E1, QX-E9)   (4)
#   EX-B  incomplete-elliptic closed form, FORCED certificate + falsifiers(8)
#   EX-C  non-elementarity ingredients (E2)                               (5)
#   EX-G  certified-interval / numeric guards (counted separately)        (5)

import sys
import sympy
import mpmath
from sympy import (symbols, sqrt, Rational, I, pi, log, exp, sin, cos, atan,
                   simplify, expand, radsimp, together, trigsimp, diff,
                   discriminant, minimal_polynomial, expand_log,
                   elliptic_f, elliptic_e, elliptic_pi)

print("ex_softcheck: sympy %s, mpmath %s" % (sympy.__version__, mpmath.__version__))

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

def nonzero_algebraic(e, var):
    """Exact nonzero decision for an algebraic number: minpoly == var iff e == 0."""
    return minimal_polynomial(e, var) != var

x = symbols('x')
z = symbols('z', positive=True)
a = symbols('a', positive=True)          # the reciprocal-factorization root, a>0
w = symbols('w')                          # w = a + 1/a = K^2 kappa^2 - 2 (free)
Ks, kaps = symbols('K_s kappa_s', positive=True)   # generic positive K, kappa
Kk = symbols('K_k', positive=True)        # the product K*kappa, kept symbolic

s5    = sqrt(5)
phi   = (1 + s5)/2
tau   = (s5 - 1)/2
gap   = phi**-4
K2    = 1 - gap
K     = sqrt(K2)
kappa = pi/(2*log(phi))
theta = symbols('theta', positive=True)

# ---------------------------------------------------------------- EX-A
print("== EX-A: rail integrand in the z-chart (from prop:trichotomy + QX-E1/E9) ==")

zfun = sqrt(1 - exp(-2*theta/kappa))
ck("EX-A1", zero(diff(zfun, theta) - (1 - zfun**2)/(kappa*zfun)),
   "dz/dtheta = (1-z^2)/(kappa z) from D1+D3 (replicates QX-E1)")

LHS = sqrt(Ks**2 + ((1 - z**2)/(kaps*z))**2)*kaps*z/(1 - z**2)
RHS = sqrt(Ks**2*kaps**2*z**2 + (1 - z**2)**2)/(1 - z**2)
ck("EX-A2", zero(LHS - RHS) and zero(LHS**2 - RHS**2),
   "change of variables: sqrt(K^2+(dz/dtheta)^2) dtheta/dz = "
   "sqrt(K^2 kap^2 z^2 + (1-z^2)^2)/(1-z^2), decided directly and at the squared level")

fact_ok = expand((z**2 + a)*(z**2 + 1/a)*a - a*(z**4 + (a + 1/a)*z**2 + 1)) == 0
pert = expand((z**2 + a)*(z**2 + 1/a) - (z**4 + (a + 1/a + Rational(1, 1000))*z**2 + 1))
ck("EX-A3", fact_ok and pert == -z**2*Rational(1, 1000) and pert != 0,
   "reciprocal factorization (z^2+a)(z^2+1/a) = z^4+(a+1/a)z^2+1 (QX-E9); "
   "FALSIFIER: perturbed middle coefficient leaves exact residue -z^2/1000 (rejected)")

ck("EX-A4", all(zero(exp(-2*(n*pi/2)/kappa) - tau**(2*n)) for n in range(2, 7))
        and all(zero(zfun.subs(theta, n*pi/2)**2 - (1 - tau**(2*n))) for n in range(2, 7)),
   "endpoint map (D2+D3): t_n = e^{-2 theta_n/kappa} = tau^{2n}, "
   "z_n^2 = 1 - tau^{2n}, n = 2..6 (the QX-G3 grid parametrization)")

# ---------------------------------------------------------------- EX-B
print("== EX-B: incomplete-elliptic closed form (E1) ==")

P  = (z**2 + a)*(z**2 + 1/a)
m  = 1 - 1/a**2                 # elliptic parameter (m = k^2 convention)
nc = 1 + 1/a                    # Pi characteristic
th = atan(z*sqrt(a))
integrand = sqrt(P)/(1 - z**2)

ck("EX-B1", expand((z**2 - 1)*(z**2 + w + 1) + (w + 2) - (z**4 + w*z**2 + 1)) == 0,
   "partial-fraction split: P/(1-z^2) = -(z^2+w+1) + (w+2)/(1-z^2), w = a+1/a")

L1 = elliptic_f(th, m)/sqrt(a)
ck("EX-B2", zero(diff(L1, z) - 1/sqrt(P)),
   "lemma 1 (first kind): d/dz [ a^{-1/2} F(arctan(z sqrt(a)) | 1-1/a^2) ] = P^{-1/2}")

L2 = z*sqrt((z**2 + a)/(z**2 + 1/a)) - sqrt(a)*elliptic_e(th, m)
ck("EX-B3", zero(diff(L2, z) - z**2/sqrt(P)),
   "lemma 2 (second kind): d/dz [ z sqrt((z^2+a)/(z^2+1/a)) - sqrt(a) E(theta|m) ] = z^2 P^{-1/2}")

L3 = (a/(a + 1)*elliptic_f(th, m) + elliptic_pi(nc, th, m)/(a + 1))/sqrt(a)
ck("EX-B4", zero(diff(L3, z) - 1/((1 - z**2)*sqrt(P))),
   "lemma 3 (third kind): d/dz [ a^{-1/2} ( a/(a+1) F + 1/(a+1) Pi(1+1/a; theta|m) ) ] "
   "= ((1-z^2) sqrt(P))^{-1}")

wa = a + 1/a
cF_asm  = -(wa + 1)/sqrt(a) + (wa + 2)/sqrt(a)*a/(a + 1)
cPi_asm = (wa + 2)/(sqrt(a)*(a + 1))
ck("EX-B5", zero(cF_asm + a**Rational(-3, 2)) and zero(cPi_asm - (a + 1)*a**Rational(-3, 2)),
   "coefficient assembly: c_E = sqrt(a), c_F = -a^{-3/2}, c_Pi = (a+1) a^{-3/2} "
   "from G = -L2 - (w+1) L1 + (w+2) L3")

G = (-z*sqrt((z**2 + a)/(z**2 + 1/a))
     + sqrt(a)*elliptic_e(th, m)
     - a**Rational(-3, 2)*elliptic_f(th, m)
     + (a + 1)*a**Rational(-3, 2)*elliptic_pi(nc, th, m))
ck("EX-B6", zero(diff(G, z) - integrand),
   "FORCED certificate: dG/dz - sqrt((z^2+a)(z^2+1/a))/(1-z^2) == 0 exactly, "
   "symbolic in (a, z), a > 0, on 0 < z < 1")

Gbad1 = G - (a + 1)*a**Rational(-3, 2)*elliptic_pi(nc, th, m) \
          + (a + 1)*a**Rational(-3, 2)*elliptic_pi(1 - 1/a, th, m)
Gbad2 = G + 2*a**Rational(-3, 2)*elliptic_f(th, m)
pt = [(a, 4), (z, Rational(1, 2))]
d_bad1 = radsimp(together((diff(Gbad1, z) - integrand).subs(pt)))
d_bad2 = radsimp(together((diff(Gbad2, z) - integrand).subs(pt)))
d_good = radsimp(together((diff(G, z) - integrand).subs(pt)))
ck("EX-B7", nonzero_algebraic(d_bad1, x) and nonzero_algebraic(d_bad2, x)
        and minimal_polynomial(d_good, x) == x,
   "FALSIFIERS (exact, via minimal_polynomial at a=4, z=1/2): wrong Pi characteristic "
   "1-1/a rejected; flipped F sign rejected; true form's defect exactly 0")

ck("EX-B8", zero(1 - nc*sin(th)**2 - (1 - z**2)/(1 + a*z**2))
        and zero(1 - m*sin(th)**2 - (a + z**2)/(a*(1 + a*z**2))),
   "domain validity: 1 - nc sin^2(theta) = (1-z^2)/(1+a z^2) > 0 and "
   "1 - m sin^2(theta) = (a+z^2)/(a(1+a z^2)) > 0 on 0<z<1: FTC applies on [z_n, z_{n+1}]")

# ---------------------------------------------------------------- EX-C
print("== EX-C: non-elementarity ingredients (E2) ==")

ap = (w + sqrt(w**2 - 4))/2
am = (w - sqrt(w**2 - 4))/2
ck("EX-C1", zero(ap*am - 1) and zero(ap + am - w)
        and zero((1 + ap)*(1 + am) - (2 + w)),
   "a = (w+sqrt(w^2-4))/2 with a * (1/a) = 1, a + 1/a = w exactly; (1+a)(1+1/a) = 2+w")

dsc = discriminant(z**4 + w*z**2 + 1, z)
ck("EX-C2", expand(dsc - 16*(w - 2)**2*(w + 2)**2) == 0,
   "disc(z^4+wz^2+1) = 16 (w-2)^2 (w+2)^2: 4 distinct roots iff w != +-2; "
   "w = K^2 kappa^2 - 2 > 2 is EX-G1/G2, so the kernel curve y^2 = P has genus 1")

ck("EX-C3", expand(-a**Rational(-3, 2)*a**Rational(3, 2) + 1) == 0
        and zero(sqrt(a)**2 - a)
        and zero((a + 1)*a**Rational(-3, 2)*a**Rational(3, 2) - (a + 1)),
   "coefficients: c_F * a^{3/2} = -1 (nonzero for ALL a>0, unconditionally); "
   "c_E^2 = a and c_Pi * a^{3/2} = a+1, both nonzero given a > 0 (EX-G2)")

ck("EX-C4", zero((nc - 1) - 1/a) and zero((nc - m) - (a + 1)/a**2)
        and zero((1 - m) - 1/a**2) and zero(nc*(2 - nc) - m),
   "Pi non-degeneracy: nc-1 = 1/a != 0, nc-m = (a+1)/a^2 != 0, 1-m = 1/a^2 != 0 "
   "(nc notin {0,1,m}, m != 1); paired-characteristic tie nc(2-nc) = m")

apk = ap.subs(w, Kk**2 - 2)
amk = am.subs(w, Kk**2 - 2)
res_sq = (1 + apk)*(1 + amk)          # y(+-1)^2
ck("EX-C5", zero(res_sq - Kk**2) and zero(sqrt(res_sq) - Kk)
        and zero((-sqrt((z**2 + a)*(z**2 + 1/a))/(1 + z)).subs(z, 1)
                 + sqrt((1 + a)*(1 + 1/a))/2),
   "residue ledger: y(+-1)^2 = (1+a)(1+1/a) = K^2 kappa^2; residue of ds/dz at "
   "z=+-1 is -+ K kappa/2 != 0 given K kappa > 2 (EX-G1) -- the pole pair is genuine")

# ---------------------------------------------------------------- EX-G guards
print("== EX-G: certified-interval and numeric guards (counted separately) ==")
from mpmath import iv, mp, mpf, quad as mquad

iv.dps = 60
iphi = (1 + iv.sqrt(5))/2
igap = iphi**-4
iK2  = 1 - igap
ikap = iv.pi/(2*iv.log(iphi))
gk("EX-G1", iv.sqrt(iK2)*ikap > iv.mpf(2),
   "interval-certified (60 dps): K kappa > 2 (re-certifies QX-G2; a real, quartic elliptic)")

icc = iK2*ikap**2 - 2
ia  = (icc + iv.sqrt(icc**2 - 4))/2
im_ = 1 - 1/ia**2
inc = 1 + 1/ia
gk("EX-G2", (ia > iv.mpf(1)) and (iv.mpf(0) < im_ < iv.mpf(1))
        and (iv.mpf(1) < inc < iv.mpf(2)) and (icc > iv.mpf(2)),
   "interval-certified (60 dps): a > 1 strictly, m in (0,1), nc in (1,2), w > 2 "
   "(distinct roots; proper modulus; Pi characteristic in the non-degenerate window)")

mp.dps = 50
fphi = (1 + mp.sqrt(5))/2
ftau = 1/fphi
fgap = fphi**-4
fK   = mp.sqrt(1 - fgap)
fkap = mp.pi/(2*mp.log(fphi))
fcc  = (fK*fkap)**2 - 2
fa   = (fcc + mp.sqrt(fcc**2 - 4))/2
fm   = 1 - 1/fa**2
fnc  = 1 + 1/fa

def Gnum(zz, char=None):
    ch = fnc if char is None else char
    tt = mp.atan(zz*mp.sqrt(fa))
    alg = -zz*mp.sqrt((zz**2 + fa)/(zz**2 + 1/fa))
    return (alg + mp.sqrt(fa)*mp.ellipe(tt, fm)
            - fa**mpf(-1.5)*mp.ellipf(tt, fm)
            + (fa + 1)*fa**mpf(-1.5)*mp.ellippi(ch, tt, fm))

zf = lambda t_: mp.sqrt(1 - mp.e**(-2*t_/fkap))
integ = lambda t_: mp.sqrt(fK**2 + ((1 - zf(t_)**2)/(fkap*zf(t_)))**2)

ok_cf, ok_br, tbl = True, True, []
for n in range(2, 6):
    tn  = ftau**(2*n)
    zn  = mp.sqrt(1 - tn)
    zn1 = mp.sqrt(1 - ftau**(2*(n + 1)))
    ds  = mquad(integ, [n*mp.pi/2, (n + 1)*mp.pi/2])
    cf  = Gnum(zn1) - Gnum(zn)
    ok_cf = ok_cf and (abs(ds - cf) < mpf(10)**-40)
    I1n = (mp.log((1 - ftau**2*tn)/(1 - tn)) - ftau*tn)/(4*fK*fkap)
    Fp  = lambda tt: 1/(1 - tt) + 3*mp.log(1 - tt) - 3*(1 - tt) + (1 - tt)**2/2
    I2n = (Fp(tn) - Fp(ftau**2*tn))/(16*fK**3*fkap**3)
    lo, hi = mp.pi*fK/2 + I1n - I2n, mp.pi*fK/2 + I1n
    ok_br = ok_br and (lo - mpf(10)**-40 <= cf <= hi + mpf(10)**-40)
    tbl.append((n, ds, cf))
gk("EX-G3", ok_cf,
   "numeric corroboration (50 dps): closed form G(z_{n+1})-G(z_n) matches quadrature "
   "of Delta s_n to <1e-40 on the QX-G3 grid n=2..5")
gk("EX-G4", ok_br,
   "closed form sits inside the corpus bracket piK/2 + I1 - I2 <= Delta s_n <= piK/2 + I1 "
   "(prop:railbounds), n=2..5")

n2  = 2
z2  = mp.sqrt(1 - ftau**4)
z3  = mp.sqrt(1 - ftau**6)
ds2 = mquad(integ, [n2*mp.pi/2, (n2 + 1)*mp.pi/2])
bad = Gnum(z3, char=1 - 1/fa) - Gnum(z2, char=1 - 1/fa)
gk("EX-G5", abs(ds2 - bad) > mpf(10)**-6,
   "numeric FALSIFIER: wrong Pi characteristic (1-1/a) misses the n=2 quadrature "
   "by more than 1e-6 -- the match in EX-G3 is not tolerance slack")

print("   n |        Delta s_n quadrature (25 sig)      |        closed form (25 sig)")
for n, ds, cf in tbl:
    print("   %d | %s | %s" % (n, mp.nstr(ds, 25), mp.nstr(cf, 25)))
print("   a          =", mp.nstr(fa, 25))
print("   m = k^2    =", mp.nstr(fm, 25))
print("   nc         =", mp.nstr(fnc, 25))

# ---------------------------------------------------------------- summary
print("=" * 64)
n_exact_fail = len([f for f in FAIL if not f[0].startswith('EX-G')])
print("EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
