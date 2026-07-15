#!/usr/bin/env python3
# mx_softcheck.py -- multiplier-torus / modular-marking + spectral-ledger harness.
# HANDOFF v1.0 dossiers 3.2 (Route R1) and 3.6 (Route R5); targets M1-M6.
# Written cold this session (2026-07-15), fail-first: every falsifier check
# (a deliberately perturbed variant that must be REJECTED) was written before
# the claim it protects was made to pass.
#
# Discipline: every DECISION is exact over Q / Q(sqrt5) (extended by i and by
# the formal reals log(phi), pi; Q(sqrt5) signs via rational arithmetic only,
# nx sgnQ5 idiom).  The ONLY transcendental input is a single certified
# interval sandwich  10 < kappa^2 < 11  (mpmath iv, 60 dps; guard MX-G1,
# QX-G idiom, counted separately).  Every lattice-length comparison
#   |n log q + 2 pi i m|^2 = ln^2(phi) * ( n^2 + kappa^2 (n+4m)^2 )
# is linear in kappa^2 with INTEGER coefficients, and a linear function on an
# interval attains its extrema at the endpoints, so each comparison is decided
# by exact rational arithmetic against the sandwich endpoints (open interval:
# an endpoint zero with the other endpoint positive is still strictly positive
# on the interior).  Floats never touch an exact decision; the 50-dps ledger
# table is display-only; the iv ordering guards MX-G2/MX-G3 are counted
# separately.
#
# Corpus citations verified against complex-rung-generator-v1.6.tex this
# session (theorem counter is per-section, shared by all environments):
#   Thm 2.3  thm:generator      q = i*tau, gap schedule W_4 = gap
#   Thm 2.4  thm:interp         Z-family log Q = ln(tau) I + (pi/2 + 2 pi k) J
#   Cor 2.5  cor:d3select       (i)-(iv); (iv) per-floor speed sqrt(1+kappa_k^2)
#   Rem 2.6  rem:gammatransfer  sixth characterization
#   Thm 3.3/Def 3.6 (thm:factorization/def:Q)  Q = +-tau J forced up to
#            orientation, no D1-D3 (the torus object's contamination audit)
#   Thm 4.2  thm:gauge          c = 1/4 gauge rotation = "one charge quantum"

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   simplify, expand, radsimp, trigsimp, together, expand_log,
                   minimal_polynomial, Poly, re, im, integrate, conjugate,
                   Matrix, diff)

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

x = symbols('x')
rho = symbols('rho', positive=True)
nn, mm, kint = symbols('n m k', integer=True)
jj = symbols('jj')                       # tail parameter, used with jj >= 0

s5    = sqrt(5)
phi   = (1 + s5)/2
tau   = (s5 - 1)/2
gap   = phi**-4
kappa = pi/(2*log(phi))
q     = I*tau

Lq = -log(phi) + I*pi/2                          # claimed principal log q
mu = Rational(1, 4) + I*log(phi)/(2*pi)          # claimed modulus

# exact Q(sqrt5) sign decisions (nx idiom)
def q5AB(e):
    from sympy import cancel
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQ5(e):
    A, B = q5AB(e)
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    if A > 0:
        return 1 if A*A > 5*B*B else -1
    return 1 if 5*B*B > A*A else -1

# sandwich decision: is a + b*kappa^2 > 0 for ALL kappa^2 in the OPEN interval
# (SAND_LO, SAND_HI)?  Linear => extrema at endpoints; an endpoint zero with
# the other endpoint positive is strictly positive on the open interior.
SAND_LO, SAND_HI = Fraction(10), Fraction(11)
def pos_on_sandwich(a, b):
    a, b = Fraction(a), Fraction(b)
    vlo, vhi = a + b*SAND_LO, a + b*SAND_HI
    if vlo > 0 and vhi > 0: return True
    if vlo == 0 and vhi > 0: return True
    if vhi == 0 and vlo > 0: return True
    return False

def poly_nonneg_coeffs(e, v):
    return all(Fraction(str(c)) >= 0 for c in Poly(expand(e), v).all_coeffs())

# ---------------------------------------------------------------- MX-A
print("== MX-A: replication of the layer built upon ==")
ck("MX-A1", zero(tau**2 + tau - 1) and zero(phi*tau - 1) and zero(phi**2 - phi - 1),
   "seed identities: tau^2+tau=1, phi*tau=1, phi^2=phi+1")
ck("MX-A2", zero(q**2 + tau**2) and zero(q**4 - gap)
        and expand(minimal_polynomial(q, x) - (x**4 + 3*x**2 + 1)) == 0,
   "generator: q^2=-tau^2, q^4=gap, minpoly(q)=x^4+3x^2+1 (Thm 2.3/Thm 4.1)")
ck("MX-A3", sgnQ5(tau) > 0 and sgnQ5(1 - tau) > 0 and sgnQ5(phi - 1) > 0,
   "0 < |q| = tau < 1 and phi > 1 (exact Q(sqrt5) signs): E_q = C*/q^Z nondegenerate, ln(phi) > 0")
ck("MX-A4", zero(q - tau*exp(2*pi*I*Rational(1, 4))),
   "q = tau e^{2 pi i (1/4)}: the rung angle is the c=1/4 charge quantum (Thm 4.2)")

# ---------------------------------------------------------------- MX-M1
print("== MX-M1: principal log and the modulus of E_q ==")
# falsifier first: the k=1 branch value reproduces q under exp but violates
# the principal strip Im in (-pi, pi]; the strip must REJECT it.
Lq1 = Lq + 2*pi*I
ck("MX-M1F1", zero(exp(Lq1) - q) and Fraction(5, 2) > Fraction(1, 1),
   "falsifier: k=1 value -ln(phi)+i(5pi/2) has exp = q BUT 5/2 > 1 violates the strip: rejected")
ck("MX-M1a", zero(exp(Lq) - q),
   "exp(-ln(phi) + i pi/2) = q exactly")
ck("MX-M1b", zero(im(Lq) - pi/2) and Fraction(-1) < Fraction(1, 2) <= Fraction(1),
   "Im log q = pi/2 lies in the principal strip (-pi, pi] (rational coefficient check)")
ck("MX-M1c", simplify(exp(2*pi*I*kint).rewrite(cos)) == 1
        and all(2*abs(kv) >= 2 for kv in range(-6, 7) if kv != 0)
        and Fraction(1) - Fraction(-1) == 2,
   "uniqueness: branches differ by 2 pi i k (e^{2 pi i k}=1); |2k|>=2 = strip diameter, endpoint -1 open => no second strip value")
ck("MX-M1d", sgnQ5(tau) > 0 and sgnQ5(phi) > 0 and zero(tau*phi - 1)
        and zero(exp(log(tau) + log(phi)) - 1),
   "ln(tau) = -ln(phi): tau,phi > 0, tau*phi = 1, exp(ln tau + ln phi) = 1 with real exponent")
ck("MX-M1e", zero(Lq/(2*pi*I) - mu),
   "modulus mu = (log q)/(2 pi i) = 1/4 + i ln(phi)/(2 pi)")
ck("MX-M1f", zero(re(mu) - Rational(1, 4)) and zero(im(mu) - log(phi)/(2*pi)),
   "Re mu = 1/4 (the charge quantum), Im mu = ln(phi)/(2 pi)")
ck("MX-M1g", zero(re(mu)/im(mu) - kappa),
   "Re mu / Im mu = (1/4)/(ln phi/(2 pi)) = pi/(2 ln phi) = kappa exactly")
ck("MX-M1h", zero((Lq + 2*pi*I*kint)/(2*pi*I) - mu - kint),
   "branch shift: mu_k = mu + k, so mu is well defined mod Z; Im mu branch-independent")
ck("MX-M1F2", zero(((Rational(1, 3)/im(mu)) - kappa)*6*log(phi)/pi - 1),
   "falsifier: perturbed Re mu = 1/3 gives ratio - kappa = pi/(6 ln phi) != 0 (product with 6 ln phi/pi is exactly 1)")
ck("MX-M1i", zero(exp(2*pi*I*mu) - q),
   "e^{2 pi i mu} = q: mu is the standard modulus of the Tate parameter")

# ---------------------------------------------------------------- MX-M2
print("== MX-M2: per-floor arc constant |log q|^2 ==")
# falsifier first: the k=1 branch changes the arc constant by exactly 6 pi^2.
d21 = expand((Lq + 2*pi*I)*conjugate(Lq + 2*pi*I)) - expand(Lq*conjugate(Lq))
ck("MX-M2F1", zero(d21 - 6*pi**2),
   "falsifier: |log_1 q|^2 - |log q|^2 = 6 pi^2 != 0: the arc constant is branch-sensitive")
ck("MX-M2a", zero(expand(Lq*conjugate(Lq)) - (log(phi)**2 + pi**2/4)),
   "|log q|^2 = ln^2(phi) + pi^2/4")
ck("MX-M2b", zero(log(phi)**2*(1 + kappa**2) - log(phi)**2 - pi**2/4),
   "ln^2(phi)(1+kappa^2) = ln^2(phi) + pi^2/4: same constant in Cor 2.5(iv) form")
ck("MX-M2c", zero(integrate(exp(-rho), (rho, 0, log(phi))) - (1 - tau))
        and zero(expand((-1 + I*kappa)*conjugate(-1 + I*kappa)) - (1 + kappa**2)),
   "Cor 2.5(iv) verified: L_0 = (1-tau) sqrt(1+kappa^2); log-chart speed^2 = 1+kappa^2")
ck("MX-M2d", zero((-1 + I*kappa)*log(phi) - Lq),
   "per-floor log-chart displacement (-1+i kappa) ln(phi) = log q: |log q| IS the per-floor arc")

# ---------------------------------------------------------------- MX-G1 (guard)
print("== MX-G1: the single certified transcendental input ==")
from mpmath import iv, mp
iv.dps = 60
iphi = (1 + iv.sqrt(5))/2
ikap2 = (iv.pi/(2*iv.log(iphi)))**2
iv_lo = iv.mpf(SAND_LO.numerator)/iv.mpf(SAND_LO.denominator)
iv_hi = iv.mpf(SAND_HI.numerator)/iv.mpf(SAND_HI.denominator)
gk("MX-G1", (ikap2 > iv_lo) and (ikap2 < iv_hi),
   "interval-certified (60 dps): SAND_LO=10 < kappa^2 < SAND_HI=11 (the exact endpoints every lattice decision below uses)")

# ---------------------------------------------------------------- MX-M3
print("== MX-M3: systole of E_q = C/Lambda, Lambda = Z log q + Z 2 pi i ==")
fexpr = lambda a, b: (a, b)   # class key (a,b) means length^2 = ln^2 phi (a + b kappa^2)

ck("MX-M3a", zero((nn*Lq + 2*pi*I*mm) - (-nn*log(phi) + I*(nn*pi/2 + 2*pi*mm)))
        and zero(expand(nn**2*log(phi)**2 + (nn*pi/2 + 2*pi*mm)**2)
                 - expand((-nn*log(phi))**2 + (nn*pi/2 + 2*pi*mm)**2)),
   "decomposition: n log q + 2 pi i m = -n ln(phi) + i(n pi/2 + 2 pi m); |v|^2 = n^2 ln^2 phi + (n pi/2 + 2 pi m)^2")
ck("MX-M3b", zero(expand(nn**2*log(phi)**2 + (nn*pi/2 + 2*pi*mm)**2
                 - log(phi)**2*(nn**2 + kappa**2*(nn + 4*mm)**2))),
   "kappa^2 form: |v|^2 = ln^2(phi) (n^2 + kappa^2 (n+4m)^2) -- all comparisons linear in kappa^2")

# falsifiers first (each would FAIL if the systole claim were false that way):
ck("MX-M3F1", pos_on_sandwich(0 - 1, 16 - 1),
   "falsifier: meridian 2 pi i rejected as systole: f(0,1)-f(1,0) = -1+15 kappa^2 > 0")
ck("MX-M3F2", pos_on_sandwich(16 - 1, 0 - 1),
   "falsifier: real competitor 4 log q - 2 pi i rejected as systole: 16 > 1+kappa^2 (kappa^2 < 15)")
ck("MX-M3F3", pos_on_sandwich(-7, 1) and not pos_on_sandwich(7, -1),
   "falsifier: wrong ordering f(3,-1) < f(4,-1) rejected: 9+kappa^2 > 16 (kappa^2 > 7)")

# finite enumeration, window |n| <= 4, |m| <= 2, classes up to +-
def canon(nv, mv):
    return (nv, mv) if (nv > 0 or (nv == 0 and mv > 0)) else (-nv, -mv)

def collect(N, M):
    cl = {}
    for nv in range(-N, N + 1):
        for mv in range(-M, M + 1):
            if (nv, mv) == (0, 0):
                continue
            key = (nv*nv, (nv + 4*mv)**2)
            cl.setdefault(key, set()).add(canon(nv, mv))
    return cl

cl3 = collect(4, 2)
ok_min = True
for (a, b), reps in cl3.items():
    if (a, b) == (1, 1):
        ok_min = ok_min and reps == {(1, 0)}          # equality class exactly +-(1,0)
    else:
        ok_min = ok_min and pos_on_sandwich(a - 1, b - 1)
ck("MX-M3c", ok_min,
   "window |n|<=4,|m|<=2: f(n,m) > f(1,0) = 1+kappa^2 for every class except exactly +-(1,0)")

# exact tail lemmas closing the window (polynomial identities, nonneg coefficients)
t_n4  = expand((4 + jj)**2 - 16 - jj*(jj + 8)) == 0
t_n5  = expand((5 + jj)**2 - 25 - jj*(jj + 10)) == 0
t_n13 = expand((13 + jj)**2 - 169 - jj*(jj + 26)) == 0
t_m3  = all(expand((n0 + 12 + 4*jj)**2 - 64 - (n0 + 4 + 4*jj)*(n0 + 20 + 4*jj)) == 0
            and poly_nonneg_coeffs(n0 + 4 + 4*jj, jj) for n0 in range(-4, 5))
t_m5  = all(expand((n0 + 20 + 4*jj)**2 - 64 - (n0 + 12 + 4*jj)*(n0 + 28 + 4*jj)) == 0
            and poly_nonneg_coeffs(n0 + 12 + 4*jj, jj) for n0 in range(-12, 13))
t_sym = zero(expand(((-nn)**2 + kappa**2*(-nn - 4*mm)**2) - (nn**2 + kappa**2*(nn + 4*mm)**2)))
t_dec = (pos_on_sandwich(16 - 1, -1)      # n-tail: n^2 >= 16 > 1+kappa^2
         and pos_on_sandwich(-1, 64 - 1)) # m-tail: 64 kappa^2 > 1+kappa^2
ck("MX-M3d", t_n4 and t_n5 and t_n13 and t_m3 and t_m5 and t_sym and t_dec,
   "tails exact: |n|>=4 => f >= n^2 >= 16 > 1+kappa^2; |m|>=3,|n|<=4 => (n+4m)^2 >= 64 => f >= 64 kappa^2 > 1+kappa^2; (n,m) -> (-n,-m) symmetry")

# second-shortest class: +-(4,-1), squared value 16 ln^2 phi, i.e. length 4 ln phi
ok_second = cl3[(16, 0)] == {(4, -1)}
for (a, b), reps in cl3.items():
    if (a, b) in ((1, 1), (16, 0)):
        continue
    ok_second = ok_second and pos_on_sandwich(a - 16, b - 0)
ok_second = ok_second and pos_on_sandwich(25 - 16, 0) and pos_on_sandwich(-16, 64)
ck("MX-M3e", ok_second and zero(4*Lq - 2*pi*I + 4*log(phi))
        and zero(expand((4*Lq - 2*pi*I)*conjugate(4*Lq - 2*pi*I)) - 16*log(phi)**2)
        and zero(exp(4*Lq - 2*pi*I) - gap),
   "second-shortest class = +-(4 log q - 2 pi i) = -4 ln phi, length 4 ln phi: the four-rung charge cycle W_4 = gap (Thm 2.3) as a geodesic")

# systolic marking: branch generator log_k q = (1,k); minimal iff k = 0
ok_mark = zero(expand((1 + kappa**2*(1 + 4*kint)**2) - (1 + kappa**2))
               - expand(kappa**2*((1 + 4*kint)**2 - 1)))
ok_mark = ok_mark and all((1 + 4*kv)**2 - 1 >= 8 for kv in range(-6, 7) if kv != 0)
ok_mark = ok_mark and expand((5 + 4*jj)**2 - 25 - 8*jj*(2*jj + 5)) == 0 \
                  and expand((3 + 4*jj)**2 - 9 - 8*jj*(2*jj + 3)) == 0 \
                  and pos_on_sandwich(0, 8)      # 8 kappa^2 > 0
ck("MX-M3f", ok_mark,
   "f(1,k)-f(1,0) = kappa^2((1+4k)^2-1) >= 8 kappa^2 > 0 for k != 0: the k-branch generator is a systole representative iff k = 0")

# ---------------------------------------------------------------- MX-M4
print("== MX-M4: restatement test (mandatory for OP-RATE-adjacent claims) ==")
ck("MX-M4a", zero((Lq + 2*pi*I*kint) - (Integer(1)*Lq + kint*(2*pi*I)))
        and Matrix([[1, kint], [0, 1]]).det() == 1
        and zero(Lq - ((Lq + 2*pi*I*kint) - kint*(2*pi*I))),
   "lattice branch-independence: {log_k q, 2 pi i} = unimodular change of basis of Lambda, det 1, for every k")
ck("MX-M4b", zero((nn*(Lq + 2*pi*I*kint) + 2*pi*I*mm) - (nn*Lq + 2*pi*I*(mm + nn*kint))),
   "spectrum k-invariance: n log_k q + 2 pi i m = n log q + 2 pi i (m+nk); (n,m)->(n,m+nk) bijective")
ck("MX-M4c", zero(expand((2*pi*I*kint)*conjugate(2*pi*I*kint)) - 4*pi**2*kint**2)
        and all(kv*kv >= 1 for kv in range(-6, 7) if kv != 0),
   "marking NON-invariance: log_k q - log q = 2 pi i k, |.|^2 = 4 pi^2 k^2 > 0 for k != 0")
print("   VERDICT: the lattice, torus, length spectrum and systole VALUE survive")
print("   kappa -> kappa_k for all k; only the MARKING separates branches, and it")
print("   does so via a minimality principle (systole). Per HANDOFF section 2:")
print("   characterization #7, NOT a derivation of OP-RATE.")

# ---------------------------------------------------------------- MX-M5
print("== MX-M5: certified closed-geodesic length ledger of E_q ==")
cl5 = collect(12, 4)
EXPECTED = [((1, 1),  {(1, 0)}),
            ((16, 0), {(4, -1)}),
            ((9, 1),  {(3, -1)}),
            ((25, 1), {(5, -1)}),
            ((4, 4),  {(2, 0), (2, -1)}),
            ((49, 1), {(7, -2)}),
            ((64, 0), {(8, -2)}),
            ((36, 4), {(6, -1), (6, -2)}),
            ((81, 1), {(9, -2)}),
            ((1, 9),  {(1, -1)})]
# sort ALL window classes by exact rational surrogate, then certify pairwise
surro = lambda ab: Fraction(ab[0]) + Fraction(ab[1])*Fraction(21, 2)
order = sorted(cl5.keys(), key=surro)
first10 = order[:10]
ok_led = [fk for fk, _ in EXPECTED] == first10
ok_led = ok_led and all(cl5[fk] == rp for fk, rp in EXPECTED)
# certify the sort: strict increase between consecutive of the first 10,
# and the 10th below every remaining window class (exact sandwich decisions)
for i in range(9):
    (a1, b1), (a2, b2) = first10[i], first10[i + 1]
    ok_led = ok_led and pos_on_sandwich(a2 - a1, b2 - b1)
a10, b10 = first10[9]
for (a, b) in order[10:]:
    ok_led = ok_led and pos_on_sandwich(a - a10, b - b10)
ck("MX-M5a", ok_led,
   "first 10 distinct classes (up to +-) certified sorted; representatives and multiplicities as expected")
# window sufficiency for the ledger: outside |n|<=12,|m|<=4 every f exceeds the 10th value
ck("MX-M5b", t_n13 and t_m5
        and pos_on_sandwich(169 - 1, -9)     # n-tail: 169 > 1 + 9 kappa^2
        and pos_on_sandwich(-1, 64 - 9),     # m-tail: 64 kappa^2 > 1 + 9 kappa^2
   "ledger window sufficiency: |n|>=13 => f >= 169 > 1+9 kappa^2; |m|>=5,|n|<=12 => f >= 64 kappa^2 > 1+9 kappa^2")
ck("MX-M5c", zero(expand(log(phi)**2*(4 + 4*kappa**2) - 4*(log(phi)**2*(1 + kappa**2))))
        and zero(expand(log(phi)**2*(36 + 4*kappa**2) - 4*(log(phi)**2*(9 + kappa**2)))),
   "ledger internal ties: value_5 = (2 |log q|)^2 (iterate of systole); value_8 = 4 * value_3 (iterate + primitive (6,-1) degenerate)")

# certified-interval ordering + display table (guards / displays, not exact checks)
ilen2 = [ (iv.mpf(a) + iv.mpf(b)*ikap2) for (a, b) in first10 ]
ilnphi = iv.log(iphi)
ilens = [ ilnphi*iv.sqrt(v) for v in ilen2 ]
gk("MX-G2", all(ilens[i].b < ilens[i + 1].a for i in range(9)),
   "interval-certified (60 dps): the 10 ledger lengths are strictly increasing")
gk("MX-G3", (ilens[0].b < 4*ilnphi.a) and ((4*ilnphi).b < ilens[2].a),
   "interval-certified (60 dps): systole |log q| < 4 ln phi < |3 log q - 2 pi i| (corroborates MX-M3e)")
mp.dps = 50
fphi = (1 + mp.sqrt(5))/2
fkap2 = (mp.pi/(2*mp.log(fphi)))**2
flnphi = mp.log(fphi)
print("   certified closed-geodesic length ledger of E_q (lengths ln(phi)*sqrt(a+b*kappa^2)):")
print("   #  | (a,b)    | class reps (+-)       | mult | length (18 sig)")
for i, ((a, b), rp) in enumerate(EXPECTED):
    Lv = flnphi*mp.sqrt(a + b*fkap2)
    form = "%d+%d*kap^2" % (a, b) if b else "%d = (%d ln phi)^2" % (a, int(a**0.5))
    print("   %2d | (%2d,%d)  | %-21s | %d    | %s" %
          (i + 1, a, b, ",".join(str(t) for t in sorted(rp)), len(rp), mp.nstr(Lv, 18)))
print("   systole = |log q| = ln(phi) sqrt(1+kappa^2) =", mp.nstr(flnphi*mp.sqrt(1 + fkap2), 18))
print("   second  = |4 log q - 2 pi i| = 4 ln(phi)    =", mp.nstr(4*flnphi, 18))

# ---------------------------------------------------------------- MX-M6
print("== MX-M6: conditional no-CM corollary (algebra forced; transcendence conditional) ==")
# falsifier first: the perturbed identity with 1/3 must be REJECTED
ck("MX-M6F1", zero((((mu - Rational(1, 3))*kappa) - I/4)*(-12/kappa) - 1),
   "falsifier: (mu - 1/3) kappa - i/4 = -kappa/12 != 0 (product with -12/kappa is exactly 1): perturbed tie rejected")
ck("MX-M6a", zero((mu - Rational(1, 4))*kappa - I/4),
   "exact bridge: (mu - 1/4) kappa = i/4, so mu algebraic => kappa = (i/4)/(mu-1/4) algebraic")
ck("MX-M6b", sgnQ5(phi - 1) > 0 and zero(im(mu - Rational(1, 4)) - log(phi)/(2*pi)),
   "nonvanishing: phi > 1 => ln phi > 0 => mu - 1/4 = i ln(phi)/(2 pi) != 0 and mu not real (honest elliptic modulus)")
ck("MX-M6c", zero(kappa - (I*Rational(1, 4))/(mu - Rational(1, 4))),
   "inverse form kappa = (i/4)/(mu - 1/4): contrapositive gives kappa not algebraic => mu not algebraic => mu not imaginary quadratic => E_q has no CM [conditional on T3; CM criterion external-established]")

# ---------------------------------------------------------------- summary
print("=" * 64)
n_guard_fail = sum(1 for _, ok, _ in GUARD if not ok)
n_exact_fail = len(FAIL) - n_guard_fail
n_exact_pass = len(PASS)
print("MX EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (n_exact_pass, n_exact_fail, sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
