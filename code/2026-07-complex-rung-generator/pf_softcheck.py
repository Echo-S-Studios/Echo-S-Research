#!/usr/bin/env python3
# pf_softcheck.py -- the PERIOD FRONTIER, mapped precisely, and the corpus
# constants' non-relations CERTIFIED as a validated-detector bounded-negative.
# Dossier PF (v1.9).  This harness RESOLVES NOTHING external-open.  It delivers,
# in strict discipline order:
#   PF-3  the FORCED sub-relations ledger -- the exact algebraic identities the
#         frontier SITS ON (decided over Q / Q(sqrt5) / Q(sqrt13), symbolic in k
#         / theta / K where the statement is a family) -- INCLUDING one new
#         forced geometric identity (the residual ellipse's focal parameter is
#         exactly kappa; its eccentricity^2 is exactly the elliptic parameter m).
#   PF-2  the precise CONDITIONAL / EXTERNAL-OPEN tool-gap map -- for each of
#         Gelfond-Schneider, Baker, Schneider/Chudnovsky, Schanuel,
#         Kontsevich-Zagier, EXACTLY what it settles and EXACTLY the gap it
#         leaves.  Sharper than PM-2: it separates Baker's REACH (the Q-linear
#         independence of {1,ln2,ln3,ln phi}, FORCED-GIVEN-ESTABLISHED) from
#         Baker's LIMIT, and gives the Schanuel-CONDITIONAL for (pi,ln phi) and
#         (kappa,ln phi) with the exact Q-linear-independence premise re-encoded.
#   PF-1  a VALIDATED-DETECTOR bounded-negative (COMPUTED): a high-precision
#         (230 dps) PSLQ integer-relation search over the dossier's bases, with
#         LOGGED precision + coefficient-height H + monomial-degree bounds.
#         FAIL-FIRST (mandatory): the detector must FIRST FIND several KNOWN
#         relations -- including a planted relation in a SIZE-6 basis, the same
#         size as the flagship negative -- before any negative is trusted; and
#         it must correctly report NOTHING on two TRUE-NEGATIVE controls (kappa,
#         pi -- known transcendental).  Every negative is COMPUTED corroboration,
#         NEVER a transcendence / algebraic-independence proof.
#
# Written cold this session (2026-07-15) from the dossier + papers; machinery
# re-encoded from scratch (no import of any corpus/companion code).
#
# DISCIPLINE (pm/tx idiom):
#  * every DECISION is exact over Q / Q(sqrt5) / Q(sqrt13) (extended by i where
#    needed) and, for algebraicity, at the minimal_polynomial / resultant level;
#    Q(sqrt D) real signs via rational arithmetic only; floats NEVER touch an
#    exact decision.
#  * the ONLY transcendental inputs are (a) the certified mpmath-iv interval
#    guards (group PF-G, 60 dps) and (b) the PSLQ integer-relation search
#    (group PF-1, 230 dps) -- BOTH counted SEPARATELY from the exact ledger.
#  * fail-first everywhere: PF-3 identities carry perturbation falsifiers;
#    PF-1's detector is validated on 6 KNOWN positives (incl. a size-6 planted
#    relation) and 2 true-negative controls before any open-basis negative.
#
# HONESTY (the deliverable this session): NOTHING here forces L_res
# transcendence or (pi,ln phi) / (kappa,L_res) algebraic independence.  Those
# stay EXTERNAL-OPEN / CONDITIONAL.  A "forced" claim about either would be a
# discipline violation; there is none.
#
# RESTATEMENT TEST (mandatory, PF-2E): kappa -> kappa_k=(4k+1)kappa; (4k+1) is a
# nonzero (odd) rational for every k, so every PF statement is verbatim per
# branch.  PF is arithmetic ABOUT the constants -- BRANCH-BLIND: a
# characterization / bounded-negative, NOT an OP-RATE derivation (HANDOFF sec 2;
# OP-RATE closed AS A CLASSIFICATION in v1.8 -- not reopened here).

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   simplify, expand, radsimp, trigsimp, together, cancel,
                   expand_log, Poly, minimal_polynomial, factor_list, resultant,
                   discriminant, CRootOf, re, im)
from sympy.polys.polyerrors import NotAlgebraic

PASS, FAIL, GUARD, COMP = [], [], [], []

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

def pk(cid, cond, desc):
    # COMPUTED (PSLQ) detector check -- counted separately; a FALSE here is a
    # real failure (detector broken on a positive, or a genuine relation found
    # where the claim says none), so it also fails the harness.
    ok = bool(cond)
    COMP.append((cid, ok, desc))
    print(("COMP-PASS" if ok else "COMP-FAIL"), cid, "-", desc)
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

x, Y, X, t, th, Kc = symbols('x Y X t theta K')
kk = symbols('k', integer=True)

s5    = sqrt(5)
s13   = sqrt(13)
phi   = (1 + s5)/2
tau   = (s5 - 1)/2
kappa = pi/(2*log(phi))
gap   = phi**-4
K2    = 1 - gap
tau0  = (1 + s13)/2
beta4 = CRootOf(x**4 - x**3 - x**2 - x + 1, 1)

def qdAB(e, D):
    r = sqrt(D)
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, r)
    assert p.degree() <= 1, "not linear in sqrt(%s): %s" % (D, e)
    B = p.coeff_monomial(r) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQd(e, D):
    A, B = qdAB(e, D)
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    if A > 0:
        return 1 if A*A > D*B*B else -1
    return 1 if D*B*B > A*A else -1

def irred(P, v):
    fl = factor_list(P, v)[1]
    return len(fl) == 1 and fl[0][1] == 1

def is_alg_integer(P, v):
    pp = Poly(P, v)
    cs = pp.all_coeffs()
    return cs[0] == 1 and all(c == int(c) for c in cs)

print("=" * 70)
print("PF -- PERIOD FRONTIER: forced ledger + tool-gap map + bounded-negative.")
print("CITE [established, NOT harness-checkable] -- the deciding tools:")
print("  Gelfond-Schneider (1934): log(alpha)/log(beta) rational or transc. for")
print("    nonzero algebraic alpha,beta.  (-1,phi): i pi/ln phi transc => kappa")
print("    transcendental [inherited tx/T3, thm:ratetrans]. REACH: alpha^beta,")
print("    log-ratios.  GAP: no products/elliptic values/alg-independence.")
print("  Baker (1966-68, linear forms in logs): logs of multiplicatively-indep")
print("    algebraics are lin. indep. over Qbar.  REACH: Q-linear independence")
print("    of {1,ln2,ln3,ln phi} (FORCED-GIVEN-ESTABLISHED here). GAP: NONlinear")
print("    (algebraic) relations, elliptic-integral values.")
print("  Schneider 1937 / Chudnovsky: elliptic periods transcendental / alg-")
print("    independent -- REQUIRE an ALGEBRAIC modulus or CM lattice. GAP: our")
print("    modulus k=sqrt m is TRANSCENDENTAL (m=kappa^2/(1+kappa^2), PF-3E),")
print("    so INAPPLICABLE to L_res.")
print("  Schanuel's conjecture: with {ln phi, i pi} Q-linearly independent and")
print("    e^{ln phi}=phi, e^{i pi}=-1 algebraic => trdeg >= 2 => (pi,ln phi)")
print("    algebraically independent.  CONDITIONAL (PF-2C).")
print("  Kontsevich-Zagier 2001 (Periods): L_res is a period-FUNCTION value at a")
print("    transcendental argument -- one step beyond the pure period ring;")
print("    relations are period-conjecture territory.  EXTERNAL-OPEN.")
print("CITE [corpus companion, QUOTED] Echo-S-Research")
print("  papers/2026-06-salem-slot/salem_slot.tex Def 4 (def:redirect):")
print("    'the redirection of a Salem number beta is tau_0 := tr(beta) =")
print("     beta + beta^{-1}, the dominant root of its trace-down T ... a")
print("     totally real integer polynomial'; Thm (thm:occupant): 'tau_0 > 2'.")
print("  papers/2026-07-complex-rung-generator v1.8 Rem (rem:periodmap),")
print("    Prop (prop:tau0alg), Prop (prop:sphlength): L_res = sqrt(1+kappa^2)")
print("    E(kappa^2/(1+kappa^2)), the exact reduction, the external-open ledger.")
print("PF makes that ledger a validated-detector bounded-negative + tool-gap map.")
print("=" * 70)

# ================================================================ PF-3
# The FORCED sub-relations ledger: the exact structure the frontier sits on.
print("== PF-3: FORCED sub-relations ledger (exact identities) ==")

ck("PF-3A", zero(phi**2 - phi - 1) and zero(phi*tau - 1)
        and zero(tau**2 + tau - 1) and zero(phi + 1/phi - s5),
   "[forced] seed + golden-trace: phi^2=phi+1, phi*tau=1, tau^2+tau=1,"
   " phi+phi^{-1}=sqrt5 (the redirection limit sqrt5, salem-slot golden floor)")

ck("PF-3B", zero(2/kappa - 4*log(phi)/pi) and zero(pi/log(phi) - 2*kappa)
        and zero(log(phi)/pi - 1/(2*kappa)) and zero(pi - 2*kappa*log(phi)),
   "[forced identity] the definitional tie: 2/kappa = 4 ln phi/pi, i.e."
   " pi/ln phi = 2 kappa, ln phi/pi = 1/(2 kappa), pi = 2 kappa ln phi"
   " (kappa := pi/(2 ln phi); tx/T4 near-miss anchor)")

ck("PF-3C", expand((pi/2 + 2*pi*kk)/log(phi) - (4*kk + 1)*kappa) == 0
        and zero((4*kk + 1)*kappa/(4*kk + 1) - kappa)
        and all((4*v + 1) != 0 and (4*v + 1) % 2 == 1 for v in range(-40, 41)),
   "[forced] branch family: kappa_k = (pi/2+2 pi k)/ln phi = (4k+1) kappa,"
   " symbolic in integer k; (4k+1) is a NONZERO ODD rational for every k")

ck("PF-3D", expand(minimal_polynomial(tau0, t) - (t**2 - t - 3)) == 0
        and zero(tau0**2 - tau0 - 3) and zero(tau0 - (1 + s13)/2)
        and Poly(t**2 - t - 3, t).degree() == 2 and irred(t**2 - t - 3, t)
        and is_alg_integer(t**2 - t - 3, t) and discriminant(Poly(t**2 - t - 3, t)) == 13,
   "[forced, inherited PM-1B] tau_0^2 = tau_0 + 3: tau_0=(1+sqrt13)/2 is the"
   " Salem-slot occupant, algebraic integer of degree 2, minpoly t^2-t-3,"
   " disc 13 (salem-slot Def 4 redirection; the algebraic anchor of the frontier)")

m_mod = Kc**2/(1 + Kc**2)
red = 1 + Kc**2*cos(th)**2 - (1 + Kc**2)*(1 - m_mod*sin(th)**2)
ck("PF-3E", zero(red) and zero(Kc**2 - m_mod/(1 - m_mod)),
   "[forced identity, prop:sphlength] the elliptic reduction, symbolic in K,th:"
   " 1 + K^2 cos^2 th = (1+K^2)(1 - m sin^2 th) with m=K^2/(1+K^2); Mobius"
   " inverse K^2 = m/(1-m).  Hence L_res = int_0^{pi/2} sqrt(1+kappa^2 cos^2)"
   " = sqrt(1+kappa^2) E(m): the FORCED defining identity of L_res")
ck("PF-3Ef", not zero(1 + Kc**2*cos(th)**2
                      - (1 + Kc**2)*(1 - (m_mod + Rational(1, 1000))*sin(th)**2)),
   "FALSIFIER: perturbing the modulus m -> m+1/1000 breaks the reduction"
   " identity (rejected) -- PF-3E has teeth")

ck("PF-3F", zero(K2 - (1 - gap)) and zero(gap - phi**-4)
        and expand(minimal_polynomial(K2, x) - (x**2 + 5*x - 5)) == 0
        and zero((x**2 + 5*x - 5).subs(x, K2)) and sgnQd(K2, 5) == 1,
   "[forced] the lens gate: K^2 = 1 - gap, gap = phi^{-4}; K^2 algebraic deg-2,"
   " minpoly x^2+5x-5, and K^2 > 0 (exact Q(sqrt5) sign) [tx ledger]")

# NEW forced geometric identity: the residual ellipse's focal parameter = kappa.
a2 = 1 + Kc**2          # semi-major^2 of the residual ellipse (axes sqrt(1+k^2),1)
c2 = a2 - 1             # linear-eccentricity^2 = a^2 - b^2, b = 1
ck("PF-3G", zero(c2 - Kc**2) and zero(c2/a2 - m_mod)
        and zero(sqrt(a2)**2 - (1 + Kc**2)),
   "[forced, NEW] geometric reading: L_res is the quarter-perimeter of the"
   " ellipse with semi-axes (sqrt(1+kappa^2), 1); its linear eccentricity c has"
   " c^2 = a^2 - b^2 = kappa^2 EXACTLY (focal parameter = kappa), and its"
   " eccentricity^2 = c^2/a^2 = m EXACTLY (the elliptic parameter). Symbolic in K")

# ================================================================ PF-2
# The precise tool-gap map: CONDITIONAL / EXTERNAL-OPEN, each gap named.
print("== PF-2: the precise CONDITIONAL / EXTERNAL-OPEN tool-gap map ==")

# PF-2A -- inherited settled layer + engine teeth (algebraicity engine not vacuous)
teeth = 0
from mpmath import iv, mp  # (mp used later; iv for the guards)
import mpmath as mnum
for tv in (kappa, pi, log(phi)):
    try:
        minimal_polynomial(tv, X)
    except NotAlgebraic:
        teeth += 1
ck("PF-2A", zero(I*pi/log(phi) - 2*I*kappa)
        and expand(minimal_polynomial(2*I, X) - (X**2 + 4)) == 0
        and expand(minimal_polynomial(tau0, X) - (X**2 - X - 3)) == 0
        and teeth == 3,
   "[forced-given-established / inherited] settled layer: kappa=(i pi/ln phi)/(2i),"
   " 2i algebraic (X^2+4) => kappa in Qbar IFF i pi/ln phi in Qbar; G-S kills the"
   " latter (kappa transcendental). tau_0 algebraic (X^2-X-3). ENGINE TEETH:"
   " minimal_polynomial RAISES NotAlgebraic on kappa, pi, ln phi (not vacuous)")

# PF-2B -- Baker's REACH: Q-linear independence of {1, ln2, ln3, ln phi}.
# Premise (multiplicative independence of 2,3,phi), grounded in-box exactly:
# 2^a 3^b phi^c = 1 with (a,b,c) in [-6,6]^3, (a,b,c)!=0, has NO solution.
BOX = range(-6, 7)
def phi_pow_AB(c):
    A, B = qdAB(phi**c, 5)   # phi^c = A + B sqrt5, exact rationals
    return A, B
mult_indep = True
for c in BOX:
    Ac, Bc = phi_pow_AB(c)
    if Bc != 0:
        continue            # phi^c irrational => 2^a 3^b phi^c irrational != 1
    # Bc == 0 <=> c == 0 (Fibonacci F_c = 0 only at c=0); then A_c = 1:
    if c != 0 or Ac != 1:
        mult_indep = False
    else:
        for a in BOX:
            for b in BOX:
                if a == 0 and b == 0:
                    continue
                if Fraction(2)**a * Fraction(3)**b == 1:   # 2^a 3^b = 1 nontrivially?
                    mult_indep = False
ck("PF-2B", mult_indep and phi_pow_AB(0) == (Fraction(1), Fraction(0))
        and all(phi_pow_AB(c)[1] != 0 for c in BOX if c != 0),
   "[FORCED-GIVEN-ESTABLISHED, Baker] Q-linear independence of {1,ln2,ln3,ln phi}."
   " Grounded premise (exact, box |a|,|b|,|c|<=6): 2^a 3^b phi^c = 1 only at"
   " (0,0,0) -- for c!=0 the sqrt5-part F_c/2 != 0 (phi^c irrational); for c=0,"
   " 2^a 3^b=1 forces a=b=0. So 2,3,phi are multiplicatively independent; by"
   " Baker their logs (with 1) are Qbar-linearly independent [established]")

# PF-2C -- Schanuel CONDITIONAL: (pi, ln phi) and (kappa, ln phi) alg-independent.
ck("PF-2C", simplify(im(log(phi))) == 0 and simplify(re(I*pi)) == 0
        and simplify(im(I*pi)) == pi and sgnQd(phi - 1, 5) == 1
        and expand(minimal_polynomial(phi, X) - (X**2 - X - 1)) == 0
        and expand(minimal_polynomial(-1, X) - (X + 1)) == 0
        and zero(pi - 2*kappa*log(phi)) and zero(kappa - pi/(2*log(phi))),
   "[CONDITIONAL on Schanuel] premise re-encoded exactly: {ln phi, i pi} is"
   " Q-linearly independent (ln phi real nonzero: im=0, phi>1; i pi purely"
   " imaginary nonzero: re=0, im=pi), and e^{ln phi}=phi, e^{i pi}=-1 are"
   " algebraic (X^2-X-1, X+1). Schanuel => trdeg Q(ln phi,i pi,phi,-1) >= 2 =>"
   " (pi,ln phi) algebraically independent. Since Q(kappa,ln phi)=Q(pi,ln phi)"
   " (pi=2 kappa ln phi, kappa=pi/(2 ln phi)), ALSO (kappa,ln phi) alg-indep."
   " NOT forced -- Schanuel is conjectural")

# PF-2D -- Schneider/Chudnovsky GAP: the modulus is transcendental.
ck("PF-2D", zero(Kc**2 - m_mod/(1 - m_mod))
        and zero(m_mod.subs(Kc, kappa) - kappa**2/(1 + kappa**2)),
   "[EXTERNAL-OPEN, gap named] Schneider/Chudnovsky elliptic-period theorems"
   " need an ALGEBRAIC modulus; here k=sqrt m with m=kappa^2/(1+kappa^2), and"
   " the Mobius inverse K^2=m/(1-m) gives m in Qbar IFF kappa in Qbar -- kappa"
   " transcendental (G-S) => m TRANSCENDENTAL => the theorems are INAPPLICABLE."
   " L_res transcendence is NOT settled by any shipped tool: EXTERNAL-OPEN")

# PF-2E -- restatement test / branch-blindness (mandatory).
ck("PF-2E", expand((pi/2 + 2*pi*kk)/log(phi) - (4*kk + 1)*kappa) == 0
        and zero((4*kk + 1)*kappa/(4*kk + 1) - kappa)
        and all((4*v + 1) % 2 == 1 and (4*v + 1) != 0 for v in range(-40, 41)),
   "RESTATEMENT TEST: kappa_k=(4k+1)kappa; (4k+1) nonzero odd rational for all k,"
   " so kappa_k in Qbar IFF kappa in Qbar and EVERY PF statement is verbatim per"
   " branch.  Branch-blind: PF SELECTS no k -- a characterization / bounded-"
   " negative, NOT an OP-RATE derivation (OP-RATE closed as a classification v1.8)")

print("   LEDGER (verdict | tool | gap):")
print("     kappa transcendental        | SETTLED  | Gelfond-Schneider (-1,phi)")
print("     tau_0 algebraic (t^2-t-3)   | SETTLED  | salem-slot / minpoly")
print("     (kappa,tau_0) no Q-relation | SETTLED  | resultant/tower + G-S [PM-1X]")
print("     {1,ln2,ln3,ln phi} Q-lin ind| SETTLED  | Baker [PF-2B]")
print("     L_res transcendental?       | OPEN     | modulus transc: no elliptic tool")
print("     (pi,ln phi) alg-independent | COND.    | Schanuel [PF-2C]; ratio=2kappa settled")
print("     (kappa,ln phi) alg-indep    | COND.    | Schanuel [PF-2C]")
print("     (kappa,L_res) alg-indep     | OPEN     | beyond Schanuel; period-conj (K-Z)")

# ================================================================ PF-1
# VALIDATED-DETECTOR bounded-negative (COMPUTED, counted separately).
print("== PF-1: validated-detector PSLQ bounded-negative (COMPUTED) ==")
DPS   = 230
Hbnd  = 10**12          # coefficient-height (max |coeff|) bound
MSTEP = 20000
GEN   = mnum.mpf(10)**-180   # residual below this <=> a GENUINE relation
mp.dps = DPS
print("   LOG: precision = %d dps ; coefficient-height bound H = 10^12 ;" % DPS)
print("        maxsteps = %d ; genuine-residual threshold = 1e-180." % MSTEP)
print("   A negative = PSLQ returns None (no integer relation with |coeff|<=H).")
print("   COMPUTED corroboration ONLY -- never a transcendence/independence proof.")

fphi  = (1 + mnum.sqrt(5))/2
flnphi = mnum.log(fphi)
fkappa = mnum.pi/(2*flnphi)
fm     = fkappa**2/(1 + fkappa**2)
fLresE = mnum.sqrt(1 + fkappa**2)*mnum.ellipe(fm)                     # closed form
fLresQ = mnum.quad(lambda a: mnum.sqrt(1 + fkappa**2*mnum.cos(a)**2),
                   [0, mnum.pi/2])                                    # quadrature
fLres  = fLresE
ftau0  = (1 + mnum.sqrt(13))/2
fK2    = 1 - fphi**-4
fbeta4 = mnum.findroot(lambda z: z**4 - z**3 - z**2 - z + 1, mnum.mpf('1.7220838'))
fln2   = mnum.log(2)
fln3   = mnum.log(3)
one    = mnum.mpf(1)

def search(vec):
    r = mnum.pslq(vec, maxcoeff=Hbnd, maxsteps=MSTEP)
    if r is None:
        return None, None
    resid = abs(sum(mnum.mpf(a)*b for a, b in zip(r, vec)))
    return r, resid

def same_rel(r, target):
    return r is not None and (list(r) == list(target)
                              or list(r) == [-c for c in target])

# ---- fail-first POSITIVES: the detector MUST FIND these known relations ----
r, res = search([one, fphi, fphi**2])
pk("PF-1P1", same_rel(r, [-1, -1, 1]) and res < GEN,
   "POSITIVE (deg-2, size 3): {1,phi,phi^2} -> %s (phi^2-phi-1=0), resid<1e-180"
   " -- detector finds the golden minimal relation" % (r,))
r, res = search([one, ftau0, ftau0**2])
pk("PF-1P2", same_rel(r, [-3, -1, 1]) and res < GEN,
   "POSITIVE (deg-2, size 3): {1,tau_0,tau_0^2} -> %s (tau_0^2-tau_0-3=0)"
   " -- detector finds the Salem-slot occupant's relation" % (r,))
r, res = search([one, fK2, fK2**2])
pk("PF-1P3", same_rel(r, [-5, 5, 1]) and res < GEN,
   "POSITIVE (deg-2, size 3): {1,K^2,K^4} -> %s (K^4+5K^2-5=0, minpoly x^2+5x-5)"
   % (r,))
r, res = search([one, fbeta4, fbeta4**2, fbeta4**3, fbeta4**4])
pk("PF-1P4", same_rel(r, [1, -1, -1, -1, 1]) and res < GEN,
   "POSITIVE (deg-4, size 5): {1,beta4,..,beta4^4} -> %s (least-degree Salem"
   " x^4-x^3-x^2-x+1) -- HIGHER-DEGREE detection validated" % (r,))
r, res = search([1/fkappa, flnphi/mnum.pi])
pk("PF-1P5", same_rel(r, [1, -2]) and res < GEN,
   "POSITIVE (definitional corpus relation): {1/kappa, ln phi/pi} -> %s"
   " (1/kappa = 2 ln phi/pi, i.e. 2/kappa = 4 ln phi/pi, PF-3B)" % (r,))
# size-6 PLANTED positive -- SAME basis size as the flagship negative PF-1N6:
r, res = search([one, fphi, fphi**2, ftau0, ftau0**2, fphi*ftau0])
pk("PF-1P6", r is not None and res < GEN and any(c != 0 for c in r),
   "POSITIVE (size 6, PLANTED): {1,phi,phi^2,tau_0,tau_0^2,phi*tau_0} -> %s"
   " resid<1e-180 -- a genuine relation IS found in a SIZE-6 basis at this"
   " precision/height, so the size-6 NEGATIVES are not false-negatives" % (r,))

# ---- TRUE-NEGATIVE controls: KNOWN transcendental => detector finds nothing ----
r, res = search([one, fkappa, fkappa**2, fkappa**3, fkappa**4, fkappa**5, fkappa**6])
pk("PF-1N1", r is None,
   "TRUE-NEGATIVE control: {1,kappa,..,kappa^6} -> None. kappa is transcendental"
   " (G-S, established), so NO algebraic relation exists at any degree; the"
   " detector correctly reports none (no false-positive) [degree<=6, H<=1e12]")
r, res = search([one, mnum.pi, mnum.pi**2, mnum.pi**3, mnum.pi**4, mnum.pi**5, mnum.pi**6])
pk("PF-1N2", r is None,
   "TRUE-NEGATIVE control: {1,pi,..,pi^6} -> None. pi transcendental (Lindemann,"
   " established); detector correctly reports none [degree<=6, H<=1e12]")

# ---- CORROBORATION bounded-negatives on the OPEN questions (COMPUTED only) ----
r2, _ = search([one, fLres, fLres**2])
r4, _ = search([one, fLres, fLres**2, fLres**3, fLres**4])
r6, _ = search([one, fLres, fLres**2, fLres**3, fLres**4, fLres**5, fLres**6])
pk("PF-1N3", r2 is None and r4 is None and r6 is None,
   "BOUNDED-NEGATIVE [COMPUTED]: L_res algebraicity ladder deg 2,4,6 all -> None"
   " => L_res is NOT algebraic of degree <= 6 with height <= 1e12. Corroborates"
   " (does NOT prove) L_res transcendental (EXTERNAL-OPEN, PF-2D)")
r, _ = search([one, fkappa, fkappa**2, fLres, fLres**2, fkappa*fLres])
pk("PF-1N4", r is None,
   "BOUNDED-NEGATIVE [COMPUTED] (flagship, size 6): {1,kappa,kappa^2,L_res,"
   "L_res^2,kappa*L_res} -> None => no such Q-linear relation with H<=1e12."
   " Corroborates (kappa,L_res) algebraic independence at this monomial"
   " degree/height (EXTERNAL-OPEN; validated size-6 by PF-1P6)")
r, _ = search([one, mnum.pi, flnphi, mnum.pi*flnphi])
pk("PF-1N5", r is None,
   "BOUNDED-NEGATIVE [COMPUTED]: {1,pi,ln phi,pi*ln phi} -> None. Corroborates"
   " (pi,ln phi) algebraic independence (CONDITIONAL on Schanuel, PF-2C);"
   " NOT a proof [H<=1e12]")
r, _ = search([one, flnphi, fln2, fln3])
pk("PF-1N6", r is None,
   "BOUNDED-NEGATIVE [COMPUTED]: {1,ln phi,ln2,ln3} -> None. Numeric witness for"
   " the Q-linear independence that Baker FORCES (PF-2B): consistent, H<=1e12")
r, _ = search([one, fkappa, fLres])
r2, _ = search([one, fkappa, ftau0, fLres])
pk("PF-1N7", r is None and r2 is None,
   "BOUNDED-NEGATIVE [COMPUTED]: {1,kappa,L_res} and {1,kappa,tau_0,L_res} both"
   " -> None (no linear relation, H<=1e12) -- extends PM-G6 to the L_res-joint"
   " bases at 230 dps")

# ================================================================ PF-G guards
print("== PF-G: certified-interval + numeric-corroboration guards (separate) ==")
iv.dps = 60
iphi   = (1 + iv.sqrt(5))/2
ilnphi = iv.log(iphi)
ikappa = iv.pi/(2*ilnphi)
ikap2  = ikappa**2
im_mod = ikap2/(1 + ikap2)
itau0  = (1 + iv.sqrt(13))/2

gk("PF-G1", (ikap2 > iv.mpf(10)) and (ikap2 < iv.mpf(11)),
   "interval-certified (60 dps): 10 < kappa^2 < 11 [inherits MX/PM guards]")
gk("PF-G2", (im_mod > iv.mpf(10)/11) and (im_mod < iv.mpf(11)/12),
   "interval-certified: modulus m in (10/11,11/12) subset (0,1) -- a proper"
   " (transcendental, PF-2D) elliptic modulus")
gk("PF-G3", (iv.sqrt(1 + ikap2) > iv.mpf('3.41')) and (iv.sqrt(1 + ikap2) < iv.mpf('3.42'))
        and (ikappa > itau0 + iv.mpf('0.9')),
   "interval-certified: semi-major axis sqrt(1+kappa^2) in (3.41,3.42) (brackets"
   " L_res=3.7312193125); kappa - tau_0 > 0.9 (the constants are distinct)")

# fail-first for L_res: two INDEPENDENT computations must agree (validates the
# constant AND the reduction identity PF-3E numerically).
gk("PF-G4", abs(fLresE - fLresQ) < mnum.mpf(10)**-200,
   "L_res fail-first: closed form sqrt(1+kappa^2)E(m) and DIRECT quadrature of"
   " int sqrt(1+kappa^2 cos^2) agree to < 1e-200 -- the constant and the"
   " reduction identity are numerically confirmed (the frontier is well-posed)")

# decimal audit vs paper/register literals (pm/tx idiom).
mp.dps = 60
def audit_ok(val, lit):
    dp = len(lit.split('.')[1])
    litv = mnum.mpf(lit)
    ulps = abs(val - litv)*mnum.mpf(10)**dp
    if ulps <= mnum.mpf('0.51'):
        return True
    trunc = mnum.floor(val*mnum.mpf(10)**dp)/mnum.mpf(10)**dp
    return abs(trunc - litv) < mnum.mpf(10)**(-dp - 4)
a_phi  = (1 + mnum.sqrt(5))/2
a_ln   = mnum.log(a_phi)
a_kap  = mnum.pi/(2*a_ln)
a_Lres = mnum.sqrt(1 + a_kap**2)*mnum.ellipe(a_kap**2/(1 + a_kap**2))
audits = [
    ("kappa",  a_kap,                 '3.2642513026'),
    ("L_res",  a_Lres,                '3.7312193125'),
    ("tau_0",  (1 + mnum.sqrt(13))/2, '2.3027756377'),
    ("m_mod",  a_kap**2/(1 + a_kap**2), '0.9142023918'),
]
okA = True
for name, val, lit in audits:
    hit = audit_ok(val, lit)
    okA = okA and hit
    print("   audit %-8s printed %s  recomputed %s  %s"
          % (name, lit, mnum.nstr(val, 14), "ok" if hit else "MISMATCH"))
gk("PF-G5", okA,
   "decimal audit at 60 dps: kappa, L_res, tau_0, m match the paper/register"
   " literals (0.51 ulp or floor-truncation)")

# ================================================================ summary
print("=" * 70)
n_guard_fail = sum(1 for _, ok, _ in GUARD if not ok)
n_comp_fail  = sum(1 for _, ok, _ in COMP if not ok)
n_exact_fail = len(FAIL) - n_guard_fail - n_comp_fail
n_pos = sum(1 for cid, ok, _ in COMP if ok and 'P' in cid.split('-')[1][:2])
n_neg = sum(1 for cid, ok, _ in COMP if ok and cid.split('-')[1].startswith('1N'))
print("PF EXACT: %d passed, %d failed | COMPUTED: %d/%d PSLQ checks (%d positives"
      " found, %d bounded-negatives) | GUARDS: %d/%d certified"
      % (len(PASS), n_exact_fail,
         sum(1 for _, ok, _ in COMP if ok), len(COMP), n_pos, n_neg,
         sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
