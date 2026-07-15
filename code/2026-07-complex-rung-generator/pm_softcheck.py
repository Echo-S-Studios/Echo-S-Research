#!/usr/bin/env python3
# pm_softcheck.py -- period-relation map for the corpus transcendentals.
# Dossier PM (v1.8): make the v1.7 register's EXTERNAL-OPEN item
#   "period relations among {kappa, L_res, tau_0}"
# precise and PARTIALLY FORCED.  Targets PM-1 (forced sub-ledger: the exact
# arithmetic TYPE of each constant + the forced kappa-vs-tau_0 negative
# settlement), PM-2 (the precise EXTERNAL-OPEN map, bounded + citable),
# PM-3 (numeric corroboration guards).
#
# Written cold this session (2026-07-15) from the dossier + papers; machinery
# re-encoded from scratch (no import of any corpus/companion code).
#
# DISCIPLINE (nx/qx/tx/mx idiom):
#  * every DECISION is exact over Q / Q(sqrt5) / Q(sqrt13) (extended by i where
#    needed) and, for algebraicity, at the polynomial (minimal_polynomial /
#    resultant) level;  Q(sqrt D) real signs via rational arithmetic only.
#  * the ONLY transcendental inputs are the certified mpmath-iv interval guards
#    (group PM-G, 60 dps) and the 50/140-dps numeric displays -- counted
#    SEPARATELY;  floats never touch an exact decision.
#  * fail-first: every group carries at least one falsifier -- a deliberately
#    perturbed / planted variant that MUST be rejected, or a genuine equality
#    the engine MUST detect -- written before the claim it protects:
#    PM-A2f, PM-1B5f, PM-1B6f, PM-1C1f, PM-1C2f, PM-1X2, PM-1X5, PM-2Bf, PM-G*.
#
# WHAT IS NEW (forced) vs INHERITED (cited):
#  * INHERITED [tx / T3, forced given Gelfond-Schneider]: kappa = pi/(2 ln phi)
#    is transcendental; kappa not in Qbar.  NOT re-derived here (G-S is not
#    harness-checkable); re-encoded only as the bridge kappa=(i pi/ln phi)/(2i).
#  * NEW [forced, PM-1B]: tau_0 (the Salem-slot occupant) is ALGEBRAIC -- its
#    exact minimal polynomial recovered and verified, degree 2, totally real,
#    Perron, algebraic integer, tau_0 > 2 > phi;  the whole redirection family
#    and its golden limit sqrt5 are algebraic.
#  * NEW [forced given established, PM-1X]: kappa transcendental vs tau_0
#    algebraic => NO polynomial P(kappa,tau_0)=0 over Q that constrains kappa;
#    the register's "kappa vs tau_0" period-relation worry is SETTLED
#    NEGATIVELY at the algebraic level.  Encoded as an exact resultant/tower
#    argument with teeth.
#  * NEW [forced reduction, PM-2B]: any algebraic relation between tau_0 and
#    L_res forces L_res algebraic -- the (tau_0,L_res) question COLLAPSES to
#    the single-constant "is L_res transcendental" (EXTERNAL-OPEN).
#  * MAP [bounded EXTERNAL-OPEN, PM-2]: L_res transcendence; (kappa,L_res) and
#    (pi,ln phi) algebraic independence remain open, with the deciding tools
#    stated precisely (G-S / Baker settle what; Schneider-elliptic needs an
#    algebraic modulus we do not have; Schanuel/Kontsevich-Zagier conjectural).
#
# RESTATEMENT TEST (mandatory, PM-1X6): substitute kappa -> kappa_k=(4k+1)kappa.
# (4k+1) is a nonzero rational for every k, so kappa_k is transcendental for
# every k (tx/T3) and every PM claim is verbatim identical per branch.  PM is
# arithmetic ABOUT the constants, BRANCH-BLIND: a characterization/exclusion,
# NOT an OP-RATE derivation.  Said plainly in-harness.

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   simplify, expand, radsimp, trigsimp, together, cancel,
                   expand_log, Poly, minimal_polynomial, factor_list, resultant,
                   discriminant, CRootOf, re, im, real_roots)
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

x, Y, X, t, th, Kc = symbols('x Y X t theta K')
kk = symbols('k', integer=True)

s5    = sqrt(5)
s13   = sqrt(13)
phi   = (1 + s5)/2
tau   = (s5 - 1)/2
kappa = pi/(2*log(phi))

# the least-degree Salem number (deg 4), already a corpus ledger constant (tx):
P_salem = x**4 - x**3 - x**2 - x + 1
beta4   = CRootOf(P_salem, 1)              # the real root > 1 (index 1)
# its trace redirection = the Salem-slot occupant of the companion paper:
tau0    = (1 + s13)/2                        # = beta4 + 1/beta4, root of t^2-t-3
muS     = CRootOf(x**3 - x - 1, 0)           # plastic (smallest Pisot)

# ------------------------------------------------------------- exact Q(sqrt D)
def qdAB(e, D):
    """Exact (A,B) with e = A + B sqrt(D), A,B in Q. Raises if not in Q(sqrt D)."""
    r = sqrt(D)
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, r)
    assert p.degree() <= 1, "not linear in sqrt(%s): %s" % (D, e)
    B = p.coeff_monomial(r) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQd(e, D):
    """Exact sign of A + B sqrt(D) via rational arithmetic only."""
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

def totally_real(P, v):
    """Exact: number of real roots (Sturm) equals the degree."""
    pp = Poly(P, v)
    return pp.count_roots() == pp.degree()

def is_alg_integer(P, v):
    """P monic with integer coefficients (=> its roots are algebraic integers)."""
    pp = Poly(P, v)
    cs = pp.all_coeffs()
    return cs[0] == 1 and all(c == int(c) for c in cs)

print("=" * 68)
print("CITE [established, NOT harness-checkable]")
print("  Gelfond-Schneider (A.O.Gelfond 1934; Th.Schneider 1934, indep.):")
print("    for nonzero algebraic alpha,beta (fixed log branches, log!=0),")
print("    log(alpha)/log(beta) is rational or transcendental.")
print("    Instance (alpha,beta)=(-1,phi): i*pi/ln(phi) transcendental => kappa")
print("    transcendental, kappa not in Qbar.  [inherited: tx/T3, thm:ratetrans]")
print("  Baker (1966-68, linear forms in logs): for nonzero algebraic alpha_j")
print("    with log alpha_j lin. indep. over Q, {1,log alpha_1,...} are lin.")
print("    indep. over Qbar; any nonzero Qbar-linear form in the logs is")
print("    transcendental.  (Re-proves kappa transcendental; does NOT reach")
print("    products/ratios beyond logs, elliptic-integral values, or algebraic")
print("    independence.)")
print("  Schneider 1937 / Chudnovsky: periods of elliptic curves are")
print("    transcendental / alg-independent -- but require an ALGEBRAIC modulus")
print("    or CM lattice (inapplicable: our modulus is transcendental, below).")
print("  Kontsevich-Zagier 2001 (Periods): the framework for L_res's status.")
print("  Schanuel's conjecture: would settle (pi,ln phi) alg. independence.")
print("CITE [corpus companion, QUOTED]  Echo-S-Research paper")
print("  papers/2026-06-salem-slot/salem_slot.tex:")
print("  Def 3.7 (def:redirect): 'the redirection of a Salem number beta is")
print("    tau_0 := trace(beta) = beta + beta^{-1}, the dominant root of its")
print("    trace-down T; the redirected object is T, a totally real integer")
print("    polynomial.'  Cor 2.4 (cor:nofourth): 'its trace-down T is totally")
print("    real'.  Thm 3.8 (thm:occupant): 'tau_0 = beta+beta^{-1} > 2 ... >")
print("    2 > phi'.  Benchmark table: beta_4 (deg 4) has trace-down t^2-t-3,")
print("    tau_0 = 2.302776.  Lem 6.1: phi+phi^{-1}=sqrt5. Thm 6.2: beta_n->phi")
print("    => tau_0->sqrt5. Prop 8.1: plastic mu_P (x^3-x-1) -> 2.079596.")
print("CITE [corpus, v1.7 tex L.497-500, Remark scan-retirement, QUOTED]")
print("  'period relations among the corpus transcendentals (kappa, L_res, and")
print("   the register's tau_0) and the algebraic independence of pi and ln phi")
print("   are beyond Gelfond-Schneider and remain external-open.'")
print("  This harness makes that sentence PRECISE and PARTIALLY FORCED.")
print("=" * 68)

# ================================================================ PM-A
print("== PM-A: replication of the inherited layer (no new content) ==")
ck("PM-A1", zero(phi**2 - phi - 1) and zero(phi*tau - 1) and zero(tau**2 + tau - 1)
        and zero(phi + 1/phi - s5),
   "seed + golden-trace identities: phi^2=phi+1, phi*tau=1, tau^2+tau=1,"
   " phi+phi^{-1}=sqrt5 (the redirection limit, salem-slot Lem 6.1)")
ck("PM-A2", zero(I*pi/log(phi) - 2*I*kappa)
        and zero((I*pi/log(phi))/(2*I) - kappa)
        and expand(minimal_polynomial(2*I, X) - (X**2 + 4)) == 0,
   "kappa bridge: i*pi/ln(phi)=2i*kappa, kappa=(i*pi/ln phi)/(2i), and 2i is"
   " algebraic (X^2+4) => kappa in Qbar IFF i*pi/ln phi in Qbar; T3 (G-S) kills"
   " the latter, so kappa NOT in Qbar [inherited tx/T3, forced given G-S]")
teeth = 0
for tv in (kappa, pi, log(phi)):
    try:
        minimal_polynomial(tv, X)
    except NotAlgebraic:
        teeth += 1
ck("PM-A2f", teeth == 3,
   "FALSIFIER (engine teeth): minimal_polynomial RAISES NotAlgebraic on kappa,"
   " pi, ln(phi) -- the algebraicity engine is not vacuous (this is a teeth"
   " check only; transcendence of kappa is G-S/T3, of pi is Lindemann)")

# ================================================================ PM-1B
print("== PM-1B: arithmetic TYPE of tau_0 = ALGEBRAIC (forced) ==")
Tsal = t**2 - t - 3                        # claimed trace-down / minpoly of tau0
ck("PM-1B1", zero(expand(x**2*Tsal.subs(t, x + 1/x) - P_salem)),
   "trace-down (salem-slot Def 3.1): for the deg-4 Salem P=x^4-x^3-x^2-x+1,"
   " x^2 * T(x+1/x) = P(x) with T(t)=t^2-t-3 -- exact symbolic in x")
ck("PM-1B2", irred(P_salem, x) and zero(expand(x**4*P_salem.subs(x, 1/x) - P_salem))
        and Poly(P_salem, x).subs(x, 1) < 0 and Poly(P_salem, x).subs(x, 2) > 0,
   "beta_4 is a Salem number: minpoly x^4-x^3-x^2-x+1 irreducible, reciprocal"
   " (x^4 P(1/x)=P), and P(1)=-1<0<P(2)=3 => a REAL root in (1,2) => NON-"
   "cyclotomic (cyclotomic roots lie on |z|=1) [least-degree Salem, in tx ledger]")
ck("PM-1B3", expand(minimal_polynomial(beta4 + 1/beta4, t) - Tsal) == 0
        and simplify((beta4 + 1/beta4) - tau0) == 0
        and expand(minimal_polynomial(tau0, t) - Tsal) == 0,
   "tau_0 = beta_4 + 1/beta_4 has minimal polynomial t^2-t-3 and equals"
   " (1+sqrt13)/2 exactly (salem-slot Def 3.7 redirection; table row beta_4)")
ck("PM-1B4", Poly(Tsal, t).degree() == 2 and irred(Tsal, t)
        and is_alg_integer(Tsal, t) and zero(Tsal.subs(t, tau0)),
   "tau_0 is ALGEBRAIC of degree 2: t^2-t-3 irreducible and MONIC INTEGER"
   " (=> tau_0 an algebraic integer) and tau_0 is its root")
ck("PM-1B5", discriminant(Poly(Tsal, t)) == 13 and totally_real(Tsal, t),
   "tau_0 TOTALLY REAL: disc(t^2-t-3)=13>0 and both roots real (Sturm count=2)"
   " (salem-slot Cor 2.4: trace-downs are totally real)")
ck("PM-1B5f", discriminant(Poly(t**2 + 1, t)) == -4
        and not totally_real(t**2 + 1, t) and not totally_real(t**2 + t + 1, t),
   "FALSIFIER (totally-real teeth): disc(t^2+1)=-4<0 and t^2+1, t^2+t+1 are"
   " correctly flagged NOT totally real (real-root count 0 != 2)")
conj0 = (1 - s13)/2                         # the captured conjugate
ck("PM-1B6", sgnQd(tau0 - (-conj0), 13) == 1 and sgnQd(tau0 - 1, 13) == 1,
   "tau_0 PERRON: real, tau_0>1, and tau_0 - |conjugate| = (1+sqrt13)/2 -"
   " (sqrt13-1)/2 = 1 > 0 => strictly dominant root (exact Q(sqrt13) sign)")
ck("PM-1B6f", sgnQd((-conj0) - tau0, 13) == -1 and sgnQd(conj0, 13) == -1,
   "FALSIFIER (Perron teeth): the CAPTURED conjugate (1-sqrt13)/2 < 0 is NOT"
   " dominant (|conj|-tau_0 < 0) -- the machinery rejects the wrong root")
ck("PM-1B7", sgnQd(tau0 - 2, 13) == 1 and sgnQd(2 - phi, 5) == 1,
   "tau_0 > 2 > phi (exact signs): tau_0-2=(sqrt13-3)/2>0 since 13>9;"
   " 2-phi=(3-sqrt5)/2>0 since 9>5 (salem-slot Thm 3.8)")
# general type-fact, not special to beta_4: two more Salem instances + the
# accumulation points, all ALGEBRAIC (golden limit + a Pisot accumulation)
Leh = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1   # Lehmer
TLeh = t**5 + t**4 - 5*t**3 - 5*t**2 + 4*t + 3
P6  = x**6 - x**5 - x**3 - x + 1                                 # deg-6 Salem
T6  = t**3 - t**2 - 3*t + 1
okgen = (zero(expand(x**5*TLeh.subs(t, x + 1/x) - Leh)) and irred(TLeh, t)
         and totally_real(TLeh, t) and is_alg_integer(TLeh, t))
okgen = okgen and (zero(expand(x**3*T6.subs(t, x + 1/x) - P6)) and irred(T6, t)
         and totally_real(T6, t) and is_alg_integer(T6, t))
ck("PM-1B8", okgen,
   "GENERAL type-fact (not special to beta_4): Lehmer (deg-10) trace-down"
   " t^5+t^4-5t^3-5t^2+4t+3 and the deg-6 Salem trace-down t^3-t^2-3t+1 are"
   " each irreducible, MONIC INTEGER, totally real => every Salem redirection"
   " tau_0=beta+1/beta is an algebraic integer, totally real [salem-slot Cor 2.4]")
mp_sqrt5   = minimal_polynomial(s5, t)
mp_plastic = minimal_polynomial(muS + 1/muS, t)
ck("PM-1B9", expand(mp_sqrt5 - (t**2 - 5)) == 0 and totally_real(t**2 - 5, t)
        and expand(mp_plastic - (t**3 + t**2 - 4*t - 5)) == 0
        and is_alg_integer(t**3 + t**2 - 4*t - 5, t),
   "the LIMIT and accumulation points are algebraic too: golden limit sqrt5="
   "phi+1/phi (minpoly t^2-5, totally real); plastic accumulation mu_P+1/mu_P"
   " (minpoly t^3+t^2-4t-5) [salem-slot Thm 6.2, Prop 8.1] -- the whole"
   " redirection family + accumulation set lies in Qbar")

# ================================================================ PM-1C
print("== PM-1C: arithmetic TYPE of L_res = elliptic-period value (forced red.) ==")
m_mod = Kc**2/(1 + Kc**2)                    # E-parameter m = kappa^2/(1+kappa^2)
red = 1 + Kc**2*cos(th)**2 - (1 + Kc**2)*(1 - m_mod*sin(th)**2)
ck("PM-1C1", zero(red),
   "exact reduction (symbolic in K,theta): 1 + K^2 cos^2(th) ="
   " (1+K^2)(1 - m sin^2(th)) with m=K^2/(1+K^2); hence L_res ="
   " int_0^{pi/2} sqrt(1+kappa^2 cos^2 th) dth = sqrt(1+kappa^2) E(m)"
   " (v1.7 prop:sphlength, boxed) -- quarter-perimeter of an ellipse, axes"
   " sqrt(1+kappa^2), 1")
ck("PM-1C1f", not zero(1 + Kc**2*cos(th)**2
                       - (1 + Kc**2)*(1 - (m_mod + Rational(1, 1000))*sin(th)**2)),
   "FALSIFIER: perturbing the modulus m -> m+1/1000 breaks the reduction"
   " identity (rejected) -- the identity has teeth")
ck("PM-1C2", zero(Kc**2 - m_mod/(1 - m_mod)),
   "modulus m is TRANSCENDENTAL: the Mobius inverse kappa^2 = m/(1-m) is exact,"
   " so m in Qbar IFF kappa^2 in Qbar IFF kappa in Qbar; kappa NOT in Qbar"
   " (T3) => m transcendental [forced given corpus/G-S]")
ck("PM-1C2f", not zero(Kc**2 - (m_mod + Rational(1, 1000))/(1 - (m_mod + Rational(1, 1000)))),
   "FALSIFIER: the Mobius inverse fails for a perturbed modulus (rejected)")
print("   STATUS: L_res is the value of the complete elliptic integral E (2nd")
print("   kind) at parameter m.  E at an ALGEBRAIC modulus is a Kontsevich-")
print("   Zagier period (perimeter of an algebraic ellipse); HERE the modulus")
print("   m (equiv. the axis sqrt(1+kappa^2)) is TRANSCENDENTAL, so L_res is the")
print("   value of the elliptic period FUNCTION at a transcendental argument --")
print("   one step beyond the pure period ring.  Its transcendence is NOT")
print("   settled by G-S (alpha^beta) nor by Schneider-elliptic (needs algebraic")
print("   modulus): EXTERNAL-OPEN (see PM-2).")

# ================================================================ PM-1X
print("== PM-1X: kappa vs tau_0 -- NO algebraic relation over Q (forced) ==")
mY = Y**2 - Y - 3                            # minpoly of tau_0, in Y
r0, r1 = (1 + s13)/2, (1 - s13)/2            # the two roots of mY
# PM-1X1 positive control: a GENUINE algebraic tie yields a nonzero rational
# polynomial witnessing algebraicity (here a = tau_0^2, tie P = X - Y^2).
Pdemo = X - Y**2
Rdemo = resultant(Pdemo, mY, Y)
ck("PM-1X1", expand(Rdemo - (X**2 - 7*X + 9)) == 0 and Rdemo != 0
        and zero(Rdemo.subs(X, tau0**2))
        and expand(minimal_polynomial(tau0**2, X) - (X**2 - 7*X + 9)) == 0,
   "positive control: for the tie P=X-Y^2 and target a=tau_0^2, R(X)="
   "Res_Y(P,minpoly(tau_0)) = X^2-7X+9 in Q[X], nonzero, R(a)=0, and R=minpoly(a)"
   " -- ANY genuine algebraic tie to tau_0 extracts a rational-coeff polynomial"
   " vanishing at the target => the target is algebraic")
# PM-1X2 falsifier-control: a VACUOUS tie (P divisible by the minpoly) gives a
# zero resultant -- the machinery distinguishes genuine from vacuous ties.
ck("PM-1X2", resultant(mY, mY, Y) == 0
        and resultant(3*mY, mY, Y) == 0,
   "FALSIFIER-CONTROL (teeth): a VACUOUS tie P=minpoly(tau_0) (or 3*minpoly)"
   " gives R identically 0 -- the resultant flags 'no constraint on X'; so"
   " nonzero-R <=> a genuine constraint on X (the case that forces algebraicity)")
# PM-1X3 the product-over-conjugates identity that makes 'P(kappa,tau_0)=0 =>
# R(kappa)=0' concrete:  R(X) = prod_{roots r of minpoly} P(X,r).
ck("PM-1X3", zero(expand(Rdemo - (X - r0**2)*(X - r1**2)))
        and zero(Pdemo.subs(Y, r0).subs(X, tau0**2)),
   "resultant factorization: R(X)=Res_Y(P,minpoly)=prod_r P(X,r) over the"
   " conjugates r of tau_0; here = (X-r0^2)(X-r1^2).  So if the r0=tau_0 factor"
   " P(X,tau_0) vanishes at X=kappa then R(kappa)=0 (the load-bearing step)")
# PM-1X4 the contradiction / tower step.
try:
    minimal_polynomial(kappa, X)
    kappa_alg = True
except NotAlgebraic:
    kappa_alg = False
ck("PM-1X4", Rdemo != 0 and zero(Rdemo.subs(X, tau0**2)) and (not kappa_alg),
   "tower step: SUPPOSE P in Q[X,Y] with P(kappa,tau_0)=0 and R:=Res_Y(P,minpoly)"
   " NOT identically 0 (genuine, PM-1X2).  Then R in Q[X]\\{0} and R(kappa)=0"
   " (PM-1X3) => kappa algebraic over Q (root of a nonzero rational polynomial)"
   " => kappa in Qbar.  But kappa NOT in Qbar (T3/G-S, PM-A2).  CONTRADICTION."
   " Hence no such P: NO polynomial relation over Q constrains kappa via tau_0."
   " [FORCED given established (tower law) + inherited transcendence (G-S)]")
# PM-1X5 planted false relation: 'kappa = tau_0' (P = X - Y).
Rp = resultant(X - Y, mY, Y)
ck("PM-1X5", expand(Rp - (X**2 - X - 3)) == 0 and Rp != 0,
   "PLANTED FALSE RELATION kappa=tau_0 (P=X-Y): R=X^2-X-3 nonzero; IF kappa"
   " were a root then kappa in Qbar (impossible).  kappa is NOT a root:"
   " kappa^2-kappa-3 != 0 -- certified > 3 by guard PM-G3, and structurally"
   " nonzero since kappa not in Qbar.  The planted tie is rejected")
# PM-1X6 restatement test / branch-blindness (mandatory, OP-RATE-adjacent).
ck("PM-1X6", all((4*kv + 1) != 0 and (4*kv + 1) % 2 == 1 for kv in range(-40, 41))
        and expand((pi/2 + 2*pi*kk)/log(phi) - (4*kk + 1)*kappa) == 0
        and zero((4*kk + 1)*kappa/((4*kk + 1)) - kappa),
   "RESTATEMENT TEST: kappa_k=(pi/2+2 pi k)/ln phi=(4k+1)kappa; (4k+1) is a"
   " NONZERO (odd) rational for every k, so kappa_k in Qbar IFF kappa in Qbar,"
   " hence kappa_k transcendental for all k (tx/T3) and PM-1X is VERBATIM"
   " identical per branch.  Branch-blind: this SELECTS no k -- an arithmetic"
   " CHARACTERIZATION/exclusion, NOT an OP-RATE derivation (HANDOFF sec 2)")

# ================================================================ PM-2
print("== PM-2: the precise EXTERNAL-OPEN map (bounded, citable) ==")
# PM-2A -- the SETTLED entries (each re-encoded / inherited)
ck("PM-2A", zero(I*pi/log(phi) - 2*I*kappa)          # kappa transc [PM-A2/T3]
        and irred(Tsal, t) and is_alg_integer(Tsal, t)  # tau_0 alg [PM-1B]
        and resultant(mY, mY, Y) == 0                   # kappa-vs-tau_0 [PM-1X]
        and Rp != 0,
   "SETTLED sub-ledger: (i) kappa transcendental [G-S, tx/T3]; (ii) tau_0"
   " algebraic deg-2 totally-real [PM-1B, forced]; (iii) (kappa,tau_0) satisfy"
   " NO constraining relation over Q [PM-1X, forced given established]")
# PM-2B -- forced REDUCTION: any (tau_0,L_res) relation forces L_res algebraic.
# Same machinery with tau_0 in the algebraic slot: a genuine tie between an
# unknown L and tau_0 yields a nonzero rational polynomial vanishing at L.
Ldemo = X - Y                                # stand-in tie P(L,tau_0)=L-tau_0
RL = resultant(Ldemo, mY, Y)
ck("PM-2B", expand(RL - (X**2 - X - 3)) == 0 and RL != 0
        and resultant(mY, mY, Y) == 0,
   "REDUCTION [forced]: tau_0 algebraic => any P(tau_0,L_res)=0 over Q with"
   " nonzero resultant R:=Res_Y(P,minpoly(tau_0)) gives R in Q[X]\\{0}, R(L_res)"
   "=0 => L_res algebraic.  So the (tau_0,L_res) pair-dependence question"
   " COLLAPSES to the single-constant 'is L_res algebraic' -- vacuous ties give"
   " R=0 (teeth), genuine ties force algebraicity")
ck("PM-2Bf", RL != 0 and zero(RL.subs(X, tau0)) and not zero((RL + 1).subs(X, tau0)),
   "FALSIFIER (teeth): the reduction resultant R=X^2-X-3 is nonzero and genuinely"
   " VANISHES at the algebraic target (R(tau_0)=0), while the perturbation R+1"
   " does NOT vanish there -- the vanishing is a real detection of algebraicity,"
   " not vacuous")
print("   SETTLED (forced): kappa transcendental (yes) | tau_0 algebraic (yes) |")
print("     (kappa,tau_0) no relation (no) | (tau_0,L_res) reduces to L_res")
print("     transcendence.")
print("   EXTERNAL-OPEN (bounded, citable), with the tool-gap reason each:")
print("     - L_res transcendental?  OPEN. G-S handles alpha^beta only;")
print("       Schneider/Chudnovsky elliptic-period theorems need an ALGEBRAIC")
print("       modulus/CM lattice, but m is transcendental (PM-1C2). No tool.")
print("     - (kappa,L_res) algebraically independent?  OPEN. Would follow from")
print("       Schanuel; NOT from G-S or Baker (linear forms in logs reach")
print("       neither elliptic-integral values nor algebraic independence).")
print("     - (pi,ln phi) algebraically independent?  OPEN (Schanuel-")
print("       conjectural).  NOTE the RATIO pi/ln phi=2 kappa IS transcendental")
print("       (settled), but independence of the PAIR is a strictly stronger,")
print("       open statement.")
print("   These are Kontsevich-Zagier period-relation questions, NOT reachable")
print("   by Gelfond-Schneider or Baker.  A bounded map of the open territory")
print("   is the deliverable; nothing here derives kappa.")

# ================================================================ PM-G guards
print("== PM-G: certified-interval + numeric-corroboration guards (separate) ==")
from mpmath import iv, mp
import mpmath as mnum
iv.dps = 60
iphi   = (1 + iv.sqrt(5))/2
ilnphi = iv.log(iphi)
ikappa = iv.pi/(2*ilnphi)
ikap2  = ikappa**2
itau0  = (1 + iv.sqrt(13))/2
im_mod = ikap2/(1 + ikap2)

gk("PM-G1", (ikap2 > iv.mpf(10)) and (ikap2 < iv.mpf(11)),
   "interval-certified (60 dps): 10 < kappa^2 < 11 [inherits MX-G1]")
gk("PM-G2", (ikappa - itau0) > iv.mpf('0.9'),
   "interval-certified: kappa - tau_0 > 0.9 (kappa != tau_0; the constants are"
   " distinct -- corroborates the negative settlement)")
gk("PM-G3", (ikap2 - ikappa - 3) > iv.mpf(3),
   "interval-certified: kappa^2 - kappa - 3 > 3 > 0 => kappa is NOT a root of"
   " tau_0's minimal polynomial t^2-t-3 (kappa not a Galois conjugate of tau_0)"
   " -- certifies the PM-1X5 planted-relation rejection")
gk("PM-G4", (im_mod > iv.mpf(10)/11) and (im_mod < iv.mpf(11)/12),
   "interval-certified: modulus m in (10/11, 11/12) subset (0,1) -- a proper"
   " (and transcendental, PM-1C2) elliptic modulus")
gk("PM-G5", (iv.sqrt(1 + ikap2) > iv.mpf('3.41')) and (iv.sqrt(1 + ikap2) < iv.mpf('3.42'))
        and ((ikappa + iv.pi/2) > iv.mpf('4.83')) and ((ikappa + iv.pi/2) < iv.mpf('4.84')),
   "interval-certified: the v1.7 exact bounds sqrt(1+kappa^2) in (3.41,3.42) and"
   " kappa+pi/2 in (4.83,4.84) bracket L_res=3.7312193125 (prop:sphlength"
   " bounds; corroborates the elliptic value's location)")

# integer-relation corroboration (mpmath pslq; tx/mx mp idiom, 140 dps).
# genuine relations resolve to residual < 1e-120; spurious Diophantine
# approximations do not.  Counted separately: corroboration, NOT proof.
mp.dps = 140
fphi = (1 + mnum.sqrt(5))/2
flnphi = mnum.log(fphi)
fkappa = mnum.pi/(2*flnphi)
ftau0  = (1 + mnum.sqrt(13))/2
fLres  = mnum.sqrt(1 + fkappa**2)*mnum.ellipe(fkappa**2/(1 + fkappa**2))
GEN = mnum.mpf(10)**-120
def genuine(vec):
    r = mnum.pslq(vec, maxcoeff=10**8, maxsteps=50000)
    if r is None:
        return None, False
    resid = abs(sum(mnum.mpf(a)*b for a, b in zip(r, vec)))
    return r, (resid < GEN)
one = mnum.mpf(1)
rt, gt = genuine([one, ftau0, ftau0**2])                 # teeth: tau_0 algebraic
_, g_k  = genuine([one, fkappa])                          # kappa combos: expect None
_, g_kt = genuine([one, fkappa, ftau0, fkappa*ftau0])
_, g_kl = genuine([one, fkappa, fLres, fkappa*fLres])
_, g_ktl = genuine([one, fkappa, ftau0, fLres])
gk("PM-G6", gt and (rt == [-3, -1, 1]) and (not g_k) and (not g_kt)
        and (not g_kl) and (not g_ktl),
   "integer-relation corroboration (pslq, 140 dps, height<=1e8): TEETH -- the"
   " genuine relation tau_0^2-tau_0-3=0 is detected ([-3,-1,1], resid<1e-120);"
   " every kappa-joint basis {1,kappa}, {1,k,t0,k*t0}, {1,k,L,k*L}, {1,k,t0,L}"
   " yields NO genuine relation (corroboration only, NOT a proof)")

# decimal audit vs register/paper literals (tx audit idiom).
mp.dps = 60
def audit_ok(val, lit):
    dp = len(lit.split('.')[1])
    litv = mnum.mpf(lit)
    ulps = abs(val - litv)*mnum.mpf(10)**dp
    if ulps <= mnum.mpf('0.51'):
        return True
    trunc = mnum.floor(val*mnum.mpf(10)**dp)/mnum.mpf(10)**dp
    return abs(trunc - litv) < mnum.mpf(10)**(-dp - 4)
fphi = (1 + mnum.sqrt(5))/2
flnphi = mnum.log(fphi)
fkappa = mnum.pi/(2*flnphi)
fLres = mnum.sqrt(1 + fkappa**2)*mnum.ellipe(fkappa**2/(1 + fkappa**2))
fLint = mnum.quad(lambda th: mnum.sqrt(1 + fkappa**2*mnum.cos(th)**2), [0, mnum.pi/2])
audits = [
    ("kappa",     fkappa,                 '3.2642513026'),
    ("L_res(E)",  fLres,                  '3.7312193125'),
    ("L_res(int)", fLint,                 '3.7312193125'),
    ("tau_0",     (1 + mnum.sqrt(13))/2,  '2.3027756377'),
    ("sqrt5",     mnum.sqrt(5),           '2.2360679775'),
    ("m_mod",     fkappa**2/(1 + fkappa**2), '0.9142023918'),
    ("plastic_tr", (lambda u: u + 1/u)(mnum.findroot(lambda z: z**3 - z - 1, 1.3)),
                  '2.0795956235'),
]
okA = True
for name, val, lit in audits:
    hit = audit_ok(val, lit)
    okA = okA and hit
    print("   audit %-11s printed %s  recomputed %s  %s"
          % (name, lit, mnum.nstr(val, 14), "ok" if hit else "MISMATCH"))
gk("PM-G7", okA,
   "decimal audit at 60 dps: kappa, L_res (closed-form E and direct quadrature"
   " agree), tau_0=(1+sqrt13)/2, sqrt5, modulus m, plastic-trace all match the"
   " register/paper literals (0.51 ulp or floor-truncation)")

# ================================================================ summary
print("=" * 68)
n_guard_fail = sum(1 for _, ok, _ in GUARD if not ok)
n_exact_fail = len(FAIL) - n_guard_fail
n_exact_pass = len(PASS)
print("PM EXACT: %d passed, %d failed | GUARDS: %d/%d certified"
      % (n_exact_pass, n_exact_fail,
         sum(1 for _, ok, _ in GUARD if ok), len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
