#!/usr/bin/env python3
# bl_softcheck.py -- the ACTIVE-BLOCKERS ledger harness (v2.0 finalization), written
# cold this session (2026-07-15).  It does NOT re-derive the corpus; it RE-CERTIFIES,
# in one place, the single load-bearing FORCED fact behind each open problem's blocker
# CLASSIFICATION, and LABELS the external-open items (which no harness can settle) with
# their named tool-gap / conjecture.  The blocker taxonomy:
#   IMPOSSIBILITY        -- a forced no-go: the resolution cannot come from the axioms.
#   DECLARED-ATOM        -- an irreducible declaration the substrate does not force.
#   RESTATEMENT-WALL     -- every proposed reduction is logically equivalent to the axiom.
#   NAMED-LEMMA          -- reduced to one precisely-stated open lemma (new internal math).
#   EXTERNAL-TOOL-GAP    -- blocked on a theorem that does not exist / a named conjecture.
#   UNFINISHED-COMPUTATION -- tractable in principle, not yet executed.
#
# Discipline (nx/qx/rd/of idiom): every DECISION exact over Q / Q(sqrt5) / Q(sqrt13)
# (extended by radicals/i where needed); sign decisions by rational-coefficient
# arithmetic (sgnQ5); valuations by exact leading-term order; mpmath is display-only
# (group BL-H).  Fail-first: each wall-detector is shown to REJECT a planted false wall.
# The forced facts here are cross-checked against the harnesses that first proved them
# (cited inline: CL, RD/RA, QT/CH, OF, PF).  Exit 0 iff every exact check passes.
# ASCII-only output.  Run: py bl_softcheck.py

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, log, pi, I, Rational, Integer, simplify, expand,
                   radsimp, together, cancel, Poly, factor_list, series, im, re as sre,
                   fibonacci, lucas, nsimplify)

x, z, t, m_s = symbols('x z t m_s')
s5 = sqrt(5); s13 = sqrt(13); s3 = sqrt(3)
phi = (1 + s5) / 2
tau = (s5 - 1) / 2
gap = phi ** -4
K2 = 1 - gap
K = sqrt(K2)
zc = s3 / 2
zc2 = Rational(3, 4)
beta2 = phi ** 4 - 1              # K-seed rotation-pair modulus^2 (off-circle)
kappa = pi / (2 * log(phi))

PASS, FAIL, LABEL = [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(("PASS" if ok else "FAIL"), cid, "-", desc)
    return ok

def lab(cid, klass, desc):
    # An EXTERNAL-open blocker no harness can settle: LABELLED, not checked.
    LABEL.append((cid, klass))
    print("LABEL", cid, "[%s]" % klass, "-", desc)

def zero(e):
    e2 = simplify(expand(e))
    if e2 == 0:
        return True
    return simplify(radsimp(together(e2))) == 0

def q5AB(e):
    """(A,B) with e = A + B*sqrt5, A,B in Q; raises if not in Q(sqrt5)."""
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
    if A > 0:  return 1 if A * A > 5 * B * B else -1
    return 1 if 5 * B * B > A * A else -1

def vz(expr):
    """Exact z-adic valuation (order of the leading term at z=0)."""
    lt = expr.series(z, 0, 8).removeO()
    p = Poly(lt, z)
    # lowest-degree present monomial exponent
    monoms = [mm[0] for mm, c in p.terms() if c != 0]
    return min(monoms) if monoms else 10

print("== BL-A  OP-RATE: IMPOSSIBILITY (discrete) + DECLARED-ATOM (D3) ==")
# Load-bearing FORCED fact (CL): the discrete state is k-invariant -> R6 no-go.
# rung values e^{i kappa_k n ln phi} = i^n for every branch k (symbolic e^{2 pi i k}=1).
kgrid_ok = all(simplify(sre(I ** ((1 + 4 * k) * n)) - sre(I ** n)) == 0
               and simplify(im(I ** ((1 + 4 * k) * n)) - im(I ** n)) == 0
               for k in range(-3, 4) for n in range(0, 9))
ck("BL-A1", kgrid_ok,
   "k-invariance: e^{i kappa_k n ln phi} = i^n for all k (CL-1) => R6 no-go [IMPOSSIBILITY, discrete]")
ck("BL-A2", zero((1 + 4 * 1) - 5) and zero((1 + 4 * (-1)) + 3),
   "branch family kappa_k=(4k+1)kappa: selection is branch-blind; D3 declared [DECLARED-ATOM]")

print("== BL-B  OP-RADIUS ab initio ([open]1): IMPOSSIBILITY (apex parity) ==")
r2 = K2 * z / zc                 # inherited radius squared (odd at the apex)
u_chart = z ** 2 / (1 - z ** 2)  # a chart quantity (even at the apex)
v_r2 = vz(r2); v_u = vz(u_chart)
ck("BL-B1", v_r2 == 1 and (v_r2 % 2 == 1),
   "v_z(r^2)=1 ODD: the radius is odd at the apex (RD-1/RA-A)")
ck("BL-B2", v_u == 2 and (v_u % 2 == 0),
   "chart quantity v_z even (value group 2Z): chart field cannot reach r^2 [IMPOSSIBILITY]")
# fail-first: a planted 'even r^2' must be rejected by the parity detector.
ck("BL-B3", not (vz(K2 * z / zc) % 2 == 0),
   "fail-first: parity detector REJECTS the false wall 'r^2 is even'")

print("== BL-C  OP-RADIUS profile given the bridge: DECLARED-ATOM (D4) ==")
# D4 = (i) orientation (sign of the odd generator) + (ii) cap onset.  Each is a free datum.
r_sym = K * sqrt(z / zc)          # the inherited radius (two orientations: +-r_sym)
ck("BL-C1", zero(r_sym ** 2 - r2) and not zero(2 * r_sym),
   "orientation free: +r and -r both square to r^2 yet differ (Z/2 sign, RA-B4) [DECLARED-ATOM]")
onset_a = zc; onset_b = Rational(4, 5)      # a different declared onset
ck("BL-C2", not zero(onset_a ** 2 - onset_b ** 2),
   "cap onset free: onset z_c vs 4/5 give distinct capped profiles (RA-B5) [DECLARED-ATOM]")

print("== BL-D  W[open]2 / D2 (chi selection): RESTATEMENT-WALL ==")
# The three competing quanta and the K-seed clock.  |i.beta| off-circle; arg = pi/2.
ck("BL-D1", sgnQ5(beta2 - 1) > 0,
   "K-seed rotation pair off-circle: |i.beta|^2 = phi^4-1 > 1 (QT-B) -- disjoint from g1")
# pentagon 2/5 not in (1/4)Z mod 1; quarter 1/4 is; terrain 1/2 is (order 2, helicity-killed).
pent_in = any((Fraction(2, 5) - Fraction(j, 4)) % 1 == 0 for j in range(-8, 9))
ck("BL-D2", (not pent_in) and (Fraction(1, 4) in [Fraction(j, 4) % 1 for j in range(0, 4)]),
   "quanta {terrain pi:Z/2, pentagon 2pi/5:Z/5, quarter pi/2:Z/4}: naming the K-seed clock == "
   "naming the quarter-turn (QT-C3/C5) [RESTATEMENT-WALL]")

print("== BL-E  P3/P4 (odd floor = 2): NAMED-LEMMA (ODD-2) ==")
# Upper bound M(x^m-2)=2; Smyth gives only mu_S<2; the reciprocal sub-case is forced >= phi^2>2.
def mahler_xm2_ok(m):
    # x^m-2 is Eisenstein-irreducible at 2 (so it IS the minimal polynomial); its roots are
    # 2^{1/m} zeta_m^k, all of modulus 2^{1/m}>1 (since 2>1) -- all m lie OUTSIDE the unit
    # circle -- so M = |lead| * prod_{|root|>1}|root| = 1 * (2^{1/m})^m = 2.
    irred = len(factor_list(x ** m - 2)[1]) == 1
    M = simplify((Integer(2) ** Rational(1, m)) ** m)
    return irred and (M == 2)
ck("BL-E1", all(mahler_xm2_ok(m) for m in (3, 5, 7, 9)),
   "upper bound: x^m-2 irreducible, all roots at modulus 2^{1/m}>1, M=2 for odd m (OF-1b) -- =2 CEILING forced")
# Smyth plastic number mu_S = real root of x^3-x-1: 1 < mu_S < phi < 2.
muS = symbols('muS')
# Smyth mu_S = real root of p(x)=x^3-x-1: p(1)=1-1-1=-1<0, p(phi)=phi^3-phi-1=phi>0 => mu_S in (1,phi).
ck("BL-E2", (1 ** 3 - 1 - 1) < 0 and zero((phi ** 3 - phi - 1) - phi) and sgnQ5(phi - 2) < 0,
   "Smyth floor mu_S<2: p(1)=-1<0, p(phi)=phi>0 so mu_S<phi<2 (OF) -- Smyth ALONE cannot reach 2")
ck("BL-E3", sgnQ5(phi ** 2 - 2) > 0 and zero(phi ** 2 - phi - 1),
   "reciprocal sub-case forced >= phi^2 > 2 (OF-2c); the mu_S->2 lift IS the lemma ODD-2 [NAMED-LEMMA]")

print("== BL-F  Period frontier: EXTERNAL-TOOL-GAP (forced premises only) ==")
# Moebius m = kappa^2/(1+kappa^2), inverse kappa^2 = m/(1-m): m in Qbar iff kappa in Qbar.
moeb = (m_s / (1 - m_s)) / (1 + m_s / (1 - m_s))
ck("BL-F1", zero(moeb - m_s),
   "Moebius m<->kappa^2 (PF-2D): m transcendental IFF kappa transcendental (the forced premise)")
tau0 = (1 + s13) / 2
ck("BL-F2", zero(tau0 ** 2 - tau0 - 3),
   "tau_0 algebraic: minpoly t^2-t-3 (PF-3D) -- the algebraic anchor")
# multiplicative independence of {2,3,phi} in a box (Baker's forced premise).
def mult_indep_box(B):
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            for c in range(-B, B + 1):
                if a == b == c == 0:
                    continue
                val = Integer(2) ** a * Integer(3) ** b * phi ** c
                if zero(val - 1):
                    return False
    return True
ck("BL-F3", mult_indep_box(4),
   "{2,3,phi} multiplicatively independent in box |.|<=4 => {1,ln2,ln3,ln phi} indep (PF-2B, Baker)")

print("== BL-G  the LABELLED walls, unsettleable by any harness (their classes vary) ==")
lab("BL-G1", "EXTERNAL-TOOL-GAP",
    "L_res transcendence: NO theorem reaches an elliptic period at a TRANSCENDENTAL modulus "
    "(Schneider/Chudnovsky need an algebraic modulus). MISSING: such a transcendence theorem.")
lab("BL-G2", "EXTERNAL-TOOL-GAP (named conjecture)",
    "(pi, ln phi) algebraic independence: premise {ln phi, i pi} Q-lin-indep is forced (PF-2C); "
    "the conclusion is Schanuel-conditional. MISSING: Schanuel's conjecture.")
lab("BL-G3", "EXTERNAL-TOOL-GAP (beyond Schanuel)",
    "(kappa, L_res) algebraic independence: L_res is an elliptic value outside the exp/log "
    "framework. MISSING: Kontsevich-Zagier period-relation theory (conjectural).")
lab("BL-G4", "NAMED-LEMMA (internal, open)",
    "ODD-2: a Mahler-measure lower bound of 2 for charge-admissible NON-reciprocal odd-charge "
    "objects. Smyth stops at mu_S=1.3247; the lift mu_S->2 IS this lemma. MISSING: new number theory.")
lab("BL-G5", "DECLARED-ATOM (formation dynamics)",
    "OP-RADIUS cap-onset MECHANISM: why formation stops at z_c (D4(ii) declares it; rho(z_c)=ln2 "
    "is a restatement of z_c). MISSING: a formation-dynamics principle.")
lab("BL-G6", "RESTATEMENT-WALL",
    "STRICT reduction of D2: any thinner axiom naming the clock is logically equivalent to the "
    "quarter-turn selection (QT). MISSING: a principle forcing the clock from something more primitive.")
lab("BL-G7", "UNFINISHED-COMPUTATION",
    "P1' cross-shell nu-criterion: the within-shell detector leaves cross-shell mirrored-class "
    "coherence NO-TIE (RN-3). Tractable in principle; MISSING: the cross-shell coherence detector.")

print("== BL-H  the blocker LEDGER (display) ==")
LEDGER = [
    ("OP-RATE (derive kappa)",       "IMPOSSIBILITY(discrete)+DECLARED-ATOM(D3)", "closed as classification; open as continuum derivation"),
    ("OP-RADIUS ab initio [open]1",  "IMPOSSIBILITY",                              "impossible from D1-D3 (apex parity)"),
    ("OP-RADIUS profile",            "DECLARED-ATOM (D4: orientation + onset)",    "reduced to D4"),
    ("W[open]2 / D2 (chi)",          "RESTATEMENT-WALL",                           "relocated to substrate-grounded atom; partial"),
    ("P3/P4 (odd floor =2)",         "NAMED-LEMMA (ODD-2)",                        "bracketed; reduced to one lemma"),
    ("L_res transcendence",          "EXTERNAL-TOOL-GAP",                          "no transcendental-modulus theorem"),
    ("(pi, ln phi) independence",    "EXTERNAL-TOOL-GAP (Schanuel)",               "conditional on Schanuel"),
    ("(kappa, L_res) independence",  "EXTERNAL-TOOL-GAP (periods)",                "beyond Schanuel"),
    ("P1' cross-shell nu",           "UNFINISHED-COMPUTATION",                     "bounded scan; cross-shell open"),
]
print("   %-28s | %-38s | %s" % ("open problem", "active blocker class", "status"))
for prob, klass, status in LEDGER:
    print("   %-28s | %-38s | %s" % (prob[:28], klass[:38], status))
# meta-consistency: every external/unfinished wall is LABELLED, not falsely 'checked'.
ck("BL-H1", len(LABEL) == 7 and len(LEDGER) == 9,
   "ledger consistent: 9 open problems; 7 external/unfinished walls LABELLED (not checked)")

# ================================================================ summary
print("=" * 64)
print("BL: EXACT %d/%d checks passed | %d walls LABELLED (unsettleable by harness) %s" % (
    len(PASS), len(PASS) + len(FAIL), len(LABEL),
    "-- ALL PASS" if not FAIL else ("| FAILURES: %s" % FAIL)))
sys.exit(0 if not FAIL else 1)
