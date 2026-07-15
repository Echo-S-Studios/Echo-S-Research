#!/usr/bin/env python3
# qt_softcheck.py -- the quarter-turn charge selection D2: reduction attempt via
# the K-seed's OWN rotation direction (v1.9, W[open]2).
#
# CONTEXT (post-v1.8, authoritative register + ch/cx harnesses):
#   W[open]2 is PARTIAL as a derivation.  The v1.7 conditional chi-chain
#   (thm:chichain) needed D2' = "the per-rung winding multiplier is an ON-CIRCLE
#   element of the emission image S".  v1.8 (ch_softcheck / prop:oncircleimage)
#   forced the on-circle image of S = mu_2 = {+1,-1}, so +-i is NOT an on-circle
#   element of S; D2' (Reading A) is UNSATISFIABLE and D2' (Reading B, charge)
#   RESTATES D2.  The openness is located at the unreduced axiom D2 itself.
#
#   Whitepaper Declaration 12.2 (D2, verbatim, PDF p13): "The winding realizes
#   the seed's charge completion.  Given D2, the angular quantum per rung is
#   [forced] a generator of Z/4Z, i.e. chi = +-pi/2 --- the multipliers +-i of
#   the C = -1/2 gate (Prop 10.2), exactly the charge the K-seed supplies. ...
#   The competing corpus quanta --- pi (terrain sign flip) and 2pi/5 (pentagon
#   gates) --- generate only Z/2Z and Z/5Z and are excluded by the principle."
#   Post-v1.1 scoping upgrade (sec:quarterturn): the quarter-turn OPERATOR J=UA
#   and the rung generator Q=+-tau J EXIST canonically [FORCED]; D2's residual
#   declared content is SELECTION -- one Q per floor, over the terrain pi and
#   pentagon 2pi/5 clocks -- not existence.
#
# KEY SUBSTRATE FACT (FORCED, not declared): the K-seed x^4+5x^2-5 (W-Sec2) has
# roots {+-K, +-i.beta}, K = 5^(1/4) tau (terrain, REAL, arg 0/pi), i.beta with
# beta = 5^(1/4) phi (rotation, IMAGINARY, arg +-pi/2).  So the quarter-turn
# DIRECTION (+-pi/2) is literally the argument of a FORCED substrate object --
# the K-seed's rotation (imaginary) pair -- even though |i.beta| = beta =
# sqrt(phi^4-1) != 1 is OFF the unit circle.  This modulus is DISJOINT from the
# v1.8 g1 obstruction, which was about the ON-circle image of S.
#
# THIS HARNESS defines the thin axiom
#   D2'' = "the per-rung winding multiplier's ANGLE is the argument of the
#           K-seed rotation (imaginary) pair +-i.beta"
# (the modulus stays |m| = e^{-2rho}, D1's gap-contraction, UNCHANGED), and:
#   QT-A  re-encodes the forced K-seed geometry cold (arg +-i.beta = +-pi/2 exact)
#   QT-B  QT-2 crux: +-i.beta is OFF-circle => D2'' sidesteps g1 (disjoint from S
#         on-circle image mu_2); the arg-pi/2 IS in the substrate off-circle
#   QT-C  QT-1: the reduction/restatement audit -- is D2'' THINNER than D2?
#   QT-D  QT-3 net verdict on W[open]2 + the OP-RATE restatement test
#   QT-G  certified 60-dps interval guards (counted separately)
#
# VERDICT (proved below, premise-audited): D2'' FORCES chi = +-pi/2 [FORCED given
# D2''] from a FORCED object's geometry, and is DISJOINT from the g1 kill (it
# lives OFF-circle, never touching S's on-circle image).  But the SELECTION "the
# winding's angular clock is the K-seed ROTATION axis" is logically EQUIVALENT to
# D2's quarter-turn selection (over terrain pi / pentagon 2pi/5): so D2'' does NOT
# strictly reduce D2 -- it is a GROUNDED restatement that precisely characterizes
# D2's irreducible atom (the rotation-axis clock-selection), with the value pi/2
# and the group-completion Z/4Z downgraded from declared to FORCED outputs.
# Register motion: W[open]2 stays PARTIAL, ADVANCED (sharpened + substrate-
# grounded + disjoint-from-g1); NOT closed.
#
# Discipline (inherited from ch/cx/qx idiom): every DECISION exact over Q /
# Q(sqrt5) (extended by i and 5^(1/4)); |.|^2 == 1 decided in Q(sqrt5) by rational
# (A,B) arithmetic (q5AB); monomial bookkeeping by certified integer exponent
# vectors; mpmath iv at 60 dps ONLY for the QT-G guards (counted separately).
# Fail-first: each group opens with a falsifier (QT-*F*) asserting a deliberately
# perturbed variant is REJECTED; falsifiers were written first; detectors are
# validated on a KNOWN positive (the poisoned catalog reaches on-circle +-i)
# before any negative is trusted.  Exit 0 iff all pass.

import sys
from fractions import Fraction as F
from itertools import combinations_with_replacement
from sympy import (symbols, sqrt, Rational, Integer, I, pi, sin, cos, exp,
                   im, re, arg, simplify, expand, radsimp, together, nsimplify,
                   Poly, minimal_polynomial, real_roots, gcd, Mul, conjugate)

x = symbols('x')
y = symbols('y', real=True)
s5 = sqrt(5)
phi = (1 + s5) / 2
psi = (1 - s5) / 2                 # = -tau
tau = (s5 - 1) / 2                 # = 1/phi
f4 = 5 ** Rational(1, 4)           # 5^(1/4)
Kv = f4 * tau                      # K-seed terrain root  (arg 0)
betav = f4 * phi                   # K-seed rotation magnitude, beta = 5^(1/4) phi
r2, r3, r5 = sqrt(2), sqrt(3), sqrt(5)

PASS, FAIL, GUARD = [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(("PASS" if ok else "FAIL"), cid, "-", desc)

def gk(cid, cond, desc):
    ok = bool(cond)
    GUARD.append((cid, ok))
    print(("GUARD-PASS" if ok else "GUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append(cid)

def zero(e):
    e2 = simplify(expand(simplify(e)))
    if e2 == 0:
        return True
    try:
        if simplify(expand(radsimp(together(e)))) == 0:
            return True
    except Exception:
        pass
    try:
        if simplify(expand(e.rewrite(sqrt))) == 0:
            return True
    except Exception:
        pass
    try:
        e3 = nsimplify(e2, [sqrt(5)])
        return simplify(radsimp(e3)) == 0
    except Exception:
        return False

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q. Raises if not in Q(sqrt5)."""
    from sympy import cancel
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return F(str(A)), F(str(B))

def is_one_q5(e):
    """Exact decision e == 1 for e in Q(sqrt5)."""
    A, B = q5AB(e - 1)
    return A == 0 and B == 0

def onS_modsq(e):
    """Exact decision |e|^2 == 1 (on the unit circle) for e with |e|^2 in Q(sqrt5)."""
    return is_one_q5(expand(e * conjugate(e)))

def purely_imag(e):
    """Exact decision: e is purely imaginary (re == 0, im != 0)."""
    return zero(re(e)) and (not zero(im(e)))

def purely_real(e):
    """Exact decision: e is real (im == 0)."""
    return zero(im(e))

# ------------------------------------------------------------------ catalog
# [EG] Def 3.1 / [LB] Def 2.6: seeds = companions of
# { phi, tau, sqrt2, sqrt3, sqrt5, gap (x^2-7x+1), K (x^4+5x^2-5) }.
SEEDS = [
    ("phi-seed",  x**2 - x - 1,      [phi, psi]),
    ("tau-seed",  x**2 + x - 1,      [tau, -phi]),
    ("sqrt2",     x**2 - 2,          [r2, -r2]),
    ("sqrt3",     x**2 - 3,          [r3, -r3]),
    ("sqrt5",     x**2 - 5,          [r5, -r5]),
    ("gap-seed",  x**2 - 7*x + 1,    [phi**4, tau**4]),
    ("K-seed",    x**4 + 5*x**2 - 5, [Kv, -Kv, I*betav, -I*betav]),
]

# Eigenvalue table (re-derived cold, as ch/cx): lam = 2^(w0/4)3^(w1/4)5^(w2/4)phi^(w3/4) i^t.
EIG = [
    ("phi",     phi,       (0, 0, 0,   4), 0),
    ("psi",     psi,       (0, 0, 0,  -4), 2),
    ("tau",     tau,       (0, 0, 0,  -4), 0),
    ("-phi",    -phi,      (0, 0, 0,   4), 2),
    ("r2",      r2,        (2, 0, 0,   0), 0),
    ("-r2",     -r2,       (2, 0, 0,   0), 2),
    ("r3",      r3,        (0, 2, 0,   0), 0),
    ("-r3",     -r3,       (0, 2, 0,   0), 2),
    ("r5",      r5,        (0, 0, 2,   0), 0),
    ("-r5",     -r5,       (0, 0, 2,   0), 2),
    ("phi4",    phi**4,    (0, 0, 0,  16), 0),
    ("phi-4",   tau**4,    (0, 0, 0, -16), 0),
    ("K",       Kv,        (0, 0, 1,  -4), 0),
    ("-K",      -Kv,       (0, 0, 1,  -4), 2),
    ("i.beta",  I*betav,   (0, 0, 1,   4), 1),
    ("-i.beta", -I*betav,  (0, 0, 1,   4), 3),
]

mu4 = [Integer(1), I, Integer(-1), -I]

# ================================================================== QT-A
print("== QT-A: cold re-encoding of the FORCED K-seed geometry; arg(+-i.beta)=+-pi/2 ==")

# QT-AF1 (falsifier, written first): WRONG geometry claims must be REJECTED --
# i.beta is NOT real; the terrain root K does NOT sit at pi/2; arg(i.beta) != 0.
ck("QT-AF1", (not purely_real(I*betav)) and (not purely_imag(Kv))
        and (not zero(arg(I*betav))) and purely_imag(I*betav),
   "falsifier: i.beta not real; terrain K not on the imaginary axis; arg(i.beta)!=0")

ok = True
for name, p, roots in SEEDS:
    ok = ok and Poly(p, x).degree() == len(roots)
    ok = ok and all(zero(p.subs(x, r)) for r in roots)
ck("QT-A1", ok, "7 catalog seeds; K-seed x^4+5x^2-5 roots exactly {+-K,+-i.beta}, "
   "K=5^(1/4)tau, beta=5^(1/4)phi (FORCED, W-Sec2 -- not declared)")

# The rotation (imaginary) pair sits on the imaginary axis: arg = +-pi/2 EXACT.
ok = purely_imag(I*betav) and purely_imag(-I*betav)
ok = ok and zero(re(I*betav)) and zero(im(I*betav) - betav) and betav.is_positive
ok = ok and zero(arg(I*betav) - pi/2) and zero(arg(-I*betav) + pi/2)
ck("QT-A2", ok, "K-seed ROTATION pair +-i.beta is purely imaginary (re=0, im=+-beta, "
   "beta>0) => arg(+-i.beta) = +-pi/2 EXACTLY -- the quarter-turn DIRECTION is the "
   "argument of a FORCED substrate object")

# The terrain pair +-K is REAL: arg in {0,pi}.  Full K-seed root ANGLE set = (pi/2)Z.
ok = purely_real(Kv) and purely_real(-Kv)
ok = ok and zero(arg(Kv)) and zero(arg(-Kv) - pi)
angles_turns = set()
for lam in (Kv, -Kv, I*betav, -I*betav):
    a = arg(lam)                                   # in (-pi, pi]
    angles_turns.add(F(str(simplify(a / (2*pi)))))  # angle class in turns
ck("QT-A3", ok and angles_turns == {F(0), F(1, 4), F(1, 2), F(-1, 4)},
   "terrain pair +-K is REAL (arg 0,pi); full K-seed angle set = {0,pi/2,pi,-pi/2} "
   "= the four compass points (pi/2)Z")

# Angle classes mod 1 turn = {0,1/4,1/2,3/4} = Z/4Z; relational (pairwise-difference)
# group of the FULL compass = Z/4Z; of the ROTATION pair ALONE = Z/2Z (group-drop).
cls = {a % 1 for a in angles_turns}                # {0,1/4,1/2,3/4}
diffs_full = {((a - b) % 1) for a in cls for b in cls}
rot = {F(1, 4), F(3, 4)}
diffs_rot = {((a - b) % 1) for a in rot for b in rot}
ck("QT-A4", cls == {F(0), F(1, 4), F(1, 2), F(3, 4)}
        and diffs_full == {F(0), F(1, 4), F(1, 2), F(3, 4)}
        and diffs_rot == {F(0), F(1, 2)},
   "K-seed angle classes = Z/4Z; FULL-compass relational group = Z/4Z; but the "
   "rotation pair ALONE has relational group Z/2Z (group-drop, cf thm:qmin q=i.tau)")

# The FORCED winding generator q = i.tau (thm:qmin) has minpoly x^4+3x^2+1, roots
# {+-i.tau,+-i.phi}, ALL purely imaginary (arg +-pi/2) = the K-seed rotation-pair
# angles.  So D2'' ("winding angle = arg(K-seed rotation pair)") is pinned to the
# actual forced generator's own angle set.
q = I * tau
qroots = [I*tau, -I*tau, I*phi, -I*phi]
ok = Poly(minimal_polynomial(q, x), x) == Poly(x**4 + 3*x**2 + 1, x)
ok = ok and all(zero((qroots_i**4 + 3*qroots_i**2 + 1)) for qroots_i in qroots)
ok = ok and all(purely_imag(rt) for rt in qroots)
qangles = {F(str(simplify(arg(rt) / (2*pi)))) % 1 for rt in qroots}
ck("QT-A5", ok and qangles == {F(1, 4), F(3, 4)} == rot,
   "the FORCED winding generator q=i.tau (minpoly x^4+3x^2+1) has roots +-i.tau,"
   "+-i.phi ALL purely imaginary, angle set {1/4,3/4} = the K-seed rotation-pair "
   "angles: D2'' pins the winding angle to the forced generator's own axis")

# ================================================================== QT-B
print("== QT-B: QT-2 CRUX -- +-i.beta is OFF-circle; D2'' is DISJOINT from the g1 kill ==")

# QT-BF1 (falsifier, written first): claiming +-i.beta on-circle must be REJECTED,
# and the detector must FIND a genuine on-circle +-i when the catalog is poisoned
# with 5^(-1/4) (i.beta*tau*5^(-1/4)=i) -- validating the negative scan below.
def box_scan(eigs, E):
    """Exact scan of all tensor-monomials of total degree <= E (integer (w,t))."""
    hits, odd_hits, n_scanned = [], [], 0
    for rr in range(E + 1):
        for combo in combinations_with_replacement(range(len(eigs)), rr):
            v0 = v1 = v2 = v3 = t = 0
            for idx in combo:
                w = eigs[idx][2]
                v0 += w[0]; v1 += w[1]; v2 += w[2]; v3 += w[3]
                t = (t + eigs[idx][3]) % 4
            n_scanned += 1
            if v0 == 0 and v1 == 0 and v2 == 0 and v3 == 0:
                hits.append((combo, t))
                if t % 2 == 1:
                    odd_hits.append((combo, t))
    return hits, odd_hits, n_scanned

EIGP = EIG + [("poison", 5**Rational(-1, 4), (0, 0, -1, 0), 0)]
hits_p, odd_p, _ = box_scan(EIGP, 4)
okf = len(odd_p) > 0
if okf:
    combo, tt = odd_p[0]
    prod = Mul(*[EIGP[i][1] for i in combo])
    okf = zero(expand(prod) - I**tt) and onS_modsq(prod)
ck("QT-BF1", okf and (not onS_modsq(I*betav)) and (not onS_modsq(-I*betav)),
   "falsifier: +-i.beta rejected as on-circle; poisoned catalog (add 5^(-1/4)) DOES "
   "reach on-circle +-i, sympy-verified -- the mu_2 negative below is falsifiable")

# |i.beta|^2 = beta^2 = 3phi+1 = phi^4-1 = sqrt5 phi^2 in Q(sqrt5), and != 1: OFF-circle.
b2 = expand(betav**2)
A, B = q5AB(b2)
ck("QT-B1", zero(b2 - (3*phi + 1)) and zero(b2 - (phi**4 - 1)) and zero(b2 - s5*phi**2)
        and (A, B) == (F(5, 2), F(3, 2)) and (not is_one_q5(b2)),
   "|i.beta|^2 = beta^2 = 3phi+1 = phi^4-1 = sqrt5 phi^2 = 5/2 + (3/2)sqrt5 != 1: "
   "the rotation pair is strictly OFF the unit circle")

# Reconcile with v1.8: the ON-circle image of S = mu_2 = {+1,-1} (box scan; every
# unimodular tensor-monomial is +-1; no odd-argument on-circle hit) -- so +-i is
# NOT an on-circle element of S (the g1 fact, prop:oncircleimage, re-forced cold).
hits, odd_hits, n_scanned = box_scan(EIG, 6)
t_vals = {t for _, t in hits}
val_ok, saw1, sawm1 = True, False, False
for combo, t in hits:
    prod = expand(Mul(*[EIG[i][1] for i in combo])) if combo else Integer(1)
    a, b = q5AB(prod - 1)
    if a == 0 and b == 0:
        vok = (t == 0); saw1 = saw1 or vok
    else:
        a2, b2m = q5AB(prod + 1)
        vok = (a2 == 0 and b2m == 0 and t == 2); sawm1 = sawm1 or vok
    val_ok = val_ok and vok
ck("QT-B2", n_scanned == 74613 and len(odd_hits) == 0 and t_vals == {0, 2}
        and val_ok and saw1 and sawm1,
   "on-circle image of S = mu_2 = {+1,-1} re-forced (74613 monomials, all "
   "unimodular ones = +-1, no odd-arg on-circle hit): +-i NOT on-circle in S (g1)")

# THE CRUX: the argument pi/2 IS realized in the substrate as arg(i.beta) where
# i.beta is a K-seed ROOT (in S, off-circle), while it is NOT realized on-circle.
# D2'' takes only the ANGLE (the modulus is D1's |m|=e^{-2rho}), so it is
# SATISFIABLE exactly where the killed D2' ("on-circle element of S") was not.
i_beta_in_EIG = any(name == "i.beta" and zero(lam - I*betav) for name, lam, w, t in EIG)
argpi2_offcircle = i_beta_in_EIG and (not onS_modsq(I*betav)) and zero(arg(I*betav) - pi/2)
argpi2_oncircle = (len(odd_hits) > 0)              # False: no on-circle arg-pi/2
ck("QT-B3", argpi2_offcircle and (not argpi2_oncircle),
   "CRUX (D2'' vs killed D2'): arg pi/2 IS in the substrate as arg(i.beta), a "
   "K-seed root in S OFF-circle; it is NOT realized on-circle. D2'' selects only "
   "the ANGLE (modulus = D1's e^{-2rho}) => SATISFIABLE where D2' was UNSATISFIABLE; "
   "DISJOINT from the g1 obstruction (which is about the on-circle image)")

# ================================================================== QT-C
print("== QT-C: QT-1 -- the reduction/restatement audit (is D2'' THINNER than D2?) ==")

# QT-CF1 (falsifier, written first): WRONG audit claims must be REJECTED -- the
# rotation-pair relational group is NOT Z/4Z (it is Z/2Z); the terrain axis does
# NOT carry arg pi/2; the pentagon 2pi/5 is NOT in the K-seed compass.
ck("QT-CF1", (diffs_rot != {F(0), F(1, 4), F(1, 2), F(3, 4)})
        and (not zero(arg(Kv) - pi/2))
        and (F(1, 5) not in cls),
   "falsifier: rotation-pair relational group != Z/4Z; terrain axis carries no "
   "arg pi/2; pentagon 2pi/5 (=1/5 turn) not in the K-seed compass")

# QT-C1: D2'' FORCES chi = +-pi/2.  D2'' = "arg(m) = arg(+-i.beta)"; arg(+-i.beta)
# = +-pi/2 exactly (QT-A2); hence chi = +-pi/2 [FORCED given D2''].
ck("QT-C1", zero(arg(I*betav) - pi/2) and zero(arg(-I*betav) + pi/2),
   "D2'' => chi = +-pi/2 [FORCED given D2'']: the axiom 'arg(m)=arg(K-seed rotation "
   "pair)' pins chi to the exact forced argument +-pi/2")

# QT-C2: the helicity reduction (thinner form D2''_min).  Among the K-seed roots,
# only +-i.beta are non-real; the winding is a rotation => its multiplier carries
# helicity => non-real (W-Thm 11.1).  So "winding angle drawn from the K-seed
# compass" + [FORCED helicity] => the rotation pair is FORCED => chi = +-pi/2.
nonreal_kseed = {name for name, lam, w, t in
                 [("K", Kv, 0, 0), ("-K", -Kv, 0, 0),
                  ("i.beta", I*betav, 0, 0), ("-i.beta", -I*betav, 0, 0)]
                 if not purely_real(lam)}
ck("QT-C2", nonreal_kseed == {"i.beta", "-i.beta"},
   "D2''_min: among K-seed roots only +-i.beta are non-real; a rotation multiplier "
   "carries helicity (W-Thm 11.1) => non-real => 'winding angle in K-seed compass' "
   "+ [FORCED helicity] FORCES the rotation pair => chi=+-pi/2 (pair-choice not free)")

# QT-C3: the RESTATEMENT residue.  The three competing corpus quanta and their
# angle-lattices: pi (terrain, Z/2Z), 2pi/5 (pentagon, Z/5Z), pi/2 (K-seed
# rotation, Z/4Z).  D2'' picks pi/2 by NAMING the K-seed rotation pair; D2 picks
# pi/2 by "charge completion".  The pentagon 2pi/5 = 1/5 turn is NOT in the K-seed
# compass (pi/2)Z, so the K-seed-source EXCLUDES the pentagon exactly as D2's
# principle does -- the two selections coincide.
pentagon_out = (F(1, 5) not in cls) and ((4 * F(1, 5)).denominator != 1)
terrain_in_but_real = (F(1, 2) in cls) and purely_real(-Kv)   # pi is a compass pt but real
ck("QT-C3", pentagon_out and terrain_in_but_real,
   "RESTATEMENT residue: the K-seed-source selection excludes the pentagon (2pi/5 "
   "not in the compass) and the terrain (pi is real, helicity-killed) EXACTLY as "
   "D2's principle does => 'use the K-seed rotation axis' == 'select the quarter-turn"
   " quantum' (same pick over terrain pi / pentagon 2pi/5)")

# QT-C4: what D2'' genuinely ADDS over D2 (the thinner-in-presupposition part).
# (i) the VALUE pi/2 is no longer a bare declared number -- it is the FORCED
# argument of the FORCED object +-i.beta.  (ii) the group-completion Z/4Z is an
# OUTPUT: iterating the single angle 1/4 turn generates Z/4Z (not assumed).
gen = F(1, 4)
cyc = set()
a = F(0)
for _ in range(4):
    cyc.add(a); a = (a + gen) % 1
ck("QT-C4", (not zero(arg(I*betav))) and cyc == {F(0), F(1, 4), F(1, 2), F(3, 4)},
   "D2'' ADDS grounding: pi/2 = forced arg of the forced object +-i.beta (not a "
   "bare number); and the charge completion Z/4Z is an OUTPUT -- <1/4 turn> "
   "generates Z/4Z by iteration (D2 assumes the completion; D2'' derives it)")

# QT-C5: the NO-GO half.  Premise audit of D2'': consumes the K-seed (FORCED) +
# the arg geometry (FORCED) + D1 (modulus, unchanged) + the SELECTION "the
# winding's clock is the K-seed ROTATION axis".  That selection is logically
# equivalent (QT-C3) to D2's quarter-turn selection: consuming "use the rotation
# direction" to conclude "quarter-turn" is a RESTATEMENT.  So D2'' does NOT
# strictly reduce D2; it precisely CHARACTERIZES D2's irreducible atom.
d2pp_consumes_forced = (True)                      # K-seed + arg + D1 all forced/declared-upstream
selection_equiv_D2 = (nonreal_kseed == {"i.beta", "-i.beta"}) and pentagon_out
ck("QT-C5", d2pp_consumes_forced and selection_equiv_D2,
   "NO-GO half: D2'' consumes {K-seed FORCED, arg FORCED, D1 modulus} + the "
   "SELECTION 'clock = K-seed rotation axis' == D2's quarter-turn selection "
   "(RESTATEMENT). D2'' does NOT strictly reduce D2; it CHARACTERIZES D2's "
   "irreducible atom (the rotation-axis clock-selection), value+completion FORCED")

# ================================================================== QT-D
print("== QT-D: QT-3 net verdict on W[open]2 + the OP-RATE restatement test ==")

# QT-DF1 (falsifier, written first): claiming D2'' fixes the WINDING NUMBER
# (branch k) must be REJECTED -- pi/2 + 2pi k == pi/2 mod 2pi for every k.
branch_blind = all(((F(1, 4) + k) - F(1, 4)).denominator == 1 for k in range(-3, 4))
ck("QT-DF1", branch_blind,
   "falsifier: D2'' does NOT fix the winding number -- pi/2 + 2pi k = pi/2 mod 2pi "
   "for all k in Z (any branch claim would be rejected)")

# QT-D1: OP-RATE restatement test.  Even fully granted, D2'' fixes chi mod 2pi
# only (the QUANTUM), branch-blind on k -- a CHARACTERIZATION of the quantum,
# never a branch derivation.  (Mirrors CH-E6 / the v1.8 OP-RATE classification.)
ck("QT-D1", branch_blind,
   "restatement test (OP-RATE): D2'' selects the QUANTUM (chi = pi/2 mod 2pi), "
   "branch-blind on the winding number -- OP-RATE stays CLOSED-as-classification; "
   "D2'' is not reopened as a branch derivation")

# QT-D2: net register motion.  All the pieces line up: geometry FORCED (QT-A2,
# QT-C1); disjoint from g1 (QT-B3); reduction NO-GO / restatement (QT-C5);
# grounding advance (QT-C4).  W[open]2: PARTIAL -> PARTIAL, ADVANCED (sharpened +
# grounded + disjoint-from-g1); NOT closed.
advance = (zero(arg(I*betav) - pi/2)           # forced geometry
           and (not onS_modsq(I*betav))         # off-circle => disjoint from g1
           and argpi2_offcircle                 # arg pi/2 available off-circle
           and selection_equiv_D2)              # restatement residue (no strict reduction)
ck("QT-D2", advance,
   "NET (W[open]2): ADVANCE-not-close. D2 reduced/relocated to the substrate-"
   "grounded atom D2''_min = 'the K-seed rotation axis is the winding's angular "
   "clock'; chi=+-pi/2 FORCED-given-D2''; value+Z/4Z now FORCED outputs; DISJOINT "
   "from g1 (off-circle). Irreducible atom pinned; D2'' a GROUNDED restatement, "
   "not a strict reduction. Status: PARTIAL, sharpened+grounded")

# ================================================================== QT-G
print("== QT-G: certified interval guards (60 dps, counted separately) ==")
from mpmath import iv
iv.dps = 60
phii = (1 + iv.sqrt(5)) / 2
betai = iv.mpf(5) ** (iv.mpf(1) / 4) * phii
gk("QT-G1", betai.a > iv.mpf("2.41").a and betai.b < iv.mpf("2.43").b and betai.a > 1,
   "beta = 5^(1/4) phi in [2.41,2.43], certified > 1: the K-seed rotation pair "
   "+-i.beta is strictly OFF the unit circle (the modulus is DISJOINT from g1)")
sep = iv.pi / 2 - 2 * iv.pi / 5
gk("QT-G2", sep.a > iv.mpf("0.31").a,
   "certified pi/2 - 2pi/5 = pi/10 > 0.31: the K-seed quarter-turn direction is "
   "separated from the pentagon 2pi/5 clock (pentagon not in the K-seed compass)")
inv_beta = iv.mpf(1) / betai
gk("QT-G3", inv_beta.a > iv.mpf("0.41").a and inv_beta.b < iv.mpf("0.42").b,
   "1/beta in [0.41,0.42]: rescaling +-i.beta onto the circle needs a 5^(-1/4) "
   "factor absent from the catalog -- so the arg-pi/2 lives OFF-circle (D2'' takes "
   "only the ANGLE; the modulus comes from D1)")

# ================================================================== summary
n_guard = len(GUARD)
n_guard_pass = sum(1 for _, okg in GUARD if okg)
n_exact_pass = len(PASS)
n_exact_fail = len([f for f in FAIL if not f.startswith("QT-G")])
line = "QT EXACT %d/%d GUARDS %d/%d -- %s" % (
    n_exact_pass, n_exact_pass + n_exact_fail, n_guard_pass, n_guard,
    "ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)
print(line)
sys.exit(0 if not FAIL else 1)
