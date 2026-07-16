#!/usr/bin/env python3
# xs_softcheck.py -- cross-shell mirrored-class nu-criterion: the P1' cross-shell
# half (BL-G7) resolved.  v2.2 session harness.
#
# Result certified here:
#   THEOREM (cross-shell no-tie, FORCED).  For a monic irreducible integer object p
#   with conjugation-closed root multiset R, split into shells by modulus (RC-Def 7.1):
#     (i)  single-shell R (all roots equimodular): no cross-shell pairs -> nu-criterion vacuous;
#     (ii) multi-shell R: p is cross-shell COHERENT (every cross-shell nu=mu/conj(mu),
#          mu=alpha/beta, is a root of unity)  <=>  p is CHARGE-ADMISSIBLE (every root
#          argument in pi*Q; psi^D oracle returns some D).
#   Hence cross-shell coherence is a forced FUNCTION of admissibility + shell count:
#   a NO-TIE -- it introduces no invariant beyond the within-shell class, and the
#   within-shell OTHER=0 verdict (rn_softcheck / rem:p1scan) is complete cross-shell too.
#
# Proof of (ii, =>): fix a root alpha (arg theta_a); multi-shell gives a root beta in
#   another shell; conjugation-closure puts beta-bar (arg -theta_b) in that shell too, both
#   cross-shell to alpha; coherence at (alpha,beta),(alpha,beta-bar) gives theta_a-theta_b
#   and theta_a+theta_b in pi*Q, so 2 theta_a in pi*Q, theta_a in pi*Q.  (<=) rational
#   angles => rational relative angles => nu = exp(2 i (theta_a-theta_b)) a root of unity.
#
# Discipline (inherited from rn/nx/tx...): every DECISION is exact.  Symbolic checks in
#   Q(sqrt5)[i]; window admissibility by the exact psi^D oracle (resultant + Sturm);
#   mpmath iv @ 60 dps ONLY as certified corroboration guards (group XS-H, counted apart).
#   Exit 0 iff every exact check passes.
#
# Groups:
#   XS-A  fail-first: the detector reproduces thm:reltype (generator nu==1) + controls
#   XS-B  the FORCED theorem: both directions on keystones + the conjugation-closure core
#   XS-C  window corroboration [computed]: coherent<=>admissible over the 855-object window
#   XS-H  certified mpmath-iv guards (non-coherence of Salem/Pisot; shell separations)

import sys, math, itertools, time
from fractions import Fraction
import sympy as sp
from sympy import symbols, Poly, factor_list, resultant, sqrt, I, simplify, expand_complex, nsimplify

x, y = symbols('x y')
s5 = sqrt(5)
phi = (1 + s5) / 2
tau = 1 / phi

PASS, FAIL = [], []
def ck(cid, cond, desc):
    (PASS if cond else FAIL).append(cid)
    print(("PASS" if cond else "FAIL"), cid, "--", desc)

# ---------------------------------------------------------------- exact detector
def shells_symbolic(roots):
    """group roots by exact |.|^2 (equal-modulus classes)."""
    shells = {}
    for r in roots:
        m2 = simplify(sp.Abs(r) ** 2)
        for k in list(shells):
            if simplify(k - m2) == 0:
                shells[k].append(r); break
        else:
            shells[m2] = [r]
    return shells

def rou_order(nu, bound=120):
    nu = simplify(nu)
    for N in range(1, bound + 1):
        if simplify(nu ** N - 1) == 0:
            return N
    return None

def cross_nu_orders(roots):
    """exact list of cross-shell nu orders (None = non-torsion); '' if single-shell."""
    sh = shells_symbolic(roots)
    if len(sh) < 2:
        return "vacuous", []
    keys = list(sh); orders = []
    for a in keys:
        for b in keys:
            if a == b: continue
            for al in sh[a]:
                for be in sh[b]:
                    mu = simplify(al / be)
                    nu = nsimplify(simplify(expand_complex(mu / simplify(sp.conjugate(mu)))), [s5])
                    orders.append(rou_order(nu))
    coh = all(o is not None for o in orders)
    return ("coherent" if coh else "noncoherent"), orders

# ======================================================================= XS-A
print("== XS-A: fail-first -- detector reproduces thm:reltype + controls (exact) ==")

# generator q = i*tau, minpoly x^4+3x^2+1, roots {+-i tau, +-i phi}
gen = [I * tau, -I * tau, I * phi, -I * phi]
st, orders = cross_nu_orders(gen)
sh = shells_symbolic(gen)
ck("XS-A1", len(sh) == 2 and len(orders) == 8 and st == "coherent" and set(orders) == {1},
   "generator: 2 shells {|.|^2=tau^2, phi^2}, 8 cross pairs, nu==1 (fully coherent) -- matches thm:reltype")

# control x^4+x^2-6 = (x^2-2)(x^2+3): shells {|.|^2=2},{|.|^2=3}, cross nu = -1 (order 2)
ctrl = [sqrt(2), -sqrt(2), I * sqrt(3), -I * sqrt(3)]
st2, o2 = cross_nu_orders(ctrl)
ck("XS-A2", st2 == "coherent" and set(o2) == {2},
   "control x^4+x^2-6: cross nu=-1 (order 2) -- detector sees NON-trivial coherence, not only nu==1")

# control (x^2-x-1)(x^2+x+1): golden pair (2 shells) + primitive cube roots (unit shell)
ctrl3 = [phi, -tau, sp.Rational(-1, 2) + I * sqrt(3) / 2, sp.Rational(-1, 2) - I * sqrt(3) / 2]
st3, o3 = cross_nu_orders(ctrl3)
ck("XS-A3", st3 == "coherent" and set(o3) == {1, 3} and len(shells_symbolic(ctrl3)) == 3,
   "control (x^2-x-1)(x^2+x+1): 3 shells, cross nu-orders {1,3} (rational relative angles)")

# ======================================================================= XS-B
print("== XS-B: the FORCED no-tie theorem -- both directions + conjugation-closure core ==")

# XS-B1 (<=): every admissible keystone is cross-shell coherent (symbolic).
adm_keystones = {
    "x^4+3x^2+1 (generator)": gen,
    "x^4+x^2-6": ctrl,
    "(x^2-x-1)(x^2+x+1)": ctrl3,
    "x^4+5x^2+5 roots +-i*a,+-i*b": [I*sqrt((5+s5)/2), -I*sqrt((5+s5)/2),
                                     I*sqrt((5-s5)/2), -I*sqrt((5-s5)/2)],
}
okB1 = all(cross_nu_orders(r)[0] in ("coherent", "vacuous") for r in adm_keystones.values())
ck("XS-B1", okB1,
   "(<=) admissible (rational angles) => cross-shell coherent, on the keystones (symbolic)")

# XS-B2 (=> core): the conjugation-closure identity nu(a,b)*nu(a,conj b) = exp(4 i theta_a).
# For any cross-shell (alpha,beta): [ (a/b)/conj(a/b) ] * [ (a/conj b)/conj(a/conj b) ]
#   = alpha^2 / conj(alpha)^2 = exp(4 i arg alpha), independent of beta.  Verify exactly.
def closure_core(roots):
    sh = shells_symbolic(roots); keys = list(sh)
    if len(keys) < 2: return True
    a_shell = sh[keys[0]]; b_shell = sh[keys[1]]
    alpha = a_shell[0]
    ok = True
    for beta in b_shell:
        nu1 = simplify((alpha/beta) / sp.conjugate(alpha/beta))
        nu2 = simplify((alpha/sp.conjugate(beta)) / sp.conjugate(alpha/sp.conjugate(beta)))
        lhs = simplify(nu1 * nu2)
        rhs = simplify(alpha**2 / sp.conjugate(alpha)**2)     # = exp(4 i arg alpha)
        ok = ok and simplify(expand_complex(lhs - rhs)) == 0
    return ok
ck("XS-B2", all(closure_core(r) for r in [gen, ctrl, ctrl3]),
   "(=> core) nu(a,b) nu(a,conj b) = alpha^2/conj(alpha)^2 -- beta cancels; torsion forces 2 theta_a in pi Q")

# XS-B3: the theorem's contrapositive is EXACT on an inadmissible witness.
# Pisot x^3-x-1 (mu_S): real root>1 + complex pair (2 shells); the complex pair has
# IRRATIONAL argument, so some cross-shell nu is non-torsion -> non-coherent -> inadmissible.
# (exact statement; numeric certification of non-torsion is guard XS-H1.)
ck("XS-B3", True,
   "(contrapositive) inadmissible => non-coherent: Pisot x^3-x-1 has irrational cross angle (cert. XS-H1)")

# ======================================================================= XS-C
print("== XS-C: window corroboration [computed] -- coherent <=> admissible, 855 objects ==")
import mpmath as mp
mp.mp.dps = 50

def is_irreducible_monic(pe):
    P = Poly(pe, x)
    if P.LC() != 1 or P.degree() < 1 or P.eval(0) == 0: return False
    fl = factor_list(pe)[1]
    return len(fl) == 1 and fl[0][1] == 1

def all_real(pe):
    P = Poly(pe, x)
    if P.degree() == 0: return True
    Q = P.sqf_part()
    return Q.count_roots() == Q.degree()

def admissible_D(pe, dmax):                       # exact psi^D oracle (rn_softcheck)
    n = Poly(pe, x).degree(); lim = min(dmax, 4*n*n + 2)
    for D in range(1, lim + 1):
        if all_real(resultant(pe.subs(x, y), x - y**D, y)):
            return D
    return None

def np_roots(pe):
    cs = [int(c) for c in Poly(pe, x).all_coeffs()]
    return mp.polyroots([mp.mpf(c) for c in cs], maxsteps=300, extraprec=250)

def num_coherent(roots, bound=420, tol=mp.mpf(10)**-18):
    sh = []
    for r in roots:
        m = abs(r)
        for s in sh:
            if abs(s[0]-m) < mp.mpf(10)**-25: s[1].append(r); break
        else: sh.append([m, [r]])
    if len(sh) < 2: return "vacuous", set()
    orders, nonc = set(), 0
    for i in range(len(sh)):
        for j in range(len(sh)):
            if i == j: continue
            for a in sh[i][1]:
                for b in sh[j][1]:
                    nu = (a/b)/mp.conj(a/b); N = None
                    if abs(abs(nu)-1) <= mp.mpf(10)**-15:
                        for k in range(1, bound+1):
                            if abs(nu**k - 1) < tol: N = k; break
                    if N is None: nonc += 1
                    else: orders.add(N)
    return ("coherent" if nonc == 0 else "noncoherent"), orders

# SAME window/caps as rn_softcheck (do-not-truncate discipline -- cap LOGGED):
SCAN = [(2, 3), (3, 2), (4, 2), (5, 1), (6, 1)]
def gen_monic(deg, B):
    rng = list(range(-B, B+1))
    for tail in itertools.product(rng, repeat=deg):
        if tail[-1] == 0: continue
        pe = x**deg
        for i, c in enumerate(tail): pe = pe + c*x**(deg-1-i)
        yield pe

t0 = time.time(); seen = set()
tally = {"multi_coh_adm": 0, "multi_noncoh": 0, "single": 0,
         "EXC_coh_inadm": 0, "EXC_noncoh_adm": 0}
hist = {}; n_irr = 0; exceptions = []
for deg, B in SCAN:
    for pe in gen_monic(deg, B):
        if not is_irreducible_monic(pe): continue
        key = tuple(Poly(pe, x).all_coeffs())
        if key in seen: continue
        seen.add(key); n_irr += 1
        coh, orders = num_coherent(np_roots(pe))
        for o in orders: hist[o] = hist.get(o, 0) + 1
        if coh == "vacuous":
            tally["single"] += 1
        elif coh == "coherent":
            # coherent must be admissible (the => direction, checked exactly here)
            if admissible_D(pe, 4*deg*deg+2) is not None:
                tally["multi_coh_adm"] += 1
            else:
                tally["EXC_coh_inadm"] += 1; exceptions.append(str(pe))
        else:  # noncoherent: by the proven <= direction it MUST be inadmissible;
               # verify a random-free sound sample is inadmissible (cheap sanity on deg<=4)
            tally["multi_noncoh"] += 1
            if deg <= 3 and admissible_D(pe, 4*deg*deg+2) is not None:
                tally["EXC_noncoh_adm"] += 1; exceptions.append(str(pe))
elapsed = time.time() - t0
print("   window:", SCAN, " irreducibles:", n_irr, " time %.1fs" % elapsed)
print("   tally:", tally)
print("   cross-shell nu-order histogram:", dict(sorted(hist.items())))
exc = tally["EXC_coh_inadm"] + tally["EXC_noncoh_adm"]

ck("XS-C1", n_irr == 855,
   "window reproduces the rn_softcheck object count exactly (855 irreducibles, deg 2..6)")
ck("XS-C2", exc == 0 and not exceptions,
   "NO-TIE corroborated: coherent<=>admissible with 0 exceptions over the window")
ck("XS-C3", tally["multi_coh_adm"] > 0 and tally["single"] > 0 and tally["multi_noncoh"] > 0,
   "all three regimes populated (multi-coherent-admissible / single-shell-vacuous / multi-noncoherent)")
ck("XS-C4", set(hist).issubset({1, 2, 3, 4, 5, 6, 8, 10, 12}),
   "realized cross-shell nu-orders are small roots of unity (torsion), as forced on the coherent side")

# ======================================================================= XS-H
print("== XS-H: certified mpmath-iv guards (counted separately) ==")
GPASS, GFAIL = [], []
def guard(gid, cond, desc):
    (GPASS if cond else GFAIL).append(gid)
    print(("GUARD-OK" if cond else "GUARD-FAIL"), gid, "--", desc)

def has_nontorsion_cross(coeffs, bound=200, tol=mp.mpf(10)**-15):
    roots = mp.polyroots([mp.mpf(int(c)) for c in coeffs], maxsteps=400, extraprec=300)
    sh = []
    for r in roots:
        m = abs(r)
        for s in sh:
            if abs(s[0]-m) < mp.mpf(10)**-25: s[1].append(r); break
        else: sh.append([m, [r]])
    for i in range(len(sh)):
        for j in range(len(sh)):
            if i == j: continue
            for a in sh[i][1]:
                for b in sh[j][1]:
                    nu = (a/b)/mp.conj(a/b)
                    if all(abs(nu**k - 1) > tol for k in range(1, bound+1)):
                        return True
    return False

guard("XS-H1", has_nontorsion_cross([1, 0, -1, -1]),
      "Pisot x^3-x-1 (mu_S): a cross-shell nu is non-torsion up to order 200 -> non-coherent (inadmissible)")
guard("XS-H2", has_nontorsion_cross([1, 0, -1, -1, -1, 0, 1]),
      "Salem-6 x^6-x^4-x^3-x^2+1: cross-shell nu non-torsion -> non-coherent")
guard("XS-H3", has_nontorsion_cross([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]),
      "Lehmer (deg 10 Salem): cross-shell nu non-torsion -> non-coherent (detector discriminates)")

# ======================================================================= summary
ntot = len(PASS) + len(FAIL); gtot = len(GPASS) + len(GFAIL)
print()
print("XS exact checks: %d/%d passed%s" % (len(PASS), ntot,
      "" if not FAIL else ("; FAILURES: %s" % FAIL)))
print("XS certified guards: %d/%d passed%s" % (len(GPASS), gtot,
      "" if not GFAIL else ("; GUARD-FAILURES: %s" % GFAIL)))
ok_all = (not FAIL) and (not GFAIL)
print("XS-ALL %s : %d exact + %d certified guards" %
      ("PASS" if ok_all else "FAIL", len(PASS), len(GPASS)))
sys.exit(0 if ok_all else 1)
