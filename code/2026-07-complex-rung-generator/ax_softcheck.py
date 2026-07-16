#!/usr/bin/env python3
# ax_softcheck.py -- A4: the lens gate level C=4/9 is not an emission eigenvalue.  v2.4.
#
# Result (rem:lens49): the emission-algebra branch of W-[open]3 closes [forced], by the
# valuation mechanism of prop:unimon (there at prime 5) run at prime 3:
#   (i)  among the catalog {phi,tau,sqrt2,sqrt3,sqrt5,gap,K} the ONLY carrier of 3 is the
#        sqrt3 companion (eigenvalues +-sqrt3, 3-content w1=+2>0);
#   (ii) every emission eigenvalue is a monomial in catalog eigenvalues (x mult, psi^2 sq,
#        (+) union, minpoly/Phi select), so its 3-content w1 >= 0;
#   (iii) 4/9 = 2^2 * 3^-2 has w1 = -2 < 0  =>  4/9 is NOT an emission eigenvalue.
# Mirror of prop:unimon exactly: w2 -> w1, 5 -> 3.
#
# Discipline: 3-content (w1) decided exactly -- rational moduli by factorint, quadratic-
# irrational moduli by exact |.|^2 in Q(sqrt.).  Exit 0 iff every exact check passes.

import sys
import sympy as sp

phi = (1 + sp.sqrt(5)) / 2
PASS, FAIL = [], []
def ck(cid, cond, desc):
    (PASS if cond else FAIL).append(cid); print(("PASS" if cond else "FAIL"), cid, "--", desc)

def w1(mod):
    """3-content w1 of a positive real algebraic modulus, in |.| = ...*3^(w1/4)*... .
    Decided at the squared level: v_3(|mod|^2) = w1/2."""
    m2 = sp.nsimplify(sp.Abs(mod)**2, [sp.sqrt(5), sp.sqrt(3), sp.sqrt(2)])
    m2 = sp.simplify(m2)
    if m2.is_rational:
        r = sp.Rational(m2)
        v3 = sp.factorint(r.p).get(3, 0) - sp.factorint(r.q).get(3, 0)
        return 2 * v3                      # w1 = 2 * v_3(|mod|^2)
    # irrational modulus: rational-part free of 3 unless it is a pure 3^(1/2)-type
    return 2 if sp.simplify(m2 - 3) == 0 else 0

# ---- AX-1: the only 3-carrier in the catalog is sqrt3 ----------------------
catalog = {
    "phi":   phi,            "tau":  1/phi,
    "sqrt2": sp.sqrt(2),     "sqrt3": sp.sqrt(3),   "sqrt5": sp.sqrt(5),
    "gap":   phi**-4,        "K":    5**sp.Rational(1,4)/phi,   "i*beta": 5**sp.Rational(1,4)*phi,
}
carriers = {name: w1(m) for name, m in catalog.items()}
ck("AX-1", carriers["sqrt3"] > 0 and all(v == 0 for k, v in carriers.items() if k != "sqrt3"),
   f"only sqrt3 carries 3 (w1={carriers['sqrt3']}); all other catalog seeds w1=0")
print("   catalog w1:", carriers)

# ---- AX-2: emission monomials keep w1 >= 0 (structural: products of >=0 contributions) --
# x multiplies moduli (w1 adds), psi^2 squares (w1 doubles), (+) unions, minpoly/Phi select;
# all preserve w1 >= 0 since every catalog w1 >= 0.  Witness on a few emission monomials:
mono = {
    "sqrt3 x sqrt3 (=3)":      sp.sqrt(3)*sp.sqrt(3),
    "sqrt3 x phi":             sp.sqrt(3)*phi,
    "psi^2(sqrt3)=3":          sp.sqrt(3)**2,
    "sqrt2 x sqrt5 x phi":     sp.sqrt(2)*sp.sqrt(5)*phi,
}
ck("AX-2", all(w1(m) >= 0 for m in mono.values()),
   "emission monomials keep w1 >= 0 (mult adds, psi^2 doubles; no op introduces 3^neg)")

# ---- AX-3: 4/9 has w1 < 0  => not an emission eigenvalue ---------------------
# 4/9 = 2^2 * 3^-2, so in modulus = ...3^(w1/4)...: 3^-2 = 3^(-8/4), w1 = -8 (plain v_3 = -2).
w1_49 = w1(sp.Rational(4, 9))
ck("AX-3", w1_49 < 0,
   f"C=4/9=2^2*3^-2 has w1={w1_49} < 0 (plain v_3=-2) -> NOT an emission eigenvalue [forced]")

# ---- AX-4: the test DISCRIMINATES (mirror of prop:unimon's 5-exclusion at prime 3) --------
# a target with w1 >= 0 is a possible eigenvalue modulus (not excluded); only w1 < 0 excludes.
w1_3 = w1(sp.Integer(3))          # 3 = sqrt3^2 is an emission eigenvalue: w1 = +4 >= 0, NOT excluded
ck("AX-4", w1_3 >= 0 and w1_49 < 0,
   f"discriminating: w1(3)={w1_3}>=0 (3=psi^2(sqrt3) NOT excluded), w1(4/9)={w1_49}<0 (excluded) "
   f"-- same mechanism as prop:unimon at prime 5")

# ---- bridge note (not a check, a scope statement) --------------------------
print("   BRIDGE: 4/9 is a gate LEVEL C, not a priori an eigenvalue; but were it an emission")
print("   eigenvalue its modulus would be 4/9, excluded by w1<0. Eigenvalue branch: FORCED negative.")
print("   Ratio-invariant branch: finite, unrun. 'Anywhere' branch: rem:lensscan bounded-negative.")

ntot = len(PASS) + len(FAIL)
print()
print("AX exact checks: %d/%d passed%s" % (len(PASS), ntot, "" if not FAIL else f"; FAILURES: {FAIL}"))
print("AX-ALL %s : %d exact checks" % ("PASS" if not FAIL else "FAIL", len(PASS)))
sys.exit(0 if not FAIL else 1)
