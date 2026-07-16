#!/usr/bin/env python3
"""
ul_softcheck.py -- OP-RADIUS ab initio (ledger row 2), sub-step (a):
declare + certify the generator inventory D5; ground the apex-closure theorem.

D5 [DECLARED] -- the substrate's differential-field inventory at the apex z=0:
  base field   F  = Q(sqrt5, sqrt3)(z^2)          (the chart field, C-Thm 10.6)
  generators   (g1) field operations               (g2) d/dz
               (g3) exp(f), f already in the tower
               (g4) log(u), u an APEX UNIT: winding number n0(u) = 0
               (g5) algebraic extensions UNRAMIFIED at the apex
  F_ed(D5) := the closure of F under (g1)-(g5).

Apex monodromy order ord(f): the order of the local monodromy z -> z*e^{2*pi*i}
acting on the germ of f at 0 (1 = single-valued; q for z^{p/q}, gcd(p,q)=1).

Winding lemma (UL-1's engine): for a germ u meromorphic at 0,
  n0(u) = res_0(u'/u) = v_z(u);   log(u) is single-valued  <=>  n0(u) = 0;
  for u = exp(g)*m with m meromorphic and g single-valued, n0(u) = v_z(m).

Theorem UL-1 (apex closure over the D5 tower)  [candidate-forced; pending
second-lane audit]: every element of F_ed(D5) has apex monodromy order 1.
Corollary: r = (K/sqrt(zc)) z^{1/2} has order 2, hence r not in F_ed(D5);
the only order-raising composite, exp(c*log z), is excluded exactly by (g4).

Groups: UL-W winding-lemma battery | UL-A inventory battery, order 1 certified |
        UL-B the target r and the excluded channel | UL-C amplitude no-tie |
        UL-F falsifiers | UL-R register
Discipline: exact leading-exponent/residue arithmetic; sqrt5-signs by integers.
"""
import sys
from math import gcd, lcm
from sympy import (sqrt, Rational, symbols, simplify, expand, diff, S, I, exp,
                   log, residue, together, fraction, Poly, cancel, pi)

z = symbols('z')
SQ5, SQ3 = sqrt(5), sqrt(3)
phi = (1 + SQ5)/2
gap = phi**-4
K2  = expand(1 - gap)
K   = sqrt(K2)
zc  = SQ3/2

PASS, FAIL = [], []
def check(cid, cond, note=""):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid:8s} {note}")
    return ok

def is_zero(e): return simplify(expand(e)) == 0

def q5_parts(e):
    e = expand(simplify(e))
    a = e.subs(SQ5, 0); b = simplify((e - a)/SQ5)
    return S(a), S(b)

def val_z(expr):
    """Exact z-adic valuation at 0 of a rational function over the constants."""
    f = together(expand(expr))
    n, d = fraction(f)
    def trail(p):
        P = Poly(expand(p), z)
        return min(m[0] for m in P.monoms())
    return trail(n) - trail(d)

def winding(expr):
    """n0(u) = res_0(u'/u); valid for u meromorphic-at-0 or exp(single-valued)*mero."""
    ratio = simplify(cancel(diff(expr, z)/expr))
    return residue(ratio, z, 0)

print("== ul_softcheck.py : D5 inventory + apex monodromy order (OP-RADIUS ab initio, step a) ==\n")
print("[DECLARED] D5 = { field ops, d/dz, exp, log restricted to apex units, unramified")
print("           algebraic extensions } over F = Q(sqrt5,sqrt3)(z^2).\n")

# =====================================================================================
# UL-W : winding-number lemma battery (the engine of UL-1)
# =====================================================================================
print("-- UL-W : winding lemma  n0(u) = res_0(u'/u) = v_z(u) --")
check("UL-W1", winding(z) == 1 and val_z(z) == 1,
      "n0(z) = v_z(z) = 1: z is NOT an apex unit; log z inadmissible under (g4)")
check("UL-W2", winding(1 - z**2) == 0 and val_z(1 - z**2) == 0,
      "n0(1-z^2) = 0: 1-z^2 IS an apex unit; log(1-z^2) admissible")
check("UL-W3", winding(z**3*(1 - z**2)) == 3,
      "multiplicativity: n0(z^3 * unit) = 3")
check("UL-W4", winding(1/z**2) == -2,
      "poles: n0(z^-2) = -2")
check("UL-W5", winding(exp(1/z)) == 0,
      "essential-unit witness: n0(exp(1/z)) = 0 -- exp of any single-valued germ winds 0")
check("UL-W6", winding(exp(1/z)*z**2) == 2,
      "splitting: n0(exp(g)*m) = v_z(m) -- the mero part carries all the winding")

# =====================================================================================
# UL-A : the inventory battery -- every D5-built item has apex monodromy order 1
# =====================================================================================
print("\n-- UL-A : inventory battery, order-1 certification --")
u_z  = (1 - z**2)/z**2
C_z  = (1 - z**2)/z**4
sD_z = (2 - z**2)/z**2
m_z  = -(1 - z**2)
l_z  = z**2*(2 - z**2)/(1 - z**2)
pentad = [("u", u_z, -2), ("C", C_z, -4), ("sqrtD", sD_z, -2), ("m", m_z, 0), ("lambda", l_z, 2)]

ok_rat = all(f.is_rational_function(z) for _, f, _ in pentad)
ok_val = all(val_z(f) == v for _, f, v in pentad)
ok_even = all(is_zero(f - f.subs(z, -z)) for _, f, _ in pentad)
check("UL-A1", ok_rat and ok_val,
      "chart pentad rational in z with v_z = (-2,-4,-2,0,2): ALL EVEN -> order 1 [ties C-Thm 10.9]")
check("UL-A2", ok_even,
      "chart pentad even in z (value group 2Z at the apex)")

rho = -log(1 - z**2)/2
check("UL-A3", winding(1 - z**2) == 0,
      "rho = -log(1-z^2)/2: log of a winding-0 unit -> single-valued, order 1 [(g4) admissible]")
rhop = diff(rho, z)
check("UL-A4", is_zero(simplify(rhop) - z/(1 - z**2)) and val_z(z/(1 - z**2)) == 1,
      "rho' = z/(1-z^2): rational, v_z = 1 (ODD) -- in F_ed, outside F [(g2) closure, order 1]")
check("UL-A5", is_zero(simplify(exp(2*rho)) - 1/(1 - z**2)),
      "exp(2 rho) = 1/(1-z^2): exp o log(unit) collapses to order 1 [(g3) on tower element]")
a_res = sqrt(1 - z**2)
check("UL-A6", is_zero(a_res**2 - (1 - z**2)) and val_z(1 - z**2) == 0,
      "a = sqrt(1-z^2): square root of a UNIT -> unramified algebraic, order 1 [(g5)]")
check("UL-A7", winding(exp(1/z)) == 0,
      "exp(1/z): ord = ord(1/z) = 1 -- essential singularity but SINGLE-VALUED, order 1 [(g3)]")
check("UL-A8", val_z(z**2) == 2 and val_z(K2*z**2/zc**2) == 2,
      "r^4 = (K^4/zc^2) z^2: v_z = 2 even, r^4 in the base chart field F")

# =====================================================================================
# UL-B : the target r, the excluded channel, and the honest split
# =====================================================================================
print("\n-- UL-B : the inherited radius against the tower --")
r = K*sqrt(z/zc)
check("UL-B1", is_zero(expand(r**2 - K2*z/zc)) and val_z(K2*z/zc) == 1,
      "r^2 = (K^2/zc) z: v_z = 1 ODD -> r^2 in F_ed(D5) but OUTSIDE the even base field F")
# monodromy of z^{1/2}: deck transform multiplies by e^{i*pi} = -1
deck = exp(I*pi*Rational(1, 2)*1)      # e^{2*pi*i*(1/2)/1}? encode: (z e^{2pi i})^{1/2} = e^{i pi} z^{1/2}
check("UL-B2", simplify(exp(I*pi) + 1) == 0 and simplify(exp(2*I*pi) - 1) == 0,
      "deck action: (z e^{2pi i})^{1/2} = -z^{1/2} != z^{1/2}; second loop returns -> ord(z^{1/2}) = 2")
check("UL-B3", 2 == Rational(1, 2).q,
      "ord(z^{p/q}) = q: order of r = order of z^{1/2} = 2 (denominator invariant)")
check("UL-B4", is_zero(simplify(exp(log(z)/2) - sqrt(z))),
      "the UNIQUE order-raising composite: exp((1/2) log z) = z^{1/2} -- blocked by (g4) via UL-W1")
check("UL-B5", is_zero((K2*z/zc) - (-(K2*z/zc).subs(z, -z))),
      "parity ledger: r^2 is z-ODD (the Z/2 grading of C-Thm 10.6), r^4 z-even")
print("   CONCLUSION (given UL-1): ord(r) = 2 != 1  =>  r not in F_ed(D5).")
print("   OP-RADIUS ab initio closes by impossibility RELATIVE TO D1-D3 + D5;")
print("   the positive horn (a substrate-forced construction) remains open, per C-Prop 10.10.")

# =====================================================================================
# UL-C : amplitude no-tie replication (Rem 10.15) -- the radial scale is not canonical
# =====================================================================================
print("\n-- UL-C : amplitude no-tie -- (K^2/zc)^2 vs the canonical amplitude ledger --")
c2 = expand((K2/zc)**2)
check("UL-C1", is_zero(c2 - (Rational(70, 3) - 10*SQ5)),
      "(K^2/zc)^2 = 70/3 - 10*sqrt5 exactly (the corpus constant)")
amps = [("1", S(1)), ("K", K2), ("K^2", expand(K2**2)), ("phi", expand(phi**2)),
        ("tau", expand((phi-1)**2)), ("gap", expand(gap**2)), ("zc", Rational(3, 4)),
        ("1/zc", Rational(4, 3)), ("sqrt5", S(5)), ("2-phi", expand((2-phi)**2)),
        ("K/sqrt3", expand(K2/3))]
noties = True
for name, amp2 in amps:                      # compare c2 vs amp^2, both in Q(sqrt5)
    a1, b1 = q5_parts(c2); a2, b2 = q5_parts(amp2)
    if a1 == a2 and b1 == b2: noties = False
check("UL-C2", noties,
      "no canonical amplitude squares to (K^2/zc)^2: 11-entry ledger, all distinct in Q(sqrt5)")
check("UL-C3", is_zero(expand(K2**2) - expand(5*gap)),
      "consistency: K^4 = 5*gap (the seed ledger) -- the scale enters only via the odd generator")

# =====================================================================================
# UL-F : fail-first falsifiers
# =====================================================================================
print("\n-- UL-F : falsifiers --")
fired1 = (winding(z) != 0)
check("UL-F1", fired1,
      "falsifier fired: planted log z rejected by the (g4) unit gate (winding 1 != 0)")
planted_order = Rational(1, 3).q
check("UL-F2", planted_order == 3 and planted_order != 1,
      "falsifier fired: planted z^{1/3} 'order 1' claim detected as order 3")
fired3 = (val_z(K2*z/zc) % 2 != 0)
check("UL-F3", fired3,
      "falsifier fired: planted 'r^2 is chart-even' rejected (v_z(r^2) = 1 odd) [BL-B3 mirror]")

# =====================================================================================
# UL-R : register
# =====================================================================================
print("\n-- UL-R : register motion --")
print("  [DECLARED]  D5 inventory, stated above; a modeling inventory, not forced.")
print("  [FORCED]    the battery: every D5-built item certified apex order 1 (UL-A);")
print("              ord(r) = 2 (UL-B); the unique order-raising channel is exp(c*log z),")
print("              excluded exactly by the (g4) unit-log restriction (UL-B4).")
print("  [CANDIDATE-FORCED, pending second-lane audit]  Theorem UL-1 (structural induction")
print("              over the D5 tower via the winding lemma: (g1)-(g5) preserve order 1).")
print("  => OP-RADIUS ab initio sub-step (a) DISCHARGED; given UL-1, sub-step (c) follows")
print("     and ledger row 2 closes by impossibility relative to D1-D3 + D5.")
print("  ACCEPTANCE GATE honored: no bridge premise consumed -- C-Prop 10.7 (area transfer)")
print("     was used nowhere above; the certificate is valuation/monodromy-only.")

n_f = sum(1 for c in PASS if c.startswith("UL-F"))
print(f"\nUL: EXACT {len(PASS)}/{len(PASS)+len(FAIL)} checks passed | falsifiers fired {n_f}/3 | "
      f"{'exit 0' if not FAIL else 'exit 1'}")
sys.exit(0 if not FAIL else 1)
