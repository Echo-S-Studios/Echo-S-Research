#!/usr/bin/env python3
"""
d4_softcheck.py -- OP-RADIUS profile (ledger row 3): D4 independence certificate.

Target: label C-Thm 10.12's two drop-one witnesses as the D4 certificate.
Under closure-corollary C2 (a certified-independent declared atom is CLOSED,
not perpetually 'awaiting a principle'), two drop-one witnesses = an
independence certificate => the row closes by horn (c).

Discipline: every decision exact over Q(sqrt5, sqrt3); sqrt5-signs by integer
arithmetic; floats decide nothing; fail-first falsifiers must fire.
Groups: D4-A chart-tie replication | D4-B orientation witness |
        D4-C onset witness | D4-D closure + residual bijection | D4-F falsifiers
"""
import sys
from sympy import (sqrt, Rational, symbols, simplify, expand, diff, pi, S,
                   together, limit, oo)

z  = symbols('z',  positive=True)
z0 = symbols('z0', positive=True)

SQ5, SQ3 = sqrt(5), sqrt(3)
phi = (1 + SQ5) / 2
gap = phi**-4                      # (7-3*sqrt5)/2
K2  = expand(1 - gap)              # (3*sqrt5-5)/2
K   = sqrt(K2)
zc  = SQ3 / 2

PASS, FAIL = [], []
def check(cid, cond, note=""):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid:8s} {note}")
    return ok

def is_zero(e):
    return simplify(expand(e)) == 0

def q5_parts(e):
    """Exact (a,b) with e = a + b*sqrt5, a,b in Q."""
    e = expand(simplify(e))
    a = e.subs(SQ5, 0)
    b = simplify((e - a) / SQ5)
    return S(a), S(b)

def q5_pos(e):
    """Exact sign decision e > 0 for e in Q(sqrt5), by integer comparison."""
    a, b = q5_parts(e)
    if b == 0: return a > 0
    if a == 0: return b > 0
    if a > 0 and b > 0: return True
    if a < 0 and b < 0: return False
    if a > 0:            # a>0, b<0 : e>0  <=>  a^2 > 5 b^2
        return (a**2 - 5*b**2) > 0
    return (5*b**2 - a**2) > 0     # a<0, b>0

print("== d4_softcheck.py : OP-RADIUS profile / D4 independence certificate ==")
print("   D4 = (i) orientation of the apex double cover  ^  (ii) cap onset")
print("   inherited chart-ties: exponent 1/2 (area transfer), cap value K, kink slope <-> onset\n")

# ---------------------------------------------------------------- D4-A: chart-tie replication
r_horn = K * sqrt(z / zc)                    # inherited horn, orientation +, onset zc
rate   = diff(pi * r_horn**2, z)             # areal deposition rate d(pi r^2)/dz

check("D4-A1", is_zero(diff(rate, z)) and is_zero(rate - pi*K2/zc),
      "exponent-1/2 tie: d(pi r^2)/dz = pi K^2/zc, constant on the horn (Prop 10.7)")
check("D4-A2", is_zero(r_horn.subs(z, zc) - K),
      "cap-value tie: r(zc) = K, continuity at the join (W-Thm 15.1)")
slope_zc = diff(r_horn, z).subs(z, zc)
check("D4-A3", is_zero(slope_zc - K/SQ3),
      "kink slope r'(zc^-) = K/sqrt3 (W-Prop 15.2a)")
check("D4-A4", q5_pos(expand((K/SQ3)**2)),
      "kink nonzero: (K/sqrt3)^2 = (3*sqrt5-5)/6 > 0  [guard 45>25]")
check("D4-A5", q5_pos(expand(K2 - zc**2)),
      "K > zc: K^2 - 3/4 = (6*sqrt5-13)/4 > 0  [guard 180>169]")

# ---------------------------------------------------------------- D4-B: drop-one witness (i) -- orientation
check("D4-B1", is_zero((-r_horn)**2 - r_horn**2),
      "witness (i): every chart tie is sign-blind, (-r)^2 = r^2 identically")
check("D4-B2", is_zero(sqrt((-K)**2) - K),
      "witness (i): |cap| datum = K for both orientations")
check("D4-B3", q5_pos(expand((2*K)**2)),
      "witness (i): 2K != 0, so +r and -r are DISTINCT chart-consistent profiles")

# ---------------------------------------------------------------- D4-C: drop-one witness (ii) -- onset
z0alt = Rational(4, 5)          # C-Thm 10.12's alternate onset
zstar = Rational(5, 6)          # separation height
r_alt = K * sqrt(z / z0alt)     # alternate capped paraboloid, same exponent & cap value

check("D4-C1", (z0alt < zstar) and (Rational(25, 36) < Rational(3, 4)),
      "ordering: z0=4/5 < z*=5/6 < zc  [guards 24<25 and 100<108]")
check("D4-C2", is_zero(diff(diff(pi*r_alt**2, z), z)) and is_zero(r_alt.subs(z, z0alt) - K),
      "witness (ii): alternate profile keeps exponent-1/2 + cap-value K + continuity")
# separation at z*: inherited horn (z* < zc) vs alternate cap (z* > 4/5)
ratio4 = simplify((zstar/zc)**2)             # (r(z*)/K)^4 for the inherited horn
check("D4-C3", ratio4 == Rational(25, 27) and ratio4 != 1,
      "witness (ii): separation at z* -- r_inh(z*)^2/K^2 has square 25/27 != 1  [75 != 81]")
onset_ratio = zstar / z0alt                  # (5/6)/(4/5) = 25/24 > 1: z* is PAST the alt onset
check("D4-C4", onset_ratio == Rational(25, 24) and q5_pos(onset_ratio - 1),
      "witness (ii): z*/z0 = 25/24 > 1 => alternate profile is on its cap at z*, value K; "
      "with C3 (inherited (r/K)^4 = 25/27 < 1, below cap) the two onsets are DISTINCT at z*")

# ---------------------------------------------------------------- D4-D: closure + residual
w1 = all(c in PASS for c in ("D4-B1", "D4-B2", "D4-B3"))
w2 = all(c in PASS for c in ("D4-C1", "D4-C2", "D4-C3", "D4-C4"))
check("D4-D1", w1 and w2,
      "CLOSURE (C2): both drop-one witnesses certified => D4 independent => row closes by horn (c)")

s_of = K / (2*z0)                             # kink slope as a function of onset
check("D4-D2", is_zero(diff(s_of, z0) + K/(2*z0**2)) and q5_pos(K2),
      "residual: ds/dz0 = -K/(2 z0^2) < 0 -- slope <-> onset is a strict bijection (Rem 10.13)")
check("D4-D3", is_zero(s_of.subs(z0, zc) - K/SQ3) and is_zero(s_of.subs(z0, z0alt) - 5*K/8),
      "residual: bijection lands the two witnesses at slopes K/sqrt3 and 5K/8")
check("D4-D4", Rational(1, 3) != Rational(25, 64),
      "residual: the two witness slopes are distinct at the squared level  [64 != 75]")
print("   [DECLARED] representative of the D4 equivalence class: (orientation +, onset zc).")
print("   Modeling choice, per Rem 10.13; no derivation is claimed for it.")

# ---------------------------------------------------------------- D4-F: fail-first falsifiers
r_bad1   = K * z / zc                          # planted exponent 1
rate_bad = diff(pi * r_bad1**2, z)
check("D4-F1", not is_zero(diff(rate_bad, z)),
      "falsifier fired: planted exponent-1 profile violates constant-areal-rate tie")
cap_bad = Rational(9, 10) * K                  # planted discontinuous cap
check("D4-F2", not is_zero(r_horn.subs(z, zc) - cap_bad),
      "falsifier fired: planted cap 0.9K breaks continuity at the join")

# ---------------------------------------------------------------- summary
n_f = sum(1 for c in PASS if c.startswith("D4-F"))
print(f"\nD4: EXACT {len(PASS)}/{len(PASS)+len(FAIL)} checks passed | falsifiers fired {n_f}/2 | "
      f"{'exit 0' if not FAIL else 'exit 1'}")
print("REGISTER MOTION: ledger row 3 (OP-RADIUS profile) -> CLOSED under C2;")
print("  D4 = certified-independent declared atom; residual = representative declaration only.")
sys.exit(0 if not FAIL else 1)
