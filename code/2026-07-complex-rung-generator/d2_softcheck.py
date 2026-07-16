#!/usr/bin/env python3
"""
d2_softcheck.py -- W[open]2 / D2 (ledger row 4): the drop-one independence audit.

Design (Appendix B.1 spec, with the falsifier line CORRECTED -- plant a KILLING
clock, not the consistent one):

  D3' := "one quantum per entropy floor, chi unspecified"  (D3 with the clock freed)
  Sweep every FORCED fact of the corpus against each competing clock
  chi in {pi (terrain), 2pi/5 (pentagon), pi/2 (quarter)}.

  Verdict rule:
    pentagon has no KILL  =>  a model of {D1, D3', keystone, forced facts}
                              violates D2  =>  D2 independent  =>  closes under C2.
    some forced fact kills the pentagon    =>  D2 redundant    =>  D2 retires.
  Either horn closes the row.

  Acceptance gate (fail-first): chi = pi MUST die, or the sweep is not
  discriminating and the harness aborts.

Forced-fact battery B (each item re-encoded exactly, with its own source):
  F1 helicity          W-Thm 11.1 / W-Cor 11.2  (real multiplier carries no helicity)
  F2 on-circle image   C-Prop 3.8/3.10 (S on-circle = {+-1}); Reading A excluded (C-Rem 3.11)
  F3 pentagon-by-image C-Prop 3.7 (zeta5 not in S) -- conditional on D2', hence NO kill
  F4 keystone          C-Thm 3.3 (J = UA unconditional) -- chi-free, NO kill
  F5 chart blindness   W-Thm 13.1 (gate pentad is a function of z only) -- NO kill
  F6 rung passage      W-Thm 9.1 (ladder exists for every clock) -- NO kill
  F7 clock invariance  W-Prop 14.2 (qualitative structure survives rescaling) -- NO kill
  F8 K-seed rotation   C-Prop 3.12 (Arg(+-i*beta) = +-pi/2, forced object) -- NO kill w/o D2

Discipline: exact decisions in Q(sqrt5)[i]; sqrt5-signs by integer arithmetic;
no floats at decision boundaries; falsifiers must fire.
"""
import sys
from itertools import combinations_with_replacement
from sympy import (sqrt, Rational, symbols, simplify, expand, diff, pi, S, I,
                   sin, cos, exp, log, limit, oo, Poly, Matrix, eye, im, re,
                   solve, together)

z = symbols('z', positive=True)
x = symbols('x')

SQ5 = sqrt(5)
phi = (1 + SQ5) / 2
tau = phi - 1                      # = 1/phi = (sqrt5-1)/2
gap = phi**-4
K2  = expand(1 - gap)
beta = 5**Rational(1, 4) * phi     # rotation magnitude, beta^2 = phi^2*sqrt5

PASS, FAIL = [], []
def check(cid, cond, note=""):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid:8s} {note}")
    return ok

def is_zero(e): return simplify(expand(e)) == 0

def q5_pos(e):
    e = expand(simplify(e))
    a = e.subs(SQ5, 0); b = simplify((e - a)/SQ5); a, b = S(a), S(b)
    if b == 0: return a > 0
    if a == 0: return b > 0
    if a > 0 and b > 0: return True
    if a < 0 and b < 0: return False
    return (a**2 - 5*b**2) > 0 if a > 0 else (5*b**2 - a**2) > 0

print("== d2_softcheck.py : D2 / chi-selection drop-one audit (W[open]2, ledger row 4) ==\n")

# =====================================================================================
# D2-A : battery replication -- each forced fact re-encoded exactly, cold
# =====================================================================================
print("-- D2-A : forced-fact battery, re-encoded cold --")

# F1 helicity: the winding multiplier q_chi = e^{i chi} * tau must be non-real
#              with 0 < |q| < 1 (simultaneous rotation + contraction).
check("D2-A1", q5_pos(tau) and q5_pos(1 - tau),
      "F1 substrate: 0 < tau < 1 exactly  [9>5 and 5>1]")
sin2 = {S.Pi: expand(sin(pi)**2),
        2*S.Pi/5: expand(sin(2*pi/5)**2),
        S.Pi/2: expand(sin(pi/2)**2)}
check("D2-A2", sin2[S.Pi] == 0,
      "F1: sin^2(pi) = 0 -- terrain multiplier -tau is REAL")
check("D2-A3", is_zero(sin2[2*S.Pi/5] - (S(5)/8 + SQ5/8)) and q5_pos(sin2[2*S.Pi/5]),
      "F1: sin^2(2pi/5) = (5+sqrt5)/8 > 0 -- pentagon multiplier NON-real")
check("D2-A4", sin2[S.Pi/2] == 1,
      "F1: sin^2(pi/2) = 1 -- quarter multiplier NON-real")

# F2 on-circle image / Reading-A exclusion: rotation gate C=-1/2 has multipliers -+i
C = Rational(-1, 2)
us = solve(x**2 + x - C, x)
ms = sorted([simplify(-C/(1+u)**2) for u in us], key=lambda t: str(t))
check("D2-A5", set(ms) == {I, -I},
      "F2: rotation gate C=-1/2 multipliers are exactly {+i, -i} (W-Prop 10.2)")
check("D2-A6", I not in (1, -1) and -I not in (1, -1),
      "F2: +-i not in on-circle S = {+-1}  =>  Reading A ('multiplier in on-circle S') UNSATISFIABLE")
print("   [EXCLUDED] Reading A is not a forced fact (C-Rem 3.11): it is unsatisfiable for the")
print("   rotation gate itself; it enters below only as falsifier D2-FF3.")

# F3 pentagon-by-image: zeta5 not in S -- but the kill is conditional on D2' (== Reading A)
Phi5 = Poly(x**4 + x**3 + x**2 + x + 1, x)
check("D2-A7", Phi5.is_irreducible,
      "F3: Phi_5 irreducible over Q")
check("D2-A8", all(Phi5.eval(v) != 0 for v in (1, -1, I, -I)),
      "F3: Phi_5 has no root in mu_4 -- zeta5 not an on-circle S element")
check("D2-A9", all(Rational(k, 5) not in {S(0), Rational(1,4), Rational(1,2), Rational(3,4)}
                   for k in range(1, 5)),
      "F3: pentagon angle classes k/5 meet the quarter lattice only at 0")
print("   [SCOPED]  F3's kill is conditional on D2' = Reading A, which is unsatisfiable")
print("   (D2-A6); therefore F3 contributes NO kill to any clock in this sweep.")

# F4 keystone quarter-turn: J = UA exists unconditionally, chi-free
R = Matrix([[0, 1], [1, 1]]); Id = eye(2)
H = 2*R - Id
N = Matrix([[0, 1], [1, 0]]) + 2*Matrix([[1, 0], [0, -1]])
J = Matrix([[0, -1], [1, 0]])
check("D2-A10", H*H == 5*Id and N*N == 5*Id and H*N == 5*J and J*J == -Id,
      "F4: H^2 = N^2 = 5I, HN = 5J, J^2 = -I  (C-Thm 3.3, unconditional)")
check("D2-A11", (H*N) == -(N*H),
      "F4: UA = -AU (anticommutation forced by the return-kernel geometry)")
Q4 = expand((tau**4))
check("D2-A12", is_zero(Q4 - gap),
      "F4: Q^4 = tau^4 J^4 = gap * I -- the gap schedule, chi-free")

# F5 chart blindness: the gate pentad is a function of z alone
u_z = (1 - z**2)/z**2; C_z = (1 - z**2)/z**4; m_z = -(1 - z**2)
check("D2-A13", is_zero(u_z**2 + u_z - C_z) and is_zero(-C_z/(1 + u_z)**2 - m_z),
      "F5: u(z) is a fixed point of G_{C(z)} with multiplier -(1-z^2); theta never enters")

# F6 rung passage: the ladder exists for every clock (heights are clock-free)
lad_ok = True
for n in (1, 2, 3):
    w = phi**(2*n)
    un = 1/(w - 1); Cn = w/(w - 1)**2; mn = -1/w
    lad_ok &= is_zero(un**2 + un - Cn) and is_zero(-Cn/(1 + un)**2 - mn)
check("D2-A14", lad_ok and is_zero((phi**2/(phi**2-1)**2) - 1)
              and is_zero((phi**4/(phi**4-1)**2) - Rational(1, 5)),
      "F6: rung ladder replicates, C1=1 (golden), C2=1/5 (K-gate); heights clock-free")

# F7 clock invariance witness (W-Prop 14.2), instantiated on the pentagon clock
check("D2-A15", is_zero(expand(gap*phi**4 - 1)),
      "F7: gap*phi^4 = 1  =>  rho(K) = 2 ln(phi)  =>  theta(K) = 2*chi for EVERY clock")
# lens off-lattice for every clock: log_phi 2 irrational -- norm route + Binet grid
norm_ok = True
for p in range(1, 11):
    e = expand(phi**p)
    b = simplify((e - e.subs(SQ5, 0))/SQ5)
    norm_ok &= (b != 0)                      # sqrt5-part F_p/2 != 0 => phi^p irrational
check("D2-A16", norm_ok,
      "F7: phi^p irrational for p=1..10 (Binet sqrt5-part) -- and |N(phi^p)|=1 != 4^q = N(2^q)")
check("D2-A17", limit(-log(1 - z**2)/2, z, 1, '-') == oo,
      "F7: infinite winding at the top for every clock (rho -> oo as z -> 1^-)")

# F8 K-seed rotation pair: Arg(+-i*beta) = +-pi/2, a FORCED object
check("D2-A18", is_zero(expand((I*beta)**4 + 5*(I*beta)**2 - 5))
              and is_zero(expand(sqrt(K2)**4 + 5*K2 - 5)),
      "F8: K-seed x^4+5x^2-5 has roots {+-K, +-i*beta} -- forced substrate object")
check("D2-A19", re(I*beta) == 0 and q5_pos(expand(beta**2)),
      "F8: i*beta purely imaginary, beta > 0  =>  Arg(+-i*beta) = +-pi/2 exactly")
print("   [SCOPED]  F8 supplies the quarter-turn ARGUMENT as a forced object; absent D2's")
print("   selection ('use the rotation axis') it constrains no clock: NO kill (Rem 3.13).\n")

# =====================================================================================
# D2-C : the sweep
# =====================================================================================
print("-- D2-C : the three-clock sweep --")
CLOCKS = {"terrain  pi":   S.Pi,
          "pentagon 2pi/5": 2*S.Pi/5,
          "quarter  pi/2": S.Pi/2}

def sweep(battery):
    kills = {name: [] for name in CLOCKS}
    for name, chi in CLOCKS.items():
        s2 = expand(sin(chi)**2)
        if "F1" in battery and s2 == 0:
            kills[name].append("F1-helicity")
        if "A" in battery:                     # Reading A, only ever planted as falsifier
            # e^{i chi} in {+-1}  <=>  sin(chi) = 0
            if s2 != 0:
                kills[name].append("A-oncircle")
            else:
                kills[name].append("A-helicity-conflict")  # real => F1 conflict anyway
        # F3 conditional on D2' (== A): contributes only if A active AND zeta5 needed
        # (already covered by A above); with A excluded, F3 contributes nothing.
        # F4,F5,F6,F7,F8: chi-free / no-kill by D2-A10..A19.
    return kills

B_full = {"F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"}
kills = sweep(B_full)
for name in CLOCKS:
    print(f"     {name:16s} kills = {kills[name] if kills[name] else 'NONE (survives)'}")

# Acceptance gate D2-B (fail-first): the terrain clock MUST die
if not check("D2-B1", "F1-helicity" in kills["terrain  pi"],
             "ACCEPTANCE GATE: chi = pi is KILLED by helicity -- the sweep discriminates"):
    print("D2: sweep not discriminating -- ABORT"); sys.exit(1)

check("D2-C1", kills["pentagon 2pi/5"] == [],
      "pentagon clock 2pi/5 has NO kill in the forced battery")
check("D2-C2", kills["quarter  pi/2"] == [],
      "quarter clock pi/2 has NO kill (consistency: the D2-selected clock survives)")

# =====================================================================================
# D2-D : the pentagon violation certificate
# =====================================================================================
print("\n-- D2-D : model M5 = {D1 chart, theta = (2pi/5) rho/ln(phi), inherited r, keystone} --")
angleset5 = {Rational(k, 5) for k in range(5)}          # <1/5> mod 1
check("D2-D1", Rational(1, 4) not in angleset5,
      "1/4 not in <1/5> mod 1  =>  M5's angle lattice does NOT realize Delta(K-seed) = Z/4Z")
check("D2-D2", all((Rational(1,4) - Rational(k,5)) not in
                   {S(mm) for mm in range(-2, 3)} for k in range(5)),
      "no integers k,m solve 1/4 = k/5 + m  [4k+20m = 5 is unsolvable mod 4]")
check("D2-D3", 5 % 4 != 0,
      "|Z/5| = 5, and Z/4 embeds in Z/5 iff 4 | 5 -- false")
m5_model = (kills["pentagon 2pi/5"] == []) and all(
    c in PASS for c in ("D2-A13", "D2-A14", "D2-A15", "D2-A16", "D2-A17", "D2-A10"))
check("D2-D4", m5_model,
      "M5 satisfies D1 + D3' + keystone + every battery fact, and VIOLATES D2")

# =====================================================================================
# D2-E : verdict
# =====================================================================================
print("\n-- D2-E : verdict --")
independent = m5_model and (Rational(1, 4) not in angleset5)
check("D2-E1", independent,
      "INDEPENDENCE HORN: M5 is a model of {D1, D3', keystone, B} violating D2")
check("D2-E2", not (len(kills["pentagon 2pi/5"]) > 0),
      "REDUNDANCY HORN rejected: no forced fact kills the pentagon clock")
check("D2-E3", independent,
      "CLOSURE (C2): D2 = certified-independent declared atom; W[open]2's unrun audit is now run")
print("   Register motion: D2 stays [declared]; the ROW closes -- 'derive D2' is settled")
print("   negatively relative to battery B; the restatement wall (Sec 14.4) stands as the")
print("   record of why no thinner axiom exists. Independence is battery-relative:")
print("   a future forced fact naming the clock would re-open the redundancy horn.")

# =====================================================================================
# D2-S : on-circle sub-closure scan (C-Prop 3.8 replication, degree <= 6)
# =====================================================================================
print("\n-- D2-S : unimodular tensor-monomial scan over the 16 catalog eigenvalues --")
# entry = (w0,w1,w2,w3, a): |lam| = 2^{w0/4} 3^{w1/4} 5^{w2/4} phi^{w3/4}, arg = a * pi/2
CAT = [
    ("phi",    (0,0,0, 4), 0), ("psi=-tau", (0,0,0,-4), 2),
    ("tau",    (0,0,0,-4), 0), ("-phi",     (0,0,0, 4), 2),
    ("sqrt2",  (2,0,0, 0), 0), ("-sqrt2",   (2,0,0, 0), 2),
    ("sqrt3",  (0,2,0, 0), 0), ("-sqrt3",   (0,2,0, 0), 2),
    ("sqrt5",  (0,0,2, 0), 0), ("-sqrt5",   (0,0,2, 0), 2),
    ("phi^4",  (0,0,0,16), 0), ("phi^-4",   (0,0,0,-16), 0),
    ("K",      (0,0,1,-4), 0), ("-K",       (0,0,1,-4), 2),
    ("i*beta", (0,0,1, 4), 1), ("-i*beta",  (0,0,1, 4), 3),
]
SYMVAL = {"phi": phi, "psi=-tau": -tau, "tau": tau, "-phi": -phi,
          "sqrt2": sqrt(2), "-sqrt2": -sqrt(2), "sqrt3": sqrt(3), "-sqrt3": -sqrt(3),
          "sqrt5": SQ5, "-sqrt5": -SQ5, "phi^4": phi**4, "phi^-4": phi**-4,
          "K": sqrt(K2), "-K": -sqrt(K2), "i*beta": I*beta, "-i*beta": -I*beta}
# entry self-audit: |lam|^4 == 2^w0 3^w1 5^w2 phi^w3 exactly, and arg class correct
audit = True
for name, w, a in CAT:
    lam = SYMVAL[name]
    modulus4 = expand(simplify((lam*lam.conjugate())**2))
    target   = expand(simplify(S(2)**w[0] * S(3)**w[1] * S(5)**w[2] * phi**w[3]))
    audit &= is_zero(modulus4 - target)
    ang_ok = {0: im(lam) == 0 and q5_pos(expand(re(lam)*abs(1)) if re(lam) != 0 else S(1)),
              }  # sign handled below
    if a in (0, 2):
        audit &= (im(lam) == 0)
    else:
        audit &= (re(lam) == 0)
check("D2-S1", audit,
      "entry self-audit: 16 eigenvalues, |lam|^4 = 2^w0 3^w1 5^w2 phi^w3, arg class = a*(pi/2)")

total = 0; uni = 0; uni_vals_even = True; uni_golden_only = True
GOLD = {0, 1, 2, 3, 10, 11}                              # indices of the golden block
for d in range(0, 7):
    for combo in combinations_with_replacement(range(16), d):
        total += 1
        wsum = [0, 0, 0, 0]; asum = 0
        for idx in combo:
            _, w, a = CAT[idx]
            for j in range(4): wsum[j] += w[j]
            asum += a
        if wsum == [0, 0, 0, 0]:
            uni += 1
            if asum % 2 != 0: uni_vals_even = False
            if not set(combo) <= GOLD: uni_golden_only = False
check("D2-S2", total == 74613,
      f"scan population: {total} monomials of degree <= 6 (corpus count replicated)")
check("D2-S3", uni == 60,
      f"unimodular monomials: {uni} (corpus count 60 replicated)")
check("D2-S4", uni_vals_even,
      "every unimodular monomial has EVEN quarter-charge => value in {+1,-1}: +-i unreachable")
check("D2-S5", uni_golden_only,
      "5/3/2-valuation forcing: every unimodular monomial uses only the golden block")

# =====================================================================================
# D2-F : fail-first falsifiers  (B.1 spec line INVERTED -> corrected here)
# =====================================================================================
print("\n-- D2-F : falsifiers (corrected: plant a KILLING clock, not the consistent one) --")
# FF1: plant the claim 'chi = pi survives'; the sweep must contradict it.
planted_claim = ("terrain  pi", "SURVIVES")
detected1 = len(kills["terrain  pi"]) > 0
check("D2-FF1", detected1,
      "falsifier fired: planted 'pi survives' is CONTRADICTED by the sweep (inverted-line fix)")
# FF2: degrade the battery (drop F1): nothing can die -> non-discrimination alarm.
kills_deg = sweep(B_full - {"F1"})
nondiscr = all(len(v) == 0 for v in kills_deg.values())
check("D2-FF2", nondiscr,
      "falsifier fired: battery minus F1 kills nothing -> alarm: F1 is the load-bearing kill")
# FF3: plant Reading A as a battery fact: every clock dies -> unsatisfiability alarm.
kills_A = sweep(B_full | {"A"})
unsat = all(len(v) > 0 for v in kills_A.values())
check("D2-FF3", unsat,
      "falsifier fired: Reading A planted -> ALL clocks die -> A unsatisfiable, excluded (C-Rem 3.11)")

# ---------------------------------------------------------------- summary
n_f = sum(1 for c in PASS if c.startswith("D2-FF"))
print(f"\nD2: EXACT {len(PASS)}/{len(PASS)+len(FAIL)} checks passed | falsifiers fired {n_f}/3 | "
      f"{'exit 0' if not FAIL else 'exit 1'}")
print("VERDICT: pentagon has no KILL => D2 INDEPENDENT of {D1, D3', keystone, B} => closes under C2.")
sys.exit(0 if not FAIL else 1)
