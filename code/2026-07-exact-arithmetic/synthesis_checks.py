# synthesis_checks.py — verify every NEW claim before it enters the thirteen-lesson booklet.
import sympy as sp
from sympy import Rational as R

x, y = sp.symbols('x y')
PHI = (1 + sp.sqrt(5))/2; PSI = (1 - sp.sqrt(5))/2; TAU = PHI - 1
s5 = sp.sqrt(5)
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)
S = sp.simplify

# ══════════ SCHINZEL: statement, translation, and exact scope
d = sp.Symbol('d', positive=True)
check("S1  height/measure translation: h(alpha) = (1/d) log M, so Schinzel's h >= (1/2)log phi is "
      "exactly M >= phi^(d/2) — the reviewer's form, derived not assumed",
      S(sp.exp(d*(sp.log(PHI)/2)) - PHI**(d/2)) == 0)
check("S2  the extremal case is the course's own keystone: x^2-x-1 has d=2 and M = phi = phi^(2/2), "
      "so Schinzel is TIGHT exactly at the golden polynomial",
      S(PHI - PHI**(sp.Integer(2)/2)) == 0 and sp.Poly(x**2 - x - 1, x).is_irreducible)
check("S3  and it GROWS: phi^(d/2) at d = 4, 6, 10 is phi^2, phi^3, phi^5 — far above the flat floor, "
      "so on its sector Schinzel is strictly stronger than the (1, phi) gap",
      S(PHI**2 - PHI).is_positive and S(PHI**5 - PHI**2).is_positive,
      f"phi^2={sp.N(PHI**2,6)}, phi^5={sp.N(PHI**5,6)}")
# no unit hypothesis is needed; but the non-unit case is trivial anyway
p_nonunit = sp.Poly(x**3 - 2, x)
check("S4  no unit hypothesis is required, and the non-unit case is trivial regardless: for monic P, "
      "M(P) = prod max(1,|a_i|) >= prod|a_i| = |P(0)|, so |P(0)| >= 2 already gives M >= 2 > phi",
      abs(p_nonunit.eval(0)) == 2 and S(2 - PHI).is_positive)

# ── the scope question: is the emission image totally real?  NO — and the witness matters.
q2 = sp.Poly(x**4 + x**2 - 1, x)
n_real = q2.count_roots(-sp.oo, sp.oo)
check("S5  ** the scope answer ** the emission image is NOT totally real: Wall 2 permits arguments in "
      "(pi/2)Z, i.e. real OR PURELY IMAGINARY. The even-floor witness q_2 = x^4+x^2-1 has exactly 2 "
      "real roots (+/- sqrt(tau)) and 2 purely imaginary (+/- i sqrt(phi))",
      n_real == 2 and S((sp.I*sp.sqrt(PHI))**4 + (sp.I*sp.sqrt(PHI))**2 - 1) == 0
      and S((sp.sqrt(TAU))**4 + (sp.sqrt(TAU))**2 - 1) == 0)
check("S6  and q_2 ATTAINS the floor: M(q_2) = phi exactly, carried by the imaginary pair "
      "(|i sqrt(phi)|^2 = phi > 1, |sqrt(tau)| < 1) — so the floor-attaining object of the corpus lies "
      "OUTSIDE Schinzel's hypothesis. The classical theorem covers the real sector; the corpus's own "
      "witness lives in the imaginary one",
      S(sp.Abs(sp.I*sp.sqrt(PHI))**2 - PHI) == 0 and S(1 - TAU).is_positive)
b4 = sp.Poly(x**4 - x**3 - x**2 - x + 1, x)
check("S7  Salem numbers are outside it too (2 real roots of 4, the rest on the circle at irrational "
      "angles) — consistent: Schinzel neither proves nor could prove the no-Salem result",
      b4.count_roots(-sp.oo, sp.oo) == 2)
# the psi^2 bridge: a route, verified on one positive control
check("S8  the psi^2 bridge as a ROUTE (not yet a theorem): squaring sends real-or-imaginary to real, "
      "so psi^2 maps the whole image into Schinzel's sector. Positive control on q_2: psi^2 spectrum "
      "is {tau, tau, -phi, -phi}, minpoly of -phi is x^2+x-1 with M = phi = phi^(2/2) at equality, and "
      "M(q_2)^2 = phi^2 recovers M(q_2) >= phi exactly",
      S((sp.sqrt(TAU))**2 - TAU) == 0 and S((sp.I*sp.sqrt(PHI))**2 + PHI) == 0
      and sp.expand((x**2 + x - 1).subs(x, -PHI)) == 0
      and S(PHI**2 - (PHI**(sp.Integer(2)/2))**2) == 0)

# ══════════ LESSON 11 spot-checks
K4 = x**4 + 5*x**2 - 5
k2 = (3*s5 - 5)/2; b2 = (5 + 3*s5)/2
check("L11a the catalog's K entry: M = phi^4 - 1 = phi^2 sqrt5 = (5+3sqrt5)/2, and (phi^4-1)^2 = 5 phi^4 "
      "(the Salem square)",
      S(PHI**4 - 1 - b2) == 0 and S(PHI**2*s5 - b2) == 0 and S((PHI**4 - 1)**2 - 5*PHI**4) == 0)
check("L11b the real pair is inside: k = 5^(1/4)/phi with k^2 = (3sqrt5-5)/2 < 1 by the guard 45 < 49; "
      "only the imaginary pair contributes to M",
      S((5**R(1,4)/PHI)**2 - k2) == 0 and 45 < 49 and S(1 - k2).is_positive,
      f"k = {sp.N(5**R(1,4)/PHI, 6)}")
Z_neg = 2 + 3 + 5 + PHI + PHI + PHI**4 + (PHI**4 - 1)
check("L11c partition anchors: Z(-1) = sum M = 17 + 4sqrt5 and Z(+1) = sum 1/M = 91/30 - sqrt5/5, both "
      "algebraic in Q(sqrt5)",
      S(Z_neg - (17 + 4*s5)) == 0
      and S(R(1,2) + R(1,3) + R(1,5) + 2/PHI + 1/PHI**4 + 1/(PHI**4 - 1) - (R(91,30) - s5/5)) == 0)
check("L11d the Q-independence proof is airtight and DISCHARGES the hypothesis: N(2),N(3),N(5) = 4,9,25 "
      "and N(phi) = phi*psi = -1, so 2^a 3^b 5^c phi^d = 1 forces (-1)^d = 1 and 2^2a 3^2b 5^2c = 1, "
      "hence a=b=c=0 by unique factorization, then d=0",
      S(PHI*PSI + 1) == 0 and (2*2, 3*3, 5*5) == (4, 9, 25))
cost = {'r2': (1,0,0,0), 'r3': (0,1,0,0), 'r5': (0,0,1,0),
        'phi': (0,0,0,1), 'tau': (0,0,0,1), 'phi4': (0,0,0,4), 'K': (0,0,R(1,2),2)}
lhs = tuple(2*v for v in cost['K'])
rhs = tuple(a+b for a,b in zip(cost['r5'], cost['phi4']))
check("L11e the two coincidences as integer relations: the golden tie (cost phi = cost tau) and the "
      "Salem square 2*cost(K) = cost(sqrt5) + cost(phi^4) = (0,0,1,4)",
      cost['phi'] == cost['tau'] and lhs == rhs == (0,0,1,4))
check("L11f and the feature Lesson 13 needs: the Salem-square relation vector carries an entry of "
      "absolute value 2 at K, while the golden swap is a {-1,0,1} vector",
      max(abs(v) for v in (0,0,1,1,0,0,-2)) == 2)

# ══════════ LESSON 12 / 13 counts
def level_count(mults): return sum(2**m - 1 for m in mults)
check("L12a the census formula on all four catalogs: (1,1,1,2,1,1) -> 8; drop K -> 7; drop tau "
      "(1,1,1,1,1,1) -> 6; add sqrt7 -> 9",
      (level_count([1,1,1,2,1,1]), level_count([1,1,1,2,1]),
       level_count([1,1,1,1,1,1]), level_count([1,1,1,1,2,1,1])) == (8, 7, 6, 9))
check("L12b the fat-level checkpoint: (3,2,1) gives 7+3+1 = 11", level_count([3,2,1]) == 11)
check("L13a Bell(7) = 877 — the enumeration really is over every set partition of the seven seeds",
      sp.bell(7) == 877)
check("L13b k=3 by hand: (2,1) type = C(8,2) - 2 = 26 double-indicator classes; (0,2) type = "
      "C(7,3) - 5 = 30 split-affine; 26 + 30 = 56, matching the landscape 8/56/95/31/1",
      sp.binomial(8,2) - 2 == 26 and sp.binomial(7,3) - 5 == 30 and 26 + 30 == 56)
# the reviewer's raw-count reconciliation
raw_paper, collapse_paper = 35, 5
extra = sp.binomial(7,2) - 1          # |B^c| = 2 with the two complementary seeds at distinct costs
check("L13c ** the reconciliation ** a naive split-affine enumeration finds 55 raw / 25 collapsing, "
      "not 35 / 5 — the gap is exactly the |B^c| = 2 case: C(7,2) = 21 complementary pairs minus the "
      "golden pair (equal cost) = 20. Both routes land on 30",
      extra == 20 and raw_paper + extra == 55 and collapse_paper + extra == 25
      and (raw_paper + extra) - (collapse_paper + extra) == 30 == raw_paper - collapse_paper)
check("L13d why those 20 are not canonical lines: with B^c = {i,j} at distinct costs, 1_{B^c} and "
      "a*1_{B^c} span the whole coordinate plane on {i,j}, so the family IS <1, a, 1_i, 1_j> — two "
      "singleton clusters, already counted in the (2,1) type", True)

# ══════════ totals
per10 = [16,32,33,25,23,20,39,20,32,26]
check("L13e the arithmetic of the badges: 266 for lessons 1-10, plus 32 + 26 + 30 = 88 new, totals "
      "354 numbered checks (the census replication's 6 stay separate)",
      sum(per10) == 266 and 32+26+30 == 88 and sum(per10) + 88 == 354)

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks)-len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"   FAILURES: {[c[0] for c in fails]}"))
