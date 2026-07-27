# lesson6_checks.py — exact verification for Lesson 6 (the relational layer).
# Discipline: no float crosses a decision boundary. Independent re-implementation:
# deposit scripts NOT consulted; criteria taken from the papers' stated definitions.
import sympy as sp
from sympy import Rational as R_
from fractions import Fraction
import time

x, y = sp.symbols('x y')
PHI = (1 + sp.sqrt(5)) / 2
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)

def companion(p):
    P = sp.Poly(p, x); d = P.degree(); assert P.LC() == 1
    Cm = sp.zeros(d, d)
    for i in range(d - 1): Cm[i + 1, i] = 1
    cf = P.all_coeffs()
    for i in range(d): Cm[i, d - 1] = -cf[d - i]
    return Cm

# ---------------- the probe core (papers' Appendix B, re-implemented) ----------------
def ratio_poly(p_expr):
    Rr = sp.resultant(p_expr.subs(x, y), sp.expand(p_expr.subs(x, x*y)), y)
    return sp.Poly(sp.expand(Rr), x).primitive()[1]

def mixed_ratio(p_expr, q_expr):
    Rr = sp.resultant(q_expr.subs(x, y), sp.expand(p_expr.subs(x, x*y)), y)
    return sp.Poly(sp.expand(Rr), x).primitive()[1]

def totients_upto(N):
    phi = list(range(N + 1))
    for i in range(2, N + 1):
        if phi[i] == i:
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi

def contacts(P):
    d = P.degree()
    N = 2 * d * d
    phi = totients_upto(N)
    hits = {}
    for m in range(1, N + 1):
        if phi[m] > d: continue
        cm = sp.Poly(sp.cyclotomic_poly(m, x), x)
        T, mult = P, 0
        while True:
            q, r = sp.div(T, cm)
            if r.is_zero: mult, T = mult + 1, q
            else: break
        if mult: hits[m] = mult
    return hits

# ================================================= A. gauge, rigidity, twist, golden
r3 = sp.root(2, 3)
check("A1  anchor instances: x^3-2 angles {0,1/3,2/3} (n=m=3); x^3+2 root (2^{1/3} e^{i pi/3})^3 = -2 "
      "so angles {1/6,1/2,5/6}: n=6, m=3 — anchor n=2m realized",
      sp.simplify((r3*sp.exp(sp.I*sp.pi/3))**3 + 2) == 0
      and sp.simplify((-r3)**3 + 2) == 0)
check("A2  sign twist: T(x^3+2) = (-1)^3 p(-x) = x^3 - 2  (ledger B)",
      sp.expand(-((-x)**3 + 2) - (x**3 - 2)) == 0)
check("A3  twist coset arithmetic (ledger L): odd e, odd m => (e+m)/2m in (1/m)Z, m in {3,5,7,9}",
      all((e + m) % 2 == 0 for m in (3, 5, 7, 9) for e in range(1, 2*m, 2)))
check("A4  the golden internal relation (ledger J): phi/phi' = -phi^2, so t_rel(phi, phi') = 1/2",
      sp.simplify(PHI/((1 - sp.sqrt(5))/2) + PHI**2) == 0)
# group drop E vs its sign-mirror K (the helix quartic)
E5 = x**4 + 5*x**2 + 5
check("A5  ledger E: x^4+5x^2+5 irreducible; x^2-roots (-5±sqrt5)/2 both negative (sqrt5 < 5); "
      "squared moduli (5∓sqrt5)/2 both > 1 (sqrt5 < 3) => M = |p(0)| = 5; angles {1/4,3/4}: "
      "absolute Z/4, relational Z/2 — the group drop",
      sp.Poly(E5, x).is_irreducible
      and ((-5 + sp.sqrt(5))/2).is_negative and ((5 - sp.sqrt(5))/2 - 1).is_positive)
K5 = x**4 + 5*x**2 - 5
check("A5b mirror K = x^4+5x^2-5: one x^2-root positive ((3sqrt5-5)/2 > 0) => real pair + imaginary "
      "pair: angles {0,1/2,1/4,3/4}, full Delta = Z/4 — type finer than signature",
      ((-5 + 3*sp.sqrt(5))/2).is_positive and ((-5 - 3*sp.sqrt(5))/2).is_negative)

# ================================================= B. ledger replications (contact scans)
t0 = time.time()
sigA = contacts(ratio_poly(x**3 - 2))
sigB = contacts(ratio_poly(x**3 + 2))
check("B1  ledger A/B: x^3-2 and x^3+2 both {Phi_1^3, Phi_3^3} — the probe is gauge-blind",
      sigA == {1: 3, 3: 3} and sigB == sigA)
check("B2  ledger C (sharpness of pinning): x^4-2, one shell => torsion ratios: {Phi_1^4, Phi_2^4, Phi_4^4}",
      contacts(ratio_poly(x**4 - 2)) == {1: 4, 2: 4, 4: 4})
check("B3  ledger D (their harness incident): q2 = x^4+x^2-1 -> {Phi_1^4, Phi_2^4} (ordered pairs)",
      contacts(ratio_poly(x**4 + x**2 - 1)) == {1: 4, 2: 4})
check("B4  ledger E/F: x^4+5x^2+5 and K both {Phi_1^4, Phi_2^4} — identical signatures, different types",
      contacts(ratio_poly(E5)) == {1: 4, 2: 4} and contacts(ratio_poly(K5)) == {1: 4, 2: 4})
b4 = x**4 - x**3 - x**2 - x + 1
Rb4 = ratio_poly(b4)
check("B5  ledger G: beta_4 -> {Phi_1^4} ONLY (complete, bound 512): circle block inert",
      contacts(Rb4) == {1: 4})
Cb4 = companion(b4)
tworoute = sp.Poly(sp.expand((Cb4 * Cb4.inv()).charpoly(x).as_expr()), x)  # placeholder; real check below
CC = companion(b4)
K16 = sp.Matrix(16, 16, lambda i, j: sum(CC[i//4, k]*(CC.inv())[k % 1, 0] for k in range(0)) if False else 0)
# two-route properly: charpoly of C (x) C^{-1}
Cinv = CC.inv()
def kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i//B.rows, j//B.cols]*B[i % B.rows, j % B.cols])
F2 = sp.Poly(sp.expand(kron(CC, Cinv).charpoly(x).as_expr()), x).primitive()[1]
check("B5b two-route (ledger G): primpart charpoly(C (x) C^{-1}) = Rat_{beta4} exactly",
      sp.expand(F2.as_expr() - Rb4.as_expr()) == 0 or sp.expand(F2.as_expr() + Rb4.as_expr()) == 0)
# Z* — hypothesis-necessity witness (Pisot paper Prop 3.3)
Zs = x**4 - 3*x**2 + 1
check("B6  Z* = x^4-3x^2+1 = (x^2-x-1)(x^2+x-1) exactly; torsion ratio -1 = phi/(-phi) across factors; "
      "scan {Phi_1^4, Phi_2^4} — irreducibility AND the disjoint-modulus clause are both necessary",
      sp.expand(Zs - (x**2 - x - 1)*(x**2 + x - 1)) == 0
      and contacts(ratio_poly(Zs)) == {1: 4, 2: 4})
# Salem certifications (ledger O) by trace-fold Sturm pattern (1, 0, n/2 - 1)
def salem_cert(p):
    P = sp.Poly(p, x); n = P.degree()
    if P.all_coeffs() != P.all_coeffs()[::-1] or not P.is_irreducible: return False
    tt = sp.symbols('tt')
    # fold: p(x)/x^{n/2} = T(x + 1/x)
    zk = {0: sp.Integer(2), 1: tt}
    for k in range(2, n//2 + 1):
        zk[k] = sp.expand(tt*zk[k-1] - zk[k-2])
    cf = P.all_coeffs()
    T = sp.expand(sum(cf[i]*zk[n//2 - i] for i in range(n//2)) + cf[n//2])
    TP = sp.Poly(T, tt)
    return (TP.count_roots(2, sp.oo), TP.count_roots(-sp.oo, -2), TP.count_roots(-2, 2)) == (1, 0, n//2 - 1)
S6p = x**6 - x**4 - x**3 - x**2 + 1
S8p = x**8 - x**5 - x**4 - x**3 + 1
Lp = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1
check("B7  ledger O: beta4, S6, S8, Lehmer all Salem by trace-fold Sturm (1, 0, n/2-1) — exact",
      all(salem_cert(p) for p in (b4, S6p, S8p, Lp)))
check("B8  plastic number (ledger W): x^3 - x - 1 -> {Phi_1^3}: relationally inert",
      contacts(ratio_poly(x**3 - x - 1)) == {1: 3})
sigI = contacts(mixed_ratio(Lp, b4))
check("B9  ledger I: mixed Rat_{L, beta4} (deg 40): contacts = EMPTY — the two keystones are not circle-locked",
      sigI == {} and sp.gcd(sp.Poly(Lp, x), sp.Poly(b4, x)) == sp.Poly(1, x))
sigH = contacts(ratio_poly(Lp))
check("B10 ledger H (full): Lehmer Rat_L, degree 100, complete bound 20000 -> {Phi_1^10} only",
      sigH == {1: 10}, f"elapsed so far {time.time()-t0:.1f}s")

# ================================================= C. imported rigidity & the measure identity
p41 = x**4 - x + 1
KX = kron(companion(p41), companion(p41))
FX = sp.Poly(sp.expand(KX.charpoly(x).as_expr()), x)
q2img = sp.Poly(sp.expand((companion(p41)**2).charpoly(x).as_expr()), x)
check("C1  x^4-x+1: charpoly(C (x) C) = S6^2 * (x^4+2x^2-x+1) exactly, and the quartic factor IS psi^2(p)",
      sp.expand(FX.as_expr() - sp.expand(S6p**2 * (x**4 + 2*x**2 - x + 1))) == 0
      and sp.expand(q2img.as_expr() - (x**4 + 2*x**2 - x + 1)) == 0)
check("C1b measure-identity chain links: p(0)=1, zero real roots, non-reciprocal, irreducible "
      "=> M(x^4 - x + 1) = tau_{S6} (a Mahler value in (theta0, phi) at irrational angles)",
      p41.subs(x, 0) == 1 and sp.Poly(p41, x).count_roots(-sp.oo, sp.oo) == 0
      and sp.Poly(p41, x).all_coeffs() != sp.Poly(p41, x).all_coeffs()[::-1]
      and sp.Poly(p41, x).is_irreducible)
Fb = sp.Poly(sp.expand(kron(Cb4, Cb4).charpoly(x).as_expr()), x)
gsq = sp.gcd(Fb, Fb.diff())
sqf = sp.quo(Fb, gsq)
check("C2  ledger S (beta4 (x) beta4 non-inert witness): (x-1)-multiplicity exactly 4; deg gcd(F,F') = 7; "
      "3 distinct real roots, all positive (Sturm)",
      contacts_dummy := True and
      (lambda: (lambda m1: m1 == 4)(max(k for k in range(1, 6)
            if sp.rem(Fb, sp.Poly((x-1)**k, x)).is_zero)))()
      and sp.degree(gsq) == 7
      and sqf.count_roots(-sp.oo, sp.oo) == 3 and sqf.count_roots(0, sp.oo) == 3)

# NOTE. Section D (the exhaustive [-2,2]^5 quintic census) previously lived here and was
# superseded during authoring by the standalone lesson6_census.py, which carries the corrected
# complex-interval unpacking. The dead copy is removed rather than left to crash: this file is
# checks A-C, which are the 20 the Lesson 6 badge names. The census is a separate artifact with
# its own six checks and its own manifest row.

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
