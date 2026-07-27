# revision_checks.py — verify each substantive repair BEFORE it enters the booklet.
import hashlib, os
from pathlib import Path

HERE = Path(__file__).resolve().parent   # sibling scripts resolve relative to this file
import sympy as sp
from fractions import Fraction as F

x, s, t, L, q, A = sp.symbols('x s t L q A', positive=True)
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)

# ── R1  the cost-floor row: 2c log mu_S at the frame-shift c, NOT the spectral gap
muS = sp.CRootOf(sp.Poly(sp.Symbol('y')**3 - sp.Symbol('y') - 1, sp.Symbol('y')), 0)
C = sp.Symbol('C', positive=True)
c_fs = sp.sqrt(1 + 4*C)/(2*C)
floor_fs = sp.simplify(2*c_fs*sp.log(muS))
at1 = sp.simplify(floor_fs.subs(C, 1))
check("R1  frame-shift cost floor = (sqrt(1+4C)/C)*log(mu_S) = lambda*log(mu_S); at C=1 that is "
      "sqrt5*log(mu_S), NOT sqrt5 — the two differ",
      sp.simplify(floor_fs - (sp.sqrt(1+4*C)/C)*sp.log(muS)) == 0
      and sp.simplify(at1 - sp.sqrt(5)*sp.log(muS)) == 0
      and sp.N(at1, 8) != sp.N(sp.sqrt(5), 8),
      f"sqrt5*log(muS) = {sp.N(at1, 6)}  vs  sqrt5 = {sp.N(sp.sqrt(5), 6)}")
check("R1b the other two rows are consistent with the same column rule floor = lambda*log(mu_S): "
      "2*log(muS) = 0.5624, 8*log(muS) = 2.2496, 16*log(muS) = 4.4992",
      abs(sp.N(2*sp.log(muS), 10) - sp.Float('0.5624')) < sp.Float('5e-5')
      and abs(sp.N(8*sp.log(muS), 10) - sp.Float('2.2496')) < sp.Float('5e-5')
      and abs(sp.N(16*sp.log(muS), 10) - sp.Float('4.4992')) < sp.Float('5e-5'))

# ── R2  real indefiniteness does NOT imply rational isotropy
qf = 2*sp.Symbol('a')**2 - 4*sp.Symbol('b')**2
check("R2  the reviewer's counterexample stands: 2a^2 - 4b^2 has signature (1,1) yet a^2 = 2b^2 has "
      "no nonzero rational solution (sqrt2 irrational) — so 'indefinite over R' does NOT give a "
      "rational isotropic vector",
      sp.Matrix([[2, 0], [0, -4]]).eigenvals() == {sp.Integer(2): 1, sp.Integer(-4): 1}
      and not sp.sqrt(2).is_rational)
# but our two witnesses ARE rationally isotropic — verify by exhibiting the vectors
def power_gram(p):
    P = sp.Poly(p, x); d = P.degree()
    Cm = sp.zeros(d, d)
    for i in range(d - 1): Cm[i + 1, i] = 1
    cf = P.all_coeffs()
    for i in range(d): Cm[i, d - 1] = -cf[d - i]
    tr = {0: sp.Integer(d)}; Pw = sp.eye(d)
    for k in range(1, 2*d - 1):
        Pw = Pw*Cm; tr[k] = sp.trace(Pw)
    return sp.Matrix(d, d, lambda i, j: tr[i + j])
Gi, G3 = power_gram(x**2 + 1), power_gram(x**3 - 2)
vi, v3 = sp.Matrix([1, 1]), sp.Matrix([0, 1, 0])
check("R2b both of the booklet's witnesses happen to BE rationally isotropic — Q(i): (1,1) is the "
      "element 1+i with Tr((1+i)^2) = 0; Q(cbrt2): (0,1,0) is theta with Tr(theta^2) = 0. The "
      "CONCLUSION survives; only the ARGUMENT (real => rational) was invalid",
      (vi.T*Gi*vi)[0, 0] == 0 and (v3.T*G3*v3)[0, 0] == 0)
check("R2c and the negative-gain mode needs no isotropy at all: over Q(i) with basis {1}, the "
      "residual of 1+i is (0,1) with norm -2 < 0 — real indefiniteness alone suffices for THAT",
      (sp.Matrix([0, 1]).T*Gi*sp.Matrix([0, 1]))[0, 0] == -2)
check("R2d so total reality is a SUFFICIENT conservative precondition, decided by one exact count: "
      "real roots == degree. Q(sqrt2+sqrt3) 4/4 PASS; Q(cbrt2) 1/3 and Q(i) 0/2 REJECT",
      sp.Poly(x**4 - 10*x**2 + 1, x).count_roots(-sp.oo, sp.oo) == 4
      and sp.Poly(x**3 - 2, x).count_roots(-sp.oo, sp.oo) == 1
      and sp.Poly(x**2 + 1, x).count_roots(-sp.oo, sp.oo) == 0)

# ── R3  Fisher cannot go Lorentzian: the Hermitian companion is the positive-definite continuation
Mi = sp.Matrix([[1, sp.I], [1, -sp.I]])
G2i = sp.simplify(Mi.conjugate().T*Mi)
check("R3  Q(i): the bilinear trace form is diag(2,-2) (indefinite) while the Hermitian embedding "
      "form M*M = 2I is positive definite — the latter is the legitimate Fisher continuation, the "
      "former is not a Fisher metric at all",
      Gi == sp.diag(2, -2) and G2i == 2*sp.eye(2)
      and all(ev > 0 for ev in G2i.eigenvals()) and any(ev < 0 for ev in Gi.eigenvals()))
Cs = sp.Symbol('C_r', real=True)
c_at_neg = (sp.sqrt(1 + 4*Cs)/(2*Cs)).subs(Cs, -1)
check("R3b what actually happens to c past the flip: the frame-shift FORMULA returns a non-real "
      f"value at C = -1 ({sp.simplify(c_at_neg)}), i.e. the canonicalization has no real solution "
      "there — a domain boundary of the recipe, not a complex Fisher metric",
      sp.simplify(sp.im(c_at_neg)) != 0)

# ── R4  the ad_R argument needs a scalar clause
a_, b_, c_, d_ = sp.symbols('a b c d')
R = sp.Matrix([[a_, b_], [c_, d_]])
adR = lambda X: R*X - X*R
check("R4  ad_R(I) = ad_R(R) = 0 symbolically for ANY 2x2 R; if R is nonscalar, I and R are "
      "independent so dim ker >= 2; if R = mI then ad_R is identically zero (ker = all of M_2, "
      "dim 4). Either way the 0-eigenvalue has multiplicity >= 2",
      adR(sp.eye(2)) == sp.zeros(2, 2) and sp.simplify(adR(R)) == sp.zeros(2, 2)
      and sp.simplify((R.subs({b_: 0, c_: 0, d_: a_})*sp.Matrix([[0, 1], [0, 0]])
                       - sp.Matrix([[0, 1], [0, 0]])*R.subs({b_: 0, c_: 0, d_: a_}))) == sp.zeros(2, 2))

# ── R5  Lehmer irreducibility is a separate certificate, not a corollary of the root count
Lp = sp.Poly(x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1, x)
fake = sp.Poly(sp.expand((x**2 - 3*x + 1)*(x**2 - 3*x + 1)*(x**6 + 1)), x)
check("R5  L is irreducible over Q — an independent exact factorization check, NOT implied by the "
      "root distribution. Witness that the implication fails: a reciprocal polynomial can carry "
      "Salem-like geometry and still factor",
      Lp.is_irreducible and not fake.is_irreducible
      and Lp.all_coeffs() == Lp.all_coeffs()[::-1])

# ── R6  the 1-D law, written with its limits
u = q*sp.integrate((L - s)/A, (s, 0, x))
check("R6  u(x) = q*int_0^x (L-s)/A(s) ds solves -(A u')' = q with u(0) = 0 and (A u')(L) = 0; the "
      "supply end is x = 0 and the FAR end is x = L, where u(L) = q*int_0^L (L-s)/A ds",
      sp.simplify(-sp.diff(A*sp.diff(u, x), x) - q) == 0
      and u.subs(x, 0) == 0
      and sp.simplify((A*sp.diff(u, x)).subs(x, L)) == 0
      and sp.simplify(u.subs(x, L) - q*L**2/(2*A)) == 0,
      "with A constant the far-end value is qL^2/(2A) — the sharp case")

# ── R7  the badge arithmetic
per = [16, 32, 33, 25, 23, 20, 39, 20, 32, 26]
check("R7  the ten numbered run badges total exactly 266 (the census replication's 6 checks are a "
      "separate script and are counted separately)", sum(per) == 266, f"{'+'.join(map(str,per))} = {sum(per)}")

# ── R8  the audit manifest: are the eleven scripts distinct artifacts?
files = sorted(f.name for f in HERE.iterdir() if f.name.startswith('lesson') and f.suffix == '.py')
digests = {}
for f in files:
    digests[f] = hashlib.sha256((HERE / f).read_bytes()).hexdigest()
check("R8  the lesson verification scripts sit beside this file and are pairwise distinct (no "
      "accidental duplicate despite two equal byte counts)",
      len(files) >= 11 and len(set(digests.values())) == len(files),
      f"{len(files)} scripts found in {HERE}")
print()
for f in files:
    print(f"    {f:24s} {(HERE / f).stat().st_size:6d} B   {digests[f][:16]}…")

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks)-len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"   FAILURES: {[c[0] for c in fails]}"))
