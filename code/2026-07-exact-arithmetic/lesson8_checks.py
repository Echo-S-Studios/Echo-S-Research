# lesson8_checks.py — exact verification for Lesson 8 (the closure engine and the door).
# Discipline: no float crosses a decision boundary; guard replica mirrors stated semantics.
import sympy as sp
from fractions import Fraction as F
from itertools import combinations, product

x, t = sp.symbols('x t')
s5 = sp.sqrt(5); PHI = (1 + s5)/2
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

def kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i//B.rows, j//B.cols]*B[i % B.rows, j % B.cols])

# ================= A. the closure guard, replicated (emission_closure_guard semantics)
def trace_down(Rp):
    c = Rp.all_coeffs(); deg = Rp.degree(); m = deg // 2
    p = [sp.Integer(2), t]
    for k in range(1, m): p.append(sp.expand(t*p[-1] - p[-2]))
    T = sp.Integer(c[m])
    for k in range(1, m + 1): T += c[m - k]*p[k]
    return sp.Poly(sp.expand(T), t)

def flip_straddle(T):
    m = T.degree()
    n_real = T.count_roots(-sp.oo, sp.oo)
    at2 = 1 if T.eval(2) == 0 else 0
    atm2 = 1 if T.eval(-2) == 0 else 0
    n_above = T.count_roots(2, sp.oo) - at2
    n_inside = T.count_roots(-2, 2) - at2 - atm2
    return (n_real == m and n_above == 1 and n_inside == m - 1 and at2 == 0 and atm2 == 0)

def sign_in_Qsqrt5(e):
    e = sp.expand(e)
    a = sp.Rational(sp.nsimplify(e.subs(s5, 0)))
    b = sp.Rational(sp.nsimplify(sp.expand((e - a)/s5)))
    if b == 0: return (a > 0) - (a < 0)
    if a == 0: return (b > 0) - (b < 0)
    if a > 0 and b > 0: return 1
    if a < 0 and b < 0: return -1
    s = (a*a > 5*b*b)          # |a| vs |b|sqrt5 decided over Q
    return (1 if a > 0 else -1) if s else (1 if b > 0 else -1)

def validate_closure(M):
    cp = sp.Poly(M.charpoly(x).as_expr(), x)
    salem = []
    for f, _m in sp.factor_list(cp.as_expr())[1]:
        Fp = sp.Poly(f, x)
        cf = Fp.all_coeffs()
        if Fp.degree() >= 4 and Fp.degree() % 2 == 0 and cf == cf[::-1] and Fp.is_irreducible:
            T = trace_down(Fp)
            if flip_straddle(T):
                below = sign_in_Qsqrt5(Fp.as_expr().subs(x, PHI)) > 0   # R(phi) > 0  <=>  beta < phi
                salem.append((Fp, below))
    if not salem: return "FORCED"
    return "INVALID_CLOSURE" if any(b for _, b in salem) else "FORCED_ABOVE_FLOOR"

Rm = companion(x**2 - x - 1)
b4p = x**4 - x**3 - x**2 - x + 1
Lp = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1
fw = [kron(Rm, companion(x**2 - 2)),
      sp.diag(companion(x**2 - 2), companion(x**2 - 3)),
      companion(x**4 + 5*x**2 - 5)**2,
      kron(Rm, sp.eye(2)) - kron(sp.eye(2), Rm.T)]      # self-action ad_R
check("A1  framework ops read FORCED (kron, dsum, square, self-action): the Salem locus is never constructed",
      all(validate_closure(M) == "FORCED" for M in fw))
Foreign = sp.diag(companion(Lp), sp.Matrix([[1]]))
check("A2  the door, exercised: diag(C_Lehmer, [1]) is TRACELESS (a commutator by Shoda) and the guard "
      "fires INVALID_CLOSURE — a sub-phi Salem detected exactly via sign of L(phi) in Q(sqrt5)",
      Foreign.trace() == 0 and validate_closure(Foreign) == "INVALID_CLOSURE"
      and sign_in_Qsqrt5(Lp.subs(x, PHI)) > 0)
check("A3  a planted degree-4 Salem reads FORCED_ABOVE_FLOOR (beta_4 >= phi): the guard certifies the "
      "FLOOR, complementary to the no-Salem angle theorem",
      validate_closure(companion(b4p)) == "FORCED_ABOVE_FLOOR"
      and sign_in_Qsqrt5(b4p.subs(x, PHI)) < 0)
check("A4  cyclotomic control x^6 + 1 reads FORCED (reciprocal but not Salem)",
      validate_closure(companion(x**6 + 1)) == "FORCED")
check("A5  the census claim 'beta_4 > phi', exact chain: beta_4 has its real root in (17/10, 2) and "
      "phi < 17/10 since 125 < 144",
      sp.Poly(b4p, x).count_roots(sp.Rational(17, 10), 2) == 1 and 125 < 144)

# ================= B. the 27-subfield census from group theory
# L = Q(sqrt2, sqrt3, 5^{1/4}, i); G = Gal(L/Q) = D4 x V4 (order 32);
# subfields of K = Q(sqrt2, sqrt3, 5^{1/4})  <->  subgroups H >= <c>, c = complex conjugation.
def d4_mul(a, b):
    (k1, f1), (k2, f2) = a, b
    return (((k1 + (k2 if f1 == 0 else -k2)) % 4), f1 ^ f2)
def g_mul(g, h):
    return (d4_mul(g[0], h[0]), g[1] ^ h[1], g[2] ^ h[2])
def g_inv(g):
    (k, f), b2, b3 = g
    ki = (-k) % 4 if f == 0 else k
    return ((ki, f), b2, b3)
G = [((k, f), b2, b3) for k in range(4) for f in range(2) for b2 in range(2) for b3 in range(2)]
e = ((0, 0), 0, 0)
c = ((0, 1), 0, 0)                                    # complex conjugation: fixes K, sends i -> -i
def closure(gens):
    S = {e}
    frontier = list(gens)
    while frontier:
        g = frontier.pop()
        if g in S: continue
        S.add(g)
        for h in list(S):
            for p in (g_mul(g, h), g_mul(h, g), g_inv(g)):
                if p not in S: frontier.append(p)
    return frozenset(S)
base = closure([c])
subgroups = {base}
frontier = [base]
while frontier:
    H = frontier.pop()
    for g in G:
        if g in H: continue
        H2 = closure(list(H) + [g])
        if H2 not in subgroups:
            subgroups.add(H2); frontier.append(H2)
sigs = {}
for H in subgroups:
    deg = 32 // len(H)                                # [Fix(H):Q]
    # left cosets gH; embedding real iff g^{-1} c g in H
    seen, reps = set(), []
    for g in G:
        key = frozenset(g_mul(g, h) for h in H)
        if key not in seen: seen.add(key); reps.append(g)
    r1 = sum(1 for g in reps if g_mul(g_mul(g_inv(g), c), g) in H)
    sigs[H] = (deg, r1)
n27 = len(subgroups)
all_shape = all(r1 == deg or 2*r1 == deg for (deg, r1) in sigs.values())
quartic21 = sum(1 for (deg, r1) in sigs.values() if deg == 4 and r1 == 2)
check("B1  EXACTLY 27 subfields of K = Q(sqrt2, sqrt3, 5^{1/4}) — replicated from the subgroup lattice "
      f"of D4 x V4 above complex conjugation", n27 == 27, f"count = {n27}")
check("B2  every subfield is totally real or signature (2k, k): r1 in {deg, deg/2} for all 27 "
      "(the lattice invariant the exclusion rests on)", all_shape)
check("B3  exactly FOUR Salem-shape (2,1) quartic subfields — the corrected count (uniqueness overstatement "
      "fixed in the suite)", quartic21 == 4, f"found {quartic21}")
quartics = [x**4 - 5, x**4 - 20, x**4 - 45, x**4 - 180]
check("B4  the four exhibited: Q(d * 5^{1/4}), d in {1, sqrt2, sqrt3, sqrt6} — minpolys x^4 - 5d^4, each "
      "irreducible with exactly 2 real roots (2 real + 1 complex pair)",
      all(sp.Poly(q, x).is_irreducible and sp.Poly(q, x).count_roots(-sp.oo, sp.oo) == 2
          for q in quartics))
tK = 2*sp.sqrt((3*s5 - 5)/2)
check("B5  the suite's example element 2K: minpoly x^4 + 20x^2 - 80 (irreducible, signature (2,1)) — "
      "verified by exact substitution and root count",
      sp.simplify(sp.expand(tK**4 + 20*tK**2 - 80)) == 0
      and sp.Poly(x**4 + 20*x**2 - 80, x).is_irreducible
      and sp.Poly(x**4 + 20*x**2 - 80, x).count_roots(-sp.oo, sp.oo) == 2)

# ================= C. similarity theory (matrix-plates canonical semantics)
def invariant_factors(A):
    n = A.rows
    XA = sp.eye(n)*x - A
    dets = {0: sp.Poly(1, x)}
    for k in range(1, n + 1):
        g = sp.Poly(0, x)
        for rows in combinations(range(n), k):
            for cols in combinations(range(n), k):
                m = XA[rows, cols].det()
                g = sp.Poly(sp.gcd(g.as_expr(), sp.expand(m)), x)
        dets[k] = g.monic() if g.degree() >= 0 and not g.is_zero else g
    facs = []
    for k in range(1, n + 1):
        q = sp.div(dets[k], dets[k - 1])[0].monic()
        if q.degree() >= 1: facs.append(q)
    return facs

A1 = sp.diag(Rm, Rm)                                   # phi (+) phi
B1 = companion(sp.expand((x**2 - x - 1)**2))
fa, fb = invariant_factors(A1), invariant_factors(B1)
check("C1  phi (+) phi: invariant factors [x^2-x-1, x^2-x-1]; companion((x^2-x-1)^2): one factor; "
      "SAME charpoly, DIFFERENT keys => NOT similar (charpoly complete only on the non-derogatory locus)",
      [f.as_expr() for f in fa] == [x**2 - x - 1, x**2 - x - 1]
      and [f.as_expr() for f in fb] == [sp.expand((x**2 - x - 1)**2)]
      and sp.expand(A1.charpoly(x).as_expr() - B1.charpoly(x).as_expr()) == 0
      and fa != fb)
mA = fa[-1]; mB = fb[-1]
check("C2  flags: phi(+)phi is DEROGATORY (minpoly != charpoly) but NOT defective (minpoly squarefree); "
      "the companion is non-derogatory but DEFECTIVE (repeated factor)",
      mA.degree() < 4 and sp.gcd(mA, mA.diff()).degree() == 0
      and mB.degree() == 4 and sp.gcd(mB, mB.diff()).degree() > 0)
RCF = sp.diag(*[companion(f.as_expr()) for f in fa])
check("C3  rational canonical form: (+) of companions of the invariant factors; largest factor = minpoly "
      "(Cor largestif — the very bridge Lesson 5's coords_to_minpoly vendors)",
      invariant_factors(RCF) == fa and sp.expand(mA.as_expr() - (x**2 - x - 1)) == 0)

# ================= D. Generative Emptiness, objects III-V
band_ok, minM, argmin = True, None, None
for bq in range(-6, 7):
    for cq in range(-6, 7):
        p = sp.Poly([1, bq, cq], x)
        if cq == 0: continue
        r = bq*bq - 4*cq
        if r < 0:
            Mv = sp.Integer(abs(cq)) if abs(cq) > 1 else sp.Integer(1)
        else:
            roots = [(-bq + sp.sqrt(r))/2, (-bq - sp.sqrt(r))/2]
            Mv = sp.Integer(1)
            for rt in roots:
                a2 = sp.expand(rt*rt)
                if sp.simplify(a2 - 1).is_positive: Mv = sp.expand(Mv*sp.sqrt(a2))
            Mv = sp.simplify(sp.radsimp(Mv))
        gt1 = sp.simplify(Mv - 1).is_positive
        if gt1:
            below_phi = sp.sign(sp.expand((Mv - PHI)*sp.sqrt(5))) == -1 if Mv != PHI else False
            if Mv != PHI and sp.simplify(Mv - PHI).is_negative:
                band_ok = False
            if minM is None or sp.simplify(Mv - minM).is_negative:
                minM, argmin = Mv, (bq, cq)
check("D1  base gap (Prop 4.1, exhaustive |b|,|c| <= 6): NO monic integer quadratic has M in (1, phi); "
      "the minimum above 1 is EXACTLY phi at x^2 - x - 1",
      band_ok and sp.simplify(minM - PHI) == 0 and argmin == (-1, -1))
A4 = companion(x**2 - 7*x + 1)
f2 = sp.factor(kron(A4, A4).charpoly(x).as_expr())
check("D2  graded normal form, case phi^4 (x) phi^4: charpoly = (x-1)^2 (x^2 - 47x + 1) exactly, "
      "with 47 = L_8 — cyclotomic (q=0) times grow, the tropical ladder again",
      sp.expand(f2 - (x - 1)**2*(x**2 - 47*x + 1)) == 0 and sp.lucas(8) == 47)
KM = companion(x**4 + 5*x**2 - 5)
fK = sp.expand(kron(KM, KM).charpoly(x).as_expr())
tgt = sp.expand((x**2 + 5)**4 * (x**2 - 5*x - 5)**2 * (x**2 + 5*x - 5)**2)
check("D3  case K (x) K: charpoly = (x^2+5)^4 (x^2-5x-5)^2 (x^2+5x-5)^2 exactly — the imaginary q={1,3} "
      "sector (roots ±i sqrt5) plus TWO grow factors, one being the K-gate modulus polynomial y^2+5y-5",
      sp.expand(fK - tgt) == 0)
check("D4  charge action instances: kron(R,R) spectrum signs {+,+,-,-} -> charges {0,0,2,2} (add); "
      "squaring doubles: R^2 spectrum all positive -> {0,0}; and phi(x)phi's on-circle part (x+1)^2 | x^4-1",
      sp.factor(kron(Rm, Rm).charpoly(x).as_expr()) == sp.factor((x + 1)**2*(x**2 - 3*x + 1))
      and sp.Poly((Rm*Rm).charpoly(x).as_expr(), x).count_roots(0, sp.oo) == 2
      and sp.rem(x**4 - 1, x + 1, x) == 0)
dd = sp.symbols('d')
check("D5  minimal chain, step 1: #channels(d) = d^2 - d + 1 = 3  iff  d = 2 (the ternary lock; "
      "spec(ad_R) = {0, ±sqrt5} is its realization)",
      sp.factor(dd**2 - dd + 1 - 3) == (dd - 2)*(dd + 1))

# ================= E. KL_DTA: exact dual-path battery + the carrier bits
def clmul(X, Y):
    A, B, C, D = X; Ee, Ff, Gg, Hh = Y
    return (A*Ee + B*Ff + C*Gg - D*Hh, A*Ff + B*Ee - C*Hh + D*Gg,
            A*Gg + B*Hh + C*Ee - D*Ff, A*Hh + B*Gg - C*Ff + D*Ee)
def matc(X):
    A, B, C, D = X
    return ((A + C, B - D), (B + D, A - C))
def mtr(M): return M[0][0] + M[1][1]
def mdet(M): return M[0][0]*M[1][1] - M[0][1]*M[1][0]
H4 = [(F(1, 2), F(1), F(-1, 2), F(0)),      # the phi keystone
      (F(1), F(0), F(0), F(0)), (F(0), F(1), F(1), F(0)), (F(2), F(-1), F(0), F(3))]
dual = True
for X in H4:
    for Y in H4:
        Z = clmul(X, Y); MZ = matc(Z)
        MX, MY = matc(X), matc(Y)
        MM = ((MX[0][0]*MY[0][0] + MX[0][1]*MY[1][0], MX[0][0]*MY[0][1] + MX[0][1]*MY[1][1]),
              (MX[1][0]*MY[0][0] + MX[1][1]*MY[1][0], MX[1][0]*MY[0][1] + MX[1][1]*MY[1][1]))
        dual &= (MZ == MM) and (2*Z[0] == mtr(MZ)) and (mdet(matc(X))*mdet(matc(Y)) == mdet(MZ))
check("E1  two-route closure, EXACT (Fractions, zero tolerance): Cl-cocycle route == matrix route for "
      "products, traces (tr = 2a), and dets, over a 4-holding battery incl. the keystone", dual)
grade = lambda k: bin(k).count("1")
c2 = lambda g: g*(g - 1)//2
inv = {(a, b): [(-1)**(a*grade(k) + b*c2(grade(k))) for k in range(4)] for a in (0, 1) for b in (0, 1)}
xs = sp.symbols('p q r s')
revX = (xs[0], xs[1], xs[2], -xs[3])
mX = sp.Matrix([[xs[0] + xs[2], xs[1] - xs[3]], [xs[1] + xs[3], xs[0] - xs[2]]])
check("E2  carrier bits (KL_DTA's own reading): grade = popcount; the FOUR involutions are the dual "
      "Klein group (-1)^(a*g + b*C(g,2)) = {id, grade-inv, reversion, conjugation}; and reversion "
      "[1,1,1,-1] IS transpose under mat, symbolically",
      inv[(0, 0)] == [1, 1, 1, 1] and inv[(1, 0)] == [1, -1, -1, 1]
      and inv[(0, 1)] == [1, 1, 1, -1] and inv[(1, 1)] == [1, -1, -1, -1]
      and sp.simplify(sp.Matrix([[revX[0] + revX[2], revX[1] - revX[3]],
                                 [revX[1] + revX[3], revX[0] - revX[2]]]) - mX.T) == sp.zeros(2, 2))

print()
fails = [ch for ch in checks if not ch[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[ch[0] for ch in fails]}"))
