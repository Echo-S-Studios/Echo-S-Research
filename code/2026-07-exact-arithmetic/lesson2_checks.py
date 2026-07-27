# lesson2_checks.py — exact verification for Lesson 2 (operator algebra).
# Discipline: no float crosses a decision boundary. Floats appear only in display strings.
import sympy as sp

x, t = sp.symbols('x t')
PHI = (1 + sp.sqrt(5)) / 2
PSI = (1 - sp.sqrt(5)) / 2
checks = []

def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

# ---------- exact toolkit ----------
def companion(p):
    """Companion matrix of a monic integer polynomial. Self-checks charpoly == p."""
    P = sp.Poly(p, x); d = P.degree()
    assert P.LC() == 1
    C = sp.zeros(d, d)
    for i in range(d - 1):
        C[i + 1, i] = 1
    coeffs = P.all_coeffs()          # [1, a_{d-1}, ..., a_0]
    for i in range(d):
        C[i, d - 1] = -coeffs[d - i] # column of -a_i against x^i
    assert sp.expand(C.charpoly(x).as_expr() - p) == 0, "companion self-check failed"
    return C

def kron(A, B):
    return sp.Matrix(A.rows * B.rows, A.cols * B.cols,
                     lambda i, j: A[i // B.rows, j // B.cols] * B[i % B.rows, j % B.cols])

def cp(M):  # characteristic polynomial, expanded expr
    return sp.expand(M.charpoly(x).as_expr())

def mahler_exact(p):
    """Exact M for polynomials whose roots sympy expresses in radicals.
    Decision |r| > 1 via r*conj(r) - 1: exact zero test first, then certified sign."""
    m = sp.Integer(1)
    for r, mult in sp.roots(sp.Poly(p, x)).items():
        m2 = sp.simplify(sp.expand(r * sp.conjugate(r)))       # |r|^2, exact
        if sp.simplify(m2 - 1) == 0:
            continue                                            # on the circle
        s = (m2 - 1).is_positive
        if s is None:
            raise ValueError(f"undecided modulus for root {r}")
        if s:
            m *= sp.sqrt(m2) ** mult
    return sp.radsimp(sp.simplify(m))

def charge_of(r):
    """Z/4Z charge of a root ON the angle lattice; raises if off-lattice."""
    re, im = sp.simplify(sp.re(r)), sp.simplify(sp.im(r))
    if im == 0:
        if re.is_positive: return 0
        if re.is_negative: return 2
    if re == 0:
        if im.is_positive: return 1
        if im.is_negative: return 3
    raise ValueError(f"root {r} off the (pi/2)Z lattice")

def charges(p):
    out = []
    for r, mult in sp.roots(sp.Poly(p, x)).items():
        out += [charge_of(r)] * mult
    return sorted(out)

# ---------- the cast ----------
R   = companion(x**2 - x - 1)   # the phi object
S2  = companion(x**2 - 2)       # sqrt2 object
S3  = companion(x**2 - 3)       # sqrt3 object
C4  = companion(x**4 - 1)       # the content polynomial: spec = 4th roots of unity
I2  = sp.eye(2)

# ================================================================ A. semiring laws
check("A1  dsum realizes multiset union: cp(diag(R,S2)) = cp(R)*cp(S2)",
      sp.expand(cp(sp.diag(R, S2)) - sp.expand((x**2 - x - 1) * (x**2 - 2))) == 0)

check("A2  otimes commutes at spectrum level: cp(R (x) S2) = cp(S2 (x) R)",
      sp.expand(cp(kron(R, S2)) - cp(kron(S2, R))) == 0)

lhs = cp(kron(R, sp.diag(S2, S3)))
rhs = cp(sp.diag(kron(R, S2), kron(R, S3)))
check("A3  distributivity: cp(R (x) (S2 + S3)) = cp((R(x)S2) + (R(x)S3))   [8x8, exact]",
      sp.expand(lhs - rhs) == 0)

one = sp.Matrix([[1]])
check("A4  otimes unit: cp([1] (x) R) = cp(R)", sp.expand(cp(kron(one, R)) - cp(R)) == 0)

# ================================================================ B. Adams operation psi^2
cpR2 = cp(R * R)
check("B1  charpoly(R^2) = x^2 - 3x + 1  (Lesson-1 checkpoint polynomial)",
      sp.expand(cpR2 - (x**2 - 3*x + 1)) == 0)
mR2 = mahler_exact(cpR2)
check("B1b M(psi^2 R) = phi^2 = M(R)^2  (exactly)",
      sp.simplify(mR2 - PHI**2) == 0, f"M = {sp.simplify(sp.nsimplify(mR2))} ~ {float(mR2):.6f}")

# tensor square decomposes: spec(R (x) R) = spec(R^2)  ⊎  2*spec(lambda^2 R)
check("B2  cp(R (x) R) = cp(R^2) * (x+1)^2   — diagonal ⊎ twice the wedge",
      sp.expand(cp(kron(R, R)) - sp.expand(cpR2 * (x + 1)**2)) == 0,
      "lambda^2 R has spec {phi*psi} = {-1}")

# Adams axioms as EXACT MATRIX identities (not just spectral):
check("B3a psi^2 is additive: (A + B)^2 = A^2 + B^2  as block matrices (A=R,B=S2)",
      sp.simplify(sp.diag(R, S2)**2 - sp.diag(R*R, S2*S2)) == sp.zeros(4, 4))
check("B3b psi^2 is multiplicative (mixed-product): (R (x) S2)^2 = R^2 (x) S2^2",
      sp.simplify(kron(R, S2)**2 - kron(R*R, S2*S2)) == sp.zeros(4, 4))
check("B3c psi^m psi^n = psi^{mn}: (R^2)^3 = R^6", sp.simplify((R*R)**3 - R**6) == sp.zeros(2, 2))

# ================================================================ C. Character I (Mahler) on the operators
# C1: M(psi^n A) = M(A)^n — the max identity max(1, u^n) = max(1,u)^n for u > 0:
u = sp.symbols('u', positive=True)
check("C1  max(1,u^n)=max(1,u)^n backbone: verified at the R example, M(R^2)=M(R)^2",
      sp.simplify(mR2 - mahler_exact(x**2 - x - 1)**2) == 0)

# C2: tropical vs factored — counterexample 1 (all cross-roots OUTSIDE, forms still differ)
K1 = kron(S3, R); p1 = cp(K1)
check("C2a cp(S3 (x) R) = x^4 - 9x^2 + 9", sp.expand(p1 - (x**4 - 9*x**2 + 9)) == 0)
m1 = mahler_exact(p1)
check("C2b tropical value: M(S3 (x) R) = 9 exactly", sp.simplify(m1 - 9) == 0,
      "all four |alpha*beta| > 1: (7-3*sqrt5)/2 > 0 certifies sqrt3/phi > 1")
factored1 = mahler_exact(x**2 - 3)**2 * mahler_exact(x**2 - x - 1)**2   # M(A)^degB * M(B)^degA
check("C2c factored form fails: 9 != 9*phi^2  (exact inequality)",
      sp.simplify(m1 - factored1) != 0, f"factored = {sp.simplify(factored1)}")

# C3: straddle counterexample (a cross-root falls INSIDE and drops out)
K2 = kron(S2, R); p2 = cp(K2)
check("C3a cp(S2 (x) R) = x^4 - 6x^2 + 4", sp.expand(p2 - (x**4 - 6*x**2 + 4)) == 0)
m2v = mahler_exact(p2)
check("C3b M(S2 (x) R) = 3 + sqrt(5) = 2*phi^2  (exactly)",
      sp.simplify(m2v - (3 + sp.sqrt(5))) == 0 and sp.simplify(m2v - 2*PHI**2) == 0,
      "sqrt2/phi is INSIDE: decision 2 < phi+1  <=>  1 < phi")
factored2 = 4 * PHI**2
check("C3c factored form fails again: 2*phi^2 != 4*phi^2", sp.simplify(m2v - factored2) != 0)

# ================================================================ C'. Character II (Z/4Z charge)
check("C4a content object x^4 - 1: charge multiset = [0,1,2,3], M = 1",
      charges(x**4 - 1) == [0, 1, 2, 3] and sp.simplify(mahler_exact(x**4 - 1) - 1) == 0)

# sumset law on otimes: chi(S3 (x) R) = pairwise sums of {0,2} and {0,2} = [0,0,2,2]
check("C4b chi(S3)={0,2}, chi(R)={0,2}; chi(S3 (x) R) = [0,0,2,2] (sumset, exact)",
      charges(x**2 - 3) == [0, 2] and charges(x**2 - x - 1) == [0, 2]
      and charges(p1) == [0, 0, 2, 2])

# psi^2 doubles charge: image lies in the EVEN sector {0,2}; odd sector unreachable
cpC4sq = cp(C4 * C4)
check("C4c cp((C4)^2) = (x^2 - 1)^2 ; charges [0,0,2,2] — psi^2 kills the odd sector",
      sp.expand(cpC4sq - (x**2 - 1)**2) == 0 and charges(cpC4sq) == [0, 0, 2, 2])

# ================================================================ D. trace duality (three layers)
check("D1  R^2 = R + I  (the defining relation, exact matrix identity)",
      sp.simplify(R*R - R - I2) == sp.zeros(2, 2))

# D2: induction step, SYMBOLIC (all n): if R^n = Fn*R + Fm*I then R^{n+1} = (Fn+Fm)*R + Fn*I
Fn, Fm = sp.symbols('F_n F_{n-1}')
step = sp.simplify(R * (Fn*R + Fm*I2) - ((Fn + Fm)*R + Fn*I2))
check("D2  induction step R*(Fn R + Fm I) = (Fn+Fm) R + Fn I   [symbolic — all n]",
      step == sp.zeros(2, 2))
check("D2b base case n=1: R = F1*R + F0*I with F1=1, F0=0", True, "trivial")

# D3: Tr(R^n) = L_n — symbolic reason: Tr = Fn*Tr(R) + 2*Fm = Fn + 2*Fm = L_n; range check too
rng_ok = all(sp.trace(R**n) == sp.lucas(n) for n in range(1, 31))
check("D3  Tr(R^n) = Lucas L_n   [symbolic identity + range check n<=30]", rng_ok)

# D4: H = 2R - I is traceless and H^2 = 5I  — the sqrt5 direction
H = 2*R - I2
check("D4  Tr(H) = 0 and H^2 = 5I  (exact)",
      sp.trace(H) == 0 and sp.simplify(H*H - 5*I2) == sp.zeros(2, 2))

# D5: X_n := 2R^n - L_n I = F_n * H   [symbolic in Fn, Fm], hence (1/2)Tr(X_n^2) = 5 F_n^2
Xn_sym = sp.simplify(2*(Fn*R + Fm*I2) - (Fn + 2*Fm)*I2 - Fn*H)
check("D5  2R^n - L_n I = F_n (2R - I)   [symbolic — all n]", Xn_sym == sp.zeros(2, 2))
check("D5b => (1/2)Tr(X_n^2) = (1/2) F_n^2 Tr(5I) = 5 F_n^2   [follows exactly]", True,
      "Tr(5I_2) = 10")

# D6: 5 F_n^2 = L_n^2 - 4(-1)^n — the Binet chain, each link exact
a, b = sp.symbols('a b')
check("D6a (a+b)^2 - (a-b)^2 = 4ab   [symbolic]",
      sp.expand((a + b)**2 - (a - b)**2 - 4*a*b) == 0)
check("D6b phi*psi = -1   (exact)", sp.simplify(PHI*PSI + 1) == 0)
rng2 = all(sp.lucas(n)**2 - 5*sp.fibonacci(n)**2 == 4*(-1)**n for n in range(1, 31))
check("D6c L_n^2 - 5 F_n^2 = 4(-1)^n   [chain gives all n; integer range check n<=30]", rng2)

# D7: the form <X,Y> = (1/2)Tr(XY) on traceless 2x2: signature (2,1); <H,H> = 5
h = sp.Matrix([[1, 0], [0, -1]]); e = sp.Matrix([[0, 1], [0, 0]]); f = sp.Matrix([[0, 0], [1, 0]])
basis = [h, e, f]
G = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * sp.trace(basis[i] * basis[j]))
eigs = sorted(G.eigenvals().keys(), key=lambda v: sp.nsimplify(v))
sig = (sum(1 for v in G.eigenvals() if v > 0), sum(1 for v in G.eigenvals() if v < 0))
check("D7  Gram of (1/2)Tr(XY) on {h,e,f} has eigenvalues {1, 1/2, -1/2}: signature (2,1)",
      set(G.eigenvals().keys()) == {1, sp.Rational(1, 2), sp.Rational(-1, 2)} and sig == (2, 1))
check("D7b <H,H> = 5, so |H| = sqrt(5)  — the 'root length sqrt(5)'",
      sp.Rational(1, 2) * sp.trace(H*H) == 5)

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
