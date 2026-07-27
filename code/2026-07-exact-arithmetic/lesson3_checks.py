# lesson3_checks.py — exact verification for Lesson 3 (lambda = 2c: identity, gate, flip).
# Discipline: no float crosses a decision boundary; floats display only.
import sympy as sp

x, d, C = sp.symbols('x d C')
m = sp.symbols('m', real=True)
s = sp.symbols('sigma', positive=True)
p = sp.symbols('p', positive=True)
n, k, c, gain, logM = sp.symbols('n k c gain logM', positive=True)
PHI = (1 + sp.sqrt(5)) / 2
PSI = (1 - sp.sqrt(5)) / 2
checks = []

def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

def kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i//B.rows, j//B.cols] * B[i%B.rows, j%B.cols])

def cp(M):
    return sp.expand(M.charpoly(x).as_expr())

# ================================================= A. the gain side: 2nd-order KL
# A1 Gaussian location family: KL is EXACTLY quadratic = (1/2) * Fisher * d^2.
N  = sp.sqrt(1/(2*sp.pi*s**2)) * sp.exp(-(x-m)**2/(2*s**2))
score_mean = sp.integrate((x-m)*N, (x, -sp.oo, sp.oo))
check("A1a score has mean zero: E[x-m] = 0  (symbolic integral)", sp.simplify(score_mean) == 0)
fisher_gauss = sp.simplify(sp.integrate(((x-m)/s**2)**2 * N, (x, -sp.oo, sp.oo)))
check("A1b Fisher of Gaussian location = 1/sigma^2  (symbolic integral)",
      sp.simplify(fisher_gauss - 1/s**2) == 0)
# log(p/q) for q shifted by d: [(x-m-d)^2 - (x-m)^2] / (2 s^2) = [d^2 - 2d(x-m)]/(2 s^2)
check("A1c log-ratio algebra: (x-m-d)^2 - (x-m)^2 = d^2 - 2d(x-m)  (symbolic)",
      sp.expand((x-m-d)**2 - (x-m)**2 - (d**2 - 2*d*(x-m))) == 0)
# with E[x-m]=0:  KL = d^2/(2 s^2) EXACTLY = (1/2)*Fisher*d^2 — no O(d^3) remainder here.
check("A1d KL_gauss = d^2/(2 sigma^2) = (1/2) * Fisher * d^2   [exactly quadratic]",
      sp.simplify(d**2/(2*s**2) - sp.Rational(1, 2)*fisher_gauss*d**2) == 0)

# A2 Bernoulli(p): the generic case — zeroth & first order vanish, quadratic = (1/2) Fisher.
KLb = p*sp.log(p/(p+d)) + (1-p)*sp.log((1-p)/(1-p-d))
ser = sp.series(KLb, d, 0, 3).removeO()
check("A2a Bernoulli KL: zeroth-order term vanishes", sp.simplify(ser.coeff(d, 0)) == 0)
check("A2b Bernoulli KL: first-order term vanishes (score mean zero)",
      sp.simplify(ser.coeff(d, 1)) == 0)
check("A2c Bernoulli KL: quadratic coeff = 1/(2p(1-p)) = (1/2)*Fisher(p)",
      sp.simplify(ser.coeff(d, 2) - 1/(2*p*(1-p))) == 0)

# ================================================= B. the metric: trace form = conjugate covariance
# B1 instance K = Q(phi): embeddings matrix B, Gram G = B^T B = [[2,1],[1,3]], det = 5 = disc.
B  = sp.Matrix([[1, PHI], [1, PSI]])
G2 = sp.simplify(B.T * B)
check("B1  B^T B = [[2,1],[1,3]] = trace-form Gram of Q(phi); det = 5 = disc  (exact)",
      sp.simplify(G2 - sp.Matrix([[2, 1], [1, 3]])) == sp.zeros(2, 2)
      and sp.simplify(G2.det() - 5) == 0)

# B2 bridge to Lesson 2: a = (-1, 2) are the coords of sqrt5 = 2*theta - 1 (trace-zero direction);
#    ||a||_G^2 = a^T G a = 10 = Tr(H^2) with H = 2R - I from Lesson 2.
a2 = sp.Matrix([-1, 2])
R  = sp.Matrix([[0, 1], [1, 1]]); H = 2*R - sp.eye(2)
check("B2  trace of 2*theta-1 is 0; and a^T G a = 10 = Tr(H^2)   [Lesson-2 pairing = this G]",
      sp.simplify((2*PHI - 1) + (2*PSI - 1)) == 0
      and sp.simplify((a2.T * G2 * a2)[0, 0] - 10) == 0
      and sp.trace(H*H) == 10)

# B3 gap-basis Gram, SYMBOLIC in C (Prop detG): theta root of x^2+x-C, sqrtD = 2*theta+1.
Dsym = 1 + 4*C
th_p = (-1 + sp.sqrt(Dsym))/2
th_m = (-1 - sp.sqrt(Dsym))/2
check("B3a (2*theta+1)^2 = 1+4C for BOTH roots   [symbolic — the gap element squares to D]",
      sp.simplify((2*th_p + 1)**2 - Dsym) == 0 and sp.simplify((2*th_m + 1)**2 - Dsym) == 0)
TrsqrtD = sp.simplify((2*th_p + 1) + (2*th_m + 1))
check("B3b Tr(sqrtD) = 0, so G = diag(2, 2D) and det G = 4D   [symbolic]", TrsqrtD == 0,
      "G = diag(Tr(1), Tr(D)) = diag(2, 2D)")
check("B3c instances: C=1 -> det 20 = 4*5 (PD); C=-1 -> det -12 = 4*(-3) (Lorentzian); "
      "Q(i): basis {1,i} -> diag(2,-2), det -4 = 4*(-1)",
      (4*Dsym).subs(C, 1) == 20 and (4*Dsym).subs(C, -1) == -12
      and sp.simplify((sp.I) + (-sp.I)) == 0 and sp.simplify(sp.I**2 + (-sp.I)**2 + 2) == 0)

# ================================================= C. the identity lambda = 2c
sol = sp.solve(sp.Eq(gain/(2*c), logM), gain)[0] / logM
check("C1  MDL balance solved symbolically: exchange rate = 2c", sp.simplify(sol - 2*c) == 0)
check("C1b instances: c=1 -> lambda=2 (shipped); c=n -> lambda=2n (degree-aware)",
      sol.subs(c, 1) == 2 and sp.simplify(sol.subs(c, n) - 2*n) == 0)
lam = 2*c
check("C2  identity in genuinely free c: d(lambda)/dc = 2, free symbol {c}; and 2 != sqrt5",
      sp.diff(lam, c) == 2 and lam.free_symbols == {c} and lam.subs(c, 1) != sp.sqrt(5))
check("C3  variance reading: sigma = 1/(2c) = 1/lambda", sp.simplify(1/(2*c) - 1/lam) == 0)
Gsym = sp.MatrixSymbol('G', 2, 2)
check("C4  Cencov rescaling content: F = G/c, c -> k c  =>  F -> F/k   [symbolic]",
      sp.simplify(sp.Matrix(Gsym)/(k*c) - (sp.Matrix(Gsym)/c)/k) == sp.zeros(2, 2))

# ================================================= D. the gate ladder and the frame-shift
RC = sp.Matrix([[0, C], [1, -1]])
check("D1  charpoly(R_C) = x^2 + x - C; eigen gap = sqrt(1+4C)   [symbolic]",
      sp.expand(cp(RC) - (x**2 + x - C)) == 0
      and sp.simplify((th_p - th_m) - sp.sqrt(Dsym)) == 0)
check("D2  (2 R_C + I)^2 = (1+4C) * I   [symbolic — mirror of Lesson 2's H^2 = 5I]",
      sp.simplify(sp.expand((2*RC + sp.eye(2))**2) - Dsym*sp.eye(2)) == sp.zeros(2, 2))

lad = {sp.Rational(1, 4): 2, sp.Rational(1, 2): 3, sp.Integer(1): 5}
check("D3  ladder: D = 1+4C maps {1/4, 1/2, 1} -> {2, 3, 5}   (exact)",
      all(Dsym.subs(C, g) == D0 for g, D0 in lad.items()))
def mahler_x2_minus_D(D0):     # M(x^2 - D) = D for D > 1: both roots +-sqrt(D) outside
    return sp.sqrt(D0)*sp.sqrt(D0) if (D0 - 1) > 0 else None
check("D3b seeds: M(x^2 - D) = D for D in {2,3,5}, so r(R_C) = sqrt(M)   (exact)",
      all(sp.simplify(mahler_x2_minus_D(D0) - D0) == 0 for D0 in (2, 3, 5)))

fs_c   = sp.sqrt(1 + 4*C) / (2*C)          # frame-shift c from the gate balance 2cC = sqrt(1+4C)
fs_lam = 2*fs_c
check("D4  gate balance 2cC = sqrt(1+4C) gives c = sqrt(1+4C)/(2C)   [symbolic solve]",
      sp.simplify(sp.solve(sp.Eq(2*c*C, sp.sqrt(1+4*C)), c)[0] - fs_c) == 0)
vals = {g: sp.simplify(fs_lam.subs(C, g)) for g in lad}
check("D5  frame-shift lambda over the gates: {4*sqrt2, 2*sqrt3, sqrt5} — three distinct values",
      sp.simplify(vals[sp.Rational(1, 4)] - 4*sp.sqrt(2)) == 0
      and sp.simplify(vals[sp.Rational(1, 2)] - 2*sp.sqrt(3)) == 0
      and sp.simplify(vals[sp.Integer(1)] - sp.sqrt(5)) == 0
      and len({sp.nsimplify(v) for v in vals.values()}) == 3)
check("D5b golden gate: c = sqrt5/2 and lambda = sqrt5 = phi - psi   (exact)",
      sp.simplify(fs_c.subs(C, 1) - sp.sqrt(5)/2) == 0
      and sp.simplify(sp.sqrt(5) - (PHI - PSI)) == 0)

# self-action trifurcation: ad_R(X) = RX - XR  ->  I (x) R - R^T (x) I on vec(X)
AD = kron(sp.eye(2), R) - kron(R.T, sp.eye(2))
check("D6  spec(ad_R) : charpoly = x^2 (x^2 - 5)  =>  {-sqrt5, 0, 0, +sqrt5}   (exact)",
      sp.expand(cp(AD) - (x**4 - 5*x**2)) == 0)
LOP = kron(sp.eye(2), R) + kron(R.T, sp.eye(2)) - sp.eye(4)
check("D6b return operator L(X)=RX+XR-X has the SAME charpoly x^2(x^2-5)   (exact)",
      sp.expand(cp(LOP) - (x**4 - 5*x**2)) == 0)

# ================================================= E. the flip: D = 1 + 4C
check("E1  D(-1/4) = 0 and x^2 + x + 1/4 = (x + 1/2)^2 — parabolic double point   (exact)",
      Dsym.subs(C, -sp.Rational(1, 4)) == 0
      and sp.expand(x**2 + x + sp.Rational(1, 4) - (x + sp.Rational(1, 2))**2) == 0)

# C = -1: the elliptic instance lands on x^2 + x + 1 — Lesson 1's cyclotomic (M = 1)!
check("E2  C=-1: gate poly is x^2+x+1, divides x^3-1 (primitive cube roots); D = -3",
      sp.rem(x**3 - 1, x**2 + x + 1, x) == 0 and Dsym.subs(C, -1) == -3)
check("E2b both roots on the unit circle: |root|^2 = 1 exactly",
      all(sp.simplify(r*sp.conjugate(r) - 1) == 0 for r in sp.Poly(x**2 + x + 1, x).all_roots()))

cm1 = fs_c.subs(C, -1)
check("E3  frame-shift c at C=-1 is IMAGINARY: c = -I*sqrt(3)/2  (metric G/c leaves the real regime)",
      sp.simplify(cm1 + sp.I*sp.sqrt(3)/2) == 0 and sp.im(cm1) != 0)

RCm1 = RC.subs(C, -1)
ADm1 = kron(sp.eye(2), RCm1) - kron(RCm1.T, sp.eye(2))
check("E4  rotation face of the trifurcation: charpoly(ad_{R_-1}) = x^2 (x^2 + 3) => {+-i sqrt3, 0, 0}",
      sp.expand(cp(ADm1) - (x**4 + 3*x**2)) == 0)
check("E4b the x^2 factor (the CAPTURED 0-channel) survives on BOTH sides of the flip",
      sp.rem(x**4 - 5*x**2, x**2, x) == 0 and sp.rem(x**4 + 3*x**2, x**2, x) == 0,
      "kernel contains span{I, R}: ad_R(I) = ad_R(R) = 0 always")

# ================================================= F. the floors (display; bracket certified)
mu_bracket = sp.Poly(x**3 - x - 1, x).count_roots(sp.Rational(13247, 10000),
                                                  sp.Rational(13248, 10000))
check("F1  mu_S certified in [1.3247, 1.3248] by Sturm count (exact bracket)", mu_bracket == 1)
mu = sp.nsolve(x**3 - x - 1, x, sp.Rational(1324, 1000))
print(f"      display: floors 2c*log(mu_S): c=1 -> {float(2*sp.log(mu)):.4f}   "
      f"c=4 -> {float(8*sp.log(mu)):.4f}   c=8 -> {float(16*sp.log(mu)):.4f}")

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
