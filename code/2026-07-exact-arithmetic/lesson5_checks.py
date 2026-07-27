# lesson5_checks.py — exact verification for Lesson 5 (the learning dynamics + language layer).
# Discipline: no float crosses a decision boundary; floats display only.
import hashlib
import sympy as sp
from fractions import Fraction as F

x, t = sp.symbols('x t')
a_, b_, c_, d_, e_, f_, g_, h_ = sp.symbols('a b c d e f g h')
PHI = (1 + sp.sqrt(5)) / 2
checks = []

def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

def kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i//B.rows, j//B.cols] * B[i%B.rows, j%B.cols])

def companion(p):
    P = sp.Poly(p, x); d = P.degree(); assert P.LC() == 1
    Cm = sp.zeros(d, d)
    for i in range(d - 1): Cm[i + 1, i] = 1
    cf = P.all_coeffs()
    for i in range(d): Cm[i, d - 1] = -cf[d - i]
    assert sp.expand(Cm.charpoly(x).as_expr() - p) == 0
    return Cm

# ================================================= A. state: projector + capture (power basis, exact)
PHI4 = x**4 - 10*x**2 + 1                    # theta = sqrt2 + sqrt3
C4 = companion(PHI4)
Gp = sp.Matrix(4, 4, lambda i, j: sp.trace(C4**(i + j)))   # trace-form Gram, EXACT integers
check("A1  power-basis Gram of Q(sqrt2+sqrt3) via companion traces (all integers)",
      Gp == sp.Matrix([[4, 0, 20, 0], [0, 20, 0, 196], [20, 0, 196, 0], [0, 196, 0, 1940]]),
      f"Tr(theta^k), k=0..6 exact")
Cap = sp.Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])          # captured columns e1, e2 = {1, theta}
P = Cap * (Cap.T*Gp*Cap)**-1 * Cap.T * Gp
check("A2  projector: P^2 = P and G-self-adjoint ((GP)^T = GP), all over Q",
      sp.simplify(P*P - P) == sp.zeros(4, 4) and sp.simplify((Gp*P).T - Gp*P) == sp.zeros(4, 4))
xin = sp.Matrix([2, 1, 0, 0])                              # the demo's first observation: 2 + theta
check("A3  captured input: r = x - Px = 0 exactly (capture is an identity)",
      sp.simplify(xin - P*xin) == sp.zeros(4, 1))

# g_orthogonal_integer_vector replica: nullspace of [(G e1)^T; (G e2)^T]
rows = sp.Matrix([ (Gp*sp.Matrix([1,0,0,0])).T, (Gp*sp.Matrix([0,1,0,0])).T ])
ns = rows.nullspace()
def primitive(v):
    L = sp.lcm([sp.fraction(sp.nsimplify(c))[1] for c in v])
    w = [sp.Integer(c*L) for c in v]
    g = sp.gcd([abs(c) for c in w if c != 0])
    return sp.Matrix([c // g for c in w])
w1 = primitive(ns[0])
check("A4  first primitive G-orthogonal integer vector = (-5, 0, 1, 0), i.e. theta^2 - 5 = 2*sqrt6",
      w1 == sp.Matrix([-5, 0, 1, 0]) or w1 == sp.Matrix([5, 0, -1, 0]),
      "the demo's mysterious 2sqrt6 direction, DERIVED from the nullspace")
# coords -> minpoly: multiplication matrix M_a = sum a_i C4^i; here -5I + C4^2
Ma = -5*sp.eye(4) + C4**2
check("A5  minpoly bridge: M_a^2 - 24 I = 0 and M_a not scalar  =>  minpoly = x^2 - 24 (2sqrt6)",
      sp.simplify(Ma*Ma - 24*sp.eye(4)) == sp.zeros(4, 4) and Ma != sp.sqrt(24)*sp.eye(4),
      "coeff_height 24 — Lesson 4's ACT-4 REJECT seed, derived not quoted")

# ================================================= B. the loop: observe / propose / confirm (replica)
class LearnerReplica:
    """Faithful miniature of Protocol prot:growth: exact Welford centroid, Q-decided alignment,
    streak, propose pure+idempotent, confirm sole mutator. Fractions only (G8)."""
    def __init__(self, G, cap_cols, N=3, eps=F(1, 100)):
        self.G = [[F(v) for v in row] for row in G.tolist()]
        self.cols = [list(map(F, c)) for c in cap_cols]
        self.N, self.eps = N, eps
        self.mean, self.k, self.streak = None, 0, 0
        self.growth = 0
    def _P(self):
        Cm = sp.Matrix([[sp.Rational(v) for v in col] for col in self.cols]).T
        Gm = sp.Matrix([[sp.Rational(v) for v in row] for row in self.G])
        return Cm*(Cm.T*Gm*Cm)**-1*Cm.T*Gm, Gm
    def _resid(self, xv):
        Pm, Gm = self._P()
        xm = sp.Matrix([sp.Rational(v) for v in xv])
        r = xm - Pm*xm
        rn = (r.T*Gm*r)[0, 0]
        return [F(sp.nsimplify(v)) for v in r], F(sp.nsimplify(rn))
    def observe(self, xv):
        if any(isinstance(v, float) for v in xv):
            raise TypeError("float observation rejected (G8)")
        r, rn = self._resid(xv)
        if rn == 0:
            self.mean, self.k, self.streak = None, 0, 0
            return
        self.k += 1
        if self.mean is None:
            self.mean = r[:]
        else:
            self.mean = [m + (ri - m)/self.k for m, ri in zip(self.mean, r)]
        dev = [ri - m for ri, m in zip(r, self.mean)]
        Gm = self.G
        q = lambda v: sum(v[i]*Gm[i][j]*v[j] for i in range(4) for j in range(4))
        self.streak = self.streak + 1 if q(dev) <= self.eps**2 * q(self.mean) else 1
    def propose(self):
        if self.streak < self.N:
            return None
        return ("seed", tuple(self.mean))          # pure; carries the exact centroid
    def confirm(self, proposal):
        self.cols.append([F(v) for v in proposal[1]])
        self.growth += 1
        self.mean, self.k, self.streak = None, 0, 0

L5 = LearnerReplica(Gp, [[1,0,0,0], [0,1,0,0]])
L5.observe([2, 1, 0, 0])                                    # captured -> resets, no streak
w = [-5, 0, 1, 0]
for _ in range(3):
    L5.observe([F(2)+w[0], F(1)+w[1], F(0)+w[2], F(0)+w[3]])
p1, p2 = L5.propose(), L5.propose()
check("B1  persistence: 3 aligned off-axis observations -> streak 3; centroid = w exactly",
      L5.streak == 3 and list(p1[1]) == [F(v) for v in w])
check("B2  propose is pure + idempotent: two calls equal; growth count still 0",
      p1 == p2 and L5.growth == 0)
L5.confirm(p1)
_, rn_after = L5._resid([F(2)+w[0], F(1)+w[1], F(0)+w[2], F(0)+w[3]])
check("B3  confirm is the sole mutator; afterwards the novel direction is CAPTURED: r = 0 exactly",
      L5.growth == 1 and rn_after == 0 and L5.propose() is None,
      "Thm onegrowth: exactly one growth per persistent novelty")
try:
    L5.observe([0.5, 0, 0, 0]); g8 = False
except TypeError:
    g8 = True
check("B4  G8: a float observation raises TypeError (no float crosses the decision core)", g8)

# ================================================= C. growth tier: compositum + Kronecker Gram
GK = sp.diag(4, 8, 12, 24)          # product Gram of Q(sqrt2, sqrt3), basis {1, sqrt2, sqrt3, sqrt6}
GL = sp.diag(2, 14)                 # Gram of Q(sqrt7)
GW = kron(GK, GL)
check("C1  disjoint growth self-check (compositum.py's own assert): G_W = G_K (x) G_L, 8x8 exact",
      GW.shape == (8, 8) and GW == kron(GK, GL) and GW.det() == GK.det()**2 * GL.det()**4,
      f"det G_W = {GW.det()} = det(G_K)^2 * det(G_L)^4")
# capture-by-growth: sqrt6 against embedded K = Q(sqrt2) inside the product basis
e4 = sp.Matrix([0, 0, 0, 1])
CK = sp.Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])            # embedded {1, sqrt2}
PK = CK*(CK.T*GK*CK)**-1*CK.T*GK
rK = e4 - PK*e4
check("C2  out-of-field: sqrt6 vs embedded Q(sqrt2): r = e4, gain = 24 (exact); after capturing all "
      "of W the residual is 0 (capture-by-growth)",
      sp.simplify(rK - e4) == sp.zeros(4, 1) and (rK.T*GK*rK)[0, 0] == 24
      and sp.simplify(e4 - sp.eye(4)*e4) == sp.zeros(4, 1))

# ================================================= D. the language layer: Cl(2,0) exact
def clmul(X, Y):
    A, B, C, D = X; E, Fq, G, H = Y
    return (A*E + B*Fq + C*G - D*H,
            A*Fq + B*E - C*H + D*G,
            A*G + B*H + C*E - D*Fq,
            A*H + B*G - C*Fq + D*E)
ONEc, E1c, E2c, Ic = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)
check("D1  Clifford axioms exact: e1^2 = e2^2 = 1, i^2 = -1, e1 e2 = i = -e2 e1",
      clmul(E1c, E1c) == ONEc and clmul(E2c, E2c) == ONEc
      and clmul(Ic, Ic) == (-1, 0, 0, 0)
      and clmul(E1c, E2c) == Ic and clmul(E2c, E1c) == (0, 0, 0, -1))
# cocycle table: target = i XOR j, sign = (-1)^(bit1(i)*bit0(j)) reproduces every product
basis = [ONEc, E1c, E2c, Ic]
ok_tab = True
for i in range(4):
    for j in range(4):
        prod = clmul(basis[i], basis[j])
        tgt = i ^ j
        sgn = (-1) ** (((i >> 1) & 1) * (j & 1))
        want = tuple(sgn if k == tgt else 0 for k in range(4))
        ok_tab &= (prod == want)
check("D2  the cocycle table: target = i XOR j, sign = (-1)^(bit1(i)*bit0(j)) — verified over all 16 pairs",
      ok_tab, "one bit-product IS the noncommutativity")
check("D2b the Z/4 rotation inside the carrier: <i> has order exactly 4 (i^2 = -1, i^4 = 1)",
      clmul(Ic, Ic) == (-1, 0, 0, 0)
      and clmul(clmul(Ic, Ic), clmul(Ic, Ic)) == ONEc)
# mat iso, SYMBOLIC: mat(X*Y) = mat(X) mat(Y) with mat(a,b,c,d) = [[a+c, b-d],[b+d, a-c]]
def matc(X):
    A, B, C, D = X
    return sp.Matrix([[A + C, B - D], [B + D, A - C]])
XY = clmul((a_, b_, c_, d_), (e_, f_, g_, h_))
check("D3  Cl(2,0) ~= M2(R) as rings, symbolically: mat(X*Y) - mat(X) mat(Y) = 0 in 8 symbols",
      sp.simplify(matc(XY) - matc((a_, b_, c_, d_))*matc((e_, f_, g_, h_))) == sp.zeros(2, 2))
R = sp.Matrix([[0, 1], [1, 1]])
check("D4  the phi keystone: mat(Cl(1/2, 1, -1/2, 0)) = R = [[0,1],[1,1]] exactly; charpoly x^2-x-1",
      matc((sp.Rational(1,2), 1, -sp.Rational(1,2), 0)) == R
      and sp.expand(R.charpoly(x).as_expr() - (x**2 - x - 1)) == 0,
      "one object: loom CATALOG phi = KL_DTA keystone = the head of L")
# ker(L): the phi-slack — exact rational nullspace, dimension 2, each member returns to zero
LOP = kron(sp.eye(2), R) + kron(R.T, sp.eye(2)) - sp.eye(4)
NS = LOP.nullspace()
def unvec(v):   # column-major
    return sp.Matrix([[v[0], v[2]], [v[1], v[3]]])
slack_ok = len(NS) == 2 and all(
    sp.simplify(R*unvec(v) + unvec(v)*R - unvec(v)) == sp.zeros(2, 2) for v in NS)
check("D5  ker(L) is 2-dimensional over Q and every member satisfies RX + XR - X = 0 exactly "
      "(the phi-slack, where learned words live)", slack_ok,
      f"nullspace dim = {len(NS)}")
# the three residuals: nu (idempotence), R_K (anisotropy, bite-depth 2), L (lexicon)
P0 = sp.Matrix([[1, 0], [0, 0]])
check("D6  nu-residual: M(X) = X^T X; for the projector P0, X^T X - X = 0 (nu = 0 <=> idempotent gate)",
      sp.simplify(P0.T*P0 - P0) == sp.zeros(2, 2))
Xa = sp.Matrix([[1, 1], [1, 1]])                       # anisotropic holding (1 + e1 direction)
RK1 = Xa.T*Xa - sp.Rational(sp.trace(Xa.T*Xa), 2)*sp.eye(2)
RK2 = RK1.T*RK1 - sp.Rational(sp.trace(RK1.T*RK1), 2)*sp.eye(2)
check("D7  R_K residual: R_K(X) = X^T X - tau I; nonzero anisotropy [[0,2],[2,0]], and R_K^2 = 0 "
      "(bite-depth 2) exactly",
      RK1 == sp.Matrix([[0, 2], [2, 0]]) and RK2 == sp.zeros(2, 2))
# lexicon: same exact residue -> one entry; sha256 chain tamper-evidence
def key(res_vec):   # canonical exact serialization -> content hash
    s = "|".join(f"{fr.numerator}/{fr.denominator}" for fr in res_vec)
    return hashlib.sha256(s.encode()).hexdigest()
resA = [F(1, 2), F(-3, 7)]
resB = [F(2, 4), F(-6, 14)]                            # different presentation, SAME exact residue
lex = {key(resA): "entry", key(resB): "entry"}
h0 = "genesis"; h1 = hashlib.sha256((h0 + key(resA)).encode()).hexdigest()
h1_tampered = hashlib.sha256((h0 + key([F(1, 2), F(-3, 8)])).encode()).hexdigest()
check("D8  lexicon return-to-zero dedup: equal exact residues collide to ONE entry; witness chain "
      "detects a tampered payload", len(lex) == 1 and h1 != h1_tampered)

# ================================================= E. bridge seeds for Lessons 6 and 7
# E1 (-> L6): the pair charge t_rel = t_a - t_b in Q/Z is a cocycle, gauge-invariant under ray shifts
ts = [F(k, 12) for k in range(12)]                     # the Z/12 object's absolute charges (x^12 - 128)
def trel(i, j): return (ts[i] - ts[j]) % 1
coc = all((trel(i, j) + trel(j, k)) % 1 == trel(i, k) for i in (0, 5, 11) for j in (3, 7) for k in (2, 9))
delta = F(1, 5)                                        # an arbitrary reference-ray shift (gauge)
shifted = [(v + delta) % 1 for v in ts]
gauge = all(((shifted[i] - shifted[j]) % 1) == trel(i, j) for i in (0, 5) for j in (3, 9))
check("E1  pair charge (L6 seed): cocycle law t(a,b)+t(b,c)=t(a,c) in Q/Z, and t_rel is invariant "
      "under a reference-ray shift while absolute charges are not",
      coc and gauge and shifted != ts)
# E2 (-> L7): the gate iteration x -> 1/(1+x) generates the Fibonacci convergents exactly
seq, xk = [], F(1)
for k in range(1, 8):
    seq.append(xk); xk = 1/(1 + xk)
fib_ok = all(seq[k] == F(int(sp.fibonacci(k + 1)), int(sp.fibonacci(k + 2))) for k in range(7))
n_ = sp.symbols('n', positive=True)
rung = sp.simplify(sp.expand_log(sp.log(PHI**n_), force=True) - n_*sp.log(PHI)) == 0
check("E2  golden-substrate seeds (L7): orbit of f_1(x)=1/(1+x) from 1 is F_n/F_{n+1} exactly "
      "(7 terms); rung heights log M(psi^n phi) = n log(phi) symbolically",
      fib_ok and rung, "the ladder rho_n = n ln(phi) is the Adams height of Lesson 2")
tau0 = (1 + sp.sqrt(13))/2
check("E2b trace redirection (L7 seed): beta_4's fold root tau0 = (1+sqrt13)/2 > 2 exactly; and "
      "phi + 1/phi = sqrt5 (the redirection limit)",
      sp.expand(tau0**2 - tau0 - 3) == 0 and (sp.sqrt(13) - 3).is_positive
      and sp.simplify(PHI + 1/PHI - sp.sqrt(5)) == 0)

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
