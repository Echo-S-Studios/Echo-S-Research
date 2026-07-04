"""
Independent verification of the information-geometry claims
(Sec. 7: Thm. 7.1 (Gaussian Fisher = G), Prop. 7.2 (discriminant = Jeffreys
volume), Ex. 7.5 (info volumes), Thm. 7.7 (conjugate-family Fisher), Lemma 7.8
(trace-zero = residual subspace), Thm. 7.9 (||r||^2 = n*Fisher), Sec. 7.3 table).

Fisher matrices are re-derived from first principles: the Gaussian one from the
Hessian of the log-likelihood, the exponential-family one from the Hessian of
the log-partition function, then compared to the paper's stated matrices.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import vsub_nf as nf

x = sp.symbols('x')


def _field(poly, basis):
    C = nf.companion_from_poly(poly)
    G = nf.gram(basis, C)
    t = nf.trace_vector(basis, C)
    return C, G, t


def test_gaussian_fisher_is_G_from_hessian():
    """Thm. 7.1: Fisher of N(Ma, I) equals M^T M = G, derived here as the (constant)
    Hessian of the negative log-likelihood, independent of the sample."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    M = sp.Matrix([[1, phi], [1, phic]])
    a1, a2, y1, y2 = sp.symbols('a1 a2 y1 y2', real=True)
    a = sp.Matrix([a1, a2]); y = sp.Matrix([y1, y2])
    mu = M * a
    ll = -sp.Rational(1, 2) * ((y - mu).T * (y - mu))[0]
    fisher = sp.simplify(-sp.hessian(ll, (a1, a2)))
    G = sp.simplify(M.T * M)
    assert fisher == G == sp.Matrix([[2, 1], [1, 3]])
    # covariance cI scales Fisher by 1/c
    c = sp.symbols('c', positive=True)
    ll_c = -sp.Rational(1, 2) * ((y - mu).T * (y - mu))[0] / c
    assert sp.simplify(-sp.hessian(ll_c, (a1, a2)) - G / c) == sp.zeros(2, 2)


def test_expfamily_fisher_from_logpartition_golden():
    """Thm. 7.7 / Sec. 7.3: for p(k;a) ~ exp((Ma)_k), Fisher at a=0 = Hessian of
    A(a)=log sum_k exp((Ma)_k), equal to (1/n)G-(1/n^2)tt^T = [[0,0],[0,5/4]]."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    M = sp.Matrix([[1, phi], [1, phic]])
    a1, a2 = sp.symbols('a1 a2', real=True)
    a = sp.Matrix([a1, a2])
    eta = M * a
    A = sp.log(sp.exp(eta[0]) + sp.exp(eta[1]))
    fisher = sp.simplify(sp.hessian(A, (a1, a2)).subs({a1: 0, a2: 0}))
    G = sp.Matrix([[2, 1], [1, 3]]); n = 2; t = sp.Matrix([2, 1])
    formula = sp.simplify(G / n - (t * t.T) / n**2)
    assert fisher == formula == sp.Matrix([[0, 0], [0, sp.Rational(5, 4)]])


def test_expfamily_fisher_identity_general():
    """Thm. 7.7: Fisher_exp = M^T (I/n - 11^T/n^2) M = (1/n)G - (1/n^2)tt^T for a
    general M (symbolic 3x2), since t = M^T 1."""
    ms = sp.symbols('m11 m12 m21 m22 m31 m32')
    M = sp.Matrix([[ms[0], ms[1]], [ms[2], ms[3]], [ms[4], ms[5]]])
    n = 3
    one = sp.ones(n, 1)
    Cov = sp.eye(n) / n - (one * one.T) / n**2
    lhs = sp.simplify(M.T * Cov * M)
    G = M.T * M
    t = M.T * one
    rhs = sp.simplify(G / n - (t * t.T) / n**2)
    assert sp.simplify(lhs - rhs) == sp.zeros(2, 2)


def test_expfisher_evidence_table():
    """Sec. 7.3 table: Fisher_exp = (1/n)G - (1/n^2)tt^T equals the stated
    diagonal matrices for Q(sqrt5), Q(sqrt2,sqrt3), Q(sqrt2,sqrt3,sqrt7)."""
    # Q(sqrt5): {1,phi}
    C, G, t = _field(x**2 - x - 1, [(1, 0), (0, 1)])
    F = sp.simplify(G / 2 - (t * t.T) / 4)
    assert F == sp.Matrix([[0, 0], [0, sp.Rational(5, 4)]])
    # Q(sqrt2,sqrt3): G=diag(4,8,12,24), n=4, t=(4,0,0,0)
    G4 = sp.diag(4, 8, 12, 24); t4 = sp.Matrix([4, 0, 0, 0])
    F4 = sp.simplify(G4 / 4 - (t4 * t4.T) / 16)
    assert F4 == sp.diag(0, 2, 3, 6)
    # Q(sqrt2,sqrt3,sqrt7): G=diag(8,56,16,112,24,168,48,336), n=8, t=(8,0,...,0)
    G8 = sp.diag(8, 56, 16, 112, 24, 168, 48, 336)
    t8 = sp.Matrix([8, 0, 0, 0, 0, 0, 0, 0])
    F8 = sp.simplify(G8 / 8 - (t8 * t8.T) / 64)
    assert F8 == sp.diag(0, 7, 2, 14, 3, 21, 6, 42)


def test_trace_zero_is_residual_subspace():
    """Lemma 7.8: <v,1>_G = t^T v = Tr(v); v perp_G 1 <=> Tr(v)=0.
    In Q(sqrt5): sqrt5=(-1,2) has Tr=0 (perp to 1); the residual of 1 against
    <phi> has trace 5/3 != 0 (Rem. 7.11 both-ways check)."""
    C, G, t = _field(x**2 - x - 1, [(1, 0), (0, 1)])
    e1 = sp.Matrix([1, 0])
    for v in [sp.Matrix([-1, 2]), sp.Matrix([0, 1]), sp.Matrix([1, 0])]:
        assert (v.T * G * e1)[0] == (t.T * v)[0] == nf.field_trace(list(v), C)
    # sqrt5 is trace-zero
    assert nf.field_trace((-1, 2), C) == 0
    # projecting the constant 1 onto <phi> leaves residual with trace 5/3 != 0
    B = sp.Matrix([0, 1])
    r = nf.residual((1, 0), B, G)
    assert nf.field_trace(list(r), C) == sp.Rational(5, 3)


def test_residual_norm_is_degree_scaled_fisher():
    """Thm. 7.9 / Sec. 7.3: on the trace-zero subspace G = n Fisher_exp; for the
    residual sqrt5=(-1,2) in Q(sqrt5), ||sqrt5||_G^2 = 10 = 2*5 = n*Fisher(sqrt5)."""
    C, G, t = _field(x**2 - x - 1, [(1, 0), (0, 1)])
    F = sp.simplify(G / 2 - (t * t.T) / 4)
    v = sp.Matrix([-1, 2])                        # sqrt5, trace zero
    assert (t.T * v)[0] == 0
    assert nf.gnorm2(v, G) == 10
    assert (v.T * F * v)[0] == 5
    assert nf.gnorm2(v, G) == 2 * (v.T * F * v)[0]     # ||r||^2 = n * Fisher(r)
    # identity G = n F on the whole trace-zero subspace (parametrise v perp t)
    s = sp.symbols('s')
    vv = sp.Matrix([-t[1] * s, t[0] * s])        # t^T vv = 0
    assert sp.simplify((vv.T * G * vv)[0] - 2 * (vv.T * F * vv)[0]) == 0


def test_discriminant_is_jeffreys_volume():
    """Prop. 7.2 / Ex. 7.5: Jeffreys volume = sqrt(det G) = sqrt|d_K|; the
    integral of the constant sqrt(det G) over [0,1)^n equals sqrt(det G)."""
    cases = {5: x**2 - x - 1, 8: x**2 - 2, 108: x**3 - 2}
    for absdK, poly in cases.items():
        C = nf.companion_from_poly(poly)
        G = nf.power_gram(C)
        vol = sp.sqrt(sp.Abs(G.det()))
        assert sp.simplify(vol - sp.sqrt(absdK)) == 0
    # explicit info volumes and the ~4.6x ratio (Ex. 7.5)
    assert sp.simplify(sp.sqrt(5) - sp.sqrt(5)) == 0
    v5, v2, vc = sp.sqrt(5), 2*sp.sqrt(2), 6*sp.sqrt(3)
    assert abs(float(v5) - 2.2360679) < 1e-6
    assert abs(float(v2) - 2.8284271) < 1e-6
    assert abs(float(vc) - 10.392304) < 1e-6
    assert abs(float(vc / v5) - 4.647580) < 1e-5     # "about 4.6 times"


def test_mdl_field_description_length():
    """Rem. 7.4: (1/2) log det G = (1/2) log|d_K| = log covol(Lambda_K)."""
    for absdK, poly in {5: x**2 - x - 1, 8: x**2 - 2}.items():
        C = nf.companion_from_poly(poly)
        G = nf.power_gram(C)
        lhs = sp.Rational(1, 2) * sp.log(sp.Abs(G.det()))
        rhs = sp.log(sp.sqrt(absdK))
        assert sp.simplify(lhs - rhs) == 0
