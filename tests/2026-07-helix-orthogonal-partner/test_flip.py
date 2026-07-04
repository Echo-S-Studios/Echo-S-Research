"""
Independent verification of Section 4 ("The flip: D = 1 + 4C") of
"The Dissolved Helix and Its Orthogonal Partner".

Covers: Definition 4.1 (gate family g_C(x)=x^2+x-C, D=1+4C) and
Theorem 4.2 (the canonical flip): root formula, the trace-form Gram eq. (5)
with det G = D, the signature flip at C=-1/4, and the D<0 rotation channel
{-i sqrt|D|, 0, +i sqrt|D|}.
"""
import sympy as sp

x, C, lam = sp.symbols('x C lambda')


# ---------------------------------------------------------------- Def. 4.1
def test_gate_discriminant():
    """Def. 4.1: g_C(x) = x^2 + x - C has discriminant D = 1 + 4C."""
    g = x**2 + x - C
    a, b, c = sp.Poly(g, x).all_coeffs()
    D = b**2 - 4 * a * c
    assert sp.simplify(D - (1 + 4 * C)) == 0


# ---------------------------------------------------------------- Thm. 4.2 (i)
def test_roots_formula():
    """Thm. 4.2(i): roots of g_C are (-1 +/- sqrt(D))/2 with D = 1+4C."""
    roots = sp.solve(x**2 + x - C, x)
    D = 1 + 4 * C
    expected = {(-1 + sp.sqrt(D)) / 2, (-1 - sp.sqrt(D)) / 2}
    assert {sp.simplify(r) for r in roots} == {sp.simplify(e) for e in expected}


def test_double_root_at_flip():
    """Thm. 4.2(i): at D=0 (C=-1/4) the double root is -1/2."""
    roots = sp.solve((x**2 + x - C).subs(C, sp.Rational(-1, 4)), x)
    assert set(roots) == {sp.Rational(-1, 2)}


# ---------------------------------------------------------------- Thm. 4.2, eq. (5)
def _field_trace(expr_of_theta, theta_sym):
    """Field trace over Q(theta)/Q with theta a root of x^2+x-C: sum over the
    two conjugate roots (-1 +/- sqrt(1+4C))/2."""
    D = 1 + 4 * C
    t1 = (-1 + sp.sqrt(D)) / 2
    t2 = (-1 - sp.sqrt(D)) / 2
    return sp.simplify(expr_of_theta.subs(theta_sym, t1) +
                       expr_of_theta.subs(theta_sym, t2))


def test_gram_matrix_and_determinant_eq5():
    """Eq. (5): trace-form Gram in basis {1, theta} is [[2,-1],[-1,1+2C]]
    and det G = 1 + 4C = D."""
    th = sp.symbols('theta')
    tr1 = _field_trace(sp.Integer(1) + 0 * th, th)
    trt = _field_trace(th, th)
    trt2 = _field_trace(th**2, th)
    assert tr1 == 2
    assert sp.simplify(trt + 1) == 0            # Tr(theta) = -1
    assert sp.simplify(trt2 - (1 + 2 * C)) == 0  # Tr(theta^2) = 1+2C
    G = sp.Matrix([[tr1, trt], [trt, trt2]])
    assert sp.simplify(G.det() - (1 + 4 * C)) == 0


def test_signature_flip():
    """Thm. 4.2: leading minor 2>0 always; the form is positive-definite
    (Riemannian, sig (2,0)) for D>0 and indefinite (Lorentzian, sig (1,1))
    for D<0 - the signature toggles with sign(det G)=sign(D)."""
    def gram(cval):
        return sp.Matrix([[2, -1], [-1, 1 + 2 * cval]])
    # leading minor always positive
    for cval in [sp.Integer(1), sp.Rational(-1, 4), sp.Integer(-1)]:
        assert gram(cval)[0, 0] == 2
    # D>0 side (C=1, D=5): both eigenvalues > 0
    evs_pos = list(gram(1).eigenvals().keys())
    assert all(sp.nsimplify(e) > 0 for e in evs_pos)
    # D<0 side (C=-1, D=-3): eigenvalues of opposite sign (signature (1,1))
    evs_ind = [sp.nsimplify(e) for e in gram(-1).eigenvals().keys()]
    assert any(e > 0 for e in evs_ind) and any(e < 0 for e in evs_ind)
    # boundary C=-1/4: det = 0 (degenerate)
    assert gram(sp.Rational(-1, 4)).det() == 0


def test_rotation_channel_D_negative():
    """Thm. 4.2: on the D<0 side the ad-channel of the companion is
    {-i sqrt|D|, 0, +i sqrt|D|}. Take C=-1 (D=-3): spectrum {0,0,+i sqrt3,-i sqrt3}."""
    cval = sp.Integer(-1)                # D = -3
    D = 1 + 4 * cval
    Mc = sp.Matrix([[0, cval], [1, -1]])  # companion of x^2 + x - C
    assert sp.simplify(Mc.charpoly(lam).as_expr() - (lam**2 + lam - cval)) == 0
    # ad_M on M_2 ~ R^4
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]),
             sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    vec = lambda M: [M[0, 0], M[0, 1], M[1, 0], M[1, 1]]
    adM = sp.Matrix.hstack(*[sp.Matrix(4, 1, vec(Mc * E - E * Mc)) for E in basis])
    eig = {sp.simplify(k): v for k, v in adM.eigenvals().items()}
    root = sp.sqrt(-D)                    # sqrt|D| = sqrt3
    assert eig.get(sp.Integer(0)) == 2
    assert eig.get(sp.I * root) == 1
    assert eig.get(-sp.I * root) == 1


def test_flip_location_and_golden_face():
    """Thm. 4.2: the flip is at C=-1/4 (D=0). The golden field lives on the
    D>0 face: at C=1, D=1+4 = 5 = (sqrt5)^2 (connective sanity: sqrt(D)=sqrt5)."""
    assert (1 + 4 * sp.Rational(-1, 4)) == 0
    D_golden = 1 + 4 * sp.Integer(1)
    assert D_golden == 5
    assert sp.simplify(sp.sqrt(D_golden) - sp.sqrt(5)) == 0
