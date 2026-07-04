"""
Independent verification of Section 5: "The gate ladder and the trifurcation".
Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - eq:gengate  : x^2 + x = C has discriminant D = 1 + 4C.
  - eq:companion: R_C = [[0,C],[1,-1]] has charpoly x^2+x-C, eigenvalues
                  (-1 +- sqrt(1+4C))/2, and gap lambda_+ - lambda_- = sqrt(1+4C)=sqrt D.
  - gate table  : C in {1/4,1/2,1} -> D in {2,3,5}; seeds sqrt2,sqrt3,sqrt5; fields.
  - golden gate : x^2+x-1 roots are {1/phi, -phi}; disc(x^2+x-1)=disc(x^2-x-1)=5.
  - prop:trifurcation : spec(ad_{R_C}) = {-sqrt(1+4C), 0, +sqrt(1+4C)}, 0 mult 2,
                  0-eigenspace = span{I, R_C}, charpoly of ad = t^2 (t^2 - (1+4C)).
  - eq:gatebalance : lambda*C = r(R_C) = sqrt(1+4C); at C=1 both sides = sqrt5.
  - rem:mahler  : Mah(x^2 - D) = D (D>1); r(R_C)=sqrt(1+4C)=sqrt(Mah(x^2-D)).
"""
import sympy as sp

x, t, C = sp.symbols('x t C')


def test_gate_discriminant():
    """eq:gengate: 'x^2 + x = C, discriminant D = 1+4C'.
    Independently compute the discriminant of x^2 + x - C.
    """
    p = sp.Poly(x**2 + x - C, x)
    D = sp.discriminant(p)
    assert sp.simplify(D - (1 + 4 * C)) == 0


def test_companion_charpoly_and_eigen_gap():
    """eq:companion: R_C=[[0,C],[1,-1]] has characteristic polynomial x^2+x-C,
    eigenvalues (-1 +- sqrt(1+4C))/2 and gap sqrt(1+4C).
    """
    R = sp.Matrix([[0, C], [1, -1]])
    cp = R.charpoly(x).as_expr()
    assert sp.simplify(cp - (x**2 + x - C)) == 0
    eigs = R.eigenvals()
    lam_plus = (-1 + sp.sqrt(1 + 4 * C)) / 2
    lam_minus = (-1 - sp.sqrt(1 + 4 * C)) / 2
    # every stated eigenvalue is an actual eigenvalue
    keys = list(eigs.keys())
    assert any(sp.simplify(k - lam_plus) == 0 for k in keys)
    assert any(sp.simplify(k - lam_minus) == 0 for k in keys)
    gap = sp.simplify(lam_plus - lam_minus)
    assert sp.simplify(gap - sp.sqrt(1 + 4 * C)) == 0


def test_gate_table_discriminants_and_seeds():
    """gate table: C in {1/4, 1/2, 1} gives D=1+4C in {2,3,5}, seeds sqrt D,
    minpolys x^2-D, and gap r(R_C)=sqrt D.
    """
    table = {sp.Rational(1, 4): 2, sp.Rational(1, 2): 3, sp.Integer(1): 5}
    for Cval, Dexpected in table.items():
        D = 1 + 4 * Cval
        assert D == Dexpected
        # companion gap equals sqrt D
        R = sp.Matrix([[0, Cval], [1, -1]])
        eigs = list(R.eigenvals().keys())
        assert len(eigs) == 2
        gap = eigs[0] - eigs[1]
        assert sp.simplify(sp.Abs(gap) - sp.sqrt(Dexpected)) == 0
        # seed minpoly x^2 - D has root sqrt(D)
        assert sp.simplify((sp.sqrt(Dexpected))**2 - Dexpected) == 0


def test_golden_gate_roots_are_reciprocal_golden_pair():
    """gate table note: golden gate C=1 is T^2+T=1 whose roots (-1+-sqrt5)/2 are
    the reciprocal-golden pair {1/phi, -phi}.
    """
    phi = (1 + sp.sqrt(5)) / 2
    roots = sp.solve(x**2 + x - 1, x)
    inv_phi = sp.nsimplify(1 / phi)
    neg_phi = -phi
    got = {sp.nsimplify(r) for r in roots}
    assert any(sp.simplify(r - inv_phi) == 0 for r in got)
    assert any(sp.simplify(r - neg_phi) == 0 for r in got)
    # 1/phi = phi - 1 = (-1+sqrt5)/2
    assert sp.simplify(1 / phi - (-1 + sp.sqrt(5)) / 2) == 0


def test_tau_and_phi_both_discriminant_five():
    """gate table note: tau = x^2+x-1 (conjugate gate) and phi's minpoly x^2-x-1
    both have discriminant 5.
    """
    assert sp.discriminant(sp.Poly(x**2 + x - 1, x)) == 5
    assert sp.discriminant(sp.Poly(x**2 - x - 1, x)) == 5


# ---------------------------------------------------------------------------
# prop:trifurcation
# ---------------------------------------------------------------------------

def _ad_matrix(R):
    """4x4 matrix of ad_R(X) = [R, X] on M_2 in basis E11,E12,E21,E22."""
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]),
             sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    cols = []
    for E in basis:
        adE = R * E - E * R
        cols.append([adE[0, 0], adE[0, 1], adE[1, 0], adE[1, 1]])
    return sp.Matrix(cols).T


def test_trifurcation_charpoly_and_spectrum_symbolic():
    """prop:trifurcation: spec(ad_{R_C}) = {-sqrt(1+4C), 0, +sqrt(1+4C)} with 0 of
    multiplicity two; charpoly of the 4x4 ad matrix is t^2 (t^2 - (1+4C)).
    """
    R = sp.Matrix([[0, C], [1, -1]])
    ad = _ad_matrix(R)
    cp = sp.factor(ad.charpoly(t).as_expr())
    assert sp.simplify(cp - t**2 * (t**2 - (1 + 4 * C))) == 0
    eigs = ad.eigenvals()
    # zero eigenvalue with algebraic multiplicity 2
    assert eigs.get(sp.Integer(0)) == 2
    # +- sqrt(1+4C) present
    assert any(sp.simplify(k - sp.sqrt(1 + 4 * C)) == 0 for k in eigs)
    assert any(sp.simplify(k + sp.sqrt(1 + 4 * C)) == 0 for k in eigs)


def test_trifurcation_zero_eigenspace_is_span_I_R():
    """prop:trifurcation: 'the 0-eigenspace = span{I, R_C}'.
    Directly: [R,I]=0 and [R,R]=0, and the kernel of ad_R is exactly 2-dimensional.
    """
    Cval = sp.Integer(1)
    R = sp.Matrix([[0, Cval], [1, -1]])
    I2 = sp.eye(2)
    assert (R * I2 - I2 * R) == sp.zeros(2)
    assert (R * R - R * R) == sp.zeros(2)
    ad = _ad_matrix(R)
    assert ad.rank() == 2          # nullity 2 -> kernel dim 2 = dim span{I,R}
    # I and R (vectorized) lie in the kernel
    vecI = sp.Matrix([1, 0, 0, 1])
    vecR = sp.Matrix([0, Cval, 1, -1])
    assert ad * vecI == sp.zeros(4, 1)
    assert ad * vecR == sp.zeros(4, 1)


def test_trifurcation_golden_gate_spectrum():
    """prop:trifurcation: 'at the golden gate C=1 this is {-sqrt5, 0, +sqrt5}'."""
    R = sp.Matrix([[0, 1], [1, -1]])
    ad = _ad_matrix(R)
    eigs = ad.eigenvals()
    assert eigs.get(sp.Integer(0)) == 2
    assert any(sp.simplify(k - sp.sqrt(5)) == 0 for k in eigs)
    assert any(sp.simplify(k + sp.sqrt(5)) == 0 for k in eigs)


def test_gate_balance_at_golden():
    """eq:gatebalance: 'lambda*C = r(R_C) = sqrt(1+4C)'.  At C=1 with the forced
    lambda=sqrt5: sqrt5 * 1 = sqrt(1+4) = sqrt5.
    """
    Cval = sp.Integer(1)
    lam = sp.sqrt(5)
    r = sp.sqrt(1 + 4 * Cval)
    assert sp.simplify(lam * Cval - r) == 0
    assert sp.simplify(r - sp.sqrt(5)) == 0


def test_mahler_of_discriminant_seed():
    """rem:mahler / eq:rmahler: 'Mah(x^2 - D) = D for D>1', so
    r(R_C) = sqrt(1+4C) = sqrt D = sqrt(Mah(x^2-D)).
    Mahler measure of x^2-D = product of |roots| outside unit circle = |sqrtD|*|−sqrtD| = D.
    """
    for Dval in [2, 3, 5, 7, 24]:
        roots = [sp.sqrt(Dval), -sp.sqrt(Dval)]
        M = 1
        for rt in roots:
            if sp.Abs(rt) > 1:
                M = M * sp.Abs(rt)
        assert sp.simplify(M - Dval) == 0
        # r(R_C) = sqrt(D) = sqrt(Mah)
        assert sp.simplify(sp.sqrt(Dval) - sp.sqrt(M)) == 0
