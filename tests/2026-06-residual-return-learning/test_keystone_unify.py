r"""
Independent verification of Section 6 ("The unifying principle and the phi
keystone"), Proposition 6.2: the phi-keystone is ONE object, verified equal three
independent ways:
  (1) the exact integer companion of x^2-x-1 (the substrate's first field route),
  (2) the symmetric matrix R of the language return operator (R^2 = R + I),
  (3) the Clifford element cl(1/2, 1, -1/2, 0),
all three the void law x^2 = x + 1, with Mahler measure = phi.
"""
import mpmath as mp
import sympy as sp
from sympy import Matrix, Rational as Q, eye, sqrt, symbols

x = symbols("x")


def companion(coeffs):
    n = len(coeffs) - 1
    C = sp.zeros(n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]
    return C


def mat(a, b, c, d):
    return Matrix([[a + c, b - d], [b + d, a - c]])


def test_route1_companion_of_golden_law():
    """Prop 6.2 route (1): the companion C(x^2-x-1) = [[0,1],[1,1]], the regular
    representation rho(phi) of the substrate's first field Q(phi)."""
    C = companion([1, -1, -1])
    assert C == Matrix([[0, 1], [1, 1]])
    assert sp.expand(C.charpoly(x).as_expr() - (x**2 - x - 1)) == 0


def test_route2_clifford_R_squared():
    """Prop 6.2 route (2): the language keystone R (mat form) satisfies R^2=R+I
    exactly (the atom of the golden law)."""
    R = mat(Q(1, 2), 1, Q(-1, 2), 0)
    assert R * R == R + eye(2)


def test_route3_clifford_element_coordinates():
    """Prop 6.2 route (3): the keystone Clifford element cl(1/2,1,-1/2,0) has
    matrix [[0,1],[1,1]] -- the SAME object as routes (1),(2)."""
    R = mat(Q(1, 2), 1, Q(-1, 2), 0)
    assert R == Matrix([[0, 1], [1, 1]])


def test_three_routes_coincide():
    """Prop 6.2: the three routes produce the identical 2x2 matrix and the
    identical minimal polynomial [1,-1,-1] = x^2-x-1."""
    C = companion([1, -1, -1])
    R = mat(Q(1, 2), 1, Q(-1, 2), 0)
    assert C == R
    assert sp.expand(C.charpoly(x).as_expr() - (x**2 - x - 1)) == 0
    assert sp.expand(R.charpoly(x).as_expr() - (x**2 - x - 1)) == 0


def test_keystone_mahler_measure_is_phi():
    """Prop 6.2: the live Mahler measure of x^2-x-1 equals phi=(1+sqrt5)/2
    (product of max(1,|root|): phi>1 outside, |psi|<1 inside)."""
    mp.mp.dps = 40
    roots = mp.polyroots([1, -1, -1])
    Mah = mp.mpf(1)
    for r in roots:
        Mah *= max(mp.mpf(1), abs(r))
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(Mah - phi) < mp.mpf("1e-30")


def test_keystone_type_gen_golden_law():
    """Prop 6.2 / Ex 5.3: TYPE(R)={gen}: R satisfies the golden law X^2=X+1
    (det R = -1, scalar part 1/2)."""
    R = mat(Q(1, 2), 1, Q(-1, 2), 0)
    assert R * R == R + eye(2)     # gen: satisfies X^2 = X + 1
    assert R.det() == -1
