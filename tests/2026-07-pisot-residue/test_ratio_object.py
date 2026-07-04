"""Section 2: the ratio object Rat_p and the scan as a decision procedure.

Independently re-derives Proposition 2.1 (ratio object), the (x-1)-multiplicity
claim, the two-route companion identity charpoly(C (x) C^{-1}) = Rat_p used to
corroborate the ledger, and the Kronecker-square factorization quoted for
x^4-x+1.
"""
import os
import sys

import sympy as sp
from sympy import symbols, Poly, resultant, prod, factor_list, expand, div, simplify, kronecker_product, Matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import rat_object, cyclotomic_scan, multiplicity_of_factor, companion

x, y = symbols('x y')

# a handful of concrete monic integer polynomials (squarefree, p(0)!=0)
CASES = {
    "x^2-x-1":            [1, -1, -1],
    "x^3-x-1":            [1, 0, -1, -1],
    "x^4-2":              [1, 0, 0, 0, -2],
    "x^4-x+1":            [1, 0, 0, -1, 1],
    "x^4-x^3-1":          [1, -1, 0, 0, -1],
    "Z*=x^4-3x^2+1":      [1, 0, -3, 0, 1],
    "x^5-2x^4-...-2":     [1, -2, -2, -2, -2, -2],
}


def test_prop21_resultant_identity_symbolic():
    """Prop. 2.1: Res_y(p(y),p(xy)) = ((-1)^n p(0))^n * prod_{i,j}(x - a_j/a_i).
    Verified symbolically on a generic monic cubic in terms of its own roots."""
    a1, a2, a3 = symbols('a1 a2 a3')
    n = 3
    p = (y - a1) * (y - a2) * (y - a3)
    p_xy = (x * y - a1) * (x * y - a2) * (x * y - a3)
    R = resultant(expand(p), expand(p_xy), y)
    p0 = (-1) ** n * a1 * a2 * a3                      # p(0) = prod(-a_i)
    const = ((-1) ** n * p0) ** n
    roots = [a1, a2, a3]
    rhs = const * prod((x - aj / ai) for ai in roots for aj in roots)
    assert simplify(R - expand(rhs)) == 0


def test_prop21_degree_is_n_squared():
    """Prop. 2.1: deg Rat_p = n^2."""
    for name, c in CASES.items():
        n = len(c) - 1
        assert rat_object(c).degree() == n * n, name


def test_prop21_integer_coefficients():
    """Prop. 2.1: Rat_p has integer coefficients."""
    for name, c in CASES.items():
        R = rat_object(c)
        assert all(coeff.is_Integer for coeff in R.all_coeffs()), name


def test_prop21_xminus1_multiplicity_equals_n_for_squarefree():
    """Prop. 2.1: mult of (x-1) in Rat_p equals sum m_alpha^2 = n when p is
    squarefree (all m_alpha = 1)."""
    for name, c in CASES.items():
        n = len(c) - 1
        assert multiplicity_of_factor(rat_object(c).as_expr(), x - 1) == n, name


def test_prop21_xminus1_multiplicity_with_repeated_root():
    """Prop. 2.1: for p = (x-1)^2 (x-3) the multiplicity of (x-1) in Rat_p is
    sum m_alpha^2 = 2^2 + 1^2 = 5 (NOT n=3)."""
    p = sp.expand((x - 1) ** 2 * (x - 3))          # x^3 -5x^2 +7x -3
    coeffs = [int(k) for k in Poly(p, x).all_coeffs()]
    assert coeffs == [1, -5, 7, -3]
    assert multiplicity_of_factor(rat_object(coeffs).as_expr(), x - 1) == 5


def test_two_route_identity_charpoly_Ckron_Cinv_equals_Rat():
    """Ledger corroboration (Sec. 2): charpoly(C (x) C^{-1}) = Rat_p, since both
    have root multiset {alpha_i/alpha_j}.  Checked on several explicit p; this is
    the mechanism used for the beta_4 instance."""
    for name, c in CASES.items():
        C = companion(c)
        Cinv = C.inv()
        K = Matrix(kronecker_product(C, Cinv))
        cp = Poly(expand(K.charpoly(x).as_expr()), x).primitive()[1]
        R = rat_object(c)
        assert cp == R or cp == Poly(-R.as_expr(), x), name


def test_kronecker_square_factorization_x4_minus_x_plus_1():
    """Sec. 2 quoted identity:
        charpoly(C_{x^4-x+1} (x) C_{x^4-x+1}) = S_6^2 * (x^4+2x^2-x+1),
    the quartic being the psi^2-image.  We recompute charpoly of the tensor
    square of the companion and check it factors as (irreducible sextic)^2 times
    the stated quartic."""
    c = [1, 0, 0, -1, 1]
    C = companion(c)
    K = Matrix(kronecker_product(C, C))
    cp = expand(K.charpoly(x).as_expr())
    quartic = x ** 4 + 2 * x ** 2 - x + 1
    q, r = div(cp, quartic, x)
    assert simplify(r) == 0                                    # quartic divides
    facs = factor_list(expand(q))[1]
    assert len(facs) == 1                                      # a single factor...
    fac, mult = facs[0]
    assert mult == 2 and Poly(fac, x).degree() == 6            # ...an irreducible sextic, squared
    assert Poly(fac, x).is_irreducible
    # the sextic is self-reciprocal (a Salem sextic S_6)
    coeffs = [int(k) for k in Poly(fac, x).all_coeffs()]
    assert coeffs == coeffs[::-1]
    assert expand(fac) == x ** 6 - x ** 4 - x ** 3 - x ** 2 + 1
