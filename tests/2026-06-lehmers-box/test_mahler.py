"""Mahler-measure preliminaries: Definition 2.1 (product form = Jensen-integral
form) and Lemma 2.3 (lower bound, multiplicativity, squaring), plus Kronecker
(Lemma 2.4).  Each identity is re-derived independently.
"""
import random

import mpmath as mp
import sympy as sp

from _helpers import mahler_product

mp.mp.dps = 40
x = sp.symbols('x')


def _mahler_integral(int_coeffs, dps=40):
    """Jensen form: exp( int_0^1 log|p(e^{2 pi i theta})| d theta ).
    Only valid (finite integrand) when p has no root on the unit circle."""
    with mp.workdps(dps):
        cs = [mp.mpf(int(c)) for c in int_coeffs]
        d = len(cs) - 1

        def integrand(theta):
            z = mp.e ** (2j * mp.pi * theta)
            val = sum(cs[i] * z ** (d - i) for i in range(len(cs)))
            return mp.log(abs(val))

        integral = mp.quad(integrand, [0, mp.mpf(1) / 4, mp.mpf(1) / 2,
                                       mp.mpf(3) / 4, 1])
        return mp.e ** integral


def test_definition_two_forms_agree():
    """Def. 2.1: product form and Jensen-integral form of Mah agree (verified
    on polynomials with no on-circle roots, where the integral converges)."""
    # NB: these have NO root on the unit circle (so the Jensen integrand is
    # finite).  Salem polynomials such as beta_4 DO have on-circle roots and are
    # excluded here; their measure is checked by the product form elsewhere.
    cases = [
        [1, -1, -1],          # x^2 - x - 1  (phi)
        [1, 0, -1, -1],       # x^3 - x - 1  (mu_S)
        [1, 0, -2],           # x^2 - 2
        [1, 0, -5],           # x^2 - 5
        [1, -7, 1],           # gap
    ]
    for coeffs in cases:
        prod = mahler_product(coeffs)
        integ = _mahler_integral(coeffs)
        assert abs(prod - integ) < mp.mpf(10) ** -12, coeffs


def test_lemma_lower_bound_monic():
    """Lemma 2.3(i): Mah(p) >= 1 for monic integer p."""
    random.seed(1)
    for _ in range(40):
        d = random.randint(1, 6)
        coeffs = [1] + [random.randint(-4, 4) for _ in range(d)]
        assert mahler_product(coeffs) >= 1 - mp.mpf(10) ** -20


def test_lemma_multiplicativity():
    """Lemma 2.3(ii): Mah(pq) = Mah(p) Mah(q)."""
    pairs = [([1, -1, -1], [1, 0, -2]),
             ([1, -7, 1], [1, -1, -1, -1, 1]),
             ([1, 0, -1, -1], [1, 1, -1])]
    for p, q in pairs:
        pq = sp.Poly(p, x) * sp.Poly(q, x)
        pq_coeffs = [int(c) for c in pq.all_coeffs()]
        lhs = mahler_product(pq_coeffs)
        rhs = mahler_product(p) * mahler_product(q)
        assert abs(lhs - rhs) < mp.mpf(10) ** -12


def test_lemma_squaring_roots():
    """Lemma 2.3(iii): if p^[2] has the squared roots of p then
    Mah(p^[2]) = Mah(p)^2.  Build p^[2] from p(x)p(-x) = +/- p^[2](x^2)."""
    cases = [[1, -1, -1], [1, 0, -2], [1, -1, -1, -1, 1], [1, 0, -1, -1]]
    for coeffs in cases:
        p = sp.Poly(coeffs, x)
        d = p.degree()
        prod = sp.expand((-1) ** d * p.as_expr() * p.as_expr().subs(x, -x))
        # prod is even in x; substitute x^2 -> y to get p^[2](y)
        y = sp.symbols('y')
        prod_poly = sp.Poly(prod, x)
        sq_coeffs = {}
        for monom, coeff in prod_poly.terms():
            deg = monom[0]
            assert deg % 2 == 0
            sq_coeffs[deg // 2] = int(coeff)
        maxd = max(sq_coeffs)
        squared = [sq_coeffs.get(maxd - i, 0) for i in range(maxd + 1)]
        lhs = mahler_product(squared)
        rhs = mahler_product(coeffs) ** 2
        assert abs(lhs - rhs) < mp.mpf(10) ** -12, coeffs


def test_kronecker_cyclotomic_measure_one():
    """Lemma 2.4 (Kronecker): cyclotomic polynomials have Mahler measure 1,
    while a non-cyclotomic integer polynomial (x^2-x-1) has measure > 1."""
    for n in range(1, 16):
        phi_n = sp.Poly(sp.cyclotomic_poly(n, x), x)
        coeffs = [int(c) for c in phi_n.all_coeffs()]
        assert abs(mahler_product(coeffs) - 1) < mp.mpf(10) ** -12, n
    assert mahler_product([1, -1, -1]) > 1 + mp.mpf(10) ** -6
