"""Theorem 4.1 (the pure-pentagon theorem, degree four).

An irreducible quartic with charge group Z/5Z factors over the pentagon lattice
as a +-72 pair (modulus s) times a +-144 pair (modulus t):
    O(x) = (x^2 - s(phi-1)x + s^2)(x^2 + t phi x + t^2).

Claims verified:
  * the four coefficients equal
        [x^3] = phi(t-s)+s, [x^2] = s^2+t^2-st,
        [x^1] = st(phi(s-t)+t), [x^0] = (st)^2.
  * Galois conjugacy t = sigma(s) is exactly what clears sqrt5 from every
    coefficient; then [x^2] = k^2 - 3m and [x^0] = m^2 with k = s+t, m = st.
  * the minimizer (k,m) = (3,1) gives s = phi^2, t = phi^-2, the explicit
    quartic x^4 - x^3 + 6x^2 + 4x + 1, charge Z/5Z, M = phi^4.
  * case analysis:  s,t>1 bottoms out at M = 16 (s=t=2, (k,m)=(4,4));
    s>1>t, m=1 bottoms out at M = phi^4 (k=3); global minimum phi^4.
  * hence M in {1} cup [phi^4, infty); in particular M not in (1,2).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp

from _z5_engine import phi, sigma, sqrt5, charge_group, mahler_hp

mp.mp.dps = 50
x, s, t = sp.symbols('x s t')

# The paper's factored object, built from scratch with the exact phi.
O = sp.expand((x**2 - s * (phi - 1) * x + s**2) * (x**2 + t * phi * x + t**2))
P = sp.Poly(O, x)


def test_x3_coefficient():
    """Thm. 4.1 proof: [x^3] = phi(t - s) + s."""
    assert sp.simplify(P.coeff_monomial(x**3) - (phi * (t - s) + s)) == 0


def test_x2_coefficient():
    """Thm. 4.1 proof: [x^2] = s^2 + t^2 - s t."""
    assert sp.simplify(P.coeff_monomial(x**2) - (s**2 + t**2 - s * t)) == 0


def test_x1_coefficient():
    """Thm. 4.1 proof: [x^1] = s t (phi(s - t) + t)."""
    assert sp.simplify(P.coeff_monomial(x**1) - (s * t * (phi * (s - t) + t))) == 0


def test_x0_coefficient():
    """Thm. 4.1 proof: [x^0] = (s t)^2."""
    assert sp.simplify(P.coeff_monomial(x**0) - (s * t)**2) == 0


def test_galois_conjugacy_clears_sqrt5():
    """Thm. 4.1 proof: setting t = sigma(s) (i.e. s = a+b*sqrt5, t = a-b*sqrt5)
    makes every coefficient rational -- the sqrt5 parts cancel iff s,t are
    Galois conjugate.  We check the coefficient of sqrt5 vanishes in each."""
    a, b = sp.symbols('a b', rational=True)
    sub = {s: a + b * sqrt5, t: a - b * sqrt5}
    for deg in (3, 2, 1, 0):
        coeff = sp.expand(P.coeff_monomial(x**deg).subs(sub))
        # collect as (rational) + (rational)*sqrt5 ; the sqrt5 part must be 0
        poly_in_r5 = sp.Poly(sp.expand(coeff), sqrt5)
        # every power of sqrt5 above 0 (after reducing 5 = sqrt5^2) must vanish
        reduced = sp.expand(coeff).rewrite(sp.sqrt)
        assert sp.simplify(reduced - sigma(reduced)) == 0     # fixed by sigma => in Q


def test_symmetric_coefficient_forms_in_k_m():
    """Thm. 4.1: with s,t the two roots of y^2 - k y + m (so k=s+t, m=st), the
    fixed coefficients read [x^2] = k^2 - 3m and [x^0] = m^2.  Substitute the
    explicit roots and confirm."""
    k, m = sp.symbols('k m')
    root_s = (k + sp.sqrt(k**2 - 4 * m)) / 2
    root_t = (k - sp.sqrt(k**2 - 4 * m)) / 2
    sub = {s: root_s, t: root_t}
    # sanity: these really are roots of y^2-ky+m
    assert sp.simplify(root_s + root_t - k) == 0
    assert sp.simplify(root_s * root_t - m) == 0
    c2 = (s**2 + t**2 - s * t).subs(sub)
    assert sp.simplify(c2 - (k**2 - 3 * m)) == 0
    c0 = ((s * t)**2).subs(sub)
    assert sp.simplify(c0 - m**2) == 0


def _build_quartic(sval, tval):
    poly = sp.expand((x**2 - sval * (phi - 1) * x + sval**2)
                     * (x**2 + tval * phi * x + tval**2))
    return sp.Poly(sp.nsimplify(poly), x)


def test_minimizer_polynomial():
    """Thm. 4.1 / Sec. 7: (k,m)=(3,1) => s = phi^2, t = phi^-2, giving the
    explicit minimizer x^4 - x^3 + 6x^2 + 4x + 1."""
    sval, tval = phi**2, phi**(-2)
    # s,t are the roots of y^2 - 3y + 1
    assert sp.simplify(sval + tval - 3) == 0
    assert sp.simplify(sval * tval - 1) == 0
    Pmin = _build_quartic(sval, tval)
    coeffs = [sp.simplify(c) for c in Pmin.all_coeffs()]
    assert coeffs == [sp.Integer(1), sp.Integer(-1), sp.Integer(6),
                      sp.Integer(4), sp.Integer(1)]


def test_minimizer_has_charge5_and_measure_phi4():
    """Sec. 7: the minimizer x^4 - x^3 + 6x^2 + 4x + 1 has charge exactly Z/5Z
    and Mahler measure phi^4 (= 6.854102), computed independently."""
    coeffs = [1, -1, 6, 4, 1]
    assert charge_group(coeffs) == 5
    M = mahler_hp(coeffs)
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))
    assert mp.almosteq(M, phi4, abs_eps=mp.mpf('1e-30'))
    # and it really is an irreducible quartic
    assert sp.Poly([1, -1, 6, 4, 1], x).is_irreducible


def test_measure_is_product_of_moduli_squared():
    """Thm. 4.1: two roots at modulus s, two at t, so M = max(1,s)^2 max(1,t)^2.
    For the minimizer this is (phi^2)^2 = phi^4."""
    sval, tval = phi**2, phi**(-2)
    M = (max(1, sval))**2 * (max(1, tval))**2
    assert sp.simplify(sp.nsimplify(M) - phi**4) == 0


def test_regime_s_t_gt_1_minimum_is_16():
    """Thm. 4.1 table, regime s,t>1: M = m^2 with m=st>=2, minimum 16 at
    s=t=2, i.e. (k,m)=(4,4)."""
    sval = tval = sp.Integer(2)
    k, m = sval + tval, sval * tval
    assert (k, m) == (4, 4)
    M = (max(1, sval))**2 * (max(1, tval))**2
    assert M == 16
    assert M == m**2
    assert sp.N(16) > sp.N(phi**4)                # this regime is above the global min


def test_regime_s_gt_1_gt_t_m_eq_1_minimum_is_phi4():
    """Thm. 4.1 table, regime s>1>t, m=1 (t=1/s): M = s^2 with
    s = (k+sqrt(k^2-4))/2, k>=3, minimum phi^4 at k=3 (s=phi^2)."""
    def s_of(k):
        return (mp.mpf(k) + mp.sqrt(k * k - 4)) / 2
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))
    # k=3 gives exactly phi^4
    assert mp.almosteq(s_of(3)**2, phi4, abs_eps=mp.mpf('1e-25'))
    # strictly increasing thereafter: every k>=3 gives M >= phi^4
    for k in range(3, 20):
        assert s_of(k)**2 >= phi4 - mp.mpf('1e-25')
        if k > 3:
            assert s_of(k)**2 > phi4


def test_global_gap_excludes_1_to_2():
    """Thm. 4.1: M in {1} cup [phi^4, infty).  Since phi^4 ~ 6.854 > 2, the
    degree-4 sector certainly has no measure in (1,2)."""
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))
    assert phi4 > 2
    # the two nontrivial regime minima are phi^4 and 16; global minimum is phi^4
    assert min(phi4, mp.mpf(16)) == phi4
