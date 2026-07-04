"""The Mahler gap, directly (Sec. 6: Lemmas 6.1-6.2; App. A).

Independent brute-force / structural re-derivations that the interval (1, phi)
is empty for degree-2 integer polynomials, and that the spectral operators keep
the image in {1} U [phi, infinity).
"""
import mpmath as mp

import emgap_util as U

mp.mp.dps = 40

PHI = (1 + mp.sqrt(5)) / 2


def test_degree2_gap_is_empty_and_min_is_phi():
    """Lemma 6.1 / App. A: scanning all 625 integer quadratics x^2+bx+c with
    |b|,|c| <= 12 finds NO Mahler measure in (1, phi); the smallest measure
    strictly above 1 is exactly phi (golden seed)."""
    min_above_one = None
    count = 0
    for b in range(-12, 13):
        for c in range(-12, 13):
            count += 1
            m = U.mahler([1, b, c])
            # nothing strictly inside the open band (1, phi)
            assert not (mp.mpf("1") + mp.mpf("1e-9") < m < PHI - mp.mpf("1e-9")), (b, c, m)
            if m > 1 + mp.mpf("1e-9"):
                if min_above_one is None or m < min_above_one:
                    min_above_one = m
    assert count == 625
    assert abs(min_above_one - PHI) < mp.mpf(10) ** (-30)


def test_reciprocal_quadratic_measure_jump():
    """Lemma 6.1: x^2-bx+1 has Mahler 1 for |b|<=2 (roots on the circle) and
    Mahler >= phi^2 = 2.618 for |b|>=3, the value phi^2 hit at |b|=3."""
    for b in (-2, -1, 0, 1, 2):
        assert abs(U.mahler([1, -b, 1]) - 1) < mp.mpf(10) ** (-25)
    assert abs(U.mahler([1, -3, 1]) - PHI**2) < mp.mpf(10) ** (-25)
    for b in range(3, 13):
        assert U.mahler([1, -b, 1]) >= PHI**2 - mp.mpf(10) ** (-20)


def test_oplus_multiplicativity():
    """Lemma 6.2 / App. A: M(p (+) q) = M(p) M(q). Example
    M(phi (+) sqrt2) = phi * 2 = 3.236."""
    # direct sum spectrum = union of roots => product polynomial (x^2-x-1)(x^2-2)
    prod = [1, -1, -3, 2, 2]                    # expand (x^2-x-1)(x^2-2)
    lhs = U.mahler(prod)
    rhs = U.mahler([1, -1, -1]) * U.mahler([1, 0, -2])
    assert abs(lhs - rhs) < mp.mpf(10) ** (-25)
    assert abs(lhs - PHI * 2) < mp.mpf(10) ** (-25)


def test_squaring_squares_the_measure():
    """Lemma 6.2 / App. A: squaring eigenvalues squares the measure. The square
    of the golden value phi has minimal polynomial x^2-3x+1 and measure phi^2,
    equal to M(x^2-x-1)^2."""
    m_phi = U.mahler([1, -1, -1])
    m_phi_sq = U.mahler([1, -3, 1])            # minpoly of phi^2
    assert abs(m_phi_sq - m_phi**2) < mp.mpf(10) ** (-25)
    assert abs(m_phi_sq - PHI**2) < mp.mpf(10) ** (-25)


def test_kron_sample_no_value_in_gap_min_phi2():
    """Lemma 6.2 / App. A: over Kronecker products of catalog seeds, no Mahler
    value lands in (1, phi); the sampled minimum (nontrivial) is phi^2."""
    catalog = [[1, -1, -1], [1, 1, -1], [1, 0, -2], [1, 0, -3], [1, 0, -5]]

    def kron_spec_poly(cA, cB):
        rA = U.mp_roots(cA)
        rB = U.mp_roots(cB)
        prod = mp.mpf(1)
        for a in rA:
            for b in rB:
                v = abs(a * b)
                if v > 1:
                    prod *= v
        return prod                             # Mahler of kron spectrum

    vals = []
    for cA in catalog:
        for cB in catalog:
            m = kron_spec_poly(cA, cB)
            assert not (1 + mp.mpf("1e-9") < m < PHI - mp.mpf("1e-9"))
            if m > 1 + mp.mpf("1e-6"):
                vals.append(m)
    assert abs(min(vals) - PHI**2) < mp.mpf(10) ** (-20)


def test_cyclotomic_quadratics_measure_one():
    """App. A (Kronecker): x^2-x+1, x^2+1, x-1 are cyclotomic, Mahler = 1."""
    for coeffs in ([1, -1, 1], [1, 0, 1], [1, -1]):
        assert abs(U.mahler(coeffs) - 1) < mp.mpf(10) ** (-25)
