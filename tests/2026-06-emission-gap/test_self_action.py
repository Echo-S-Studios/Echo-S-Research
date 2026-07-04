"""The self-action ad_M as a derivation and its difference spectrum
(Sec. 10: Lemma 10.1, Lemma 10.2; App. A).
"""
import mpmath as mp
import numpy as np
import sympy as sp

import emgap_util as U

mp.mp.dps = 40

x = U.x


def test_golden_discriminant_is_five():
    """Lemma 10.1: sqrt5 = sqrt(disc(x^2-x-1)); disc = b^2-4ac = 1+4 = 5."""
    b, a, c = -1, 1, -1
    assert b * b - 4 * a * c == 5


def test_adR_spectrum_is_zero_pm_sqrt5():
    """Lemma 10.1 / App. A: ad_R (on 2x2 matrix space, R = C(x^2-x-1)) has
    spectrum {-sqrt5, 0, 0, sqrt5} -- the channel gap sqrt5, totally real,
    off the unit circle."""
    R = U.companion([1, -1, -1])
    adR = U.ad_operator(R)
    ev = sorted(np.linalg.eigvals(adR).real)
    s5 = float(mp.sqrt(5))
    expected = sorted([-s5, 0.0, 0.0, s5])
    for got, exp in zip(ev, expected):
        assert abs(got - exp) < 1e-9
    # off the unit circle: nonzero eigenvalues have |.| = sqrt5 != 1
    for e in ev:
        assert abs(e) < 1e-9 or abs(abs(e) - 1) > 0.5


def test_difference_spectrum_general():
    """Lemma 10.2: for any M with eigenvalues {mu_i}, ad_M has spectrum
    {mu_i - mu_j}. Checked against a generic small integer matrix."""
    M = np.array([[3.0, 1.0, 0.0],
                  [0.0, -2.0, 4.0],
                  [1.0, 0.0, 5.0]])
    mu = np.linalg.eigvals(M)
    diffs = sorted(round(float((a - b).real), 6) for a in mu for b in mu)
    ad_ev = sorted(round(float(e.real), 6) for e in np.linalg.eigvals(U.ad_operator(M)))
    # imaginary parts negligible for this real-spectrum matrix
    assert max(abs(e.imag) for e in np.linalg.eigvals(M)) < 1e-9
    for a, b in zip(ad_ev, diffs):
        assert abs(a - b) < 1e-6


def test_K_difference_2K_minpoly_and_signature():
    """Lemma 10.2 / App. A: the K-formation's real eigenvalue difference 2K has
    minimal polynomial x^4+20x^2-80, signature (2,1), with complex places at
    modulus 2*beta ~ 4.84 (off the circle) -- real, in K, not Salem."""
    # 2K where K is the real root of x^4+5x^2-5; substitute x -> u/2
    u = sp.symbols("u")
    Kpoly = sp.Poly([1, 0, 5, 0, -5], x)
    sub = sp.expand(Kpoly.as_expr().subs(x, u / 2) * 16)
    assert sp.expand(sub - (u**4 + 20 * u**2 - 80)) == 0
    assert U.signature_from_minpoly([1, 0, 20, 0, -80]) == (2, 1)
    rts = U.mp_roots([1, 0, 20, 0, -80])
    cplx = [r for r in rts if abs(r.real) < 1e-25]
    beta = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)
    assert abs(abs(cplx[0]) - 2 * beta) < 1e-6
    assert abs(2 * beta - mp.mpf("4.84")) < 1e-2
    # 2K is real and NOT a Salem number (no on-circle conjugate)
    assert not U.has_salem_factor(sp.Poly([1, 0, 20, 0, -80], x))
