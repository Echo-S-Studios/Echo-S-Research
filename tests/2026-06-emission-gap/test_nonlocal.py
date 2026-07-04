"""The non-local identity: one gap, four domains (Sec. 7, Remark 7.2; App. A).

The four "gaps" share emptiness, not a number. Here we check the numeric
content: the entropy/height relation h = log M (Lind-Schmidt-Ward) and the
three DISTINCT right endpoints sqrt5, phi, log phi.
"""
import mpmath as mp

import emgap_util as U

mp.mp.dps = 40

PHI = (1 + mp.sqrt(5)) / 2


def test_entropy_equals_log_mahler_maps_gap():
    """Table in Sec. 7 (LSW): entropy h = log M, so the height gap (1, phi)
    is the image of the entropy gap (0, log phi) under exp; exp(log phi) = phi."""
    h_end = mp.log(PHI)                       # entropy end
    assert abs(mp.exp(h_end) - PHI) < mp.mpf(10) ** (-30)
    # a concrete companion measure and its entropy: golden seed
    M_phi = U.mahler([1, -1, -1])
    assert abs(mp.log(M_phi) - mp.log(PHI)) < mp.mpf(10) ** (-30)


def test_three_gap_endpoints_are_distinct_numbers():
    """Remark 7.2 / App. A: the channel end sqrt5 = 2.23607, the height end
    phi = 1.61803, and the entropy end log phi = 0.48121 are three different
    quantities -- 'the correspondence is of emptiness, not a shared endpoint'."""
    s5 = mp.sqrt(5)
    logphi = mp.log(PHI)
    assert abs(s5 - mp.mpf("2.23607")) < 1e-5
    assert abs(PHI - mp.mpf("1.61803")) < 1e-5
    assert abs(logphi - mp.mpf("0.48121")) < 1e-5
    # pairwise distinct, well separated
    assert abs(s5 - PHI) > mp.mpf("0.5")
    assert abs(PHI - logphi) > mp.mpf("1.0")
    assert abs(s5 - logphi) > mp.mpf("1.5")


def test_channel_threshold_is_sqrt_discriminant():
    """Sec. 7 / Lemma 10.1: the channel threshold sqrt(D) = sqrt5 is the square
    root of the golden discriminant disc(x^2-x-1) = 5, an ad_R eigenvalue."""
    assert (-1) ** 2 - 4 * 1 * (-1) == 5
    assert abs(mp.sqrt(5) ** 2 - 5) < mp.mpf(10) ** (-30)
