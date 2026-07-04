"""Numeric constants stated in the Emission-Gap paper (abstract, Sec. 1, 7, 11,
Cor. 6.1, Cor. 10.4, Appendix A).

Each constant is rebuilt from its DEFINING polynomial / expression at >=40 dps,
then compared to the paper's printed digits.
"""
import mpmath as mp
import emgap_util as U

mp.mp.dps = 50


def test_plastic_number_muS():
    """Sec. 1 / App. A: mu_S = 1.3247179572... is the plastic number,
    the real root of x^3 - x - 1 (Smyth's bound)."""
    root = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    assert abs(root - mp.mpf("1.3247179572")) < 1e-10
    # it is a root of x^3-x-1 (independent sanity)
    assert abs(root**3 - root - 1) < mp.mpf(10) ** (-40)


def test_log_muS():
    """Cor. 6.1 / App. A: log mu_S = 0.281200... nats."""
    root = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    assert abs(mp.log(root) - mp.mpf("0.281200")) < 5e-7


def test_golden_ratio_and_log():
    """Cor. 6.1 / App. A: log phi = 0.481212..., phi = 1.618034."""
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(phi - mp.mpf("1.618034")) < 5e-7
    assert abs(mp.log(phi) - mp.mpf("0.481212")) < 5e-7


def test_phi_powers():
    """App. A: phi^2 = 2.618 and phi^4 = 6.854... = (7+3 sqrt5)/2."""
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(phi**2 - mp.mpf("2.618")) < 1e-3
    assert abs(phi**2 - (phi + 1)) < mp.mpf(10) ** (-40)          # phi^2 = phi+1
    assert abs(phi**4 - (7 + 3 * mp.sqrt(5)) / 2) < mp.mpf(10) ** (-40)


def test_muS_between_lehmer_and_golden():
    """Sec. 1: Lehmer's number 1.1762808 < mu_S < ... ; the Salem band is (1, mu_S)."""
    muS = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    lehmer = mp.mpf("1.1762808")
    assert lehmer < muS < 2


def test_lehmer_number_value():
    """Sec. 1 / Lemma 5.1: Lehmer's number M(L) = 1.1762808..., the largest
    root of the degree-10 Lehmer polynomial, and it sits below mu_S."""
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]           # x^10+x^9-x^7-...-x^3+x+1
    rts = U.mp_roots(L)
    beta = max(rts, key=lambda r: abs(r))
    assert abs(abs(beta) - mp.mpf("1.1762808")) < 1e-7
    # Mahler measure of Lehmer's polynomial equals Lehmer's number
    assert abs(U.mahler(L) - mp.mpf("1.1762808")) < 1e-6
    muS = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    assert abs(beta) < muS


def test_min_degree4_salem_beta4():
    """Cor. 10.4 / App. A: beta_4 = 1.722084..., real root of
    x^4 - x^3 - x^2 - x + 1, exceeds phi."""
    b4 = mp.findroot(lambda z: z**4 - z**3 - z**2 - z + 1, 1.72)
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(b4 - mp.mpf("1.722084")) < 1e-6
    assert b4 > phi
