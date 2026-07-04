"""
Section 6 -- Self-action: the three powers.

Proposition 6.2 (prop:rates) [OA-IT-01..04; GR-table]:
  From the golden seed (degree 2, M = phi):
    additive        (+)^k        : degree 2k , M = phi^k ,   log M / deg = (1/2) log phi
    Adams        (psi^2)^{o k}    : degree 2  , M = phi^{2^k}
    multiplicative  (x)^k         : degree 2^k, M = phi^{S_k}
  with  S_k = k * C(k-1, floor((k-1)/2))  = 1,2,6,12,30,60,140,...
        S_k ~ 2^k sqrt(k / 2 pi) ,  and S_k = the mean-absolute-deviation sum
        of the symmetric binomial.
"""

from math import comb

import mpmath as mp
import sympy as sp

from _opalg_ops import (
    golden_seed,
    mahler_exact,
    oplus,
    otimes,
    phi,
    psi,
)


# --------------------------------------------------------------------------
# S_k : three independent derivations must agree
# --------------------------------------------------------------------------
def Sk_closed(k):
    """Paper's closed form S_k = k * C(k-1, floor((k-1)/2))."""
    return k * comb(k - 1, (k - 1) // 2)


def Sk_from_tensor_spectrum(k):
    """Independent: exponent of phi in M((x)^k golden).

    golden = {phi, -1/phi}; a tensor word choosing phi j times contributes
    phi^{2j-k} with multiplicity C(k,j); only |.|>1 (i.e. 2j-k>0) enters M.
    """
    return sum(comb(k, j) * (2 * j - k) for j in range(k + 1) if 2 * j - k > 0)


def Sk_mad_sum(k):
    """Independent: total absolute-deviation sum of Binomial(k,1/2)."""
    return sum(sp.Rational(comb(k, j)) * abs(sp.Rational(2 * j - k, 2)) for j in range(k + 1))


def test_Sk_three_ways_and_sequence():
    """All three formulas agree and give 1,2,6,12,30,60,140 for k=1..7."""
    expected = [1, 2, 6, 12, 30, 60, 140]
    for i, k in enumerate(range(1, 8)):
        assert Sk_closed(k) == expected[i]
        assert Sk_from_tensor_spectrum(k) == expected[i]
        assert Sk_mad_sum(k) == expected[i]


# --------------------------------------------------------------------------
# additive power : M((+)^k) = phi^k , degree 2k
# --------------------------------------------------------------------------
def test_additive_power_measure_and_degree():
    """(+)^k golden has degree 2k and M = phi^k; entropy density = (1/2) log phi."""
    A = golden_seed()
    cur = A
    for k in range(1, 6):
        if k > 1:
            cur = oplus(cur, A)
        assert len(cur) == 2 * k
        assert sp.simplify(mahler_exact(cur) - phi**k) == 0
    # constant entropy density (1/2) log phi per dimension
    mp.mp.dps = 40
    val_phi = (1 + mp.sqrt(5)) / 2
    dens = mp.log(val_phi**5) / (2 * 5)
    assert mp.almosteq(dens, mp.log(val_phi) / 2, rel_eps=mp.mpf(10) ** -30)


# --------------------------------------------------------------------------
# Adams power : M((psi^2)^{o k}) = phi^{2^k} , degree 2
# --------------------------------------------------------------------------
def test_adams_power_measure_and_degree():
    """(psi^2)^{o k} golden keeps degree 2 and has M = phi^{2^k}."""
    A = golden_seed()
    cur = A
    for k in range(1, 6):
        cur = psi(2, cur)
        assert len(cur) == 2
        assert sp.simplify(mahler_exact(cur) - phi ** (2**k)) == 0


# --------------------------------------------------------------------------
# multiplicative power : M((x)^k) = phi^{S_k} , degree 2^k
# --------------------------------------------------------------------------
def test_tensor_power_measure_and_degree():
    """(x)^k golden has degree 2^k and M = phi^{S_k}."""
    A = golden_seed()
    cur = A
    for k in range(1, 6):
        if k > 1:
            cur = otimes(cur, A)
        assert len(cur) == 2**k
        assert sp.simplify(mahler_exact(cur) - phi ** Sk_closed(k)) == 0


# --------------------------------------------------------------------------
# asymptotic  S_k ~ 2^k sqrt(k / 2 pi)
# --------------------------------------------------------------------------
def test_Sk_asymptotic():
    """S_k / (2^k sqrt(k/2pi)) -> 1 (ratio within 3% by k=200, monotone-ish)."""
    mp.mp.dps = 50
    ratios = []
    for k in (50, 100, 200):
        approx = mp.mpf(2) ** k * mp.sqrt(mp.mpf(k) / (2 * mp.pi))
        ratios.append(mp.mpf(Sk_closed(k)) / approx)
    assert ratios[-1] > mp.mpf("0.995")
    assert ratios[-1] < mp.mpf("1.005")
    # convergence toward 1 as k grows
    assert abs(ratios[-1] - 1) < abs(ratios[0] - 1)
