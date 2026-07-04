"""Named numeric constants of 'Lehmer's Box'.

Each test independently *derives* the constant from the paper's stated premise
(a polynomial, a closed form) with mpmath/sympy, then compares to the digits the
paper prints. Transcendental checks use mp.dps >= 40.
"""
import mpmath as mp
import sympy as sp

from _helpers import (mahler_product, root_classification, is_palindromic,
                      is_salem)

mp.mp.dps = 50


# --- golden ratio and its inverse (Def. 2.9 / abstract) ---------------------
def test_phi_value():
    """Abstract & Sec 1.3: phi = (1+sqrt5)/2 = 1.6180339887..."""
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(phi - mp.mpf("1.6180339887")) < mp.mpf(10) ** -9


def test_tau_is_phi_inverse():
    """Def. 2.9: tau = (-1+sqrt5)/2 = phi^{-1}."""
    phi = (1 + mp.sqrt(5)) / 2
    tau = (-1 + mp.sqrt(5)) / 2
    assert abs(tau - 1 / phi) < mp.mpf(10) ** -40
    # tau is the smaller root of x^2 + x - 1
    assert abs(tau ** 2 + tau - 1) < mp.mpf(10) ** -40


# --- Smyth / plastic number mu_S (Sec 1.2) ----------------------------------
def test_muS_plastic_number():
    """Sec 1.2 & Lemma 2.7: mu_S = 1.3247179572..., the real root of x^3-x-1,
    the smallest Pisot number."""
    muS = mp.findroot(lambda t: t ** 3 - t - 1, mp.mpf("1.3"))
    assert abs(muS - mp.mpf("1.3247179572")) < mp.mpf(10) ** -9
    # it is the *real* root; the other two are a complex conjugate pair
    roots = mp.polyroots([1, 0, -1, -1])
    reals = [r for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -20]
    assert len(reals) == 1
    assert abs(mp.re(reals[0]) - muS) < mp.mpf(10) ** -20


def test_log_floor_values_nats():
    """Corollary (eq. floorval): log mu_S = 0.281200..., log phi = 0.481212...
    nats (paper prints 6 places, rounded)."""
    muS = mp.findroot(lambda t: t ** 3 - t - 1, mp.mpf("1.3"))
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(mp.log(muS) - mp.mpf("0.281200")) < mp.mpf(10) ** -6
    assert abs(mp.log(phi) - mp.mpf("0.481212")) < mp.mpf(10) ** -6
    # the realised floor strictly exceeds the Smyth floor, both positive
    assert 0 < mp.log(muS) < mp.log(phi)


# --- Lehmer's number (Sec 1.2) ----------------------------------------------
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]  # x^10+x^9-x^7-...-x^3+x+1


def test_lehmer_mahler_measure():
    """Sec 1.2: Mah(L) = 1.1762808182... for
    L(x)=x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1, derived by the product formula."""
    m = mahler_product(LEHMER)
    assert abs(m - mp.mpf("1.1762808182")) < mp.mpf(10) ** -9


def test_lehmer_is_salem_below_muS():
    """Sec 1.2: Lehmer's number is a Salem number lying below mu_S."""
    assert is_palindromic(LEHMER)                     # reciprocal
    out, ins, onc = root_classification(LEHMER)
    assert (out, ins, onc) == (1, 1, 8)               # Salem root pattern, deg 10
    assert is_salem(LEHMER)
    muS = mp.findroot(lambda t: t ** 3 - t - 1, mp.mpf("1.3"))
    assert mahler_product(LEHMER) < muS


# --- minimal degree-four Salem beta_4 (Cor. 5.5) ----------------------------
BETA4 = [1, -1, -1, -1, 1]  # x^4 - x^3 - x^2 - x + 1


def test_beta4_value_and_above_floor():
    """Cor. 5.5: beta_4 = 1.7220838057..., dominant root of
    x^4-x^3-x^2-x+1, and beta_4 > phi."""
    roots = mp.polyroots([mp.mpf(c) for c in BETA4], extraprec=100)
    beta4 = max((mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -20))
    assert abs(beta4 - mp.mpf("1.7220838057")) < mp.mpf(10) ** -9
    phi = (1 + mp.sqrt(5)) / 2
    assert beta4 > phi


def test_beta4_is_salem():
    """Cor. 5.5: beta_4 is a (degree-four) Salem number.  Its Mahler measure
    equals beta_4 itself."""
    assert is_salem(BETA4)
    out, ins, onc = root_classification(BETA4)
    assert (out, ins, onc) == (1, 1, 2)               # deg 4: one pair on circle
    m = mahler_product(BETA4)
    roots = mp.polyroots([mp.mpf(c) for c in BETA4], extraprec=100)
    beta4 = max((mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -20))
    assert abs(m - beta4) < mp.mpf(10) ** -30
