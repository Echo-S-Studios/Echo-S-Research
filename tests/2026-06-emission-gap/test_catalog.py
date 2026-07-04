"""The emission catalog (Def. 2.1) and its eigenvalue arguments / Mahler measures
(Lemma 4.1, Cor. 6.1, App. A).
"""
import cmath
import math

import mpmath as mp
import numpy as np
import sympy as sp

import emgap_util as U

mp.mp.dps = 50

# catalog seeds -> stated minimal polynomial (highest-first coeffs)  (Def. 2.1)
CATALOG = {
    "phi":   [1, -1, -1],      # x^2 - x - 1
    "tau":   [1, 1, -1],       # x^2 + x - 1
    "sqrt2": [1, 0, -2],       # x^2 - 2
    "sqrt3": [1, 0, -3],       # x^2 - 3
    "sqrt5": [1, 0, -5],       # x^2 - 5
    "gap":   [1, -7, 1],       # x^2 - 7x + 1
    "K":     [1, 0, 5, 0, -5], # x^4 + 5x^2 - 5
}


def test_seed_values_are_roots_of_stated_minpolys():
    """Def. 2.1: the named seed value is a root of its stated minimal polynomial."""
    phi = (1 + mp.sqrt(5)) / 2
    checks = {
        "phi": phi,                                   # golden
        "sqrt2": mp.sqrt(2),
        "sqrt3": mp.sqrt(3),
        "sqrt5": mp.sqrt(5),
        "gap": (7 + 3 * mp.sqrt(5)) / 2,              # = phi^4
        "K": mp.sqrt((-5 + 3 * mp.sqrt(5)) / 2),      # real root of x^4+5x^2-5
    }
    for name, val in checks.items():
        coeffs = CATALOG[name]
        p = sum(mp.mpf(str(c)) * val ** (len(coeffs) - 1 - i)
                for i, c in enumerate(coeffs))
        assert abs(p) < mp.mpf(10) ** (-40), name
    # tau: positive root of x^2+x-1 is (sqrt5-1)/2
    tau = (mp.sqrt(5) - 1) / 2
    assert abs(tau**2 + tau - 1) < mp.mpf(10) ** (-40)


def test_catalog_arguments_in_half_pi_Z():
    """Lemma 4.1 / App. A: every catalog eigenvalue has argument in (pi/2)Z,
    i.e. an integer multiple of 90 degrees."""
    for name, coeffs in CATALOG.items():
        for r in U.mp_roots(coeffs):
            if abs(r) == 0:
                continue
            ang = mp.degrees(mp.atan2(r.imag, r.real)) % 90
            ang = min(ang, 90 - ang)
            assert ang < mp.mpf("1e-30"), f"{name}: {r}"


def test_K_seed_spectrum_shape():
    """Lemma 4.1: x^4+5x^2-5 has spectrum {+/-K real, +/- i*beta}, with the
    real pair at 0/180 deg and the imaginary pair at 90/270 deg;
    complex-place modulus beta = 2.4195 (Prop. 9.1 / App. A)."""
    rts = U.mp_roots(CATALOG["K"])
    reals = [r for r in rts if abs(r.imag) < 1e-30]
    imags = [r for r in rts if abs(r.real) < 1e-30]
    assert len(reals) == 2 and len(imags) == 2
    assert abs(abs(imags[0]) - mp.mpf("2.4195")) < 1e-4      # beta
    assert abs(abs(reals[0]) - mp.sqrt((-5 + 3 * mp.sqrt(5)) / 2)) < 1e-30


def test_fifth_root_modulus_on_circle_claim():
    """Prop. 9.1 / App. A: |5^{1/4} * i| = 1.4953 (5^{1/4} = 1.4953...)."""
    assert abs(mp.mpf(5) ** mp.mpf("0.25") - mp.mpf("1.4953")) < 1e-4


def test_gap_seed_is_phi_fourth():
    """App. A: the 'gap' seed x^2-7x+1 realises Mahler measure phi^4 (= 6.854),
    its dominant root equalling phi^4."""
    phi = (1 + mp.sqrt(5)) / 2
    rts = U.mp_roots(CATALOG["gap"])
    dom = max(rts, key=lambda r: abs(r))
    assert abs(dom - phi**4) < mp.mpf(10) ** (-30)
    assert abs(U.mahler(CATALOG["gap"]) - phi**4) < mp.mpf(10) ** (-30)


def test_catalog_mahler_measures():
    """Cor. 6.1 / App. A: catalog Mahler measures are {phi, phi, 2, 3, 5,
    phi^4, beta^2}, minimum phi. beta^2 = (5+3 sqrt5)/2 = 5.854."""
    phi = (1 + mp.sqrt(5)) / 2
    expected = {
        "phi": phi,
        "tau": phi,
        "sqrt2": mp.mpf(2),
        "sqrt3": mp.mpf(3),
        "sqrt5": mp.mpf(5),
        "gap": phi**4,
        "K": (5 + 3 * mp.sqrt(5)) / 2,           # beta^2
    }
    measures = {}
    for name, coeffs in CATALOG.items():
        m = U.mahler(coeffs)
        measures[name] = m
        assert abs(m - expected[name]) < mp.mpf(10) ** (-30), name
    assert abs(min(measures.values()) - phi) < mp.mpf(10) ** (-30)
    # none of the catalog measures lies in the forbidden band (1, mu_S)
    for m in measures.values():
        assert not (1 + 1e-9 < m < U.MU_S - 1e-9)
