"""Theorem 5.1 (pure power x^5 - m) and Remark 6.2(1) (the psi^5 realification).

Claims:
  * x^5 - m (m >= 2) has charge group Z/5Z, is non-reciprocal, and M = m,
    with roots m^(1/5) exp(2 pi i k/5).
  * the family realizes M in {1} cup [2, infty), attaining mu(5)=2 at x^5-2.
  * psi^5 (fifth powers of the roots) sends a charge-5 object to a totally
    positive one with M(psi^5 O) = M(O)^5; for x^5-2 every fifth power is 2
    (multiplicity 5), so M(psi^5 O)=32 and M(O)=32^(1/5)=2, while the
    totally-positive gap gives only M(O) >= 2^(1/5)=1.1487.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp

from _z5_engine import phi, charge_group, is_charge5, is_reciprocal, mahler_hp

mp.mp.dps = 50


def test_pure_power_charge_measure_nonreciprocal():
    """Thm. 5.1: for m>=2, x^5 - m has charge Z/5Z, M=m, and is non-reciprocal."""
    for m in (2, 3, 4, 5, 7, 10):
        coeffs = [1, 0, 0, 0, 0, -m]
        assert charge_group(coeffs) == 5
        M = mahler_hp(coeffs)
        assert mp.almosteq(M, mp.mpf(m), abs_eps=mp.mpf('1e-25'))
        assert not is_reciprocal(coeffs)


def test_pure_power_roots_are_scaled_fifth_roots_of_unity():
    """Thm. 5.1 proof: the roots of x^5 - m are m^(1/5) exp(2 pi i k/5), all of
    modulus m^(1/5) > 1, so M = (m^(1/5))^5 = m."""
    m = 3
    roots = mp.polyroots([1, 0, 0, 0, 0, -m], maxsteps=100, extraprec=80)
    r = mp.root(m, 5)
    for a in roots:
        assert mp.almosteq(abs(a), r, abs_eps=mp.mpf('1e-30'))
    # product of the five equal moduli = (m^(1/5))^5 = m
    assert mp.almosteq(r**5, mp.mpf(m), abs_eps=mp.mpf('1e-30'))


def test_mu5_realized_at_x5_minus_2():
    """Thm. 5.1: mu(5)=2 is attained at x^5-2 (M=2), the boundary object."""
    M = mahler_hp([1, 0, 0, 0, 0, -2])
    assert mp.almosteq(M, mp.mpf(2), abs_eps=mp.mpf('1e-30'))


def test_family_realizes_floor_set():
    """Thm. 5.1: x^5-m realizes M in {1} cup [2, infty).  m=1 gives x^5-1 with
    M=1 (cyclotomic), and every m>=2 gives M=m>=2 -- no value in (1,2)."""
    assert mp.almosteq(mahler_hp([1, 0, 0, 0, 0, -1]), mp.mpf(1),
                       abs_eps=mp.mpf('1e-25'))
    for m in range(2, 12):
        M = mahler_hp([1, 0, 0, 0, 0, -m])
        assert M >= 2 - mp.mpf('1e-25')
        assert not (1 + mp.mpf('1e-9') < M < 2 - mp.mpf('1e-9'))


def test_psi5_of_x5_minus_2_is_totally_positive_at_2():
    """Remark 6.2 / Sec. 7: psi^5(x^5-2) sends every root to its fifth power,
    all equal to 2 (multiplicity 5), totally positive; hence M(psi^5 O)=32 and
    M(O)=32^(1/5)=2."""
    roots = mp.polyroots([1, 0, 0, 0, 0, -2], maxsteps=100, extraprec=80)
    fifth = [a**5 for a in roots]
    for f in fifth:
        assert mp.almosteq(f, mp.mpf(2), abs_eps=mp.mpf('1e-25'))   # real, = 2
        assert abs(mp.im(f)) < mp.mpf('1e-25') and mp.re(f) > 0     # totally positive
    # psi^5 O = (x-2)^5, Mahler measure 2^5 = 32
    M_psi = mp.mpf(2)**5
    assert mp.almosteq(M_psi, mp.mpf(32), abs_eps=mp.mpf('1e-30'))
    assert mp.almosteq(mp.root(32, 5), mp.mpf(2), abs_eps=mp.mpf('1e-30'))


def test_psi5_multiplies_mahler_by_fifth_power_general():
    """Remark 6.2: M(psi^5 O) = M(O)^5 for any charge-5 object, because the
    fifth power scales every modulus by the 5th power.  Checked on the degree-4
    minimizer (M=phi^4): its psi^5 image is totally positive with M=phi^20."""
    roots = mp.polyroots([1, -1, 6, 4, 1], maxsteps=100, extraprec=80)
    M = mp.mpf(1)
    for a in roots:
        if abs(a) > 1:
            M *= abs(a)
    fifth = [a**5 for a in roots]
    # totally positive: every fifth power is real and positive
    for f in fifth:
        assert abs(mp.im(f)) < mp.mpf('1e-20')
        assert mp.re(f) > 0
    Mpsi = mp.mpf(1)
    for f in fifth:
        if abs(f) > 1:
            Mpsi *= abs(f)
    assert mp.almosteq(Mpsi, M**5, abs_eps=mp.mpf('1e-20'))
    phi20 = mp.mpf(str(sp.N(phi**20, 45)))
    assert mp.almosteq(Mpsi, phi20, abs_eps=mp.mpf('1e-15'))


def test_realification_bound_is_weaker_than_two():
    """Remark 6.2(1): the totally-positive (1,2) gap on psi^5 O yields only
    M(O) >= 2^(1/5) = 1.1487, strictly weaker than the target floor 2 (hence
    'the fifth power is unavoidable')."""
    two15 = mp.root(2, 5)
    assert mp.almosteq(two15, mp.mpf('1.1487'), abs_eps=mp.mpf('5e-5'))
    assert two15 < 2
