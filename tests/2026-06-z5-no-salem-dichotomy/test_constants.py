"""Named constants of the paper, each re-derived from its defining property
before being compared to the paper's decimal.

Covered claims:
  * mu_S = 1.3247179... is the real root of x^3 - x - 1 (Lemma 2.4), the
    plastic number, and is a Pisot number (all conjugates inside the unit disk).
  * realification value 2^(1/5) = 1.1487  (Remark 6.2 (1)).
  * phi = (1+sqrt5)/2 with phi^2 - phi - 1 = 0  (Prop. 3.1).
  * phi^2 = (3+sqrt5)/2 ~ 2.618, phi^4 = (7+3 sqrt5)/2 ~ 6.854102 = 3 phi + 2.
  * the forced-floor improvement ordering 2^(1/5) < mu_S < 2  (Thm. 3.1 / Table 1).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp

from _z5_engine import phi

mp.mp.dps = 50
x = sp.symbols('x')


def test_mu_S_is_real_root_of_x3_minus_x_minus_1():
    """Lemma 2.4: mu_S = 1.3247179... is 'the real root of x^3 - x - 1'."""
    r = mp.findroot(lambda t: t**3 - t - 1, mp.mpf('1.3'))
    # independently a genuine root
    assert abs(r**3 - r - 1) < mp.mpf(10) ** (-45)
    # matches the paper's quoted digits 1.3247179...
    assert mp.almosteq(r, mp.mpf('1.3247179'), abs_eps=mp.mpf('5e-8'))


def test_mu_S_is_pisot_number():
    """Lemma 2.4 calls mu_S 'the plastic number, the smallest Pisot number.'
    Verify the Pisot property directly: the real root exceeds 1 while both
    conjugates lie strictly inside the unit circle.  ('Smallest' is Siegel's
    theorem, cited, not re-proved here.)"""
    roots = mp.polyroots([1, 0, -1, -1], maxsteps=100, extraprec=80)
    reals = [r for r in roots if abs(mp.im(r)) < mp.mpf('1e-30')]
    conj = [r for r in roots if abs(mp.im(r)) >= mp.mpf('1e-30')]
    assert len(reals) == 1 and mp.re(reals[0]) > 1
    assert len(conj) == 2
    for c in conj:
        assert abs(c) < 1                      # strictly inside the unit disk
    # product of all three moduli = |constant term| = 1  (algebraic-integer sanity)
    prod = mp.mpf(1)
    for r in roots:
        prod *= abs(r)
    assert mp.almosteq(prod, mp.mpf(1), abs_eps=mp.mpf('1e-30'))


def test_realification_value_2_pow_1_over_5():
    """Remark 6.2(1): the realification bound is 2^(1/5) = 1.1487."""
    val = mp.root(2, 5)
    assert mp.almosteq(val, mp.mpf('1.1487'), abs_eps=mp.mpf('5e-5'))
    # its fifth power is exactly 2
    assert mp.almosteq(val**5, mp.mpf(2), abs_eps=mp.mpf('1e-40'))


def test_phi_defining_relation():
    """Prop. 3.1: phi = (1+sqrt5)/2 satisfies phi^2 - phi - 1 = 0."""
    assert sp.simplify(phi**2 - phi - 1) == 0


def test_phi_squared_closed_form_and_decimal():
    """Lemma 2.5 / Thm. 3.1(b): phi^2 = (3+sqrt5)/2 ~ 2.618."""
    assert sp.simplify(phi**2 - (3 + sp.sqrt(5)) / 2) == 0
    assert mp.almosteq(mp.mpf(str(sp.N(phi**2, 40))), mp.mpf('2.618033988749895'),
                       abs_eps=mp.mpf('1e-15'))


def test_phi_fourth_closed_form_and_decimal():
    """Thm. 4.1 / Sec. 7: phi^4 = 6.854102, and phi^4 = 3 phi + 2 = (7+3 sqrt5)/2."""
    assert sp.simplify(phi**4 - (3 * phi + 2)) == 0
    assert sp.simplify(phi**4 - (7 + 3 * sp.sqrt(5)) / 2) == 0
    assert mp.almosteq(mp.mpf(str(sp.N(phi**4, 40))), mp.mpf('6.854102'),
                       abs_eps=mp.mpf('5e-7'))


def test_floor_improvement_ordering():
    """Table 1 / Thm. 3.1: the forced floor rises from the realification value
    2^(1/5)=1.1487 to mu_S=1.3247, still strictly below the realized floor 2."""
    two15 = mp.root(2, 5)
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf('1.3'))
    assert two15 < muS < 2
    assert mp.almosteq(muS - two15, mp.mpf('0.17601'), abs_eps=mp.mpf('1e-3'))
