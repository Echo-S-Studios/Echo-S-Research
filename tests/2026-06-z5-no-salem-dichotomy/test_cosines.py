"""Proposition 3.1 (the irrational, Galois-coupled pentagon cosines) and the
symbolic rows of the Sec. 7 verification table.

Claims:
  * 2 cos 72 deg = phi - 1.
  * 2 cos 144 deg = -phi.
  * {2 cos 72, 2 cos 144} are exactly the two roots of x^2 + x - 1.
  * they are Galois-conjugate over Q under sqrt5 -> -sqrt5.
  * both are irrational.
  * cross-term collapse identity phi^2 - phi = 1.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp

from _z5_engine import phi, sigma, sqrt5

x = sp.symbols('x')
c72 = 2 * sp.cos(sp.rad(72))     # 2 cos 72 degrees, kept symbolic/exact
c144 = 2 * sp.cos(sp.rad(144))   # 2 cos 144 degrees


def test_two_cos_72_equals_phi_minus_1():
    """Prop. 3.1: 2 cos 72 deg = phi - 1."""
    assert sp.simplify(c72 - (phi - 1)) == 0


def test_two_cos_144_equals_minus_phi():
    """Prop. 3.1: 2 cos 144 deg = -phi."""
    assert sp.simplify(c144 - (-phi)) == 0


def test_cosines_are_the_two_roots_of_x2_plus_x_minus_1():
    """Prop. 3.1 / Sec. 7: {2cos72, 2cos144} are the roots of x^2 + x - 1."""
    roots = sp.solve(x**2 + x - 1, x)
    got = {sp.nsimplify(sp.simplify(c72)), sp.nsimplify(sp.simplify(c144))}
    want = {sp.nsimplify(r) for r in roots}
    assert got == want
    # each individually annihilates the quadratic
    for c in (c72, c144):
        assert sp.simplify(c**2 + c - 1) == 0


def test_cosines_are_galois_conjugate():
    """Prop. 3.1: the two cosine values are conjugate under sqrt5 -> -sqrt5.
    Write each in Q(sqrt5) and check sigma swaps them."""
    a = sp.nsimplify(sp.simplify(c72), [sp.sqrt(5)])    # = -1/2 + sqrt5/2
    b = sp.nsimplify(sp.simplify(c144), [sp.sqrt(5)])   # = -1/2 - sqrt5/2
    assert sp.simplify(sigma(a) - b) == 0
    assert sp.simplify(sigma(b) - a) == 0
    # their sum (= -coeff of x) is rational, product (= const) is rational
    assert sp.simplify(a + b + 1) == 0        # sum = -1
    assert sp.simplify(a * b + 1) == 0        # product = -1


def test_cosines_are_irrational():
    """Prop. 3.1: both cosines are irrational (this is why the Z/3Z proof,
    which relies on 2cos120 = -1 being rational, breaks at Z/5Z)."""
    assert not sp.nsimplify(sp.simplify(c72)).is_rational
    assert not sp.nsimplify(sp.simplify(c144)).is_rational
    # contrast: the Z/3Z value 2 cos 120 deg IS rational (= -1)
    assert sp.simplify(2 * sp.cos(sp.rad(120)) + 1) == 0


def test_cross_term_collapse_identity():
    """Sec. 7 / Thm. 4.1 proof: phi^2 - phi = 1 collapses the quartic cross term."""
    assert sp.simplify(phi**2 - phi - 1) == 0
    assert sp.simplify((phi**2 - phi) - 1) == 0
