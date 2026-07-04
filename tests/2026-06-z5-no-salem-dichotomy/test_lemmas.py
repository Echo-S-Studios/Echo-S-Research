"""Lemma 2.2 (odd modulus forbids negative reals; on-unit roots are fifth roots
of unity) and Lemma 2.5 (real reciprocal units).

Claims:
  * the lattice (2 pi/5) Z = {0, 72, 144, 216, 288} deg omits 180 deg, so a
    charge-5 object has no negative real root; a negative real root forces an
    even charge (contradiction).
  * on-unit charge-5 roots are exactly the fifth roots of unity.
  * a reciprocal real pair {r, 1/r}, r > 1, has integer trace r + 1/r >= 3
    (trace 2 forces r = 1), hence r >= (3+sqrt5)/2 = phi^2, giving
    M in {1} cup [phi^2, infty).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp

from _z5_engine import phi, is_charge5, charge_group, mahler_hp

mp.mp.dps = 50


# --------------------------- Lemma 2.2 -------------------------------------
def test_pentagon_lattice_omits_180_degrees():
    """Lemma 2.2 proof: (2pi/5)Z = {0,72,144,216,288} deg does NOT contain 180."""
    lattice = {sp.Rational(k, 5) for k in range(5)}          # fractions of a turn
    assert lattice == {sp.Rational(0), sp.Rational(1, 5), sp.Rational(2, 5),
                       sp.Rational(3, 5), sp.Rational(4, 5)}
    assert sp.Rational(1, 2) not in lattice                  # 180 deg absent
    # equivalently, no lattice angle equals pi
    assert all((sp.Rational(k, 5) * 2 - 1) != 0 for k in range(5))


def test_fifth_roots_of_unity_sit_on_lattice():
    """Lemma 2.2: on-unit charge-5 roots are the fifth roots of unity."""
    for k in range(5):
        z = mp.expjpi(mp.mpf(2 * k) / 5)                     # exp(2 pi i k/5)
        assert mp.almosteq(abs(z), mp.mpf(1), abs_eps=mp.mpf('1e-40'))
        theta = (mp.arg(z) / (2 * mp.pi)) % 1
        assert mp.almosteq((5 * theta) % 1, 0, abs_eps=mp.mpf('1e-30')) or \
               mp.almosteq((5 * theta) % 1, 1, abs_eps=mp.mpf('1e-30'))
    # Phi_5 (the fifth-root-of-unity minimal poly) is charge exactly 5, M=1
    ok, M = is_charge5([1, 1, 1, 1, 1])
    assert ok and mp.almosteq(M, 1, abs_eps=mp.mpf('1e-9'))


def test_negative_real_root_forces_even_charge():
    """Lemma 2.2: an object with a negative real root cannot have charge 5.
    (x+1)*Phi_5 carries the root -1 (angle 180 deg, denominator 2), so its
    charge is lcm(2,5)=10, not 5."""
    # (x+1)(x^4+x^3+x^2+x+1) = x^5 + 2x^4 + 2x^3 + 2x^2 + 2x + 1
    p = [1, 2, 2, 2, 2, 1]
    ok, _ = is_charge5(p)
    assert not ok
    assert charge_group(p) == 10
    # x^5 - 2 (charge 5) has exactly one real root and it is positive
    roots = mp.polyroots([1, 0, 0, 0, 0, -2], maxsteps=100, extraprec=80)
    reals = [r for r in roots if abs(mp.im(r)) < mp.mpf('1e-30')]
    assert len(reals) == 1 and mp.re(reals[0]) > 0


# --------------------------- Lemma 2.5 -------------------------------------
def test_reciprocal_unit_trace_two_forces_r_equal_one():
    """Lemma 2.5: r + 1/r = 2 forces r = 1."""
    r = sp.symbols('r', positive=True)
    sol = sp.solve(sp.Eq(r + 1 / r, 2), r)
    assert sol == [sp.Integer(1)]


def test_reciprocal_unit_trace_three_gives_phi_squared():
    """Lemma 2.5: the smallest integer trace > 2 is 3, giving
    r = (3+sqrt5)/2 = phi^2 ~ 2.618, the forced Mahler floor for a real
    reciprocal unit."""
    r = sp.symbols('r', positive=True)
    # larger root of r + 1/r = 3  <=>  r^2 - 3r + 1 = 0
    sol = [s for s in sp.solve(r**2 - 3 * r + 1, r) if s > 1]
    assert len(sol) == 1
    rval = sol[0]
    assert sp.simplify(rval - (3 + sp.sqrt(5)) / 2) == 0
    assert sp.simplify(rval - phi**2) == 0
    # its Mahler contribution is r itself
    assert mp.almosteq(mp.mpf(str(sp.N(rval, 40))), mp.mpf('2.618033988749895'),
                       abs_eps=mp.mpf('1e-15'))


def test_reciprocal_unit_general_trace_formula_and_monotonicity():
    """Lemma 2.5: for integer trace k >= 3 the unit is r = (k+sqrt(k^2-4))/2,
    strictly increasing in k; k=3 -> phi^2, k=4 -> 2+sqrt3.  Hence any real
    reciprocal unit off the circle has M in {1} cup [phi^2, infty)."""
    def r_of(k):
        return (mp.mpf(k) + mp.sqrt(k * k - 4)) / 2

    r3, r4 = r_of(3), r_of(4)
    assert mp.almosteq(r3, mp.mpf(str(sp.N(phi**2, 40))), abs_eps=mp.mpf('1e-20'))
    assert mp.almosteq(r4, 2 + mp.sqrt(3), abs_eps=mp.mpf('1e-30'))
    assert r3 < r4                                            # monotone up in k
    # each satisfies r + 1/r = k exactly
    for k in (3, 4, 5, 6):
        rk = r_of(k)
        assert mp.almosteq(rk + 1 / rk, k, abs_eps=mp.mpf('1e-30'))
        assert rk >= r3                                       # floor is phi^2


def test_real_reciprocal_unit_gap_is_empty_below_phi_squared():
    """Lemma 2.5: no real reciprocal unit has M in (1, phi^2).  Scan every
    integer trace and confirm the induced Mahler measure is 1 or >= phi^2."""
    phi2 = mp.mpf(str(sp.N(phi**2, 45)))
    for k in range(2, 25):
        disc = k * k - 4
        if disc < 0:
            continue
        r = (mp.mpf(k) + mp.sqrt(disc)) / 2      # dominant root, = M contribution
        assert (mp.almosteq(r, 1, abs_eps=mp.mpf('1e-30'))
                or r >= phi2 - mp.mpf('1e-30'))
