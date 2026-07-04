"""Section 7.2: the Pisot-quartic sweep and the emission-gap probe.

Box [-3,3]^4 (2401 monic quartics): 103 certified Pisot quartics, of which 102
are complex-pair and 1 is totally real.  The totally-real hit is
x^4-3x^3-2x^2+2x+1 (theta ~ 3.390 > phi), the smallest complex-pair instance is
x^4-x^3-1 (theta ~ 1.3803), and P4's emission-gap prediction theta >= phi holds
for the totally-real subfamily.
"""
import os
import sys

from sympy import symbols, Poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import is_pisot, dominant_root, phi

x = symbols('x')


def _run_sweep():
    pis = 0
    totally_real = 0
    complex_pair = 0
    tr_hit = None
    min_cp = None            # (theta, coeffs) of smallest-theta complex-pair
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                for d in range(-3, 4):
                    coeffs = [1, a, b, c, d]
                    if d == 0:
                        continue
                    if not Poly(coeffs, x).is_irreducible:
                        continue
                    ok, n_real_in, n_pairs = is_pisot(coeffs, 40)
                    if not ok:
                        continue
                    pis += 1
                    if n_pairs == 0:
                        totally_real += 1
                        tr_hit = coeffs
                    elif n_pairs == 1:
                        complex_pair += 1
                        th = dominant_root(coeffs, 40)
                        if min_cp is None or th < min_cp[0]:
                            min_cp = (th, coeffs)
    return pis, totally_real, complex_pair, tr_hit, min_cp


_CACHE = None


def sweep():
    global _CACHE
    if _CACHE is None:
        _CACHE = _run_sweep()
    return _CACHE


def test_box_size_is_7_pow_4():
    assert 7 ** 4 == 2401


def test_pisot_quartic_counts():
    """103 certified Pisot quartics = 102 complex-pair + 1 totally real."""
    pis, tr, cp, _, _ = sweep()
    assert pis == 103
    assert tr == 1
    assert cp == 102
    assert tr + cp == pis


def test_unique_totally_real_hit_and_gap():
    """The unique totally-real Pisot quartic in the box is x^4-3x^3-2x^2+2x+1,
    with theta ~ 3.390 satisfying theta > phi (emission-gap not falsified)."""
    import mpmath as mp
    _, _, _, tr_hit, _ = sweep()
    assert tr_hit == [1, -3, -2, 2, 1]
    th = dominant_root([1, -3, -2, 2, 1], 50)
    with mp.workdps(50):
        assert abs(th - mp.mpf('3.3903')) < mp.mpf('1e-3')
        assert th > phi(50)


def test_smallest_complex_pair_instance():
    """The smallest (minimal-theta) complex-pair Pisot quartic in the box is
    x^4-x^3-1 with theta ~ 1.3803."""
    import mpmath as mp
    _, _, _, _, min_cp = sweep()
    assert min_cp[1] == [1, -1, 0, 0, -1]
    ok, n_real_in, n_pairs = is_pisot([1, -1, 0, 0, -1], 50)
    assert ok and n_pairs == 1
    with mp.workdps(50):
        assert abs(min_cp[0] - mp.mpf('1.3803')) < mp.mpf('1e-3')


def test_golden_ratio_is_the_degree2_totally_real_extreme():
    """P4: totally-real Pisot numbers satisfy theta >= phi with equality only at
    theta = phi (degree 2).  Sanity: phi is a Pisot number, its conjugate -1/phi
    is inside the unit disk, and it is a root of x^2-x-1."""
    import mpmath as mp
    with mp.workdps(50):
        p = phi(50)
        assert abs(p * p - p - 1) < mp.mpf('1e-40')      # golden ratio identity
        assert abs(-1 / p) < 1                            # conjugate inside disk
    ok, n_real_in, n_pairs = is_pisot([1, -1, -1], 50)    # x^2-x-1
    assert ok and n_pairs == 0
