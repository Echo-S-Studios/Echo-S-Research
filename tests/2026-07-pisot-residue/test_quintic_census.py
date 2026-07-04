"""Section 6: the exhaustive quintic execution.

Independently re-runs the certification cascade over the 3125 monic quintics
x^5+ax^4+bx^3+cx^2+dx+e with (a,...,e) in [-2,2]^5:

  reject e=0 (625) | reject +-reciprocal (50) | reject reducible (638) |
  reject real pattern (1318) | reject disk certificate (411) | Pisot (83),

then the pattern split real5=0 / mixed=16 / two-pair=67, and the smallest
two-pair instance x^5-2x^4-2x^3-2x^2-2x-2 (Theorem 6.1).

Pisot classification uses exact irreducibility plus high-precision root moduli;
this is sound because the certified set is irreducible and non-reciprocal, hence
has NO root on the unit circle (a genuine gap to classify against).
"""
import os
import sys

from sympy import symbols, Poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import is_pm_reciprocal, is_pisot, hp_roots

x = symbols('x')


def _run_census():
    import mpmath as mp
    tol = mp.mpf(10) ** -12
    rej = dict(e0=0, recip=0, reducible=0, realpat=0, disk=0)
    pis = 0
    patt = dict(real5=0, mixed=0, twopair=0)
    first_twopair = None
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                for d in range(-2, 3):
                    for e in range(-2, 3):
                        coeffs = [1, a, b, c, d, e]
                        if e == 0:
                            rej['e0'] += 1
                            continue
                        if is_pm_reciprocal(coeffs):
                            rej['recip'] += 1
                            continue
                        if not Poly(coeffs, x).is_irreducible:
                            rej['reducible'] += 1
                            continue
                        # single high-precision root computation, reused
                        rts = hp_roots(coeffs, 40)
                        reals = [r.real for r in rts if abs(r.imag) < tol]
                        gt1 = [r for r in reals if r > 1]
                        lem1 = [r for r in reals if r <= -1]
                        if len(gt1) != 1 or len(lem1) > 0:
                            rej['realpat'] += 1
                            continue
                        # disk certificate: every non-dominant root strictly inside
                        oncirc = any(abs(abs(r) - 1) < tol for r in rts)
                        outside = [r for r in rts if abs(r) > 1]
                        if oncirc or len(outside) != 1 or abs(outside[0].imag) > tol:
                            rej['disk'] += 1
                            continue
                        pis += 1
                        nonreal = [r for r in rts if abs(r.imag) > tol]
                        npairs = len(nonreal) // 2
                        if npairs == 0:
                            patt['real5'] += 1
                        elif npairs == 1:
                            patt['mixed'] += 1
                        elif npairs == 2:
                            patt['twopair'] += 1
                            if first_twopair is None:
                                first_twopair = coeffs
    return rej, pis, patt, first_twopair


# run once; cache on the module
_CACHE = None


def census():
    global _CACHE
    if _CACHE is None:
        _CACHE = _run_census()
    return _CACHE


def test_box_size_is_5_pow_5():
    assert 5 ** 5 == 3125


def test_reject_tally_and_certified_count():
    """Stage-1 tallies: 625 / 50 / 638 / 1318 / 411 rejects, 83 certified Pisot,
    and the whole cascade partitions the 3125-box."""
    rej, pis, _, _ = census()
    assert rej['e0'] == 625
    assert rej['recip'] == 50
    assert rej['reducible'] == 638
    assert rej['realpat'] == 1318
    assert rej['disk'] == 411
    assert pis == 83
    assert sum(rej.values()) + pis == 3125


def test_pattern_split():
    """Patterns among the 83: real5 = 0, mixed = 16, two-pair = 67; and
    16 + 67 = 83."""
    _, pis, patt, _ = census()
    assert patt['real5'] == 0
    assert patt['mixed'] == 16
    assert patt['twopair'] == 67
    assert patt['mixed'] + patt['twopair'] == pis == 83


def test_reciprocal_reject_count_is_exactly_50():
    """The +-reciprocal reject count is a clean combinatorial 25+25=50
    (palindromic: e=1, d=a, c=b -> 25; anti: e=-1, d=-a, c=-b -> 25)."""
    n = 0
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                for d in range(-2, 3):
                    for e in range(-2, 3):
                        if e == 0:
                            continue
                        if is_pm_reciprocal([1, a, b, c, d, e]):
                            n += 1
    assert n == 50


def test_smallest_two_pair_instance():
    """Theorem 6.1 / Sec. 6: the first two-pair instance in enumeration order is
    x^5-2x^4-2x^3-2x^2-2x-2, and it is a genuine two-pair Pisot quintic."""
    _, _, _, first = census()
    assert first == [1, -2, -2, -2, -2, -2]
    ok, n_real_in, n_pairs = is_pisot([1, -2, -2, -2, -2, -2])
    assert ok and n_pairs == 2 and n_real_in == 0
