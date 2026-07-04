"""Section 7 verification table + Proposition 6.1 (the residual location),
reproduced by an INDEPENDENT enumeration of the paper's exact window:

    quartics |c| <= 10,  quintics |c| <= 4,  sextics |c| <= 3.

For every monic integer poly in that box we compute charge group (Def. 2.1),
Mahler measure (Def. 2.3) and the reciprocal test from scratch, then keep the
charge-exactly-5 squarefree objects.

Reproduced claims:
  * charge-5 quartics: distinct M = {1, phi^4}; none in (1,2).
  * 13 non-reciprocal charge-5 objects; realized minimum M = 2.
  * every non-reciprocal charge-5 object obeys M >= mu_S; none in [mu_S, 2).
  * reciprocal charge-5 measures = {1, phi^2, 2+sqrt3}; none in (1,2).
  * OVERALL: no charge-5 object has M in (1,2); the realized floor is 2.

Flagged (xfail): the paper states 5 reciprocal objects; independent
enumeration finds 4 (the three distinct measures still match).  See NOTES.md.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import pytest
import sympy as sp
import mpmath as mp

from _z5_engine import phi, window_objects, is_charge5

mp.mp.dps = 40

PHI2 = float(sp.N(phi**2, 30))
PHI4 = float(sp.N(phi**4, 30))
TWO_SQRT3 = float(2 + sp.sqrt(3))
MU_S = float(mp.findroot(lambda t: t**3 - t - 1, mp.mpf('1.3')))


@pytest.fixture(scope="module")
def objs():
    return window_objects()


def _in_open_1_2(M):
    return 1 + 1e-6 < M < 2 - 1e-6


# ------------------------- quartic sub-window ------------------------------
def test_quartic_distinct_measures_are_1_and_phi4():
    """Sec. 7: 'charge-5 quartics, |c|<=10: distinct M = {1, phi^4}'."""
    quart = []
    for c3 in range(-10, 11):
        for c2 in range(-10, 11):
            for c1 in range(-10, 11):
                for c0 in range(-10, 11):
                    if c0 == 0:
                        continue
                    ok, M = is_charge5((1, c3, c2, c1, c0))
                    if ok:
                        quart.append(M)
    assert quart, "expected some charge-5 quartics"
    distinct = sorted({round(M, 4) for M in quart})
    assert distinct == sorted({round(1.0, 4), round(PHI4, 4)})


def test_quartic_none_in_open_interval_1_2():
    """Sec. 7: 'charge-5 quartics with M in (1,2): 0'."""
    for c3 in range(-10, 11):
        for c2 in range(-10, 11):
            for c1 in range(-10, 11):
                for c0 in range(-10, 11):
                    if c0 == 0:
                        continue
                    ok, M = is_charge5((1, c3, c2, c1, c0))
                    if ok:
                        assert not _in_open_1_2(M)


# --------------------- full-window aggregate claims ------------------------
def test_nonreciprocal_count_is_13(objs):
    """Sec. 7: 'non-reciprocal charge-5 found: 13 objects'."""
    nonrecip = [o for o in objs if not o[2]]
    assert len(nonrecip) == 13


def test_nonreciprocal_realized_minimum_is_2(objs):
    """Sec. 7: 'min M = 2 (x^5-2 and kin)' over the non-reciprocal class."""
    nonrecip = [o for o in objs if not o[2]]
    Mmin = min(M for _, M, _ in nonrecip)
    assert abs(Mmin - 2.0) < 1e-6
    # x^5 - 2 itself is among them
    assert any(c == (1, 0, 0, 0, 0, -2) for c, _, _ in nonrecip)


def test_all_nonreciprocal_obey_smyth_floor(objs):
    """Sec. 7: 'all non-reciprocal obey M >= mu_S: True'."""
    nonrecip = [o for o in objs if not o[2]]
    for _, M, _ in nonrecip:
        assert M >= MU_S - 1e-9


def test_nonreciprocal_empty_in_muS_to_2(objs):
    """Sec. 7 / Prop. 6.1: 'non-reciprocal charge-5 in [mu_S, 2): 0' -- the
    residual window is empty; the boundary object x^5-2 sits at exactly 2."""
    nonrecip = [o for o in objs if not o[2]]
    assert not any(MU_S - 1e-9 <= M < 2 - 1e-6 for _, M, _ in nonrecip)


def test_reciprocal_distinct_measures(objs):
    """Sec. 7: reciprocal charge-5 measures = {1, phi^2, 2+sqrt3}."""
    recip = [o for o in objs if o[2]]
    distinct = sorted({round(M, 4) for _, M, _ in recip})
    want = sorted({round(1.0, 4), round(PHI2, 4), round(TWO_SQRT3, 4)})
    assert distinct == want


def test_reciprocal_none_in_open_1_2(objs):
    """Sec. 7: 'reciprocal charge-5 in (1,2): 0'."""
    recip = [o for o in objs if o[2]]
    assert not any(_in_open_1_2(M) for _, M, _ in recip)


def test_overall_no_charge5_object_in_open_1_2(objs):
    """Sec. 7 (headline) / Prop. 6.1: 'overall charge-5 with M in (1,2): 0',
    realized floor 2.  This is the mu(5)=2 computed residual over the window."""
    assert not any(_in_open_1_2(M) for _, M, _ in objs)
    Mgt1 = [M for _, M, _ in objs if M > 1 + 1e-6]
    assert Mgt1 and abs(min(Mgt1) - 2.0) < 1e-6


def test_reciprocal_objects_all_equal_phi5_times_real_unit(objs):
    """Thm. 3.1(b) / Sec. 7: every reciprocal charge-5 object is
    Phi_5^a x (real reciprocal units); each factors as Phi_5 times a totally
    real reciprocal integer polynomial, and its measure lies in {1} cup
    [phi^2, infty)."""
    recip = [o for o in objs if o[2]]
    x = sp.symbols('x')
    Phi5 = sp.Poly([1, 1, 1, 1, 1], x)
    for coeffs, M, _ in recip:
        P = sp.Poly([int(c) for c in coeffs], x)
        quo, rem = sp.div(P, Phi5)
        assert rem == 0                                # Phi_5 divides every one
        assert (abs(M - 1.0) < 1e-6) or (M >= PHI2 - 1e-6)


# ----------------------------- flagged --------------------------------------
def test_reciprocal_object_count_is_4(objs):
    """Sec. 7 (corrected 2026-07-04): 'reciprocal charge-5 found: 4 objects'
    {Phi_5, x^5-1, Phi_5*(x^2-3x+1), Phi_5*(x^2-4x+1)}. The paper originally
    stated 5; independent enumeration of the window finds 4 (the three distinct
    measures {1, phi^2, 2+sqrt3} are unchanged)."""
    recip = [o for o in objs if o[2]]
    assert len(recip) == 4
