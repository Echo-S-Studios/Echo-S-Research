"""Sanity / boundary checks for the argumentative claims that have no single
closed-form number: the Kronecker and Smyth lemmas, multiplicativity of the
Mahler measure, the Z/3Z contrast that motivates the paper, and the
totally-positive (1,2) gap that Remark 6.2's realification leans on.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import itertools

import numpy as np
import sympy as sp
import mpmath as mp

from _z5_engine import phi, mahler_hp, mahler_np, is_reciprocal, charge_group

mp.mp.dps = 45
MU_S = mp.findroot(lambda t: t**3 - t - 1, mp.mpf('1.3'))


def test_mahler_multiplicativity():
    """Def. 2.3: M(pq) = M(p) M(q).  Checked on a few integer polynomials."""
    x = sp.symbols('x')
    pairs = [([1, 0, -2], [1, -3, 1]),          # (x^2-2)(x^2-3x+1)
             ([1, 1, 1, 1, 1], [1, 0, 0, 0, 0, -2]),   # Phi_5 * (x^5-2)
             ([1, -1, -1], [1, 0, 0, -2])]      # (x^2-x-1)(x^3-2)
    for p, q in pairs:
        Pp, Pq = sp.Poly(p, x), sp.Poly(q, x)
        prod = (Pp * Pq).all_coeffs()
        lhs = mahler_hp([int(c) for c in prod])
        rhs = mahler_hp(p) * mahler_hp(q)
        assert mp.almosteq(lhs, rhs, abs_eps=mp.mpf('1e-20'))


def test_kronecker_cyclotomics_have_measure_one():
    """Def. 2.3 (Kronecker): M = 1 exactly on products of cyclotomics / powers
    of x.  Every cyclotomic Phi_n has M=1; a non-cyclotomic integer poly does
    not."""
    x = sp.symbols('x')
    for n in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12):
        coeffs = [int(c) for c in sp.Poly(sp.cyclotomic_poly(n, x), x).all_coeffs()]
        assert mp.almosteq(mahler_hp(coeffs), mp.mpf(1), abs_eps=mp.mpf('1e-20'))
    # a non-cyclotomic example exceeds 1
    assert mahler_hp([1, 0, -2]) > 1                      # x^2 - 2, M = sqrt2


def test_smyth_bound_on_nonreciprocal_examples():
    """Lemma 2.4 (Smyth): every non-reciprocal integer poly with M>1 has
    M >= mu_S.  mu_S itself is realized by x^3-x-1 (non-reciprocal)."""
    muS_poly = [1, 0, -1, -1]                             # x^3 - x - 1
    assert not is_reciprocal(muS_poly)
    assert mp.almosteq(mahler_hp(muS_poly), MU_S, abs_eps=mp.mpf('1e-25'))
    # a spread of non-reciprocal polys, all with M > 1, must clear mu_S
    examples = [[1, -1, -1],            # x^2-x-1, M = phi
                [1, 0, 0, -2],          # x^3-2,   M = 2
                [1, 0, 0, 0, 0, -2],    # x^5-2,   M = 2
                [1, -1, -1, 0, 0],      # (x^2-x-1) x^2 -> handled below w/o x-root
                [1, 0, -3]]             # x^2-3,   M = sqrt3
    for p in examples:
        if p[-1] == 0:
            continue
        assert is_reciprocal(p) or mahler_hp(p) >= MU_S - mp.mpf('1e-20')


def test_z3_contrast_collapses_to_x3_minus_2():
    """Sec. 1.2: the Z/3Z proof closes because 2 cos 120 deg = -1 is rational,
    collapsing the search to x^3-2 with M=2.  Verify the rational cosine, the
    charge, and the measure."""
    assert sp.simplify(2 * sp.cos(sp.rad(120)) + 1) == 0     # 2cos120 = -1 in Z
    assert charge_group([1, 0, 0, -2]) == 3                  # x^3-2 has charge 3
    assert mp.almosteq(mahler_hp([1, 0, 0, -2]), mp.mpf(2), abs_eps=mp.mpf('1e-25'))


def test_totally_positive_gap_floor_is_two():
    """Remark 6.2(1) premise: a totally positive algebraic integer (!= 1) has
    M >= 2 (the totally-positive (1,2) gap).  Sanity-scan irreducible degree
    <=3 monic integer polys with all-positive real roots: none lands in (1,2),
    and the floor above 1 is exactly 2 (e.g. x-2)."""
    x = sp.symbols('x')
    found_min = None
    # degree 1: x - a
    for a in range(1, 6):
        M = float(a)
        if M > 1 + 1e-9:
            found_min = M if found_min is None else min(found_min, M)
        assert not (1 + 1e-6 < M < 2 - 1e-6)
    # degree 2 and 3: enumerate small coeffs, keep totally-positive squarefree
    # (numpy screening -- margins to the (1,2) boundary are large)
    for deg, b in [(2, 6), (3, 6)]:
        for tail in itertools.product(range(-b, b + 1), repeat=deg):
            if tail[-1] == 0:
                continue
            coeffs = [1] + list(tail)
            roots = np.roots(coeffs)
            if any(abs(r.imag) > 1e-9 or r.real <= 1e-9 for r in roots):
                continue                       # not totally positive
            rr = sorted(roots, key=lambda z: z.real)
            if any(abs(rr[i] - rr[i + 1]) < 1e-7 for i in range(len(rr) - 1)):
                continue                       # not squarefree
            M = mahler_np(roots)
            assert not (1 + 1e-6 < M < 2 - 1e-6)
            if M > 1 + 1e-9:
                found_min = M if found_min is None else min(found_min, M)
    assert found_min is not None and abs(found_min - 2.0) < 1e-6
