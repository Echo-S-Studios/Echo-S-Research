"""Angle confinement: the (pi/2)Z invariant and its consequences
(Sec. 4: Lemmas 4.1-4.4; Lemma 5.1; App. A).
"""
import mpmath as mp
import numpy as np
import sympy as sp

import emgap_util as U

mp.mp.dps = 40

x = U.x


def test_subgroup_closed_under_add_and_double():
    """Lemma 4.2 proof: (pi/2)Z = {0, pi/2, pi, 3pi/2} is closed under addition
    (kron) and under doubling (squaring); modelled as Z/4 = {0,1,2,3}."""
    S = {0, 1, 2, 3}
    for a in S:
        for b in S:
            assert (a + b) % 4 in S          # kron: args add
        assert (2 * a) % 4 in S              # squaring: args double


def _companion(coeffs):
    return U.companion(coeffs)


def test_kron_phi_phi_yields_minus_one_on_circle():
    """App. A: for A = C(x^2-x-1), A (x) A has the on-circle eigenvalue -1
    (arg 180 deg), namely the product phi*psi = -1 of the two golden roots."""
    A = _companion([1, -1, -1])
    K = np.kron(A, A)
    ev = np.linalg.eigvals(K)
    on_circle = [e for e in ev if abs(abs(e) - 1) < 1e-9]
    assert any(abs(e - (-1)) < 1e-9 for e in on_circle)
    # every on-circle eigenvalue is a 4th root of unity
    for e in on_circle:
        assert min(abs(e - z) for z in (1, -1, 1j, -1j)) < 1e-9


def test_spectral_products_keep_on_circle_eigs_fourth_roots():
    """Lemma 4.3: on-circle eigenvalues of kron / direct-sum / squares of
    catalog companions are all in {1, i, -1, -i}."""
    catalog = [[1, -1, -1], [1, 1, -1], [1, 0, -2], [1, 0, -3],
               [1, 0, -5], [1, -7, 1], [1, 0, 5, 0, -5]]
    comps = [_companion(c) for c in catalog]
    fourth_roots = (1, -1, 1j, -1j)
    for A in comps:
        for B in comps:
            for M in (np.kron(A, B),
                      np.block([[A, np.zeros((A.shape[0], B.shape[1]))],
                                [np.zeros((B.shape[0], A.shape[1])), B]]),
                      A @ A):
                for e in np.linalg.eigvals(M):
                    if abs(abs(e) - 1) < 1e-9:
                        assert min(abs(e - z) for z in fourth_roots) < 1e-8


def test_lehmer_oncircle_conjugates_not_roots_of_unity():
    """Lemma 5.1 / App. A: the Lehmer polynomial is irreducible with a root > 1
    (so not cyclotomic); its 8 on-circle conjugates therefore are NOT roots of
    unity, and their arguments (~62.81, 107.00, 137.26, 160.61 deg and mirrors)
    are not multiples of 90 deg."""
    Lc = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    Lp = sp.Poly(Lc, x)
    # irreducible over Q
    assert len(sp.factor_list(Lp.as_expr(), x)[1]) == 1
    rts = U.mp_roots(Lc)
    # a real root > 1 => not a product of cyclotomics => on-circle conj not roots of unity
    assert any(r.imag == 0 and r.real > 1 for r in rts) or \
        any(abs(r.imag) < 1e-30 and r.real > 1 for r in rts)
    on_circle = sorted(
        [mp.degrees(mp.atan2(r.imag, r.real)) % 360 for r in rts
         if abs(abs(r) - 1) < 1e-20]
    )
    assert len(on_circle) == 8
    expected = [62.81, 107.00, 137.26, 160.61, 199.39, 222.74, 253.00, 297.19]
    for got, exp in zip(on_circle, expected):
        assert abs(float(got) - exp) < 0.02
    # none is a multiple of 90 deg (not a 4th root of unity)
    for ang in on_circle:
        m = float(ang) % 90
        assert min(m, 90 - m) > 1.0
