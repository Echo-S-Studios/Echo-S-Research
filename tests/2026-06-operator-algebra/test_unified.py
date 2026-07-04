"""
Section 8 -- The unified picture.

Theorem 8.1 (thm:unified): the spectral operators form a commutative lambda-ring
with two characters (Mahler measure, angle charge).  The computationally testable
content beyond Thms 2.2-5.2 is the meeting point:

    "sqrt5 = phi + phi^{-1} is where the two characters meet -- the (x)-logarithm
     of the floor and the grow generator of the self-action."

We verify the exact identity sqrt5 = phi + 1/phi and the companion identity
phi - 1/phi = 1, plus the lambda-ring endomorphism axiom on a concrete object.
"""

import sympy as sp

from _opalg_ops import (
    golden_seed,
    K_seed,
    mahler_exact,
    ms_equal,
    oplus,
    otimes,
    phi,
    psi,
    seed_from_poly,
    x,
)

sqrt5 = sp.sqrt(5)


def test_sqrt5_is_phi_plus_inverse():
    """sqrt5 = phi + phi^{-1}  (exact)."""
    assert sp.simplify(sqrt5 - (phi + 1 / phi)) == 0


def test_phi_minus_inverse_is_one():
    """phi - phi^{-1} = 1  (companion identity; phi^{-1} = phi - 1)."""
    assert sp.simplify((phi - 1 / phi) - 1) == 0
    assert sp.simplify((1 / phi) - (phi - 1)) == 0


def test_phi_is_measure_floor_and_grow_generator():
    """phi = M(golden) is both the least Mahler generator and the additive
    self-action rate base (M((+)^k) = phi^k)."""
    A = golden_seed()
    assert sp.simplify(mahler_exact(A) - phi) == 0
    assert sp.simplify(mahler_exact(oplus(A, A)) - phi**2) == 0


def test_lambda_ring_endomorphism_axiom_concrete():
    """psi^2 respects BOTH operations on a concrete mixed object (lambda-ring):
    psi^2(A(+)B)=psi^2A(+)psi^2B and psi^2(A(x)B)=psi^2A(x)psi^2B."""
    A = golden_seed()
    B = K_seed()
    assert ms_equal(psi(2, oplus(A, B)), oplus(psi(2, A), psi(2, B)))
    assert ms_equal(psi(2, otimes(A, B)), otimes(psi(2, A), psi(2, B)))


def test_two_characters_read_off_same_object():
    """On K both characters are simultaneously defined and consistent:
    Character I gives M(K)=beta^2=phi^2 sqrt5; Character II gives full Z/4Z."""
    from _opalg_ops import charge

    K = K_seed()
    assert sp.simplify(mahler_exact(K) - phi**2 * sqrt5) == 0
    assert set(charge(K).keys()) == {0, 1, 2, 3}
