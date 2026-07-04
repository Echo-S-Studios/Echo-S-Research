"""
Theorem 2.2 (thm:semiring) -- "The emission semiring" [OA-SR-01..08].

    "(S, (+), (x), 0, 1) is a commutative semiring: (+) and (x) are commutative
     and associative, (x) distributes over (+), A (x) (B (+) C) = (A(x)B)(+)(A(x)C),
     0 is the additive identity and a multiplicative annihilator, and 1 is the
     multiplicative identity."

Objects here are multisets of *symbols* (generic eigenvalues); the semiring
laws are multiset identities, so we verify them exactly over the free commutative
setting -- if they hold for independent symbols they hold for every substitution.
"""

import sympy as sp

from _opalg_ops import ONE, ZERO, ms_equal, oplus, otimes

# Generic eigenvalue multisets (independent symbols).
a1, a2, a3 = sp.symbols("a1 a2 a3")
b1, b2 = sp.symbols("b1 b2")
c1, c2 = sp.symbols("c1 c2")
A = (a1, a2, a3)
B = (b1, b2)
C = (c1, c2)


def test_oplus_commutative():
    """(+) is commutative:  A (+) B = B (+) A  (multiset union)."""
    assert ms_equal(oplus(A, B), oplus(B, A))


def test_oplus_associative():
    """(+) is associative:  (A (+) B) (+) C = A (+) (B (+) C)."""
    assert ms_equal(oplus(oplus(A, B), C), oplus(A, oplus(B, C)))


def test_otimes_commutative():
    """(x) is commutative:  A (x) B = B (x) A."""
    assert ms_equal(otimes(A, B), otimes(B, A))


def test_otimes_associative():
    """(x) is associative:  (A (x) B) (x) C = A (x) (B (x) C)."""
    assert ms_equal(otimes(otimes(A, B), C), otimes(A, otimes(B, C)))


def test_distributivity():
    """(x) distributes over (+):  A (x) (B (+) C) = (A(x)B) (+) (A(x)C)."""
    lhs = otimes(A, oplus(B, C))
    rhs = oplus(otimes(A, B), otimes(A, C))
    assert ms_equal(lhs, rhs)


def test_right_distributivity():
    """Right distributivity too (commutative semiring): (B (+) C) (x) A = ..."""
    lhs = otimes(oplus(B, C), A)
    rhs = oplus(otimes(B, A), otimes(C, A))
    assert ms_equal(lhs, rhs)


def test_zero_additive_identity():
    """0 = empty multiset is the additive identity:  A (+) 0 = A."""
    assert ms_equal(oplus(A, ZERO), A)
    assert ms_equal(oplus(ZERO, A), A)


def test_zero_multiplicative_annihilator():
    """0 annihilates under (x):  A (x) 0 = 0."""
    assert ms_equal(otimes(A, ZERO), ZERO)
    assert ms_equal(otimes(ZERO, A), ZERO)


def test_one_multiplicative_identity():
    """1 = {1} is the multiplicative identity:  A (x) 1 = A."""
    assert ms_equal(otimes(A, ONE), A)
    assert ms_equal(otimes(ONE, A), A)


def test_cardinality_bookkeeping():
    """deg(A(+)B) = deg A + deg B and deg(A(x)B) = deg A * deg B."""
    assert len(oplus(A, B)) == len(A) + len(B)
    assert len(otimes(A, B)) == len(A) * len(B)
