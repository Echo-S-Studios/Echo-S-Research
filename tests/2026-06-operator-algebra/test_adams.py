"""
Theorem 3.2 (thm:adams) -- "Squaring is the Adams operation" [OA-PSI-01..05]
Proposition 3.3 (prop:plethysm) -- "the square is a diagonal" [OA-PLE-01..06].

Claims verified:
 * psi^n is a semiring endomorphism:
       psi^n(A(+)B) = psi^n A (+) psi^n B,
       psi^n(A(x)B) = psi^n A (x) psi^n B,
       psi^n(1) = 1,   psi^m . psi^n = psi^{mn}
 * Newton's identity  p_2 = e_1^2 - 2 e_2   (and general Newton, k<=4)
 * Plethysm  A(x)A = Sym^2 A (+) Lambda^2 A ,  psi^2 A = diagonal of Sym^2 A
 * traces  tr Sym^2 A = (p_1^2 + p_2)/2 ,  tr Lambda^2 A = (p_1^2 - p_2)/2 = e_2
 * psi^{2^k} = (psi^2)^{o k} sends lambda -> lambda^{2^k}
"""

import sympy as sp

from _opalg_ops import (
    ONE,
    elem_sym,
    ms_equal,
    oplus,
    otimes,
    power_sum,
    psi,
    sym2,
    wedge2,
)

a1, a2, a3, a4, a5 = sp.symbols("a1 a2 a3 a4 a5")
b1, b2 = sp.symbols("b1 b2")


# --------------------------------------------------------------------------
# psi^n as a semiring endomorphism
# --------------------------------------------------------------------------
def test_psi_additive():
    """psi^n(A (+) B) = psi^n A (+) psi^n B  for n = 2,3,5."""
    A, B = (a1, a2, a3), (b1, b2)
    for n in (2, 3, 5):
        assert ms_equal(psi(n, oplus(A, B)), oplus(psi(n, A), psi(n, B)))


def test_psi_multiplicative():
    """psi^n(A (x) B) = psi^n A (x) psi^n B  (uses (l*m)^n = l^n * m^n)."""
    A, B = (a1, a2), (b1, b2)
    for n in (2, 3, 4):
        assert ms_equal(psi(n, otimes(A, B)), otimes(psi(n, A), psi(n, B)))


def test_psi_unit_fixed():
    """psi^n(1) = 1  for every n."""
    for n in (2, 3, 7):
        assert ms_equal(psi(n, ONE), ONE)


def test_psi_composition():
    """psi^m . psi^n = psi^{mn}."""
    A = (a1, a2, a3)
    for m in (2, 3):
        for n in (2, 5):
            assert ms_equal(psi(m, psi(n, A)), psi(m * n, A))


# --------------------------------------------------------------------------
# Newton's identities
# --------------------------------------------------------------------------
def test_newton_p2():
    """p_2 = e_1^2 - 2 e_2  (the identity the paper singles out)."""
    for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4), (a1, a2, a3, a4, a5)]:
        p2 = power_sum(A, 2)
        e1 = elem_sym(A, 1)
        e2 = elem_sym(A, 2)
        assert sp.expand(p2 - (e1**2 - 2 * e2)) == 0


def test_newton_general_recursion():
    """General Newton p_k = sum_{i=1}^{k} (-1)^{i-1} e_i p_{k-i}, checked k<=4."""
    A = (a1, a2, a3, a4, a5)
    for k in (1, 2, 3, 4):
        pk = power_sum(A, k)
        rhs = sp.Integer(0)
        for i in range(1, k):
            rhs += (-1) ** (i - 1) * elem_sym(A, i) * power_sum(A, k - i)
        rhs += (-1) ** (k - 1) * k * elem_sym(A, k)
        assert sp.expand(pk - rhs) == 0


# --------------------------------------------------------------------------
# Plethysm  A (x) A = Sym^2 A (+) Lambda^2 A
# --------------------------------------------------------------------------
def test_tensor_square_splits():
    """A (x) A = Sym^2 A (+) Lambda^2 A  as multisets."""
    for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4)]:
        assert ms_equal(otimes(A, A), oplus(sym2(A), wedge2(A)))


def test_psi2_is_diagonal_of_sym2():
    """psi^2 A = { lambda_i^2 } is exactly the i=j diagonal inside Sym^2 A."""
    for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4)]:
        diag = tuple(sp.expand(a**2) for a in A)  # the i==j part
        # psi^2 A equals that diagonal ...
        assert ms_equal(psi(2, A), diag)
        # ... and removing the diagonal from Sym^2 leaves the strict upper part.
        L = list(A)
        offdiag = tuple(
            sp.expand(L[i] * L[j]) for i in range(len(L)) for j in range(i + 1, len(L))
        )
        assert ms_equal(oplus(diag, offdiag), sym2(A))


def test_sym2_trace():
    """tr Sym^2 A = (p_1^2 + p_2)/2."""
    for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4)]:
        lhs = sp.expand(sum(sym2(A)))
        p1, p2 = power_sum(A, 1), power_sum(A, 2)
        assert sp.expand(lhs - (p1**2 + p2) / 2) == 0


def test_wedge2_trace_equals_e2():
    """tr Lambda^2 A = (p_1^2 - p_2)/2 = e_2."""
    for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4)]:
        lhs = sp.expand(sum(wedge2(A)))
        p1, p2 = power_sum(A, 1), power_sum(A, 2)
        assert sp.expand(lhs - (p1**2 - p2) / 2) == 0
        assert sp.expand(lhs - elem_sym(A, 2)) == 0


# --------------------------------------------------------------------------
# doubling tower  psi^{2^k} = (psi^2)^{o k}
# --------------------------------------------------------------------------
def test_psi2_iterate_is_power_2k():
    """(psi^2)^{o k} sends lambda -> lambda^{2^k}  [OA-PSI-05]."""
    A = (a1, a2, a3)
    cur = A
    for k in range(1, 5):
        cur = psi(2, cur)  # apply squaring once more
        target = tuple(sp.expand(a ** (2**k)) for a in A)
        assert ms_equal(cur, target)
