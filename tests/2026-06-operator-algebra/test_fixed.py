"""
Section 7 -- Fixed points and the spectrum-preserving operators.

Proposition 7.1 (prop:fixed) [OA-FX-01]:
    M(psi^2 A) = M(A)  iff  M(A) = 1, i.e. every eigenvalue is a root of unity
    or 0.  Proof: M(A)^2 = M(A) forces M(A) in {0,1}; as M >= 1, M(A) = 1.

[OA-FX-02]:  minpoly and Phi re-present an object without changing its
    eigenvalue multiset -- idempotent, spectrum- and character-preserving.
"""

import mpmath as mp
import sympy as sp

from _opalg_ops import (
    charge,
    golden_seed,
    mahler_exact,
    psi,
    seed_from_poly,
    x,
)


# --------------------------------------------------------------------------
# the fixed-point algebra:  M^2 = M  <=>  M in {0,1}
# --------------------------------------------------------------------------
def test_only_fixed_measure_is_one():
    """M(A)^2 = M(A) has the sole solution M(A) = 1 in the range M >= 1."""
    M = sp.Symbol("M")  # no positivity assumption, so both roots appear
    sols = sp.solve(sp.Eq(M**2, M), M)  # {0, 1}
    assert set(sols) == {0, 1}
    feasible = [s for s in sols if s >= 1]  # M >= 1 always for a Mahler measure
    assert feasible == [1]


def test_mahler_at_least_one():
    """M(A) >= 1 always (product of factors max(1,|.|) >= 1)."""
    for A in [golden_seed(), seed_from_poly(x**2 - 2), seed_from_poly(x**4 + 5 * x**2 - 5)]:
        assert mahler_exact(A) >= 1


# --------------------------------------------------------------------------
# cyclotomic / Kronecker locus has M = 1 and is Adams-fixed
# --------------------------------------------------------------------------
def test_cyclotomic_objects_are_adams_fixed():
    """Roots of unity (M=1) satisfy M(psi^2 A) = M(A); non-cyclotomic do not."""
    # roots of unity: cyclotomic Phi_n
    for poly in [x**2 + 1, x**2 + x + 1, x**4 + 1, x**6 + x**3 + 1, x - 1, x + 1]:
        A = seed_from_poly(poly)
        assert sp.simplify(mahler_exact(A) - 1) == 0
        assert sp.simplify(mahler_exact(psi(2, A)) - mahler_exact(A)) == 0


def test_noncyclotomic_not_adams_fixed():
    """A non-Kronecker object (golden, M=phi>1) is NOT fixed by squaring."""
    A = golden_seed()
    assert sp.simplify(mahler_exact(A) - 1) != 0
    assert sp.simplify(mahler_exact(psi(2, A)) - mahler_exact(A)) != 0
    # instead it squares:
    assert sp.simplify(mahler_exact(psi(2, A)) - mahler_exact(A) ** 2) == 0


def test_zero_eigenvalue_excluded_but_M_one_boundary():
    """The M=1 locus is exactly {all |lambda|<=1}; a unit-circle object has M=1."""
    # eighth roots of unity via Phi_8 = x^4 + 1: all on unit circle
    A = seed_from_poly(x**4 + 1)
    mp.mp.dps = 40
    from _opalg_ops import modulus_mp

    assert all(abs(modulus_mp(a, 40) - 1) < mp.mpf(10) ** -30 for a in A)
    assert sp.simplify(mahler_exact(A) - 1) == 0


# --------------------------------------------------------------------------
# spectrum-preserving idempotents: minpoly / Phi
# --------------------------------------------------------------------------
def test_minpoly_reconstruction_preserves_spectrum():
    """Rebuilding roots from an object's minimal polynomial returns the same
    eigenvalue multiset (idempotent on the spectrum) -> both characters fixed."""
    # golden object: minimal polynomial of phi is x^2 - x - 1; its roots are the
    # golden seed again.
    A = golden_seed()
    reconstructed = seed_from_poly(x**2 - x - 1)
    from _opalg_ops import ms_equal

    assert ms_equal(A, reconstructed)
    # characters unchanged
    assert sp.simplify(mahler_exact(A) - mahler_exact(reconstructed)) == 0
    assert charge(A) == charge(reconstructed)
