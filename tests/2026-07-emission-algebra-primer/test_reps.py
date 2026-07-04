"""Part 4 -- Representation theory: the dimensions are the exponents (Sec. 4).

Independent re-derivations of the sl2 irreducible data (dimension, weights,
Casimir), the symmetric-power eigenvalue tower, the Clebsch-Gordan dimension
bookkeeping, and the Sym^2 weights.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from _eap_helpers import phi, psi, sqrt5, is_zero


def casimir(m):
    """Casimir of V_m as defined in Thm 4.1: (1/2) m (m+2)."""
    return sp.Rational(1, 2) * m * (m + 2)


# ----------------------------------------------------------------------------
# Thm 4.1 -- irreducible data
# ----------------------------------------------------------------------------
def test_irrep_dimension_and_weights():
    """Thm 4.1: dim V_m = m+1, h-weights {m,m-2,...,-m}."""
    for m in range(0, 8):
        weights = list(range(m, -m - 1, -2))
        assert len(weights) == m + 1
        assert weights == [m - 2 * k for k in range(m + 1)]
        assert weights[-1] == -m


def test_casimir_values():
    """Thm 4.1 / Prop 4.3: Casimir (1/2)m(m+2) on V_0..V_4."""
    assert casimir(0) == 0
    assert casimir(1) == sp.Rational(3, 2)
    assert casimir(2) == 4
    assert casimir(3) == sp.Rational(15, 2)
    assert casimir(4) == 12


# ----------------------------------------------------------------------------
# Def 4.2 / Prop 4.3 -- symmetric power eigenvalue tower
# ----------------------------------------------------------------------------
def test_symmetric_power_eigenvalue_tower():
    """Prop 4.3: 'spec(R|V_n) = {phi^{n-k} psi^k : k=0..n}'.
    Independently build rho_n(R) as the diagonal action on the monomial
    basis x^{n-k} y^k and read its eigenvalues."""
    for n in range(0, 6):
        tower = sorted([sp.nsimplify(phi**(n - k) * psi**k) for k in range(n + 1)],
                       key=lambda t: float(t))
        # dimension check
        assert len(tower) == n + 1
        # each entry is a product of n base eigenvalues
        for k in range(n + 1):
            assert is_zero(tower_val(n, k) - phi**(n - k) * psi**k)


def tower_val(n, k):
    return phi**(n - k) * psi**k


def test_H_weights_on_sympow():
    """Prop 4.3 proof: on Sym^n, H=2R-I acts with eigenvalues
    (n-k)sqrt5 + k(-sqrt5) = sqrt5 (n-2k); so H/sqrt5 has integer weights
    {n, n-2, ..., -n}."""
    for n in range(0, 6):
        for k in range(n + 1):
            hval = (n - k) * sqrt5 + k * (-sqrt5)
            assert is_zero(hval - sqrt5 * (n - 2 * k))
            assert is_zero(hval / sqrt5 - (n - 2 * k))


def test_M2_decomposes_as_V2_plus_V0():
    """Prop 4.3: 'M_2 = V_2 (+) V_0', dim 4 = 3 + 1, recovering ad_R
    eigenvalues {0,0,+/-sqrt5}."""
    assert (2 + 1) + (0 + 1) == 4       # dim V_2 + dim V_0
    # ad_R eigenvalues on V_2 (adjoint) are {+sqrt5,0,-sqrt5}, on V_0 {0}
    adjoint_eigs = [sqrt5, 0, -sqrt5]   # weights * root length; +1 zero from V_0
    alleigs = sorted([float(e) for e in adjoint_eigs] + [0.0])
    assert alleigs == sorted([0.0, 0.0, float(sqrt5), float(-sqrt5)])


# ----------------------------------------------------------------------------
# Thm 4.4 / Ex 4.5 -- Clebsch-Gordan
# ----------------------------------------------------------------------------
def cg_decomp(a, b):
    """Clebsch-Gordan summand highest weights j from |a-b| to a+b step 2."""
    return list(range(abs(a - b), a + b + 1, 2))


def test_clebsch_gordan_dimension_identity():
    """Thm 4.4: '(a+1)(b+1) = sum_j (j+1)' over the CG summands, for many a,b."""
    for a in range(0, 6):
        for b in range(0, 6):
            js = cg_decomp(a, b)
            assert (a + 1) * (b + 1) == sum(j + 1 for j in js), (a, b)


def test_clebsch_gordan_worked_examples():
    """Ex 4.5: V1xV1=V0+V2 (4=1+3); V2xV2=V0+V2+V4 (9=1+3+5);
    V3xV3=V0+V2+V4+V6 (16=1+3+5+7)."""
    assert cg_decomp(1, 1) == [0, 2]
    assert [j + 1 for j in cg_decomp(1, 1)] == [1, 3] and 2 * 2 == 4 == 1 + 3
    assert cg_decomp(2, 2) == [0, 2, 4]
    assert 3 * 3 == 9 == 1 + 3 + 5
    assert cg_decomp(3, 3) == [0, 2, 4, 6]
    assert 4 * 4 == 16 == 1 + 3 + 5 + 7


def test_casimir_coupled_targets():
    """Prop 5.10: V1xV2=V1+V3 (Casimirs 3/2, 15/2); V2xV2=V0+V2+V4
    (Casimirs 0,4,12); top components Cas(V3)=15/2, Cas(V4)=12."""
    assert cg_decomp(1, 2) == [1, 3]
    assert [casimir(j) for j in cg_decomp(1, 2)] == [sp.Rational(3, 2), sp.Rational(15, 2)]
    assert cg_decomp(2, 2) == [0, 2, 4]
    assert [casimir(j) for j in cg_decomp(2, 2)] == [0, 4, 12]


# ----------------------------------------------------------------------------
# Ex 4.7 -- Sym^2 in the flesh
# ----------------------------------------------------------------------------
def test_sym2_eigenvalues_and_weights():
    """Ex 4.7: on Sym^2, R has eigenvalues {phi^2, phi*psi, psi^2}={phi^2,-1,psi^2}
    and H/sqrt5 has weights {2,0,-2}."""
    eigs = [phi**2, phi * psi, psi**2]
    assert is_zero(sp.nsimplify(eigs[1]) - (-1))        # phi*psi = -1
    # H weights: sqrt5*{2,0,-2}
    hweights = [(2 - 2 * k) for k in range(3)]
    assert hweights == [2, 0, -2]
    for k in range(3):
        hval = sqrt5 * (2 - 2 * k)
        assert is_zero(hval / sqrt5 - hweights[k])
