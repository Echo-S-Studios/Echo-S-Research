"""
Independent verification of the Salem exclusion (Lem 8.1) and the commutator
escape (Prop 8.2 / ledger F, L).

Computational content only: the cited matrix theorems (Kronecker, Smyth,
Shoda, Albert-Muckenhoupt, Laffey-Reams) are taken as external inputs; here we
check the concrete witnesses -- that beta_4 and Lehmer's number are Salem-type
and charge-inadmissible, and that the Lehmer commutator is a trace-zero integer
matrix carrying Lehmer's number below the emission floor.
"""

import sys
import pathlib
import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _cmc_helpers as H

mp.mp.dps = 50
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
x = sp.symbols("x")


def _salem_structure(coeffs):
    """Return (n_real_outside, n_on_circle, n_real_inside) for the roots."""
    rs = H.roots_mp(coeffs)
    outside = [r for r in rs if abs(r) > 1 + mp.mpf(10) ** (-20)]
    inside = [r for r in rs if abs(r) < 1 - mp.mpf(10) ** (-20)]
    oncirc = [r for r in rs if abs(abs(r) - 1) <= mp.mpf(10) ** (-20)]
    return outside, oncirc, inside


# --------------------------------------------------------------------------
# Salem exclusion  (Lem 8.1 / ledger F)
# --------------------------------------------------------------------------
def test_beta4_is_salem_and_inadmissible():
    """Lem 8.1 / ledger F: beta_4 = x^4-x^3-x^2-x+1 is a Salem number
    (reciprocal; one real root >1, its inverse, two conjugates on |z|=1 at an
    irrational angle) and therefore charge-inadmissible."""
    c = [1, -1, -1, -1, 1]
    assert H.is_reciprocal(c)
    outside, oncirc, inside = _salem_structure(c)
    assert len(outside) == 1                 # the Salem number tau
    assert len(inside) == 1                  # its inverse 1/tau (real, in (0,1))
    assert len(oncirc) == 2                  # a conjugate pair on the unit circle
    assert mp.im(inside[0]) == 0 and 0 < mp.re(inside[0]) < 1
    assert H.charge_group(c, Nmax=400) is None       # irrational angle => bottom
    # Mahler equals the Salem number itself (only one root outside)
    assert abs(H.mahler(c) - abs(outside[0])) < mp.mpf(10) ** (-30)


def test_lehmer_is_salem_and_inadmissible():
    """Lem 8.1 / Prop 8.2: Lehmer's number is a Salem number of degree 10
    (one root >1, one inside, eight on the unit circle) and inadmissible."""
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    assert H.is_reciprocal(L)
    outside, oncirc, inside = _salem_structure(L)
    assert len(outside) == 1                 # Lehmer's number tau
    assert len(inside) == 1                  # 1/tau
    assert len(oncirc) == 8                  # eight conjugates on the circle
    assert H.charge_group(L, Nmax=400) is None


# --------------------------------------------------------------------------
# Commutator escape  (Prop 8.2 / ledger L)
# --------------------------------------------------------------------------
def test_lehmer_trace_minus_one():
    """Prop 8.2: Lehmer's polynomial L has trace -1 (coeff of x^9 is +1)."""
    L = sp.Poly(x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1, x)
    # sum of roots = -(coeff x^9)/(coeff x^10)
    assert -L.all_coeffs()[1] == -1


def test_commutator_polynomial_is_trace_zero():
    """Prop 8.2 / ledger L: L(x)*(x-1) = x^11-x^9-x^8+x^3+x^2-1 has trace 0
    (coeff of x^10 vanishes), so its companion matrix is a trace-zero integer
    matrix -- a commutator by Albert-Muckenhoupt / Laffey-Reams."""
    L = sp.Poly(x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1, x)
    prod = sp.expand(L.as_expr() * (x - 1))
    target = x**11 - x**9 - x**8 + x**3 + x**2 - 1
    assert sp.simplify(prod - target) == 0
    coeffs = sp.Poly(prod, x).all_coeffs()
    assert coeffs[1] == 0                    # coeff of x^10 == 0  => trace 0


def test_commutator_companion_matrix_trace_zero_carries_lehmer():
    """Prop 8.2: the companion matrix of L(x)(x-1) is an integer matrix with
    trace 0 whose spectrum contains Lehmer's number tau (and 1)."""
    # companion of monic x^11 - x^9 - x^8 + x^3 + x^2 - 1
    p = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]
    n = len(p) - 1
    C = np.zeros((n, n))
    C[0, :] = [-p[k] for k in range(1, n + 1)]        # top-row companion form
    for i in range(1, n):
        C[i, i - 1] = 1
    assert abs(np.trace(C)) < 1e-9                     # trace zero (integer)
    assert np.allclose(C, np.round(C))                 # integer matrix
    tau = float(H.mahler([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]))
    eig = np.linalg.eigvals(C)
    assert min(abs(eig - tau)) < 1e-6                  # Lehmer's number is eigenvalue
    assert min(abs(eig - 1.0)) < 1e-6                  # and 1 is an eigenvalue


def test_commutator_below_floor_and_uncharged():
    """Prop 8.2 (i),(iii): the commutator spectrum has M = tau = 1.17628 in
    (1, phi) -- below the emission floor -- and charge group bottom."""
    Lx = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]
    m = H.mahler(Lx)
    assert 1 < m < PHI                                 # strictly below the floor phi
    assert abs(m - mp.mpf("1.17628")) < 5e-6
    assert H.charge_group(Lx, Nmax=400) is None        # off every angle lattice
