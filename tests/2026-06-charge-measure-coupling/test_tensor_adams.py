"""
Independent verification of the operator table (Sec 2.1 backbone) and the
tensor / Adams computations: the tropical tensor-Mahler law (Prop 3.3), the
Adams measure law M(psi^k A)=M(A)^k, the sumset charge rule, and the CRT /
Adams primary decomposition (Thm 3.4 / ledger C, M, D).
"""

import sys
import pathlib
import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _cmc_helpers as H

mp.mp.dps = 50
TOL = mp.mpf(10) ** (-25)
x = sp.symbols("x")


def _monic_from_roots(rs):
    """Build monic polynomial coeffs (highest-first) from numeric roots."""
    p = np.poly([complex(mp.re(r), mp.im(r)) for r in rs])
    return p


# --------------------------------------------------------------------------
# Tensor eigenvalues are pairwise products; tropical (max,+) Mahler law
# (Prop 3.3 / ledger M)
# --------------------------------------------------------------------------
def test_phi_tensor_phi_charpoly_symbolic():
    """Ledger M: (x^2-x-1) (x) (x^2-x-1) has charpoly (x+1)^2(x^2-3x+1).
    Eigenvalues are the products {phi^2, -1, -1, phi'^2}."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    prods = [phi * phi, phi * phic, phic * phi, phic * phic]  # phi^2,-1,-1,phi'^2
    poly = sp.prod([x - p for p in prods])
    expected = (x + 1) ** 2 * (x**2 - 3 * x + 1)
    assert sp.simplify(sp.expand(poly) - sp.expand(expected)) == 0


def test_phi_tensor_phi_tropical_not_factored():
    """Prop 3.3: M(phi (x) phi) = phi^2 (tropical), NOT M(phi)^2 M(phi)^2 = phi^4.
    The pair -1,-1 lies ON the unit circle, so the factored law overcounts."""
    PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
    rs = H.tensor_charpoly_roots([1, -1, -1], [1, -1, -1])
    m = H.mahler_from_roots(rs)
    assert abs(m - PHI**2) < TOL                     # tropical value
    assert abs(m - PHI**4) > mp.mpf("0.1")           # NOT the factored value phi^4
    # numeric charpoly root-set matches (x+1)^2(x^2-3x+1)
    got_roots = np.sort_complex(np.roots(_monic_from_roots(rs)))
    exp_roots = np.sort_complex(np.roots(np.polymul([1, 2, 1], [1, -3, 1])))
    assert np.allclose(got_roots, exp_roots, atol=1e-9)


def test_tensor_offcircle_factored_agrees():
    """Thm 3.7 note: for (x^3-2)(x)(x^4-2) every product is off-circle, so the
    tropical law reduces to (2^{7/12})^{12} = 2^7 = 128 (factored agrees)."""
    rs = H.tensor_charpoly_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    assert abs(H.mahler_from_roots(rs) - 128) < mp.mpf(10) ** (-18)
    # all twelve products have modulus 2^{7/12} > 1
    for r in rs:
        assert abs(r) > 1


def test_tensor_charpoly_is_x12_minus_128():
    """Ledger D: charpoly of (x^3-2)(x)(x^4-2) is x^12 - 128."""
    rs = H.tensor_charpoly_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    coeffs = _monic_from_roots(rs)
    target = np.array([1] + [0] * 11 + [-128], dtype=complex)
    assert np.allclose(coeffs, target, atol=1e-6)


# --------------------------------------------------------------------------
# Adams measure law  M(psi^k A) = M(A)^k   (Sec 2.1 table)
# --------------------------------------------------------------------------
def test_adams_measure_law_power():
    """Sec 2.1 table: psi^k raises every modulus to the k-th power, so
    M(psi^k A) = M(A)^k.  Checked on x^n-2 (M=2 -> 2^k)."""
    for n in (3, 5):
        base = [1] + [0] * (n - 1) + [-2]
        rs = H.roots_mp(base)
        for k in (2, 3):
            powered = [r**k for r in rs]
            assert abs(H.mahler_from_roots(powered) - 2**k) < TOL


# --------------------------------------------------------------------------
# CRT / Adams primary decomposition realised on x^6-2  (Thm 3.4 / ledger C)
# --------------------------------------------------------------------------
def test_adams_primary_projection_on_x6_minus_2():
    """Thm 3.4 / ledger C: x^6-2 has charge Z/6; psi^3 leaves residual charge
    group Z/2, psi^2 leaves Z/3 (Adams operations realise the CRT projectors)."""
    rs = H.roots_mp([1, 0, 0, 0, 0, 0, -2])           # x^6 - 2, charge Z/6
    assert H.charge_group_from_roots(rs) == 6
    # psi^3 : cube the roots -> residual charge group Z/2
    cubed = [r**3 for r in rs]
    assert H.charge_group_from_roots(cubed) == 2
    # psi^2 : square the roots -> residual charge group Z/3
    squared = [r**2 for r in rs]
    assert H.charge_group_from_roots(squared) == 3


# --------------------------------------------------------------------------
# Sumset charge rule under tensor  (Sec 2.1 table / Thm 3.5)
# --------------------------------------------------------------------------
def test_tensor_charges_are_sumset():
    """Sec 2.1 table: charges add under tensor.  For (x^3-2)(x)(x^4-2) the
    composite charges in Z/12 equal { (4a+3b) mod 12 } = sumset(Z/3, Z/4)."""
    rs = H.tensor_charpoly_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    got = set(_charges_from_roots(rs, 12))
    # independent sumset: root of A at 2pi a/3, root of B at 2pi b/4
    #   -> t = a/3 + b/4 = (4a+3b)/12
    expected = {(4 * a + 3 * b) % 12 for a in range(3) for b in range(4)}
    assert got == expected == set(range(12))          # coprime -> full Z/12


def _charges_from_roots(rs, n):
    two_pi = 2 * mp.pi
    return [int(mp.nint(n * (mp.arg(r) / two_pi))) % n for r in rs]
