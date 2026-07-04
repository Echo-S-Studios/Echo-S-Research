"""
Independent verification of the Mahler-measure (Character I) computations.

Each test rebuilds M(p) = |lead| * prod_{|root|>1}|root| by finding roots at 50
digits (helper `mahler`) and compares to the value the paper states, computed
independently (e.g. phi from sqrt5, 2^7, etc.).
"""

import sys
import pathlib
import mpmath as mp
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _cmc_helpers as H

mp.mp.dps = 50

PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
TOL = mp.mpf(10) ** (-30)


# --------------------------------------------------------------------------
# Realizability family x^n - 2 : M = 2   (Thm 3.2 / ledger A)
# --------------------------------------------------------------------------
def test_xn_minus_2_measure_is_two():
    """Thm 3.2 / ledger A: x^n-2 has M = 2 for n = 3..7."""
    for n in range(3, 8):
        c = [1] + [0] * (n - 1) + [-2]
        assert abs(H.mahler(c) - 2) < TOL


# --------------------------------------------------------------------------
# Even-floor construction q_k = x^{2k}+x^k-1 : M = phi  (Thm 4.1 / ledger E)
# --------------------------------------------------------------------------
def test_qk_measure_is_phi():
    """Thm 4.1 / ledger E: q_k = x^{2k}+x^k-1 has M = phi for k = 2..5."""
    for k in range(2, 6):
        c = [0] * (2 * k + 1)
        c[0] = 1
        c[k] = 1
        c[2 * k] = -1
        assert abs(H.mahler(c) - PHI) < TOL


def test_qk_quadratic_factorization_k2():
    """Thm 4.1: for k=2, x^4+x^2-1 = (x^2+phi)(x^2+phi'); the measure-bearing
    pair is +- i sqrt(phi) (modulus phi^{1/2}), the inert pair +- sqrt(1/phi)."""
    x = sp.symbols("x")
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    lhs = sp.expand(x**4 + x**2 - 1)
    rhs = sp.expand((x**2 + phi) * (x**2 + phic))
    assert sp.simplify(lhs - rhs) == 0
    # x^2 = -phi  -> modulus sqrt(phi) > 1 (measure-bearing, on imaginary axis)
    assert abs(mp.sqrt(PHI) ** 2 - PHI) < TOL
    # measure = (sqrt(phi))^2 from the two conjugate imaginary roots
    assert abs(H.mahler([1, 0, 1, 0, -1]) - PHI) < TOL


# --------------------------------------------------------------------------
# Salem / Lehmer measures  (ledger F, L; Prop 8.2)
# --------------------------------------------------------------------------
def test_beta4_measure():
    """Ledger F: beta_4 = x^4-x^3-x^2-x+1 has M = 1.72208..."""
    assert abs(H.mahler([1, -1, -1, -1, 1]) - mp.mpf("1.72208380573904")) < 1e-12


def test_lehmer_measure():
    """Prop 8.2 / ledger L: Lehmer polynomial has M = tau = 1.17628..."""
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    assert abs(H.mahler(L) - mp.mpf("1.17628081825992")) < 1e-12


def test_commutator_charpoly_measure():
    """Prop 8.2 / ledger L: M(L(x)(x-1)) = M(L)*M(x-1) = tau, since M(x-1)=1."""
    # L(x)(x-1) = x^11 - x^9 - x^8 + x^3 + x^2 - 1
    Lx = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]
    assert abs(H.mahler(Lx) - H.mahler([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])) < TOL
    assert abs(H.mahler([1, -1]) - 1) < TOL            # M(x-1) = 1


# --------------------------------------------------------------------------
# Totally-positive least value  (Thm 4.6 / ledger G)
# --------------------------------------------------------------------------
def test_totally_positive_x1_x2_measure_two():
    """Ledger G: (x-1)(x-2) is totally positive with M = 2."""
    assert abs(H.mahler([1, -3, 2]) - 2) < TOL         # x^2-3x+2


# --------------------------------------------------------------------------
# Pentagon quartic minimiser  (Thm 6.5 / ledger I)
# --------------------------------------------------------------------------
def test_pentagon_quartic_measure_phi4():
    """Ledger I: x^4-x^3+6x^2+4x+1 has M = phi^4 = 6.854..."""
    assert abs(H.mahler([1, -1, 6, 4, 1]) - PHI**4) < TOL


# --------------------------------------------------------------------------
# Reciprocal Z/5 witness  (ledger K)
# --------------------------------------------------------------------------
def test_recip_z5_measure_phi2():
    """Ledger K: Phi_5 * (x^2-3x+1) has M = phi^2 = 2.618..."""
    prod = sp.Poly((sp.Poly([1, 1, 1, 1, 1], sp.symbols("x"))
                    * sp.Poly([1, -3, 1], sp.symbols("x")))).all_coeffs()
    prod = [int(a) for a in prod]
    assert abs(H.mahler(prod) - PHI**2) < TOL


# --------------------------------------------------------------------------
# Pure-power Z/5 family x^5 - m : M = m  (Thm 6.7 / ledger K)
# --------------------------------------------------------------------------
def test_x5_minus_m_measure_is_m():
    """Thm 6.7: x^5 - m has M = m (m>=2); attains mu(5)=2 at x^5-2."""
    for m in range(2, 7):
        c = [1, 0, 0, 0, 0, -m]
        assert abs(H.mahler(c) - m) < TOL


# --------------------------------------------------------------------------
# Tensor / lcm example  (Thm 3.7 / ledger D)
# --------------------------------------------------------------------------
def test_tensor_x3m2_x4m2_measure_128():
    """Ledger D: (x^3-2) (x) (x^4-2) has M = 128 = 2^7."""
    rs = H.tensor_charpoly_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    assert abs(H.mahler_from_roots(rs) - 128) < mp.mpf(10) ** (-20)
    assert 2**7 == 128
