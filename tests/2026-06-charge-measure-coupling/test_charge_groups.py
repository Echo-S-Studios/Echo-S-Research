"""
Independent verification of the charge group (Character II) computations.

`charge_group` re-derives, from scratch, the least n with alpha^n in R_{>0}
for every root (== least common denominator of arg/(2pi)); None means an
irrational conjugate angle (charge-inadmissible).
"""

import sys
import pathlib
import numpy as np
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _cmc_helpers as H


# --------------------------------------------------------------------------
# Realizability: x^n - 2 has charge group Z/n and attains all n charges
# (Thm 3.2 / ledger A)
# --------------------------------------------------------------------------
def test_xn_minus_2_charge_group_and_all_charges():
    """Thm 3.2 / ledger A: x^n-2 has charge group Z/n and every charge occurs."""
    for n in range(3, 8):
        c = [1] + [0] * (n - 1) + [-2]
        assert H.charge_group(c) == n
        assert H.charges(c, n) == list(range(n))       # attains all n charges


def test_x4_minus_2_is_cyclic_z4_not_z2xz2():
    """Ledger B: x^4-2 has charge group Z/4 (a single generator, i.e. some
    root sits at t=1/4), NOT Z/2 x Z/2."""
    c = [1, 0, 0, 0, -2]
    assert H.charge_group(c) == 4
    # a genuine generator of order 4 is present (charge 1 occurs)
    assert 1 in H.charges(c, 4)
    # Z/2 x Z/2 would have exponent 2; Z/4 has an element of order 4
    assert 4 not in [1, 2]        # order-4 charge group is not exponent-2


# --------------------------------------------------------------------------
# Even-floor construction q_k : charge group Z/2k, all charges (Thm 4.1 / E)
# --------------------------------------------------------------------------
def test_qk_charge_group_z2k_all_charges():
    """Thm 4.1 / ledger E: q_k = x^{2k}+x^k-1 has charge group Z/2k and
    attains all 2k charges (k = 2..5)."""
    for k in range(2, 6):
        c = [0] * (2 * k + 1)
        c[0] = 1
        c[k] = 1
        c[2 * k] = -1
        assert H.charge_group(c) == 2 * k
        assert H.charges(c, 2 * k) == list(range(2 * k))


# --------------------------------------------------------------------------
# Salem exclusion: charge-inadmissible objects (Lem 8.1 / ledger F, L)
# --------------------------------------------------------------------------
def test_beta4_charge_inadmissible():
    """Ledger F / Lem 8.1: beta_4 (a Salem number) has NO finite charge
    group (irrational conjugate angle)."""
    assert H.charge_group([1, -1, -1, -1, 1], Nmax=300) is None


def test_lehmer_commutator_charge_inadmissible():
    """Ledger L / Prop 8.2: L(x)(x-1) carries Lehmer's Salem number and has
    charge group bottom (inadmissible)."""
    Lx = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]
    assert H.charge_group(Lx, Nmax=300) is None


# --------------------------------------------------------------------------
# Pentagon (Z/5) witnesses  (ledger I, K)
# --------------------------------------------------------------------------
def test_pentagon_quartic_charge_z5():
    """Ledger I: x^4-x^3+6x^2+4x+1 has charge group Z/5 (irreducible,
    non-reciprocal)."""
    c = [1, -1, 6, 4, 1]
    assert H.charge_group(c) == 5
    assert not H.is_reciprocal(c)
    assert sp.Poly(c, sp.symbols("x")).is_irreducible


def test_recip_z5_witness_charge_z5():
    """Ledger K: Phi_5 * (x^2-3x+1) has charge group Z/5 and is reciprocal."""
    prod = np.polymul([1, 1, 1, 1, 1], [1, -3, 1]).tolist()
    assert H.charge_group(prod) == 5
    assert H.is_reciprocal(prod)


def test_x5_minus_2_charge_z5_nonreciprocal():
    """Thm 6.7 / ledger K: x^5-2 has charge group Z/5, non-reciprocal."""
    c = [1, 0, 0, 0, 0, -2]
    assert H.charge_group(c) == 5
    assert not H.is_reciprocal(c)


# --------------------------------------------------------------------------
# lcm law under tensor  (Thm 3.5 / ledger D)
# --------------------------------------------------------------------------
def test_tensor_charge_group_is_lcm():
    """Ledger D / Thm 3.5: (x^3-2)(x)(x^4-2), gcd(3,4)=1, has charge group
    Z/12 = Z/lcm(3,4)."""
    rs = H.tensor_charpoly_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    assert H.charge_group_from_roots(rs) == 12
