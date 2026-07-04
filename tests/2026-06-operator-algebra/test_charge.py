"""
Character II -- the angle charge.

Theorem 5.2 (thm:charge) [OA-CH-01..05]:
    charge q(lambda) = round(2 arg lambda / pi) mod 4  in  Z/4Z, and
        charge(A (x) B) = charge(A) + charge(B)   (sumset mod 4)
        charge(psi^2 A) = 2 charge(A)
        charge(A (+) B) = charge(A) union charge(B)
    Full group Z/4Z realised: real seeds carry {0,2}; the Lorentzian seed
    K = x^4 + 5x^2 - 5 carries {0,1,2,3}, its imaginary pair +- i beta at +- pi/2.
"""

import sympy as sp

from _opalg_ops import (
    charge,
    charge_one,
    golden_seed,
    K_seed,
    oplus,
    otimes,
    psi,
    seed_from_poly,
    sumset_mod4,
    x,
)

I = sp.I


# --------------------------------------------------------------------------
# the charge of individual roots-of-i*R
# --------------------------------------------------------------------------
def test_charge_of_axis_points():
    """q(+1)=0, q(-1)=2, q(+i)=1, q(-i)=3."""
    assert charge_one(sp.Integer(1)) == 0
    assert charge_one(sp.Integer(-1)) == 2
    assert charge_one(I) == 1
    assert charge_one(-I) == 3


# --------------------------------------------------------------------------
# transformation laws
# --------------------------------------------------------------------------
def test_charge_sumset_on_tensor():
    """charge(A (x) B) = charge(A) + charge(B)  (sumset mod 4)."""
    seeds = [
        golden_seed(),
        seed_from_poly(x**2 - 2),
        K_seed(),
        (I, -I, sp.Integer(2), sp.Integer(-3)),
    ]
    for A in seeds:
        for B in seeds:
            assert charge(otimes(A, B)) == sumset_mod4(charge(A), charge(B))


def test_charge_doubles_under_psi2():
    """charge(psi^2 A) = 2 * charge(A)  (elementwise, mod 4)."""
    seeds = [golden_seed(), K_seed(), (I, -I, sp.Integer(5), sp.Integer(-7))]
    for A in seeds:
        doubled = {}
        cA = charge(A)
        expected = {}
        for q, n in cA.items():
            expected[(2 * q) % 4] = expected.get((2 * q) % 4, 0) + n
        got = dict(charge(psi(2, A)))
        assert got == expected


def test_charge_union_on_oplus():
    """charge(A (+) B) = charge(A) union charge(B)  (multiset union)."""
    A, B = golden_seed(), K_seed()
    assert charge(oplus(A, B)) == charge(A) + charge(B)


# --------------------------------------------------------------------------
# realisation of the full group
# --------------------------------------------------------------------------
def test_real_seeds_carry_0_and_2():
    """Real seeds (real +- eigenvalues) carry exactly charges {0, 2}."""
    for A in [golden_seed(), seed_from_poly(x**2 - 2), seed_from_poly(x**2 - 3)]:
        assert set(charge(A).keys()) == {0, 2}


def test_K_realises_full_Z4():
    """K = x^4 + 5x^2 - 5 carries the full group {0,1,2,3}."""
    assert set(charge(K_seed()).keys()) == {0, 1, 2, 3}


def test_K_imaginary_pair_at_half_pi():
    """K's imaginary pair +- i beta sits at +- pi/2 (charges 1 and 3)."""
    roots = K_seed()
    imag = [r for r in roots if sp.simplify(sp.re(r)) == 0]
    assert len(imag) == 2
    assert sorted(charge_one(r) for r in imag) == [1, 3]


def test_charge_is_Z4_homomorphism_closure():
    """The doubling law is the charge-shadow of M(psi^2 A)=M(A)^2: 2*Z/4Z = {0,2}."""
    # Applying psi^2 twice collapses every charge to 0 (since 4q = 0 mod 4).
    A = K_seed()
    assert set(charge(psi(2, psi(2, A))).keys()) == {0}
