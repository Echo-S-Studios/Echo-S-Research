"""Part 7 -- The cyclotomic floor: where magnitude becomes charge (Sec. 7).

Independent re-derivations of the Kronecker floor behaviour, the five boundary
conditions (B1-B5), the mu_4 = Z/4Z collapse, the golden isolation gap, and
Lehmer's number (~1.17628) computed from its polynomial at high precision.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp
from _eap_helpers import (phi, psi, sqrt5, is_zero, mahler, charge_one, charge_set)


# ----------------------------------------------------------------------------
# Thm 7.1 (Kronecker) -- illustrated on cyclotomic vs non-cyclotomic
# ----------------------------------------------------------------------------
def test_kronecker_cyclotomic_have_mahler_one():
    """Thm 7.1: for a monic integer polynomial M=1 iff all roots are roots of
    unity.  Cyclotomic polynomials -> M=1."""
    x = sp.symbols('x')
    for poly in [x**2 + 1, x**2 + x + 1, x - 1, x**4 + x**3 + x**2 + x + 1,
                 x**4 - 1, x**6 + x**5 + x**4 + x**3 + x**2 + x + 1]:
        roots = sp.Poly(poly, x).all_roots()
        assert mahler(roots) == 1


def test_kronecker_noncyclotomic_exceed_one():
    """Thm 7.1: a non-cyclotomic monic integer polynomial has M>1
    (e.g. x^2-x-1 -> phi, and Lehmer's polynomial)."""
    x = sp.symbols('x')
    roots = sp.Poly(x**2 - x - 1, x).all_roots()
    assert mahler(roots) > 1


# ----------------------------------------------------------------------------
# Thm 7.2 -- the five boundary conditions
# ----------------------------------------------------------------------------
def test_B1_floor_and_identity():
    """(B1): M>=1 with equality on the cyclotomic locus; M=1 is the identity
    of the magnitude monoid ([1,inf), x)."""
    # M(A_phi)=phi>=1, M(mu_4 obj)=1
    assert mahler([phi, psi]) >= 1
    assert mahler([sp.I, -sp.I]) == 1
    # identity: 1 * m = m
    m = mahler([phi, psi])
    assert m * mp.mpf(1) == m


def test_B2_sub_semiring_closed():
    """(B2): A_1 closed under +, x, psi^n -- products, unions and Adams powers
    of roots of unity stay roots of unity (M stays 1)."""
    units = [sp.I, -sp.I, sp.Integer(-1), sp.Integer(1),
             sp.exp(2 * sp.pi * sp.I / 3)]
    # union
    assert mahler(units) == 1
    # products
    prods = [a * b for a in units for b in units]
    assert mahler(prods) == 1
    # Adams powers
    for n in range(1, 6):
        pows = [u**n for u in units]
        assert mahler(pows) == 1


def test_B3_collapse_to_charge_group():
    """(B3): on A_1, chi:mu_4 -> Z/4Z with x -> + and psi^n -> xn (mod4);
    M is constant 1."""
    mu4 = [sp.Integer(1), sp.I, sp.Integer(-1), -sp.I]
    # isomorphism chi(i^k)=k
    for k in range(4):
        assert charge_one(sp.I**k) == k % 4
    # x -> + : chi(i^a i^b) = (a+b) mod 4
    for a in range(4):
        for b in range(4):
            assert charge_one(sp.I**a * sp.I**b) == (a + b) % 4
    # psi^n -> xn : chi((i^a)^n) = n a mod 4
    for a in range(4):
        for n in range(1, 5):
            assert charge_one((sp.I**a)**n) == (n * a) % 4
    assert mahler(mu4) == 1


def test_B4_golden_isolation_gap():
    """(B4): in A's emission class M in {1} U [phi, inf); the floor is isolated
    by a gap of width >= phi-1 (~0.618), with A_phi realising M=phi (the next
    value).  Lehmer's number ~1.176 sits inside (1,phi) but is OUTSIDE A."""
    assert is_zero((phi - 1) - 1 / phi)                 # gap width phi-1 = 1/phi
    assert abs(float(phi - 1) - 0.618) < 0.001
    # golden object realises phi, the next emitted value above 1
    assert abs(float(mahler([phi, psi])) - float(phi)) < 1e-12
    # Lehmer number is in the open gap (1,phi) -- hence not an A-emission value
    lehmer = _lehmer_mahler()
    assert 1 < lehmer < float(phi)


def test_B5_one_way_nondecreasing():
    """(B5): x and + are magnitude-nondecreasing (each factor >=1), so the
    floor is a closed superselection sector -- never left downward, never
    reached from M>1."""
    A = [phi, psi]                   # M=phi
    B = [phi**2, psi**2]             # M=phi^2
    # union: M(A+B) = M(A)M(B) >= max(M(A),M(B))
    assert mahler(A + B) >= max(mahler(A), mahler(B))
    # coupling stays >= 1 and cannot fall below a factor
    coupled = [a * b for a in A for b in B]
    assert mahler(coupled) >= 1


# ----------------------------------------------------------------------------
# Ex 7.4 -- the floor, the gap, and beyond
# ----------------------------------------------------------------------------
def test_floor_examples_phi4_phi3():
    """Ex 7.4: Phi_4=x^2+1 -> {i,-i}, M=1, chi={1,3}; Phi_3=x^2+x+1 ->
    e^{+-2pi i/3}, M=1, chi={1,3}; A_phi has M=phi (first rung)."""
    x = sp.symbols('x')
    r4 = sp.Poly(x**2 + 1, x).all_roots()
    assert mahler(r4) == 1
    assert charge_set(r4) == {1, 3}
    r3 = sp.Poly(x**2 + x + 1, x).all_roots()
    assert mahler(r3) == 1
    assert charge_set(r3) == {1, 3}
    assert abs(float(mahler([phi, psi])) - float(phi)) < 1e-12


def _lehmer_mahler():
    """Independently compute the Mahler measure of Lehmer's degree-10
    polynomial x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1 at high precision."""
    mp.mp.dps = 50
    coeffs = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    roots = mp.polyroots(coeffs, maxsteps=500, extraprec=300)
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    return float(M)


def test_lehmer_number_value():
    """Ex 7.4: Lehmer's polynomial has M ~ 1.17628, inside the gap (1,phi)."""
    L = _lehmer_mahler()
    assert abs(L - 1.17628) < 1e-4
    assert 1 < L < float(phi)


def test_lehmer_polynomial_not_cyclotomic():
    """Ex 7.4: Lehmer's polynomial is non-cyclotomic (M>1), so its roots are
    not all roots of unity (Kronecker), consistent with it lying above 1."""
    assert _lehmer_mahler() > 1


# ----------------------------------------------------------------------------
# Ex 7.5 -- the floor-collapse orbit
# ----------------------------------------------------------------------------
def test_exercise_floor_collapse_orbit():
    """Ex 7.5: from {i,-i} (charge {1,3}), psi^2 -> {-1,-1} (M=1, charge {2}),
    matching 2*{1,3}={2} mod4; psi^2 again -> {1} (charge {0}); M never moves
    off 1."""
    start = [sp.I, -sp.I]
    assert charge_set(start) == {1, 3}
    assert mahler(start) == 1
    step1 = [sp.I**2, (-sp.I)**2]     # {-1,-1}
    assert all(is_zero(v - (-1)) for v in step1)
    assert mahler(step1) == 1
    assert charge_set(step1) == {2}
    assert {(2 * c) % 4 for c in charge_set(start)} == {2}
    step2 = [v**2 for v in step1]     # {1,1}
    assert all(is_zero(v - 1) for v in step2)
    assert mahler(step2) == 1
    assert charge_set(step2) == {0}
