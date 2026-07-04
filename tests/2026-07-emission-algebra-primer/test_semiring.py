"""Part 3 -- The semiring layer: magnitude and phase (Sec. 3).

Independent re-derivations of the two characters (Mahler measure & charge) on
the worked objects, the grading laws, the tropical coupling subtlety, and the
Adams operator's pullback to the core.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
import mpmath as mp
from _eap_helpers import (phi, psi, sqrt5, is_zero, mahler, mahler_exact,
                      charge_one, charge_set, outside_unit)


# ----------------------------------------------------------------------------
# Ex 3.3 -- the golden object
# ----------------------------------------------------------------------------
def test_golden_object_magnitude():
    """Ex 3.3: 'M(A_phi)=phi' (only phi>1 contributes; |psi|=phi^{-1}<1)."""
    out = outside_unit([phi, psi])
    # exactly one eigenvalue outside the unit circle, and it is phi
    assert len(out) == 1
    assert is_zero(mahler_exact([phi, psi]) - phi)


def test_golden_object_charge():
    """Ex 3.3: 'chi(A_phi) = {0,2}' -- phi at angle 0, psi at angle pi."""
    assert charge_set([phi, psi]) == {0, 2}


# ----------------------------------------------------------------------------
# Ex 3.4 -- roots of unity and mu_4
# ----------------------------------------------------------------------------
def test_roots_of_unity_weightless():
    """Ex 3.4: 'M({zeta})=1' for any root of unity (here mu_4 and a cube root)."""
    for z in [sp.Integer(1), sp.I, sp.Integer(-1), -sp.I,
              sp.exp(2 * sp.pi * sp.I / 3)]:
        assert mahler([z]) == mp.mpf(1)


def test_mu4_charges_exact():
    """Ex 3.4: 'chi(1)=0, chi(i)=1, chi(-1)=2, chi(-i)=3' -> mu_4 = Z/4Z."""
    assert charge_one(sp.Integer(1)) == 0
    assert charge_one(sp.I) == 1
    assert charge_one(sp.Integer(-1)) == 2
    assert charge_one(-sp.I) == 3
    assert charge_set([sp.Integer(1), sp.I, sp.Integer(-1), -sp.I]) == {0, 1, 2, 3}


def test_cube_root_rounds_to_charge_1():
    """Ex 3.4: 'a primitive cube root e^{2pi i/3} ... rounds to charge 1'."""
    assert charge_one(sp.exp(2 * sp.pi * sp.I / 3)) == 1


# ----------------------------------------------------------------------------
# Prop 3.6 -- the grading laws (magnitude & charge under +, x, psi^n)
# ----------------------------------------------------------------------------
def test_superposition_multiplies_magnitude():
    """Prop 3.6: 'M(A + B) = M(A) M(B)' (union multiplies the measure)."""
    A = [phi**2, psi**2]
    B = [phi, psi]
    MA = mahler_exact(A)            # phi^2
    MB = mahler_exact(B)            # phi
    MAB = mahler_exact(A + B)       # multiset union -> phi^3
    assert is_zero(MAB - MA * MB)


def test_superposition_unions_charge():
    """Prop 3.6: 'chi(A + B) = chi(A) U chi(B)'."""
    A = [phi, psi]                   # charges {0,2}
    B = [sp.I, -sp.I]                # charges {1,3}
    assert charge_set(A + B) == charge_set(A) | charge_set(B) == {0, 1, 2, 3}


def test_coupling_charge_adds_mod4():
    """Prop 3.6: 'chi(A x B) = chi(A) + chi(B) (mod 4)' (sumset mod 4)."""
    A = [sp.I]                       # charge {1}
    B = [sp.I]                       # charge {1}
    prod = [a * b for a in A for b in B]     # {i^2 = -1} -> charge {2}
    got = charge_set(prod)
    want = {(a + b) % 4 for a in charge_set(A) for b in charge_set(B)}
    assert got == want == {2}


def test_adams_magnitude_and_charge():
    """Prop 3.6: 'M(psi^n A)=M(A)^n', 'chi(psi^n A)=n chi(A) (mod 4)'."""
    A = [phi, psi]                   # M=phi, charge {0,2}
    for n in range(1, 6):
        An = [phi**n, psi**n]
        # M(psi^n A) = phi^n = M(A)^n, checked exactly
        assert is_zero(mahler_exact(An) - phi**n)
        assert is_zero(mahler_exact(An) - mahler_exact(A)**n)
        want = {(n * c) % 4 for c in charge_set(A)}
        assert charge_set(An) == want


# ----------------------------------------------------------------------------
# Rem 3.7 / Ex 3.8 -- the tropical coupling law bites
# ----------------------------------------------------------------------------
def test_tropical_coupling_products_match_paper_forms():
    """Ex 3.8: A(x)B with A=spec(R^2)={phi^2,psi^2}, B=spec(R)={phi,psi} has
    coupled spectrum '{phi^3, -phi, phi^{-1}, -phi^{-3}}'.  Compare the raw
    pairwise products to the paper's simplified forms symbolically."""
    raw = [phi**2 * phi, phi**2 * psi, psi**2 * phi, psi**2 * psi]
    claimed = [phi**3, -phi, 1 / phi, -phi**(-3)]
    for r, c in zip(raw, claimed):
        assert is_zero(sp.nsimplify(r) - sp.nsimplify(c))


def test_tropical_magnitude_is_phi4_not_phi3():
    """Ex 3.8: 'M(A x B) = phi^3 * phi = phi^4' by the tropical law, whereas
    the naive product M(A)M(B) = phi^2 * phi = phi^3 is wrong."""
    A = [phi**2, psi**2]
    B = [phi, psi]
    coupled = [a * b for a in A for b in B]
    # exact: which products are outside the unit circle?
    out = outside_unit(coupled)
    assert len(out) == 2
    exact_prod = sp.nsimplify(sp.prod([sp.Abs(v) for v in out]))
    assert is_zero(exact_prod - phi**4)                 # tropical answer is phi^4
    # the naive product M(A)M(B) = phi^2 * phi = phi^3, which is wrong: it
    # differs from the true tropical magnitude by a factor of phi.
    MA_MB = phi**2 * phi                                # = phi^3
    assert is_zero(MA_MB - phi**3)
    assert not is_zero(exact_prod - MA_MB)              # phi^4 != phi^3
    assert is_zero(exact_prod / MA_MB - phi)            # off by exactly phi


def test_superposition_does_factor():
    """Ex 3.8: 'M(A + B) = phi^3 = M(A) M(B)' -- union DOES factor."""
    A = [phi**2, psi**2]
    B = [phi, psi]
    assert is_zero(mahler_exact(A + B) - phi**3)
    assert is_zero(mahler_exact(A + B) - mahler_exact(A) * mahler_exact(B))


# ----------------------------------------------------------------------------
# Ex 3.9 -- Adams grading on the golden object
# ----------------------------------------------------------------------------
def test_exercise_adams_psi3():
    """Ex 3.9: psi^3(A_phi)={phi^3,psi^3}: M=phi^3=M(A_phi)^3, chi=3{0,2}={0,2};
    and psi^3(psi^2(A_phi))=psi^6(A_phi)={phi^6,psi^6}, M=phi^6."""
    A = [phi, psi]
    A3 = [phi**3, psi**3]
    assert is_zero(sp.nsimplify(sp.prod([sp.Abs(v) for v in outside_unit(A3)])) - phi**3)
    assert charge_set(A3) == {(3 * c) % 4 for c in charge_set(A)} == {0, 2}
    A6 = [phi**6, psi**6]
    assert is_zero(sp.nsimplify(sp.prod([sp.Abs(v) for v in outside_unit(A6)])) - phi**6)


# ----------------------------------------------------------------------------
# Prop 3.10 -- Adams pulls back to the core
# ----------------------------------------------------------------------------
def test_adams_pullback_binet():
    """Prop 3.10: on A_phi each eigenvalue phi^n = F_n phi + F_{n-1} (Binet)."""
    from _eap_helpers import fib
    for n in range(0, 9):
        assert is_zero(phi**n - (fib(n) * phi + fib(n - 1)))


def test_psi2_is_the_keystone():
    """Prop 3.10: 'psi^2 IS the keystone': it sends phi -> phi^2 = phi + 1."""
    assert is_zero(phi**2 - (phi + 1))


def test_adams_multiplicative_not_additive():
    """Prop 3.10: psi^n respects x (multiplicative) but psi^2 is NOT additive:
    'psi^2(phi+1) - [psi^2(phi)+psi^2(1)] = 2phi' (differs by 2phi != 0)."""
    # multiplicativity on a product: (a*b)^n = a^n b^n
    a, b, n = phi, psi, 3
    assert is_zero((a * b)**n - a**n * b**n)
    # non-additivity witnessed exactly
    lhs = (phi + 1)**2               # psi^2 applied to (phi+1)
    rhs = phi**2 + 1**2              # psi^2(phi) + psi^2(1)
    assert is_zero(sp.expand(lhs - rhs) - 2 * phi)
