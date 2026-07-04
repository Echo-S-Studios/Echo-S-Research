"""Section 3: the hypothesis-necessity witness Z* (Prop. 3.4) and the parent
paper's sharpness example x^4-2.

Z* = f1 f2 = (x^2-x-1)(x^2+x-1) = x^4-3x^2+1 carries the torsion ratio -1 between
phi (of f1) and -phi (of f2) at the shared modulus phi, so its scan hits Phi_2;
x^4-2 (all roots on one shell, ratios = 4th roots of unity) hits Phi_2 and Phi_4.
Both witness that a hypothesis of the pinning theorem cannot be dropped.
"""
import os
import sys

import sympy as sp
from sympy import symbols, expand, simplify, Poly, sqrt, Rational, nsimplify, Abs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import rat_object, cyclotomic_scan

x = symbols('x')
PHI = (1 + sqrt(5)) / 2


def test_zstar_factorization():
    """Prop. 3.4: Z* = (x^2-x-1)(x^2+x-1) = x^4-3x^2+1 (exact)."""
    f1 = x ** 2 - x - 1
    f2 = x ** 2 + x - 1
    assert expand(f1 * f2) == x ** 4 - 3 * x ** 2 + 1


def test_f1_f2_irreducible():
    """Prop. 3.4: each f_i is irreducible over Q."""
    assert Poly(x ** 2 - x - 1, x).is_irreducible
    assert Poly(x ** 2 + x - 1, x).is_irreducible


def _same_set(actual, expected):
    """Set equality of algebraic numbers up to symbolic simplification."""
    if len(actual) != len(expected):
        return False
    return all(any(simplify(e - a) == 0 for a in actual) for e in expected)


def test_f1_f2_roots_and_moduli():
    """Prop. 3.4 proof: f1 has roots phi, -1/phi ; f2 has roots 1/phi, -phi.
    Per-factor moduli are {phi, phi^{-1}}, each attained once."""
    r1 = list(sp.roots(x ** 2 - x - 1, x).keys())
    r2 = list(sp.roots(x ** 2 + x - 1, x).keys())
    assert _same_set(r1, [PHI, -1 / PHI])
    assert _same_set(r2, [1 / PHI, -PHI])
    # moduli per factor (all roots real): |phi| = phi, |-1/phi| = 1/phi
    mod1 = sorted(simplify(Abs(r)) for r in r1)
    mod2 = sorted(simplify(Abs(r)) for r in r2)
    assert all(simplify(a - b) == 0 for a, b in zip(mod1, mod2))
    assert simplify(mod1[0] - 1 / PHI) == 0 and simplify(mod1[1] - PHI) == 0


def test_zstar_torsion_ratio_minus_one_at_shared_modulus():
    """Prop. 3.4: phi (root of f1) and -phi (root of f2) share modulus phi and
    their ratio phi/(-phi) = -1 is a primitive 2nd root of unity (torsion)."""
    assert simplify(Abs(PHI) - Abs(-PHI)) == 0            # shared modulus
    ratio = simplify(PHI / (-PHI))
    assert ratio == -1
    assert simplify(ratio ** 2 - 1) == 0                  # order 2 -> Phi_2


def test_zstar_scan_is_phi1_4_phi2_4():
    """Prop. 3.4: complete scan gives Rat_{Z*} ~> {Phi_1^4, Phi_2^4}
    (Phi_1 from the diagonal n=4, Phi_2 from the two antipodal pairs)."""
    R = rat_object([1, 0, -3, 0, 1])
    assert R.degree() == 16
    assert cyclotomic_scan(R) == {1: 4, 2: 4}


def test_x4_minus_2_sharpness_scan():
    """Remark 3.5 / parent sharpness example: x^4-2 is irreducible with all four
    roots on one shell 2^{1/4}, ratios = powers of i, so the scan carries
    Phi_2 and Phi_4 -- exhibiting torsion contacts and thus that the pinning
    hypothesis (uniquely-attained modulus) is necessary."""
    assert Poly(x ** 4 - 2, x).is_irreducible
    R = rat_object([1, 0, 0, 0, -2])
    assert R.degree() == 16
    assert cyclotomic_scan(R) == {1: 4, 2: 4, 4: 4}
