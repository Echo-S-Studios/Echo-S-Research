"""Canonical corroboration (Sec. 2): the Lehmer scan.

Lehmer's degree-10 Salem polynomial L has deg Rat_L = 100, scan bound 20000,
and the complete contact scan returns exactly {Phi_1^{10}} -- i.e. Lehmer's
number is relationally inert (no non-trivial root-of-unity ratio among its
conjugates).  Re-derived here by building Rat_L and factoring it.
"""
import os
import sys

from sympy import symbols, Poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import rat_object, cyclotomic_scan

x = symbols('x')

# Lehmer polynomial x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]


def test_lehmer_is_the_lehmer_polynomial():
    """Sanity: L is reciprocal (Salem) of degree 10 and irreducible."""
    assert LEHMER == LEHMER[::-1]
    assert Poly(LEHMER, x).is_irreducible


def test_lehmer_rat_degree_and_bound():
    """deg Rat_L = n^2 = 100 and scan bound 2 n^4 = 20000."""
    R = rat_object(LEHMER)
    assert R.degree() == 100
    n = len(LEHMER) - 1
    assert 2 * n ** 4 == 20000


def test_lehmer_scan_is_phi1_10():
    """Complete contact scan of Rat_L returns exactly {Phi_1^{10}}: Lehmer's
    number is relationally inert."""
    R = rat_object(LEHMER)
    assert cyclotomic_scan(R) == {1: 10}
