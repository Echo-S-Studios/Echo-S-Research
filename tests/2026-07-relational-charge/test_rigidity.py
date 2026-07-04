r"""Rigidity and anchors: Theorem 4.1 (thm:rigidity), Example 4.6 (ex:drop),
Lemma 4.4 (lem:oddfull), Example 3.5 (ex:firsttype).

Paper (Theorem 4.1): a relationally admissible object with relational group
Z/mZ has all angles rational, is charge-admissible, and its absolute charge
order n satisfies n in {m, 2m}.  "Both anchors occur: x^m-2 has angles {j/m},
so n=m=Delta; x^3+2 has angles {1/6,1/2,5/6}, so n=6 while the differences
generate 1/3 Z/Z, i.e. m=3, n=2m."

Paper (Example 4.6): g=x^4+5x^2+5 irreducible (Eisenstein at 5), all roots
purely imaginary, angles {1/4,3/4}, absolute Z/4, relational Z/2, M(g)=5;
K=x^4+5x^2-5 keeps the full Delta=Z/4 (real pair at {0,1/2}).

Paper (Lemma 4.4): odd absolute order => relational and absolute groups
coincide (n odd => Delta = Z/nZ).

Angles are recovered independently from high-precision numeric roots.
"""
from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6, 7])
def test_anchor_xm_minus_2_has_n_equals_m(m):
    """x^m - 2 realizes the anchor n = m = |Delta|: its angles are exactly
    {j/m}, so both the absolute order (lcd of angles) and the relational order
    equal m."""
    p = x**m - 2
    ang = set(R.angles_of(p))
    assert ang == {Fraction(j, m) for j in range(m)}
    assert R.absolute_lcd(p) == m
    assert R.relational_order(p) == m


def test_anchor_x3_plus_2_has_n_equals_2m():
    """x^3 + 2: angles {1/6, 1/2, 5/6}; absolute order n = 6, relational order
    m = 3, so n = 2m (the second anchor, ledger B)."""
    p = x**3 + 2
    assert set(R.angles_of(p)) == {Fraction(1, 6), Fraction(1, 2), Fraction(5, 6)}
    assert R.absolute_lcd(p) == 6
    assert R.relational_order(p) == 3


def test_rigidity_n_in_m_or_2m_across_admissible_objects():
    """Theorem 4.1: for every charge-admissible object, n in {m, 2m}."""
    admissible = [
        x**2 - 2, x**3 - 2, x**4 - 2, x**5 - 2, x**6 - 2,
        x**3 + 2, x**4 + 2,
        x**4 + 5 * x**2 + 5, x**4 + 5 * x**2 - 5,
        x**4 + x**2 - 1,           # q2
        x**2 - x - 1,              # golden
        x**6 + x**3 - 1,           # q3
    ]
    for p in admissible:
        n = R.absolute_lcd(p)
        m = R.relational_order(p)
        assert n in (m, 2 * m), (sp.srepr(p), n, m)


def test_group_drop_x4_plus_5x2_plus_5():
    """Example 4.6 (ledger E): g=x^4+5x^2+5 is Eisenstein-at-5 irreducible,
    all four roots purely imaginary (angles {1/4,3/4}), absolute Z/4 but
    relational Z/2 (the group drop), and M(g)=5."""
    g = x**4 + 5 * x**2 + 5
    assert sp.Poly(g, x).is_irreducible
    # Eisenstein at 5: 5 | lower coeffs, 5 !| lead, 25 !| constant
    coeffs = sp.Poly(g, x).all_coeffs()  # [1,0,5,0,5]
    assert coeffs[0] % 5 != 0
    assert all(c % 5 == 0 for c in coeffs[1:])
    assert coeffs[-1] % 25 != 0
    # x^2-roots (-5 +- sqrt5)/2 are both negative -> purely imaginary roots
    r1 = (-5 + sp.sqrt(5)) / 2
    r2 = (-5 - sp.sqrt(5)) / 2
    assert r1 < 0 and r2 < 0
    assert set(R.angles_of(g)) == {Fraction(1, 4), Fraction(3, 4)}
    assert R.absolute_lcd(g) == 4           # absolute Z/4
    assert R.relational_order(g) == 2       # relational Z/2 (drop)
    assert abs(R.mahler_measure(g) - 5) < mp.mpf(10) ** -30


def test_sign_mirror_K_keeps_full_group():
    """Example 4.6 (ledger F): K=x^4+5x^2-5 (sign-mirror of g) keeps the full
    Delta = Z/4 because its real pair sits at angles {0, 1/2}."""
    K = x**4 + 5 * x**2 - 5
    ang = set(R.angles_of(K))
    assert ang == {Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)}
    assert R.relational_order(K) == 4
    assert R.absolute_lcd(K) == 4


def test_squared_moduli_of_g_exceed_one():
    """Example 4.6 detail: the squared moduli (5 -+ sqrt5)/2 are ~1.38 and
    ~3.62, both > 1 iff 5 - sqrt5 > 2 iff sqrt5 < 3; hence M = |g(0)| = 5."""
    a = (5 - sp.sqrt(5)) / 2
    b = (5 + sp.sqrt(5)) / 2
    assert sp.simplify(a - 1) > 0
    assert sp.simplify(b - 1) > 0
    assert sp.sqrt(5) < 3
    assert sp.nsimplify(a * b) == 5  # product of squared moduli = g(0)


def test_first_type_examples():
    """Example 3.5 (ex:firsttype): x^3-2 is a single class with Delta=Z/3
    (angles {0,1/3,2/3}); x^4+5x^2+5 is a single class with Delta=Z/2."""
    assert set(R.angles_of(x**3 - 2)) == {Fraction(0), Fraction(1, 3), Fraction(2, 3)}
    assert R.relational_order(x**3 - 2) == 3
    assert R.relational_order(x**4 + 5 * x**2 + 5) == 2


@pytest.mark.parametrize(
    "name,p",
    [
        ("x^3-2 (n=3)", x**3 - 2),
        ("x^5-2 (n=5)", x**5 - 2),
        ("x^3+2 (n=6)", x**3 + 2),
        ("x^6+x^3-1 (n=6)", x**6 + x**3 - 1),
        ("x^7-2 (n=7)", x**7 - 2),
    ],
)
def test_odd_anchors_are_full(name, p):
    """Lemma 4.4 (lem:oddfull): if the ABSOLUTE order n is odd then
    Delta = Z/nZ (relational order m equals n; no 2m doubling can occur at odd
    n).  For objects with even n (e.g. x^3+2 has n=6) the lemma is silent."""
    m = R.relational_order(p)
    n = R.absolute_lcd(p)
    if n % 2 == 1:
        assert m == n
