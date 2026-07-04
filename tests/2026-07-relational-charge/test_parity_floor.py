r"""The reference-free parity floor: Lemmas 5.1 (golden), 5.2 (q_k), 5.3
(twist), Definition/Theorem 5.5 (parityrel), Corollary 5.7 (parity bit), and
ledger entries J, K, L, R.

Paper (Lemma 5.1, lem:golden): phi/phi' = -phi^2, hence
t_rel(phi, phi') = 1/2.
Paper (Lemma 5.2, lem:qk): q_k = x^{2k}+x^k-1 has a real root of each sign,
angle set all of (1/2k)Z/Z, Delta(q_k)=Z/2kZ, anchor n=2k, and M(q_k)=phi.
Paper (Lemma 5.3, lem:twist): T(p)=(-1)^{deg p} p(-x) multiplies roots by -1,
shifts angles by 1/2, fixes differences, preserves Mahler measure and monic
integrality; T(x^3+2)=x^3-2; odd-m coset arithmetic (ledger L).
Paper (Theorem 5.5): mu_rel(m)=phi (even, attained by q_{m/2}); =mu(m) (odd).
Paper (Corollary 5.7): 1/2 in Delta(O) iff |Delta(O)| even.
Paper (ledger R): (1+sqrt5)/2=phi and (3+sqrt5)/2=phi^2.
"""
from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x
phi = (1 + sp.sqrt(5)) / 2
phip = (1 - sp.sqrt(5)) / 2  # the conjugate phi'


def test_golden_ratio_and_square_identities():
    """Ledger R: (1+sqrt5)/2 = phi and (3+sqrt5)/2 = phi^2."""
    assert sp.simplify(phi - (1 + sp.sqrt(5)) / 2) == 0
    assert sp.simplify(phi**2 - (3 + sp.sqrt(5)) / 2) == 0
    assert sp.simplify(phi**2 - phi - 1) == 0  # phi^2 = phi + 1


def test_golden_internal_relation_is_order_two():
    """Lemma 5.1 (ledger J): phi' = -1/phi, so phi/phi' = -phi^2 (negative
    real), giving relative angle 1/2, i.e. t_rel(phi, phi') = 1/2."""
    assert sp.simplify(phip - (-1 / phi)) == 0
    ratio = sp.simplify(phi / phip)
    assert sp.simplify(ratio - (-phi**2)) == 0
    assert ratio < 0            # negative real => angle 1/2
    # relative charge t_rel = Arg(ratio)/2pi = 1/2
    trel = (mp.arg(mp.mpf(str(sp.N(ratio, 40)))) / (2 * mp.pi)) % 1
    assert abs(trel - mp.mpf(1) / 2) < mp.mpf(10) ** -25


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
def test_qk_substitution_roots_have_opposite_signs(k):
    """Lemma 5.2 (ledger K): the x^k-substitution of q_k gives y^2+y-1 with
    roots 1/phi > 0 and -phi < 0."""
    y = sp.symbols("y")
    sub = y**2 + y - 1
    roots = sp.solve(sub, y)
    assert sp.simplify(max(roots) - (1 / phi)) == 0
    assert sp.simplify(min(roots) - (-phi)) == 0
    assert max(roots) > 0 and min(roots) < 0


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5, 6])
def test_qk_mahler_measure_is_phi(k):
    """Lemma 5.2: M(q_k) = phi for every k (the floor value phi is attained)."""
    q = x ** (2 * k) + x**k - 1
    assert abs(R.mahler_measure(q) - mp.mpf(str(sp.N(phi, 45)))) < mp.mpf(10) ** -25


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
def test_qk_angle_set_and_relational_group(k):
    """Lemma 5.2: q_k has angle set all of (1/2k)Z/Z and relational group
    Z/2kZ (anchor n = 2k)."""
    q = x ** (2 * k) + x**k - 1
    ang = set(R.angles_of(q))
    expected = {Fraction(j, 2 * k) for j in range(2 * k)}
    assert ang == expected
    assert R.relational_order(q) == 2 * k
    assert R.absolute_lcd(q) == 2 * k  # anchor n = 2k


def test_qk_realizes_the_half_difference():
    """Lemma 5.2: q_k realizes t_rel = 1/2 (a real root of each sign gives a
    relative angle of 1/2).  Check 1/2 is among the pairwise angle
    differences."""
    for k in [1, 2, 3]:
        q = x ** (2 * k) + x**k - 1
        ang = R.angles_of(q)
        diffs = {((a - b) % 1) for a in ang for b in ang}
        assert Fraction(1, 2) in diffs


@pytest.mark.parametrize(
    "p,expected",
    [
        (x**3 + 2, x**3 - 2),           # ledger B / L
        (x**2 - x - 1, x**2 + x - 1),   # golden <-> its twist
        (x**4 + x**2 - 1, x**4 + x**2 - 1),  # even-degree, even powers: fixed
    ],
)
def test_sign_twist_action(p, expected):
    """Lemma 5.3: T(p)(x) = (-1)^{deg p} p(-x)."""
    T = sp.expand((-1) ** sp.Poly(p, x).degree() * p.subs(x, -x))
    assert sp.expand(T - expected) == 0


@pytest.mark.parametrize(
    "p",
    [x**3 + 2, x**2 - x - 1, x**4 + 5 * x**2 - 5, x**5 - 2, x**6 + x**3 - 1],
)
def test_twist_preserves_mahler_and_relational_data(p):
    """Lemma 5.3: the sign twist preserves the Mahler measure and every
    pairwise angle difference (hence the relational group Delta)."""
    T = sp.expand((-1) ** sp.Poly(p, x).degree() * p.subs(x, -x))
    assert abs(R.mahler_measure(p) - R.mahler_measure(T)) < mp.mpf(10) ** -25
    assert R.relational_order(p) == R.relational_order(T)
    assert sp.Poly(T, x).is_monic  # preserves monic integrality


def test_twist_exchanges_anchor_sectors_at_odd_m():
    """Lemma 5.3: for odd m the twist exchanges the n=m and n=2m sectors.
    x^3+2 (m=3, n=6) maps to x^3-2 (m=3, n=3)."""
    p = x**3 + 2
    T = sp.expand((-1) ** 3 * p.subs(x, -x))
    assert R.relational_order(p) == 3 and R.absolute_lcd(p) == 6
    assert R.relational_order(T) == 3 and R.absolute_lcd(T) == 3


@pytest.mark.parametrize("m", [3, 5, 7, 9])
def test_coset_arithmetic_ledger_L(m):
    """Ledger L: for odd e and odd m, (e+m)/2 is an integer and
    (e+m)/(2m) lies in (1/m)Z/Z (the one-line step of Lemma 5.3)."""
    for e in range(1, 2 * m, 2):  # odd e < 2m
        assert (e + m) % 2 == 0
        # (e+m)/(2m) = ((e+m)/2)/m has denominator dividing m
        val = Fraction(e + m, 2 * m)
        assert m % val.denominator == 0


def test_parity_bit_is_evenness_of_relational_order():
    """Corollary 5.7: 1/2 in Delta(O) iff |Delta(O)| is even (Z/mZ contains an
    order-2 element iff m is even).  Equivalently, 1/2 in (1/m)Z/Z iff m even."""
    for m in range(1, 40):
        # 1/2 in (1/m)Z/Z iff m even:
        in_group = any(Fraction(j, m) == Fraction(1, 2) for j in range(m))
        assert in_group == (m % 2 == 0)


def test_even_floor_attained_by_q_half_m():
    """Theorem 5.5: mu_rel(m) = phi for even m, attained by q_{m/2}
    (M = phi, relational group Z/mZ)."""
    for m in [2, 4, 6, 8]:
        k = m // 2
        q = x ** (2 * k) + x**k - 1
        assert R.relational_order(q) == m
        assert abs(R.mahler_measure(q) - mp.mpf(str(sp.N(phi, 45)))) < mp.mpf(10) ** -25


def test_odd_floor_value_two_realized_by_xm_minus_2():
    """Theorem 5.5 / Remark 5.6: the odd attainment at 2 (by x^m-2):
    mu_rel(3) = 2 realized by x^3-2 (relational Z/3, M=2)."""
    for m in [3, 5, 7]:
        p = x**m - 2
        assert R.relational_order(p) == m
        assert abs(R.mahler_measure(p) - 2) < mp.mpf(10) ** -30
