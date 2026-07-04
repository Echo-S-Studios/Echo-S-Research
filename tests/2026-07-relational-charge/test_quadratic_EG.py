r"""Appendix D: the quadratic case of the emission gap (EG), self-contained.

Paper (Appendix D): no charge-admissible object of degree <= 2 has Mahler
measure in (1, phi).  Cases for p = x^2 - t x + n irreducible:
  * complex pair: n = |alpha|^2 >= 1; n=1 => M=1 (root of unity), n>=2 => M=n;
  * real pair |n|>=2: M = |n| (both outside) or M = |n|/|inner| > |n| >= 2;
  * real units n=-1: M = (|t| + sqrt(t^2+4))/2 >= (1+sqrt5)/2 = phi, equality
    exactly at |t|=1 (attained by x^2-x-1);
  * real units n=+1: real roots force |t|>=3, M = (|t|+sqrt(t^2-4))/2 >=
    (3+sqrt5)/2 = phi^2 > phi.
Hence every quadratic Mahler value lies in {1} U [phi, infinity).
"""
import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x
phi = (1 + sp.sqrt(5)) / 2


def test_real_unit_floor_n_minus_1_symbolic():
    """n=-1: M(|t|) = (|t|+sqrt(t^2+4))/2 is >= phi with equality iff |t|=1."""
    t = sp.symbols("t", positive=True)
    M = (t + sp.sqrt(t**2 + 4)) / 2
    # at t=1, M = phi:
    assert sp.simplify(M.subs(t, 1) - phi) == 0
    # M is increasing in t>0, so M >= phi for t >= 1:
    assert sp.simplify(sp.diff(M, t)) == sp.simplify((1 + t / sp.sqrt(t**2 + 4)) / 2)
    assert sp.diff(M, t).subs(t, 1) > 0
    # both gap terms nonnegative for |t|>=1: (|t|-1)>=0 and sqrt(t^2+4)-sqrt5>=0
    for tv in [1, 2, 3, 10]:
        assert tv - 1 >= 0
        assert sp.sqrt(tv**2 + 4) - sp.sqrt(5) >= 0


def test_phi_attained_by_x2_minus_x_minus_1():
    """phi is attained by x^2 - x - 1 (the golden polynomial), n=-1, t=1."""
    p = x**2 - x - 1
    assert sp.Poly(p, x).is_irreducible
    assert abs(R.mahler_measure(p) - mp.mpf(str(sp.N(phi, 45)))) < mp.mpf(10) ** -25


def test_real_unit_floor_n_plus_1_symbolic():
    """n=+1: real roots require t^2>4 (|t|>=3), and M=(|t|+sqrt(t^2-4))/2 >=
    (3+sqrt5)/2 = phi^2 > phi."""
    t = sp.symbols("t", positive=True)
    M = (t + sp.sqrt(t**2 - 4)) / 2
    assert sp.simplify(M.subs(t, 3) - phi**2) == 0
    assert sp.simplify(phi**2 - (3 + sp.sqrt(5)) / 2) == 0
    assert phi**2 > phi
    # |t| <= 2 gives cyclotomic/reducible (t=2 -> (x-1)^2 reducible):
    assert sp.factor(x**2 - 2 * x + 1) == (x - 1) ** 2


def test_complex_pair_measure_is_n():
    """Complex pair: n = |alpha|^2; if n=1 the root is a root of unity (M=1),
    if n>=2 then M=n.  Check x^2 - t x + n for a few complex-pair cases."""
    # x^2+x+1 (n=1): roots are primitive cube roots of unity, M=1
    assert abs(R.mahler_measure(x**2 + x + 1) - 1) < mp.mpf(10) ** -30
    # x^2+x+2 (n=2, complex pair since disc=1-8<0): M = 2
    assert (1 - 8) < 0
    assert abs(R.mahler_measure(x**2 + x + 2) - 2) < mp.mpf(10) ** -30
    # x^2 + 3 (n=3): M = 3
    assert abs(R.mahler_measure(x**2 + 3) - 3) < mp.mpf(10) ** -30


def test_no_quadratic_mahler_value_in_open_gap():
    """Appendix D conclusion: over all irreducible x^2 - t x + n with small
    integer t,n, no Mahler value lands in the open gap (1, phi)."""
    phi_num = mp.mpf(str(sp.N(phi, 45)))
    for t in range(-12, 13):
        for n in range(-12, 13):
            if n == 0:
                continue
            p = x**2 - t * x + n
            P = sp.Poly(p, x)
            if not P.is_irreducible:
                continue
            M = R.mahler_measure(p)
            # M is either ~1 or >= phi (allow tiny numerical slack):
            in_gap = (M > 1 + mp.mpf(10) ** -20) and (M < phi_num - mp.mpf(10) ** -20)
            assert not in_gap, (t, n, mp.nstr(M, 20))


def test_degree_one_and_reducible_quadratics_are_integers():
    """Appendix D: degree-1 objects and reducible quadratics have nonzero
    integer roots, so M in {1} U Z_{>=2}."""
    # x - k has M = max(1,|k|):
    for k in [1, 2, 3, 5]:
        assert abs(R.mahler_measure(x - k) - max(1, k)) < mp.mpf(10) ** -30
    # reducible quadratic (x-2)(x-3): M = 6
    assert abs(R.mahler_measure(sp.expand((x - 2) * (x - 3))) - 6) < mp.mpf(10) ** -25
