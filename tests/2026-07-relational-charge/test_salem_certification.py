r"""Salem certificates and Mahler values: Theorem 7.9 (thm:inert),
Theorem 7.10 (thm:inertmore), ledger O, and Example 7.19 (plastic number).

Paper (Theorem 7.10, ledger O): beta_4, S_6, S_8, L are each reciprocal and
irreducible, and the trace polynomial T (with p(x)=x^{deg/2} T(x+1/x)) has, by
Sturm count, exactly one real root in (2, inf), none in (-inf, -2], and the
remaining deg/2 - 1 in (-2, 2).

Paper (Theorem 7.9): beta_4 has M = 1.72208..., Lehmer L has M = 1.17628...

Paper (Example 7.19, ledger W): the plastic number (root of x^3-x-1) is
irreducible with one real root > 1 (Sturm) and conjugate pair at modulus
theta_0^{-1/2}; contacts {Phi_1^3} -- inert.
"""
import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x

SALEMS = [
    ("beta4", R.B4, 4),
    ("S6", R.S6, 6),
    ("S8", R.S8, 8),
    ("Lehmer", R.LEHMER, 10),
]


@pytest.mark.parametrize("name,p,deg", SALEMS)
def test_reciprocal_and_irreducible(name, p, deg):
    """Theorem 7.10: each certified object is reciprocal and irreducible."""
    P = sp.Poly(p, x)
    assert P.degree() == deg
    assert R.is_reciprocal(p)
    assert P.is_irreducible


@pytest.mark.parametrize("name,p,deg", SALEMS)
def test_trace_sturm_pattern(name, p, deg):
    """Theorem 7.10 / ledger O: trace polynomial has Sturm pattern
    (1, 0, deg/2 - 1): one root in (2,inf), none in (-inf,-2], rest in (-2,2),
    and no root at +-2."""
    a, b, mid, at2, atm2 = R.trace_sturm_pattern(p)
    d = deg // 2
    assert (a, b, mid) == (1, 0, d - 1)
    assert at2 == 0 and atm2 == 0


@pytest.mark.parametrize("name,p,deg", SALEMS)
def test_is_salem_by_full_certificate(name, p, deg):
    """The composite exact Salem certificate accepts all four."""
    assert R.is_salem_polynomial(p)


def test_beta4_mahler_measure():
    """Theorem 7.9: M(beta_4) = 1.72208... (the Salem number itself)."""
    val = R.mahler_measure(R.B4)
    assert abs(val - mp.mpf("1.72208")) < mp.mpf(10) ** -5  # paper's stated digits
    # M(beta4) equals its largest real root (the Salem number), and is the
    # dominant root of x^4-x^3-x^2-x+1:
    largest = max((r.real for r in R.roots_mp(R.B4) if abs(r.imag) < mp.mpf(10) ** -20))
    assert abs(val - largest) < mp.mpf(10) ** -25
    assert abs(sp.Poly(R.B4, x).eval(sp.Float(mp.nstr(val, 40), 40))) < mp.mpf(10) ** -30


def test_lehmer_mahler_measure():
    """Theorem 7.9: M(L) = 1.17628... (Lehmer's number, the smallest known
    Salem number)."""
    val = R.mahler_measure(R.LEHMER)
    assert abs(val - mp.mpf("1.17628")) < mp.mpf(10) ** -5  # paper's stated digits
    # it is a genuine root of Lehmer's polynomial:
    assert abs(sp.Poly(R.LEHMER, x).eval(sp.Float(mp.nstr(val, 40), 40))) < mp.mpf(10) ** -28


def test_salem_number_is_unique_root_outside_unit_circle():
    """Definitional sanity for a Salem polynomial: exactly one root of modulus
    > 1, one of modulus < 1, the rest exactly on the unit circle."""
    for name, p, deg in SALEMS:
        mods = sorted(abs(r) for r in R.roots_mp(p))
        outside = [m for m in mods if m > 1 + mp.mpf(10) ** -18]
        inside = [m for m in mods if m < 1 - mp.mpf(10) ** -18]
        oncirc = [m for m in mods if abs(m - 1) <= mp.mpf(10) ** -18]
        assert len(outside) == 1
        assert len(inside) == 1
        assert len(oncirc) == deg - 2


def test_plastic_number_certificate():
    """Example 7.19 / ledger W: x^3-x-1 is irreducible with exactly one real
    root > 1, and a complex pair at modulus theta_0^{-1/2} (root product = 1),
    with contact signature {Phi_1^3} -- relationally inert."""
    p = R.PLASTIC
    P = sp.Poly(p, x)
    assert P.is_irreducible
    # one real root, and it is > 1:
    assert P.count_roots(1, sp.oo) == 1
    assert P.count_roots(-sp.oo, sp.oo) == 1  # exactly one real root total
    theta0 = max((r.real for r in R.roots_mp(p) if abs(r.imag) < mp.mpf(10) ** -20))
    assert theta0 > 1
    # complex pair modulus = theta0^{-1/2}  (product of all three roots = -p(0)=1)
    cmods = [abs(r) for r in R.roots_mp(p) if abs(r.imag) > mp.mpf(10) ** -20]
    assert len(cmods) == 2
    assert abs(cmods[0] - theta0 ** mp.mpf("-0.5")) < mp.mpf(10) ** -20
    # inert:
    assert R.contact_signature(R.ratio_poly(p)) == {1: 3}


def test_smyth_constant_is_plastic_number():
    """Remark 2.3 context: Smyth's constant theta_0 = 1.32471... is the real
    root of x^3-x-1 (the smallest Pisot / plastic number)."""
    theta0 = max((r.real for r in R.roots_mp(R.PLASTIC) if abs(r.imag) < mp.mpf(10) ** -20))
    assert abs(theta0 - mp.mpf("1.3247179572447460")) < mp.mpf(10) ** -14
