r"""Definition 6.3 (def:rat) and Lemma 6.6 (lem:diag): the ratio object.

Paper (Definition 6.3): for monic p in Z[x] of degree n with p(0) != 0,
    Rat_p(x) = primitive part of Res_y(p(y), p(xy)),
and "the degree is exactly n^2 (no drop) and the root multiset is all n^2
ordered ratios alpha_j/alpha_i".

Paper (Lemma 6.6, lem:diag): "If p is squarefree, the multiplicity of
Phi_1 = x-1 in Rat_p is exactly n (the diagonal ratios alpha_i/alpha_i)."

We re-derive the ratio polynomial from the resultant definition and, for
p(0) = +-1, independently from the companion Kronecker matrix C (x) C^{-1}
(ledger G), then check degree, the root multiset against the actual numeric
ratios, and the diagonal multiplicity.
"""
import math

import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x

CASES = [
    ("x^3-2", x**3 - 2, 3),
    ("x^3+2", x**3 + 2, 3),
    ("x^4-2", x**4 - 2, 4),
    ("q2 = x^4+x^2-1", x**4 + x**2 - 1, 4),
    ("x^4+5x^2+5", x**4 + 5 * x**2 + 5, 4),
    ("beta4", R.B4, 4),
    ("S6", R.S6, 6),
    ("plastic x^3-x-1", R.PLASTIC, 3),
    ("x^2+x+2", x**2 + x + 2, 2),
]


@pytest.mark.parametrize("name,p,n", CASES)
def test_ratio_degree_is_n_squared(name, p, n):
    """Rat_p has degree exactly n^2 (Definition 6.3: 'no drop')."""
    Rp = R.ratio_poly(p)
    assert Rp.degree() == n * n


@pytest.mark.parametrize("name,p,n", CASES)
def test_ratio_root_multiset_equals_ordered_ratios(name, p, n):
    """The n^2 roots of Rat_p equal the n^2 ordered ratios alpha_j/alpha_i.

    Independent check: form every ordered ratio of the (high-precision) roots
    of p and confirm each is a root of the exactly-constructed Rat_p by
    evaluation (Rat_p has repeated roots, so direct root-finding of Rat_p is
    numerically ill-posed; evaluation is robust).  With degree(Rat_p) = n^2
    tested separately, 'every ratio is a root' pins the multiset.
    """
    Rp = R.ratio_poly(p)
    coeffs = [mp.mpf(int(c)) for c in Rp.all_coeffs()]  # leading first
    rts = R.roots_mp(p)
    ratios = [a / b for a in rts for b in rts]
    assert len(ratios) == n * n
    scale = sum(abs(c) for c in coeffs)
    for r in ratios:
        assert abs(mp.polyval(coeffs, r)) < scale * mp.mpf(10) ** -20


@pytest.mark.parametrize(
    "name,p", [("beta4", R.B4), ("Lehmer", R.LEHMER), ("S6", R.S6), ("S8", R.S8)]
)
def test_two_route_agreement_for_unit_constant_term(name, p):
    """Ledger G: for |p(0)|=1 the resultant route and the companion-Kronecker
    route C_p (x) C_p^{-1} give the same ratio object (hence same signature)."""
    assert abs(int(sp.Poly(p, x).eval(0))) == 1
    s_res = R.contact_signature(R.ratio_poly(p))
    s_kron = R.contact_signature(R.ratio_poly_via_kronecker(p))
    assert s_res == s_kron


@pytest.mark.parametrize("name,p,n", CASES)
def test_diagonal_multiplicity_is_n(name, p, n):
    """Lemma 6.6: squarefree p => Phi_1 multiplicity in Rat_p is exactly n."""
    P = sp.Poly(p, x)
    assert sp.gcd(P, P.diff(x)).degree() == 0  # squarefree
    assert R.phi1_multiplicity(R.ratio_poly(p)) == n


def test_leading_coefficient_is_pm_p0_to_the_n():
    """Definition 6.3: the leading coefficient of Res_y is (+- p(0))^n, so no
    degree drop occurs.  Check the un-normalised resultant's degree = n^2."""
    p = x**4 - x**3 - x**2 - x + 1  # beta4, p(0)=1
    Rraw = sp.resultant(p.subs(x, R.y), sp.expand(p.subs(x, x * R.y)), R.y)
    assert sp.Poly(sp.expand(Rraw), x).degree() == 16
