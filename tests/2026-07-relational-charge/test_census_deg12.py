r"""Exhaustive degree-12 census: Theorem 6.13 (thm:census), Remark 6.14
(rem:censusscope), ledger T and W.

Paper (Theorem 6.13): all 3^6 = 729 monic reciprocal polynomials
x^12 + c1 x^11 + ... + c6 x^6 + ... + c1 x + 1 with c_i in {-1,0,1}.  The sign
twist p(x)->p(-x) preserves the family and fixes Rat_p; it fixes exactly the
3^3 = 27 vectors with vanishing odd coefficients, giving (729+27)/2 = 378
twist-classes (Burnside).  Exactly 37 twist-classes are Salem; every complete
scan (deg Rat_p = 144, bound 41472) returns {Phi_1^12}: every certified Salem
number in the family is relationally inert.

Paper (Remark 6.14 / ledger W): rejection tally of the 341 non-Salem classes:
39 vanish at x=+-1, 257 fail the trace-polynomial Sturm pattern (1,0,5), 45 are
reducible (none fail only trace-polynomial irreducibility); 39+257+45+37=378.

NOTE ON THE 257/45 SPLIT: a purely count-based Sturm pattern test gives
256 trace-failures and 46 reducibles (their sum, 302, matches 257+45).  Exactly
one reducible polynomial -- x^12-x^11-x^10-x^9-x^7-x^6-x^5-x^3-x^2-x+1, which
has the cyclotomic factor Phi_10 and thus a coincidental trace pattern (1,0,5)
-- is bucketed differently.  All robust invariants (729, 378, 27, 37 Salem, 39
vanish, sum 378) reproduce exactly; only the internal trace/reducible split
differs by one, with no clean mathematical rule selecting that single poly.
The exact 257 value is therefore recorded as an xfail (see NOTES.md).
"""
import itertools

import pytest
import sympy as sp

import _relcharge as R

x = R.x


def build(c):
    """Reciprocal monic degree-12 polynomial from (c1,...,c6)."""
    c1, c2, c3, c4, c5, c6 = c
    return sp.Poly([1, c1, c2, c3, c4, c5, c6, c5, c4, c3, c2, c1, 1], x)


def twist(c):
    """p(x)->p(-x) on coefficient vectors: odd-index coeffs flip sign."""
    c1, c2, c3, c4, c5, c6 = c
    return (-c1, c2, -c3, c4, -c5, c6)


def pattern_ok(p):
    """Count-based trace Sturm pattern (1,0,5) with no root at +-2."""
    a, b, mid, at2, atm2 = R.trace_sturm_pattern(p.as_expr())
    return a + b == 1 and mid == 5 and at2 == 0 and atm2 == 0


# ---- build twist-class representatives once -------------------------------
_ALL = list(itertools.product([-1, 0, 1], repeat=6))
_SEEN = set()
_REPS = []
_FIXED = 0
for _c in _ALL:
    if _c in _SEEN:
        continue
    _t = twist(_c)
    if _t == _c:
        _FIXED += 1
    _SEEN.add(_c)
    _SEEN.add(_t)
    _REPS.append(_c)

# ---- classify each representative ------------------------------------------
_BUCKETS = {"vanish": 0, "trace": 0, "reducible": 0, "salem": 0}
_SALEM_REPS = []
for _c in _REPS:
    _p = build(_c)
    if _p.eval(1) == 0 or _p.eval(-1) == 0:
        _BUCKETS["vanish"] += 1
    elif not pattern_ok(_p):
        _BUCKETS["trace"] += 1
    elif not _p.is_irreducible:
        _BUCKETS["reducible"] += 1
    else:
        _BUCKETS["salem"] += 1
        _SALEM_REPS.append(_c)


def test_family_size_729():
    """Theorem 6.13: the family has 3^6 = 729 members."""
    assert len(_ALL) == 729 == 3**6


def test_twist_fixed_count_27():
    """Theorem 6.13: the twist fixes exactly 3^3 = 27 vectors (vanishing odd
    coefficients c1=c3=c5=0)."""
    assert _FIXED == 27 == 3**3
    fixed_direct = sum(1 for c in _ALL if twist(c) == c)
    assert fixed_direct == 27


def test_burnside_378_orbits():
    """Theorem 6.13: Burnside gives (729 + 27)/2 = 378 twist-classes, matching
    the direct enumeration."""
    assert len(_REPS) == 378 == (729 + 27) // 2


def test_salem_count_37():
    """Theorem 6.13: exactly 37 twist-classes are certified Salem."""
    assert _BUCKETS["salem"] == 37
    assert len(_SALEM_REPS) == 37


def test_vanish_count_39():
    """Remark 6.14: 39 twist-classes vanish at x = +-1."""
    assert _BUCKETS["vanish"] == 39


def test_rejection_tally_sums_to_378():
    """Remark 6.14: all buckets account for every twist-class."""
    assert sum(_BUCKETS.values()) == 378


def test_combined_trace_and_reducible_is_302():
    """Convention-independent invariant: (trace-fail) + (reducible) = 302,
    matching the paper's 257 + 45 = 302 (our clean count-based split is
    256 + 46 = 302; see module docstring / NOTES.md)."""
    assert _BUCKETS["trace"] + _BUCKETS["reducible"] == 302
    assert _BUCKETS["trace"] == 256      # our count-based split
    assert _BUCKETS["reducible"] == 46


@pytest.mark.xfail(
    reason="under review: the paper's 257/45 internal split differs from the "
    "clean count-based split (256/46) by exactly one reducible poly "
    "(x^12-x^11-...-x+1, having factor Phi_10 and a coincidental (1,0,5) trace "
    "pattern). Robust invariants -- 729/378/27/37 Salem/39 vanish/302 combined "
    "-- all reproduce exactly. See NOTES.md 'Flagged for human review'.",
    strict=True,
)
def test_paper_exact_trace_reducible_split():
    """Remark 6.14 verbatim: 257 trace-failures and 45 reducibles."""
    assert _BUCKETS["trace"] == 257
    assert _BUCKETS["reducible"] == 45


def test_the_split_edge_polynomial_is_reducible_with_phi10_factor():
    """Document the single boundary polynomial: it is reducible (factor
    Phi_10 = x^4-x^3+x^2-x+1) yet has trace pattern (1,0,5) by root count."""
    edge = build((-1, -1, -1, 0, -1, -1))
    assert not edge.is_irreducible
    phi10 = sp.cyclotomic_poly(10, x)
    quo, rem = sp.div(edge, sp.Poly(phi10, x))
    assert rem.is_zero                     # Phi_10 divides it
    assert pattern_ok(edge)                # yet trace pattern is (1,0,5)


@pytest.mark.parametrize("c", _SALEM_REPS)
def test_every_census_salem_is_inert(c):
    """Theorem 6.13 / ledger T: each of the 37 certified Salem numbers has
    ratio object of degree 144 with complete signature {Phi_1^12}."""
    p = build(c)
    Rp = R.ratio_poly(p.as_expr())
    assert Rp.degree() == 144
    assert R.contact_signature(Rp) == {1: 12}
