r"""Cross-object locking: Definition 8.1, Theorem 8.2 (thm:notlocked),
Remark 8.4 (rem:allpairs), ledger I and Q.

Paper (Definition 8.1): Rat_{p,q}(x)=Res_y(q(y), p(xy)) has root multiset
{alpha/beta : p(alpha)=0, q(beta)=0}, of degree (deg p)(deg q).

Paper (Theorem 8.2 / ledger I): gcd(beta_4, L)=1, and the complete contact
scan of the degree-40 mixed object Rat_{L,beta_4} returns NO cyclotomic
contact whatsoever -- beta_4 and Lehmer's number are not circle-locked.

Paper (Remark 8.4 / ledger Q): the five remaining pairs among
{beta_4, S_6, S_8, L} (degrees 24, 32, 48, 60, 80) carry no cyclotomic
contact either.
"""
import pytest
import sympy as sp

import _relcharge as R

x = R.x

PAIRS = [
    ("L x beta4", R.LEHMER, R.B4, 40),
    ("beta4 x S6", R.B4, R.S6, 24),
    ("beta4 x S8", R.B4, R.S8, 32),
    ("S6 x S8", R.S6, R.S8, 48),
    ("L x S6", R.LEHMER, R.S6, 60),
    ("L x S8", R.LEHMER, R.S8, 80),
]


@pytest.mark.parametrize("name,p,q,deg", PAIRS)
def test_mixed_ratio_degree(name, p, q, deg):
    """Definition 8.1: deg Rat_{p,q} = (deg p)(deg q)."""
    Rpq = R.mixed_ratio_poly(p, q)
    assert Rpq.degree() == deg
    assert deg == sp.Poly(p, x).degree() * sp.Poly(q, x).degree()


@pytest.mark.parametrize("name,p,q,deg", PAIRS)
def test_no_circle_locking(name, p, q, deg):
    """Theorem 8.2 / Remark 8.4: every mixed ratio object among the four
    certified Salem numbers has an EMPTY cyclotomic-contact signature -- no two
    are circle-locked."""
    assert sp.gcd(sp.Poly(p, x), sp.Poly(q, x)).degree() == 0  # coprime
    assert R.contact_signature(R.mixed_ratio_poly(p, q)) == {}


def test_mixed_ratio_root_multiset():
    """Definition 8.1: the roots of Rat_{p,q} are exactly the cross ratios
    alpha/beta.  Check numerically for the (small) beta4 x S6 pair."""
    import mpmath as mp

    p, q = R.B4, R.S6
    Rpq = R.mixed_ratio_poly(p, q)
    rp = R.roots_mp(p)
    rq = R.roots_mp(q)
    cross = sorted(
        (a / b for a in rp for b in rq),
        key=lambda w: (float(mp.re(w)), float(mp.im(w))),
    )
    got = sorted(
        R.roots_mp(Rpq.as_expr()),
        key=lambda w: (float(mp.re(w)), float(mp.im(w))),
    )
    assert len(cross) == len(got) == 24
    for u, v in zip(cross, got):
        assert abs(u - v) < mp.mpf(10) ** -16


def test_real_locking_is_trivial_but_present():
    """Remark 8.3: the real conjugates tau_{beta4} and tau_L are both at angle
    0, so the (uninformative) rational-block locking is present; the mixed
    object nonetheless has no NON-trivial (Phi_M, M>=2) contact, and in fact no
    Phi_1 contact because gcd=1 (no shared root)."""
    sig = R.contact_signature(R.mixed_ratio_poly(R.LEHMER, R.B4))
    assert 1 not in sig       # no shared root (gcd=1)
    assert sig == {}          # no informative locking either
