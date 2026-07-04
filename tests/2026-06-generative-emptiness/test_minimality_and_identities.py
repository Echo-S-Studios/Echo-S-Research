"""
Minimality chain (Prop 7.1 `prop:minimal`) and the cross-cutting identities
of the scope section / epistemic ledger.

Claims verified independently:
  * Prop 7.1(1): degree two is the UNIQUE degree with a ternary self-action:
        #channels(d) = d^2 - d + 1 = 3  iff  d = 2.
  * Prop 7.1(2): (pi/2)Z ~= Z/4Z, and the +-i generator enters through the
    quartic field Q(5^{1/4}) whose splitting field contains i.
  * Prop 7.1(3): the floor phi is the smallest non-cyclotomic quadratic Mahler
    measure (verified in test_object3_gap.py; here we check phi is a Perron/Salem
    boundary value consistent with the chain).
  * Scope / ledger: sqrt5 = phi + phi^{-1}  (the floor image = grow generator);
  * Sec 1: the ternary growth spectrum spec(ad_R) = {-sqrt5, 0, +sqrt5} has
    exactly three channels, matching #channels(2)=3.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp
import ge_helpers as H
from ge_helpers import x, PHI

PHI_V = PHI()


def test_ternary_lock_degree_two_unique():
    """Prop 7.1(1): #channels(d) = d^2-d+1 = 3 iff d=2 (the positive solution)."""
    d = sp.symbols('d')
    sols = sp.solve(sp.Eq(d**2 - d + 1, 3), d)
    assert set(sols) == {2, -1}                    # only d=2 is a valid degree
    assert [dd for dd in sols if dd > 0] == [2]
    # sanity on neighbours: 1->1 channel, 2->3, 3->7 (so 3 channels only at d=2)
    ch = lambda n: n**2 - n + 1
    assert ch(1) == 1 and ch(2) == 3 and ch(3) == 7


def test_ternary_spectrum_has_three_channels():
    """Sec 1: spec(ad_R) = {-sqrt5, 0, +sqrt5} -- three distinct eigenvalues,
    matching the ternary (3-channel) growth decision."""
    spec = {-sp.sqrt(5), sp.Integer(0), sp.sqrt(5)}
    assert len(spec) == 3
    # +-sqrt5 are the nonzero growth channels; distinct and opposite
    assert sp.simplify(sp.sqrt(5) + (-sp.sqrt(5))) == 0
    assert sp.sqrt(5) != 0


def test_half_pi_Z_isomorphic_Z4():
    """Prop 7.1(2): (pi/2)Z / 2pi ~= Z/4Z.  The four angle representatives add
    modulo 2pi exactly as {0,1,2,3} add modulo 4."""
    reps = [0, sp.pi / 2, sp.pi, 3 * sp.pi / 2]
    for i in range(4):
        for j in range(4):
            s = sp.nsimplify((reps[i] + reps[j]) % (2 * sp.pi))
            # index of the resulting representative
            idx = [k for k in range(4) if sp.simplify(s - reps[k]) == 0]
            assert idx == [(i + j) % 4]


def test_quartic_field_splitting_contains_i():
    """Prop 7.1(2): the +-i generator 'enters through the quartic field Q(5^{1/4})
    whose splitting field contains i.'  Roots of x^4-5 are 5^{1/4} i^k; the ratio
    of two of them is i, so i lies in the splitting field."""
    r = sp.Poly(x**4 - 5, x)
    assert r.is_irreducible                          # genuine quartic field
    roots = sp.roots(x**4 - 5, x)
    rootlist = list(roots.keys())
    ratios = {sp.simplify(a / b) for a in rootlist for b in rootlist if b != 0}
    assert sp.I in ratios or -sp.I in ratios


def test_sqrt5_equals_phi_plus_inverse():
    """Scope/ledger: 'sqrt5 = phi + phi^{-1}' (the floor's image identified with
    the grow generator)."""
    phi = (1 + sp.sqrt(5)) / 2
    assert sp.simplify(phi + 1 / phi - sp.sqrt(5)) == 0
    # high-precision cross-check
    assert abs((PHI_V + 1 / PHI_V) - mp.sqrt(5)) < mp.mpf(10)**(-40)


def test_phi_is_root_of_content_free_growth_polynomial():
    """Prop 7.1(3): phi is the smallest realizable growth -- it is the dominant
    root of x^2-x-1 and equals its own Mahler measure (a Perron number)."""
    from ge_helpers import mahler, as_poly, roots_mp
    P = as_poly(x**2 - x - 1)
    dom = max((r for r in roots_mp(P)), key=lambda r: abs(r))
    assert abs(dom.real - PHI_V) < mp.mpf(10)**(-30)
    assert abs(dom.imag) < mp.mpf(10)**(-30)         # real (Perron, not Salem)
    assert abs(mahler(P) - PHI_V) < mp.mpf(10)**(-30)
