r"""Cyclotomic-contact signatures: Appendix A ledger entries A-H, P, U, and the
signature-coincidence remarks (rem:scope), Theorems 7.9 / 7.10 (inertness).

Each signature is the complete list {Phi_M^mult} of cyclotomic factors of the
ratio object.  We compute it by factoring Rat_p over Q and naming the
cyclotomic factors (a different algorithm from the paper's bounded scan, and
complete by construction).  For the low-degree instances we additionally
reproduce the signature with a purely numerical mpmath engine.

Ledger targets:
  A x^3-2        -> {Phi_1^3, Phi_3^3}
  B x^3+2        -> {Phi_1^3, Phi_3^3}   (identical to A: gauge-blind)
  C x^4-2        -> {Phi_1^4, Phi_2^4, Phi_4^4}
  D q2=x^4+x^2-1 -> {Phi_1^4, Phi_2^4}
  E x^4+5x^2+5   -> {Phi_1^4, Phi_2^4}
  F x^4+5x^2-5   -> {Phi_1^4, Phi_2^4}   (identical to E; full type differs)
  G beta4        -> {Phi_1^4}
  H Lehmer       -> {Phi_1^10}
  P S6,S8        -> {Phi_1^6}, {Phi_1^8}
  U x^2+x+2 -> {Phi_1^2};  x^4+x^2+2 -> {Phi_1^4, Phi_2^4}
  W plastic x^3-x-1 -> {Phi_1^3}
"""
import pytest

import _relcharge as R

x = R.x

# name, polynomial, expected exact signature, numeric-crosscheck Mmax (or None)
LEDGER = [
    ("A x^3-2", x**3 - 2, {1: 3, 3: 3}, 30),
    ("B x^3+2", x**3 + 2, {1: 3, 3: 3}, 30),
    ("C x^4-2", x**4 - 2, {1: 4, 2: 4, 4: 4}, 40),
    ("D q2=x^4+x^2-1", x**4 + x**2 - 1, {1: 4, 2: 4}, 60),
    ("E x^4+5x^2+5", x**4 + 5 * x**2 + 5, {1: 4, 2: 4}, 60),
    ("F K=x^4+5x^2-5", x**4 + 5 * x**2 - 5, {1: 4, 2: 4}, 60),
    ("G beta4", R.B4, {1: 4}, 600),
    ("P S6", R.S6, {1: 6}, None),
    ("U q=x^2+x+2", x**2 + x + 2, {1: 2}, 30),
    ("U p=x^4+x^2+2", x**4 + x**2 + 2, {1: 4, 2: 4}, 60),
    ("W plastic x^3-x-1", R.PLASTIC, {1: 3}, 30),
]


@pytest.mark.parametrize("name,p,expected,_m", LEDGER)
def test_exact_contact_signature(name, p, expected, _m):
    """Exact signature via factor-and-identify matches the ledger."""
    assert R.contact_signature(R.ratio_poly(p)) == expected


@pytest.mark.parametrize(
    "name,p,expected,Mmax", [c for c in LEDGER if c[3] is not None]
)
def test_numeric_contact_signature(name, p, expected, Mmax):
    """Independent mpmath engine reproduces the same signature."""
    assert R.contact_signature_numeric(p, Mmax=Mmax) == expected


def test_lehmer_signature_is_phi1_to_the_10():
    """Ledger H / Theorem 7.9: deg Rat_L = 100 and the only contact is
    Phi_1^10 -- all 8 circle conjugates of Lehmer's number are mutually inert."""
    Rp = R.ratio_poly(R.LEHMER)
    assert Rp.degree() == 100
    assert R.contact_signature(Rp) == {1: 10}


def test_S8_signature_is_phi1_to_the_8():
    """Ledger P: deg Rat_{S8}=64, contacts {Phi_1^8}."""
    Rp = R.ratio_poly(R.S8)
    assert Rp.degree() == 64
    assert R.contact_signature(Rp) == {1: 8}


def test_gauge_blindness_x3_pm_2_identical():
    """rem:scope (ledger A,B): x^3-2 (absolute Z/6? no, Z/3) and x^3+2
    (absolute Z/6) return the *same* signature -- the probe reads only the
    common relational group Z/3, blind to the reference ray."""
    sa = R.contact_signature(R.ratio_poly(x**3 - 2))
    sb = R.contact_signature(R.ratio_poly(x**3 + 2))
    assert sa == sb == {1: 3, 3: 3}


def test_shell_signature_coarser_than_type_x4_quartics():
    """rem:scope (ledger E,F): x^4+5x^2+5 and K=x^4+5x^2-5 have identical
    contact signatures {Phi_1^4, Phi_2^4} while their full relational groups
    differ (Z/2 vs Z/4) -- the shell signature is strictly coarser."""
    se = R.contact_signature(R.ratio_poly(x**4 + 5 * x**2 + 5))
    sf = R.contact_signature(R.ratio_poly(x**4 + 5 * x**2 - 5))
    assert se == sf == {1: 4, 2: 4}
    # but the relational orders differ:
    assert R.relational_order(x**4 + 5 * x**2 + 5) == 2
    assert R.relational_order(x**4 + 5 * x**2 - 5) == 4


def test_beta4_and_lehmer_are_relationally_inert():
    """Theorem 7.9 (thm:inert): both keystones are inert -- the only contact of
    Rat_p is the diagonal Phi_1^n."""
    assert R.contact_signature(R.ratio_poly(R.B4)) == {1: 4}
    assert R.contact_signature(R.ratio_poly(R.LEHMER)) == {1: 10}


def test_twisted_shell_witness_is_non_inert():
    """Example 7.21 / ledger U: p = x^4+x^2+2 = q(x^2), q=x^2+x+2 irreducible
    with irrational angle, so p is inadmissible yet NON-inert: Rat_p carries a
    Phi_2 contact (the +-sqrt fibres form size-2 non-real coherence classes)."""
    import sympy as sp

    p = x**4 + x**2 + 2
    q = x**2 + x + 2
    assert sp.expand(q.subs(x, x**2) - p) == 0           # p = q(x^2)
    assert sp.Poly(p, x).is_irreducible
    assert R.contact_signature(R.ratio_poly(q)) == {1: 2}       # q angle irrational
    assert R.contact_signature(R.ratio_poly(p)) == {1: 4, 2: 4}  # non-inert
    assert sp.Poly(p, x).eval(0) == 2                    # p(0)=2 => Mahler 2
