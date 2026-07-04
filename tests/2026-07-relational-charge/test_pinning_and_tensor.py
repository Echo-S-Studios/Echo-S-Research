r"""Modulus pinning and tensor constructions: Theorem 6.15 (thm:pinning),
Example 6.20 (ex:x4quartic), Example 7.21 (ex:twistshell), Example 8.5
(ex:tensorsq), ledger S and X.

Paper (Theorem 6.15): an irreducible p with a uniquely attained root modulus
has no root-of-unity ratio between distinct roots.  Sharpness (Example 7.21):
x^4-2 (all roots one shell, modulus NOT unique) IS admissible with torsion
ratios zeta_4^j.

Paper (Example 6.20 / ledger X): p=x^4-x+1 is irreducible, non-reciprocal, no
real roots; charpoly(C_p (x) C_p) = S_6(x)^2 (x^4+2x^2-x+1); Rat_p contacts
{Phi_1^4}; Rat_{Rat_p} (deg 256) contacts {Phi_1^28}; Gal(p)=S_4;
M(x^4-x+1) = tau_{S_6}.

Paper (Example 8.5 / ledger S): the charpoly F of the 16x16 Kronecker square
C_{beta4} (x) C_{beta4} has (x-1)-multiplicity exactly 4, deg gcd(F,F')=7, and
exactly 3 distinct real roots, all positive.
"""
import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x


# ---- Theorem 6.15 pinning, sharpness, and Salem/Pisot corollaries ----------
@pytest.mark.parametrize(
    "name,p",
    [("beta4", R.B4), ("S6", R.S6), ("S8", R.S8), ("Lehmer", R.LEHMER),
     ("plastic", R.PLASTIC)],
)
def test_pinned_modulus_forbids_torsion_ratios(name, p):
    """Theorem 6.15 / Cor 6.16, 6.18: Salem and Pisot numbers have a uniquely
    attained (dominant) modulus, so no two distinct roots differ by a root of
    unity: off the real axis the ratio object has NO nontrivial cyclotomic
    contact (signature is Phi_1^n only)."""
    P = sp.Poly(p, x)
    n = P.degree()
    # dominant modulus is uniquely attained:
    mods = sorted((abs(r) for r in R.roots_mp(p)), reverse=True)
    assert mods[0] - mods[1] > mp.mpf(10) ** -12
    # hence inert:
    assert R.contact_signature(R.ratio_poly(p)) == {1: n}


def test_pinning_hypothesis_is_necessary_x4_minus_2():
    """Example 7.21 / ledger C: x^4-2 has ALL four roots on one shell (no
    unique modulus), and is irreducible with torsion ratios zeta_4^j -- so its
    ratio object carries Phi_2 and Phi_4 contacts.  Pinning does not apply."""
    p = x**4 - 2
    assert sp.Poly(p, x).is_irreducible
    mods = sorted(abs(r) for r in R.roots_mp(p))
    assert max(mods) - min(mods) < mp.mpf(10) ** -20  # all equal modulus
    assert R.contact_signature(R.ratio_poly(p)) == {1: 4, 2: 4, 4: 4}


def test_mixed_pinning_no_shared_modulus_between_distinct_salems():
    """Theorem 6.15 (mixed form) / Cor 6.17: distinct Salem numbers share no
    root modulus off the real axis, so their mixed ratio object has no
    cyclotomic contact.  (Checked here for beta4, S6 as a representative of the
    'no two Salem numbers circle-lock' corollary.)"""
    assert R.contact_signature(R.mixed_ratio_poly(R.B4, R.S6)) == {}


# ---- Example 6.20: the fully-rigid quartic x^4 - x + 1 ---------------------
def test_x4_minus_x_plus_1_basic_structure():
    """Example 6.20: x^4-x+1 is irreducible, non-reciprocal, with no real
    roots (0 real roots by Sturm)."""
    p = x**4 - x + 1
    P = sp.Poly(p, x)
    assert P.is_irreducible
    assert not R.is_reciprocal(p)
    assert P.count_roots(-sp.oo, sp.oo) == 0  # no real roots


def test_x4_minus_x_plus_1_kronecker_factorization():
    """Example 6.20 / ledger X: charpoly(C_p (x) C_p) = S_6^2 (x^4+2x^2-x+1)."""
    p = x**4 - x + 1
    C = R.companion_matrix(p)
    K = sp.Matrix(sp.kronecker_product(C, C))
    F = sp.factor(K.charpoly(x).as_expr())
    expected = (x**6 - x**4 - x**3 - x**2 + 1) ** 2 * (x**4 + 2 * x**2 - x + 1)
    assert sp.expand(F - expected) == 0
    # the second factor is irreducible (the psi^2 image):
    assert sp.Poly(x**4 + 2 * x**2 - x + 1, x).is_irreducible


def test_x4_minus_x_plus_1_shell_scan():
    """Example 6.20: Rat_p contacts {Phi_1^4} -- no within-shell torsion."""
    assert R.contact_signature(R.ratio_poly(x**4 - x + 1)) == {1: 4}


def test_x4_minus_x_plus_1_nu_scan_nested_ratio():
    """Example 6.20 / ledger X: the criterion verbatim -- Rat_{Rat_p} has
    degree 256 and complete signature {Phi_1^28} (multiplicity 16+4+8=28),
    certifying full relational inertness."""
    Rp = R.ratio_poly(x**4 - x + 1)
    RRp = R.ratio_poly(Rp.as_expr())
    assert RRp.degree() == 256
    assert R.phi1_multiplicity(RRp) == 28
    assert R.contact_signature(RRp) == {1: 28}


def test_x4_minus_x_plus_1_galois_is_S4():
    """Example 6.20: Gal(x^4-x+1) = S_4 (order 24; the only transitive group of
    degree 4 with order 24)."""
    g, _ = sp.Poly(x**4 - x + 1, x).galois_group()
    assert g.order() == 24
    assert g.is_transitive()


def test_x4_minus_x_plus_1_mahler_equals_smallest_deg6_salem():
    """Example 6.20: M(x^4-x+1) = tau_{S_6}, the smallest degree-6 Salem number
    (the dominant root of S_6)."""
    m_quartic = R.mahler_measure(x**4 - x + 1)
    tau_S6 = max(
        (r.real for r in R.roots_mp(R.S6) if abs(r.imag) < mp.mpf(10) ** -20)
    )
    assert abs(m_quartic - tau_S6) < mp.mpf(10) ** -25
    # it lies strictly between Smyth's constant and phi:
    theta0 = max((r.real for r in R.roots_mp(R.PLASTIC) if abs(r.imag) < mp.mpf(10) ** -20))
    phi = (1 + mp.sqrt(5)) / 2
    assert theta0 < m_quartic < phi


# ---- Example 8.5: beta4 (x) beta4 is NON-inert -----------------------------
def test_beta4_tensor_beta4_kronecker_structure():
    """Example 8.5 / ledger S: for F = charpoly(C_{beta4} (x) C_{beta4}):
    (x-1)-multiplicity = 4; deg gcd(F, F') = 7 (pattern 1^4 . four doubles);
    exactly 3 distinct real roots, all positive."""
    C = R.companion_matrix(R.B4)
    K = sp.Matrix(sp.kronecker_product(C, C))
    F = sp.Poly(K.charpoly(x).as_expr(), x)
    assert F.degree() == 16
    # (x-1) multiplicity:
    assert R.phi1_multiplicity(F) == 4
    # deg gcd(F, F') = 7:
    Fp = F.diff(x)
    assert sp.gcd(F, Fp).degree() == 7
    # squarefree part -> distinct real roots:
    sqfree = F.quo(sp.gcd(F, Fp))
    assert sqfree.count_roots(-sp.oo, sp.oo) == 3     # 3 distinct real roots
    assert sqfree.count_roots(0, sp.oo) == 3          # all positive
    assert sqfree.count_roots(-sp.oo, 0) == 0


def test_beta4_tensor_square_manufactures_rational_block():
    """Example 8.5: the tensor square produces a rational block including the
    diagonal products alpha*conj(alpha) = |alpha|^2 at angle 0 -- offset
    cancellation (Prop 4.9(iii)) creates coherence from an inert factor.

    The six positive reals {1,1,1,1, tau^2, tau^-2} are the rational block."""
    rts = R.roots_mp(R.B4)
    prods = [a * b for a in rts for b in rts]
    positive_reals = [
        p for p in prods
        if abs(p.imag) < mp.mpf(10) ** -18 and p.real > mp.mpf(10) ** -18
    ]
    assert len(positive_reals) == 6  # {1,1,1,1, tau^2, tau^-2}
    ones = [p for p in positive_reals if abs(p - 1) < mp.mpf(10) ** -15]
    assert len(ones) == 4
