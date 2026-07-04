"""The lattice walls: the four directions (pi/2)Z.  Section 4 (Lemmas 4.1-4.4)
and the Box exclusion (Def. 4.5, Thms 4.6-4.7).  Angles are handled exactly
(a root is real or purely imaginary) and the (pi/2)Z group facts are checked in
the exact model (pi/2)Z ~ Z (units of pi/2).
"""
import mpmath as mp
import sympy as sp

from _helpers import (mahler_product, root_classification, is_salem, PHI,
                      to_ab_sqrt5, sign_ab_sqrt5)

mp.mp.dps = 50
_PHI = (1 + mp.sqrt(5)) / 2

# catalog minimal polynomials (Def. 2.9)
CATALOG = {
    "phi":   [1, -1, -1],
    "tau":   [1, 1, -1],
    "sqrt2": [1, 0, -2],
    "sqrt3": [1, 0, -3],
    "sqrt5": [1, 0, -5],
    "gap":   [1, -7, 1],
    "K":     [1, 0, 5, 0, -5],
}


def _arg_in_halfpi_lattice(root, tol=mp.mpf(10) ** -25):
    """True iff arg(root) in {0, pi/2, pi, 3pi/2}, i.e. root is real OR purely
    imaginary (exact characterisation of the (pi/2)Z directions)."""
    return abs(mp.re(root)) < tol or abs(mp.im(root)) < tol


def test_catalog_arguments_in_halfpi_Z():
    """Lemma 4.1 (P2-ANG-01): every eigenvalue of every catalog matrix has
    argument in (pi/2)Z = {0, pi/2, pi, 3pi/2}."""
    for name, coeffs in CATALOG.items():
        roots = mp.polyroots([mp.mpf(c) for c in coeffs], extraprec=120)
        for r in roots:
            assert _arg_in_halfpi_lattice(r), (name, mp.nstr(r, 10))


def test_quartic_K_x2_split():
    """Lemma 4.1 proof: x^4+5x^2-5 gives x^2 = (-5 +/- 3 sqrt5)/2; the upper
    sign is +0.854.. (real pair +/-K), the lower is -5.854.. (imag pair
    +/- i beta).  Both real |K| and |beta| lie off the unit circle."""
    yp = (-5 + 3 * mp.sqrt(5)) / 2
    ym = (-5 - 3 * mp.sqrt(5)) / 2
    # these solve y^2 + 5y - 5 = 0 (the x^2 quadratic)
    assert abs(yp ** 2 + 5 * yp - 5) < mp.mpf(10) ** -30
    assert abs(ym ** 2 + 5 * ym - 5) < mp.mpf(10) ** -30
    assert yp > 0 and abs(yp - mp.mpf("0.854101966")) < mp.mpf(10) ** -8
    assert ym < 0 and abs(ym - mp.mpf("-5.854101966")) < mp.mpf(10) ** -8


def test_halfpi_lattice_closed_under_add_and_double_not_half():
    """Lemma 4.2 (P2-ANG-02): in units where pi/2 = 1, the lattice is Z (mod 4).
    Closed under addition and doubling; NOT under halving:
    (1/2)(pi/2) = pi/4 is not a multiple of pi/2."""
    lattice = [sp.Integer(k) for k in range(4)]           # {0,1,2,3} * (pi/2)
    # closed under addition (subgroup of Z/4)
    for a in lattice:
        for b in lattice:
            assert (a + b) % 4 in [x % 4 for x in lattice]
    # closed under doubling
    for a in lattice:
        assert (2 * a) % 4 in [x % 4 for x in lattice]
    # NOT closed under halving: half of pi/2 (a=1) is 1/2, not an integer
    half = sp.Rational(1, 2)
    assert half not in lattice            # pi/4 is off the lattice


def test_oncircle_halfpi_are_fourth_roots_of_unity():
    """Lemma 4.3 (P2-ANG-03): |z|=1 and arg z in (pi/2)Z  =>  z in {1,i,-1,-i},
    and each is a fourth root of unity (z^4 = 1)."""
    mu4 = [mp.mpc(1, 0), mp.mpc(0, 1), mp.mpc(-1, 0), mp.mpc(0, -1)]
    for k in range(4):
        z = mp.e ** (1j * (mp.pi / 2) * k)
        assert min(abs(z - w) for w in mu4) < mp.mpf(10) ** -30
        assert abs(z ** 4 - 1) < mp.mpf(10) ** -30


def test_salem_oncircle_conjugates_not_roots_of_unity():
    """Lemma 4.4 (P2-ANG-04): a Salem number's on-circle conjugates are not
    roots of unity.  Instance: Lehmer L is irreducible and non-cyclotomic, so no
    root (in particular no on-circle root) is a root of unity; and those
    on-circle roots carry argument NOT in (pi/2)Z."""
    x = sp.symbols('x')
    L = sp.Poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], x)
    # irreducible over Q
    assert L.as_expr().is_polynomial()
    assert sp.factor_list(L.as_expr())[1].__len__() == 1
    assert sp.Poly(sp.factor(L.as_expr()), x).degree() == 10
    # not equal to any cyclotomic of degree 10 (phi(n)=10 -> n in {11,22})
    for n in (11, 22):
        assert sp.expand(L.as_expr() - sp.cyclotomic_poly(n, x)) != 0
    # on-circle roots have argument off (pi/2)Z
    roots = mp.polyroots([mp.mpf(c) for c in L.all_coeffs()],
                         maxsteps=400, extraprec=400)
    oncirc = [r for r in roots if abs(abs(r) - 1) < mp.mpf(10) ** -12]
    assert len(oncirc) == 8
    for r in oncirc:
        # not real and not purely imaginary  => arg not in (pi/2)Z
        assert abs(mp.re(r)) > mp.mpf(10) ** -6
        assert abs(mp.im(r)) > mp.mpf(10) ** -6


def _in_box(coeffs):
    """Box membership (Def. 4.5): Mah in {1} U [phi, inf) AND every on-circle
    root has argument in (pi/2)Z (i.e. is real-or-imaginary on the circle)."""
    m = mahler_product(coeffs)
    floor_ok = (abs(m - 1) < mp.mpf(10) ** -12) or (m >= _PHI - mp.mpf(10) ** -12)
    roots = mp.polyroots([mp.mpf(c) for c in coeffs], maxsteps=400, extraprec=400)
    oncirc = [r for r in roots if abs(abs(r) - 1) < mp.mpf(10) ** -10]
    lattice_ok = all(_arg_in_halfpi_lattice(r) for r in oncirc)
    return floor_ok and lattice_ok


def test_catalog_lies_in_box():
    """Thm 4.6 (Emission subset Box): every catalog polynomial is in the box
    (floor satisfied; catalog has no on-circle roots so lattice wall holds)."""
    for name, coeffs in CATALOG.items():
        assert _in_box(coeffs), name


def test_lehmer_outside_box_both_walls():
    """Thm 4.7: Lehmer's number is outside the box on BOTH counts -- its measure
    lies in the forbidden strip (< phi) and its on-circle conjugates are off the
    lattice."""
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    assert mahler_product(L) < _PHI          # floor wall violated
    assert not _in_box(L)


def test_beta4_above_floor_but_outside_box_by_lattice():
    """Thm 4.7 / Fig.: the minimal degree-four Salem beta_4 is ABOVE the floor
    (beta_4 > phi) yet still outside the box, excluded by the lattice wall
    (its on-circle conjugate has irrational argument)."""
    B = [1, -1, -1, -1, 1]
    assert is_salem(B)
    assert mahler_product(B) > _PHI          # floor wall satisfied
    # but on-circle pair is off-lattice -> not in box
    roots = mp.polyroots([mp.mpf(c) for c in B], extraprec=120)
    oncirc = [r for r in roots if abs(abs(r) - 1) < mp.mpf(10) ** -12]
    assert len(oncirc) == 2
    for r in oncirc:
        assert abs(mp.re(r)) > mp.mpf(10) ** -6 and abs(mp.im(r)) > mp.mpf(10) ** -6
    assert not _in_box(B)
