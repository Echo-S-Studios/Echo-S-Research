"""
Independent verification of Section 6 ("The orthogonal partner completes Z/4Z")
of "The Dissolved Helix and Its Orthogonal Partner".

Covers Proposition 6.1 (chi(Kf) = {0,1,2,3}; imaginary roots have Re=0) and
Proposition 6.2 (the parity criterion): evenness of Kf vs. the non-even registry
quartics cons, res; their Q(sqrt5) factorisations; complex-root real parts
phi^2 and -phi; the shared field Q(5^{1/4}); signature (2,1); Galois group D_4
(order 8); non-reciprocity; and the Mahler measures
Mah(Kf)=phi^2 sqrt5, Mah(cons)=2 phi^5=11+5 sqrt5, Mah(res)=12+19 phi=(43+19 sqrt5)/2.
"""
import sympy as sp
import mpmath as mp
from sympy import galois_group
from sympy.polys.numberfields.subfield import field_isomorphism

mp.mp.dps = 60

x = sp.symbols('x')
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2

Kf = x**4 + 5 * x**2 - 5
cons = x**4 - 6 * x**3 + 26 * x**2 - 16 * x - 4
res = x**4 + 2 * x**3 + 39 * x**2 - 52 * x + 11
beta = sp.sqrt((5 + 3 * sqrt5) / 2)
K = 5**sp.Rational(1, 4) / phi


def chi_of(z):
    """Def. 2.4 Character II: round(2 arg(z)/pi) mod 4."""
    a = sp.arg(z)
    return int(sp.floor(2 * a / sp.pi + sp.Rational(1, 2)) % 4)


def mahler_numeric(poly):
    """Mahler measure = |lead| * prod max(1,|root|), computed at high precision."""
    P = sp.Poly(poly, x)
    coeffs = [mp.mpf(str(c)) for c in P.all_coeffs()]
    rts = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    m = abs(coeffs[0])
    for r in rts:
        if abs(r) > 1:
            m *= abs(r)
    return m


# ---------------------------------------------------------------- Prop. 6.1
def test_chi_Kf_is_full_group():
    """Prop. 6.1: chi(Kf) = {0,1,2,3}. Real roots +/-K -> {0,2},
    imaginary roots +/- i beta -> {1,3}."""
    charges = {chi_of(K), chi_of(-K), chi_of(sp.I * beta), chi_of(-sp.I * beta)}
    assert charges == {0, 1, 2, 3}


def test_chi_rounding_at_cardinal_points():
    """Prop. 6.1 proof: round(2*0/pi)=0, round(2*pi/pi)=2,
    round(2*(pi/2)/pi)=1, round(2*(3pi/2)/pi)=3."""
    assert chi_of(sp.Integer(1)) == 0                      # arg 0
    assert chi_of(sp.Integer(-1)) == 2                     # arg pi
    assert chi_of(sp.I) == 1                               # arg pi/2
    assert chi_of(-sp.I) == 3                              # arg -pi/2 ~ 3pi/2


def test_imaginary_roots_have_zero_real_part():
    """Prop. 6.1: Kf's complex roots are purely imaginary (Re = 0)."""
    roots = sp.Poly(Kf, x).all_roots()
    complex_roots = [r for r in roots if sp.im(r) != 0]
    assert len(complex_roots) == 2
    for r in complex_roots:
        assert sp.re(r) == 0


# ---------------------------------------------------------------- Prop. 6.2
def test_Kf_is_even_registry_are_not():
    """Prop. 6.2: Kf is even (poly in x^2); cons and res are non-even."""
    def is_even(p):
        return sp.simplify(p.subs(x, -x) - p) == 0
    assert is_even(Kf)
    assert not is_even(cons)
    assert not is_even(res)


def test_cons_factorisation_over_Qsqrt5():
    """Prop. 6.2: cons = (x^2-(3+sqrt5)x+(11+5 sqrt5))(x^2-(3-sqrt5)x+(11-5 sqrt5))."""
    stated = ((x**2 - (3 + sqrt5) * x + (11 + 5 * sqrt5)) *
              (x**2 - (3 - sqrt5) * x + (11 - 5 * sqrt5)))
    assert sp.expand(stated - cons) == 0


def test_res_factorisation_over_Qsqrt5():
    """Prop. 6.2: res = (x^2+(1-sqrt5)x+(43-19 sqrt5)/2)(x^2+(1+sqrt5)x+(43+19 sqrt5)/2)."""
    stated = ((x**2 + (1 - sqrt5) * x + (43 - 19 * sqrt5) / 2) *
              (x**2 + (1 + sqrt5) * x + (43 + 19 * sqrt5) / 2))
    assert sp.expand(stated - res) == 0


def test_cons_complex_block_real_part():
    """Prop. 6.2: cons's complex block has real part phi^2 = (3+sqrt5)/2.
    (The block x^2-(3+sqrt5)x+(11+5 sqrt5) has negative discriminant.)"""
    b = 3 + sqrt5
    c = 11 + 5 * sqrt5
    disc = b**2 - 4 * c
    assert (disc < 0) == True                              # noqa: E712 (complex)
    real_part = b / 2
    assert sp.simplify(real_part - phi**2) == 0
    assert sp.simplify(phi**2 - (3 + sqrt5) / 2) == 0
    # the conjugate block is the real one
    b2, c2 = 3 - sqrt5, 11 - 5 * sqrt5
    assert (b2**2 - 4 * c2 > 0) == True                    # noqa: E712 (real)


def test_res_complex_block_real_part():
    """Prop. 6.2: res's complex block has real part -phi = -(1+sqrt5)/2."""
    a1 = 1 + sqrt5
    c = (43 + 19 * sqrt5) / 2
    disc = a1**2 - 4 * c
    assert (disc < 0) == True                              # noqa: E712 (complex)
    real_part = -a1 / 2
    assert sp.simplify(real_part - (-phi)) == 0
    # the other block is real
    a2, c2 = 1 - sqrt5, (43 - 19 * sqrt5) / 2
    assert (a2**2 - 4 * c2 > 0) == True                    # noqa: E712 (real)


def test_all_three_generate_Q_fifth_root():
    """Prop. 6.2: Kf, cons, res each generate Q(5^{1/4}).
    For Kf: 5^{1/4} = K(K^2+4)/3 lies in Q(K). For cons, res: their root fields
    are isomorphic to Q(5^{1/4}) (field_isomorphism returns an embedding)."""
    # Kf: 5^{1/4} expressible in Q(K)
    assert sp.simplify(K * (K**2 + 4) / 3 - 5**sp.Rational(1, 4)) == 0
    target = sp.AlgebraicNumber(5**sp.Rational(1, 4))
    for p in (cons, res):
        assert sp.Poly(p, x).is_irreducible
        r = sp.AlgebraicNumber(sp.CRootOf(p, 0))
        assert field_isomorphism(r, target) is not None


def test_signature_2_1():
    """Prop. 6.2: each quartic has signature (2,1): 2 real roots, 1 complex pair."""
    for p in (Kf, cons, res):
        n_real = sp.Poly(p, x).count_roots(-sp.oo, sp.oo)
        assert n_real == 2                                 # r1 = 2
        assert (4 - n_real) // 2 == 1                      # r2 = 1


def test_galois_group_D4_order8():
    """Prop. 6.2: each quartic has Galois group D_4 (order 8). Among transitive
    subgroups of S_4 the unique one of order 8 is D_4."""
    for p in (Kf, cons, res):
        G, _ = galois_group(sp.Poly(p, x))
        assert G.order() == 8
        assert G.is_transitive()


def test_non_reciprocal():
    """Prop. 6.2: none of Kf, cons, res is reciprocal (palindromic up to sign),
    hence none is a Salem number."""
    def reciprocal(p):
        c = sp.Poly(p, x).all_coeffs()
        rc = c[::-1]
        return c == rc or c == [-v for v in rc]
    for p in (Kf, cons, res):
        assert not reciprocal(p)


# ---------------------------------------------------------------- Mahler measures
def test_mahler_Kf_closed_form():
    """Prop. 6.2 / eq. (8): Mah(Kf) = phi^2 sqrt5."""
    target = phi**2 * sqrt5
    assert abs(mahler_numeric(Kf) - mp.mpf(str(sp.N(target, 50)))) < 1e-30


def test_mahler_cons_closed_form():
    """Prop. 6.2: Mah(cons) = 2 phi^5 = 11 + 5 sqrt5."""
    assert sp.simplify(2 * phi**5 - (11 + 5 * sqrt5)) == 0        # identity
    target = 11 + 5 * sqrt5
    assert abs(mahler_numeric(cons) - mp.mpf(str(sp.N(target, 50)))) < 1e-25


def test_mahler_res_closed_form():
    """Prop. 6.2: Mah(res) = 12 + 19 phi = (43 + 19 sqrt5)/2."""
    assert sp.simplify((12 + 19 * phi) - (43 + 19 * sqrt5) / 2) == 0  # identity
    target = (43 + 19 * sqrt5) / 2
    assert abs(mahler_numeric(res) - mp.mpf(str(sp.N(target, 50)))) < 1e-25


def test_golden_field_identities():
    """Prop. 6.2 supporting identities: phi^5 = (11+5 sqrt5)/2 and
    2 phi^5 = 11+5 sqrt5; 12+19 phi = (43+19 sqrt5)/2."""
    assert sp.simplify(phi**5 - (11 + 5 * sqrt5) / 2) == 0
    assert sp.simplify(2 * phi**5 - (11 + 5 * sqrt5)) == 0
    assert sp.simplify((12 + 19 * phi) - (43 + 19 * sqrt5) / 2) == 0
