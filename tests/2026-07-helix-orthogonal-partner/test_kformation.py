"""
Independent verification of Section 5 ("The K-formation astride the fold") of
"The Dissolved Helix and Its Orthogonal Partner".

Covers Proposition 5.1 and eqs. (6)-(8) [eq:ypm, eq:Kval, eq:MK]:
Kf(x)=x^4+5x^2-5; inner quadratic y^2+5y-5 with roots y_+ , y_- straddling 0;
real roots +/-K with K = 5^{1/4}/phi; imaginary roots +/- i beta;
Mah(Kf) = beta^2 = (5+3 sqrt5)/2 = phi^2 sqrt5.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

x = sp.symbols('x')
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
Kf = x**4 + 5 * x**2 - 5


# ---------------------------------------------------------------- eq. (6)
def test_inner_quadratic_roots_eq6():
    """Eq. (6): y = x^2 gives y^2+5y-5, roots y_+ = (-5+3 sqrt5)/2,
    y_- = (-5-3 sqrt5)/2  (sqrt(25+20)=sqrt45=3 sqrt5)."""
    y = sp.symbols('y')
    roots = sp.solve(y**2 + 5 * y - 5, y)
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2
    assert {sp.simplify(r) for r in roots} == {sp.simplify(yp), sp.simplify(ym)}
    assert sp.simplify(sp.sqrt(45) - 3 * sqrt5) == 0


def test_inner_roots_straddle_zero():
    """Eq. (6): y_+ > 0 > y_-  (sign of 3 sqrt5 - 5 > 0)."""
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2
    assert (yp > 0) == True         # noqa: E712
    assert (ym < 0) == True         # noqa: E712
    assert (3 * sqrt5 - 5 > 0) == True   # noqa: E712


def test_inner_roots_numeric():
    """Eq. (6): y_+ ~ 0.8541, y_- ~ -5.8541 (paper's displayed values)."""
    yp = (-5 + 3 * mp.sqrt(5)) / 2
    ym = (-5 - 3 * mp.sqrt(5)) / 2
    assert abs(yp - mp.mpf('0.8541')) < 1e-4
    assert abs(ym - mp.mpf('-5.8541')) < 1e-4


# ---------------------------------------------------------------- eq. (7)
def test_K_equals_fifth_root_over_phi_eq7():
    """Eq. (7): K = sqrt(y_+) = 5^{1/4}/phi (exact, via rationalisation
    (5^{1/4}/phi)^2 = sqrt5/phi^2 = (3 sqrt5 -5)/2 = y_+)."""
    yp = (-5 + 3 * sqrt5) / 2
    K_candidate = 5**sp.Rational(1, 4) / phi
    assert sp.simplify(K_candidate**2 - yp) == 0
    # sqrt5/phi^2 rationalises to (3 sqrt5 - 5)/2
    assert sp.simplify(sqrt5 / phi**2 - (3 * sqrt5 - 5) / 2) == 0
    assert sp.simplify(K_candidate - sp.sqrt(yp)) == 0


def test_K_is_root_of_Kf():
    """Eq. (7): the minimal polynomial of 5^{1/4}/phi is exactly Kf
    (independent confirmation that K is a real root of Kf)."""
    mp_poly = sp.minimal_polynomial(5**sp.Rational(1, 4) / phi, x)
    assert sp.simplify(mp_poly - Kf) == 0


def test_K_numeric():
    """Eq. (7): K ~ 0.9242 and K < 1."""
    yp = (-5 + 3 * mp.sqrt(5)) / 2
    K = mp.sqrt(yp)
    assert abs(K - mp.mpf('0.9242')) < 1e-4
    assert K < 1


# ---------------------------------------------------------------- beta / eq. (8)
def test_beta_value_and_numeric():
    """Prop. 5.1: beta = sqrt(|y_-|), beta ~ 2.4195, beta > 1."""
    ym = (-5 - 3 * sqrt5) / 2
    beta = sp.sqrt(sp.Abs(ym))
    assert sp.simplify(beta**2 - (5 + 3 * sqrt5) / 2) == 0
    bnum = mp.sqrt(abs((-5 - 3 * mp.sqrt(5)) / 2))
    assert abs(bnum - mp.mpf('2.4195')) < 1e-4
    assert bnum > 1


def test_mahler_Kf_eq8():
    """Eq. (8): Mah(Kf) = beta^2 = (5+3 sqrt5)/2 = phi^2 sqrt5 ~ 5.8541.
    Roots are +/-K (|K|<1) and +/- i beta (|.|=beta>1), so Mahler = beta^2."""
    beta2 = (5 + 3 * sqrt5) / 2
    assert sp.simplify(beta2 - phi**2 * sqrt5) == 0
    # independent: Mahler = product of |root|>1 over the four roots, high precision
    coeffs = [mp.mpf(str(c)) for c in sp.Poly(Kf, x).all_coeffs()]
    rts = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    mah = mp.mpf(1)
    for r in rts:
        if abs(r) > 1:
            mah *= abs(r)
    assert abs(mah - mp.mpf(str(sp.N(beta2, 40)))) < 1e-30
    assert abs(mp.mpf(str(sp.N(beta2, 40))) - mp.mpf('5.8541')) < 1e-4


def test_K_lt_one_lt_beta():
    """Prop. 5.1 proof: K < 1 < beta, which is why Mahler collapses to beta^2."""
    K = mp.sqrt((-5 + 3 * mp.sqrt(5)) / 2)
    beta = mp.sqrt(abs((-5 - 3 * mp.sqrt(5)) / 2))
    assert K < 1 < beta


def test_root_structure_real_and_imaginary():
    """Prop. 5.1: Kf has two real roots (+/-K) and two purely imaginary
    roots (+/- i beta)."""
    roots = sp.solve(Kf, x)                       # explicit radical forms
    assert len(roots) == 4
    reals = [r for r in roots if sp.im(r) == 0]
    imags = [r for r in roots if sp.re(r) == 0 and sp.im(r) != 0]
    assert len(reals) == 2
    assert len(imags) == 2
    # each pair is closed under negation (+/-K and +/- i beta)
    assert sp.simplify(sum(reals)) == 0
    assert sp.simplify(sum(imags)) == 0
    # imaginary radii satisfy beta^2 = (5+3 sqrt5)/2
    for r in imags:
        assert sp.simplify(sp.im(r)**2 - (5 + 3 * sqrt5) / 2) == 0
