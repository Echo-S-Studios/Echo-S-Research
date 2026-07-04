"""
Independent verification of Section 12: "The K-formation seed astride the fold".
Paper: "The Exchange Rate lambda = 2c ...".

Claims (prop:kform):
  f(x) = x^4+5x^2-5, as g(y)=y^2+5y-5 in y=x^2:
    y_+ = (-5+3 sqrt5)/2 ~ 0.8541 > 0,  y_- = (-5-3 sqrt5)/2 ~ -5.8541 < 0.
  real roots +-K with K = sqrt(y_+) = 5^{1/4}/phi ~ 0.9242 (K<1, inside circle);
  imaginary roots +- i beta with beta = sqrt(|y_-|) ~ 2.4195 (outside circle);
  Mah(f) = beta^2 = (5+3 sqrt5)/2 ~ 5.8541.
Also (ssec:kformface): |i beta| = 2.4195 and |5^{1/4} i| = 1.4953, both >> 1.
"""
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 45
x, y = sp.symbols('x y')
phi = (1 + sp.sqrt(5)) / 2


def test_kformation_y_roots_straddle_zero():
    """prop:kform: 'y_+ = (-5+3 sqrt5)/2 > 0 and y_- = (-5-3 sqrt5)/2 < 0'."""
    ys = sp.solve(y**2 + 5 * y - 5, y)
    yp = (-5 + 3 * sp.sqrt(5)) / 2
    ym = (-5 - 3 * sp.sqrt(5)) / 2
    got = {sp.simplify(r) for r in ys}
    assert any(sp.simplify(r - yp) == 0 for r in got)
    assert any(sp.simplify(r - ym) == 0 for r in got)
    assert sp.N(yp) > 0
    assert sp.N(ym) < 0
    # numeric values quoted by the paper
    assert abs(float(sp.N(yp)) - 0.8541019662) < 1e-9
    assert abs(float(sp.N(ym)) - (-5.8541019662)) < 1e-9


def test_K_equals_fifth_root_over_phi():
    """prop:kform: 'K = sqrt(y_+) = 5^{1/4}/phi ~ 0.9242'."""
    yp = (-5 + 3 * sp.sqrt(5)) / 2
    K = sp.sqrt(yp)
    K_alt = 5**sp.Rational(1, 4) / phi
    assert sp.simplify(K**2 - K_alt**2) == 0            # both positive -> equal
    assert abs(float(sp.N(K_alt)) - 0.9241763718) < 1e-9
    assert float(sp.N(K_alt)) < 1                        # inside the unit circle


def test_beta_and_mahler():
    """prop:kform: 'beta = sqrt(|y_-|) ~ 2.4195, Mah(f) = beta^2 = (5+3 sqrt5)/2'.
    Verify Mahler(x^4+5x^2-5) independently from its roots.
    """
    beta2 = (5 + 3 * sp.sqrt(5)) / 2
    ym = (-5 - 3 * sp.sqrt(5)) / 2
    assert sp.simplify(beta2 - (-ym)) == 0
    beta = sp.sqrt(beta2)
    assert abs(float(sp.N(beta)) - 2.4195251530) < 1e-9
    # independent Mahler measure of the quartic
    roots = mp.polyroots([1, 0, 5, 0, -5], maxsteps=200, extraprec=200)
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    assert abs(M - mp.mpf(str(sp.N(beta2, 40)))) < mp.mpf(10) ** (-20)


def test_kform_complex_place_off_circle():
    """ssec:kformface: 'the one Lorentzian field places its complex embedding
    nowhere near the unit circle': |i beta| = 2.4195 and |5^{1/4} i| = 1.4953, both >> 1.
    """
    beta = mp.sqrt(mp.mpf(str(sp.N((5 + 3 * sp.sqrt(5)) / 2, 40))))
    assert abs(abs(mp.mpc(0, 1) * beta) - mp.mpf('2.4195251530')) < mp.mpf('1e-9')
    fifth = mp.power(5, mp.mpf(1) / 4)
    assert abs(abs(mp.mpc(0, 1) * fifth) - mp.mpf('1.4953487812')) < mp.mpf('1e-9')
    assert fifth > 1


def test_kform_real_and_imaginary_root_structure():
    """prop:kform: f(x)=x^4+5x^2-5 has real roots +-K (K<1) and imaginary roots
    +- i beta (beta>1): i.e. the real pair is inside, the imaginary pair outside
    the unit circle.
    """
    roots = np.roots([1, 0, 5, 0, -5])
    reals = sorted(r.real for r in roots if abs(r.imag) < 1e-9)
    imags = sorted(r.imag for r in roots if abs(r.real) < 1e-9)
    assert len(reals) == 2 and len(imags) == 2
    assert abs(reals[1] - 0.9241763718) < 1e-8
    assert abs(reals[0] + 0.9241763718) < 1e-8
    assert abs(imags[1] - 2.4195251530) < 1e-8
    assert abs(imags[0] + 2.4195251530) < 1e-8
    assert abs(reals[1]) < 1 < abs(imags[1])
