"""
Object II -- the content polynomial x^4-1 (Prop 3.1 `prop:content`,
Rem 3.2 `rem:algebra-not-circle`).

Claims verified independently:
  * x^4-1 = Phi_1 Phi_2 Phi_4 = (x-1)(x+1)(x^2+1);
  * its roots sit exactly at the Z/4Z angle lattice {0,pi/2,pi,3pi/2};
  * realised on-circle roots are {+1,-1}: phi(x)phi contributes -1,
    phi^4(x)phi^4 contributes +1;
  * the charge-{1,3} sector Phi_4=x^2+1 is realised only OFF the circle,
    by K's place +-i*beta at modulus beta=2.4195 (Rem 3.2);
  * full Z/4Z in the charge, only Z/2Z={+-1} on the circle.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp
import ge_helpers as H
from ge_helpers import x, phi_seed, K_seed, tensor, sq, roots_mp, on_circle_roots


def test_content_polynomial_cyclotomic_factorization():
    """Prop 3.1: x^4-1 = Phi_1 Phi_2 Phi_4 = (x-1)(x+1)(x^2+1)."""
    Phi1, Phi2, Phi4 = x - 1, x + 1, x**2 + 1
    assert sp.expand((x**4 - 1) - Phi1 * Phi2 * Phi4) == 0
    # and these are the 1st, 2nd, 4th cyclotomic polynomials
    assert sp.cyclotomic_poly(1, x) == Phi1
    assert sp.cyclotomic_poly(2, x) == Phi2
    assert sp.cyclotomic_poly(4, x) == Phi4


def test_content_roots_sit_on_charge_lattice():
    """Prop 3.1 / Fig 1: the four roots of x^4-1 are the four charge directions
    {0, pi/2, pi, 3pi/2}."""
    roots = roots_mp(sp.Poly(x**4 - 1, x))
    args = sorted(float((mp.arg(r) % (2 * mp.pi)) / (mp.pi / 2)) for r in roots)
    # charges 0,1,2,3 each realised exactly once
    assert [round(a) for a in args] == [0, 1, 2, 3]
    for a in args:
        assert abs(a - round(a)) < 1e-20


def test_phi_tensor_phi_contributes_minus_one():
    """Prop 3.1: 'phi(x)phi contributes -1'.  -1 must be a root of phi(x)phi."""
    P = tensor(phi_seed(), phi_seed())
    assert P.eval(-1) == 0
    assert P.eval(1) != 0            # it does NOT contribute +1


def test_phi4_tensor_phi4_contributes_plus_one():
    """Prop 3.1: 'phi^4(x)phi^4 contributes +1'.  phi^4 = square twice."""
    phi4 = sq(sq(phi_seed()))
    P = tensor(phi4, phi4)
    assert P.eval(1) == 0
    assert P.eval(-1) != 0


def test_K_imaginary_place_modulus_beta():
    """Rem 3.2: the charge-{1,3} place of K is +-i*beta with
    beta = 2.4195 (!= 1, hence OFF the unit circle)."""
    imag_roots = [r for r in roots_mp(K_seed())
                  if abs(r.real) < mp.mpf(10)**(-20) and abs(r.imag) > 0]
    assert len(imag_roots) == 2
    beta = abs(imag_roots[0])
    beta_exact = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)     # sqrt of |(-5-3sqrt5)/2|
    assert abs(beta - beta_exact) < mp.mpf(10)**(-30)
    assert round(float(beta), 4) == 2.4195
    assert abs(beta - 1) > mp.mpf('0.1')               # genuinely off circle


def test_full_Z4_in_charge_only_Z2_on_circle():
    """Rem 3.2: 'the full Z/4Z lives in the charge; only Z/2Z={+-1} is realised
    on the circle' -- this is exactly why x^4-1 (not x^2-1) is the content poly.
    Sweep the realised objects: every modulus-1 root lies in {+1,-1}, while the
    charge group attained (via K) is the full {0,1,2,3}."""
    objects = [phi_seed(), K_seed(), tensor(phi_seed(), phi_seed()),
               tensor(sq(sq(phi_seed())), sq(sq(phi_seed()))), tensor(K_seed(), K_seed())]
    seen_on_circle = set()
    for P in objects:
        for r in on_circle_roots(P):
            # snap to nearest Gaussian-integer-ish unit root
            re, im = round(float(r.real)), round(float(r.imag))
            seen_on_circle.add((re, im))
    assert seen_on_circle <= {(1, 0), (-1, 0)}          # on-circle = Z/2Z
    assert (0, 1) not in seen_on_circle and (0, -1) not in seen_on_circle
    # but the charge lattice really is order 4 (K supplies 1 and 3)
    from ge_helpers import charges
    assert set(charges(K_seed())) == {0, 1, 2, 3}


def test_x2_minus_1_would_undercount():
    """Rem 3.2: 'x^2-1 would record only its on-circle shadow.'  x^2-1 has NO
    root at +-i, so it cannot host the charge-{1,3} sector; x^4-1 can."""
    assert sp.Poly(x**2 - 1, x).eval(sp.I) != 0
    assert sp.Poly(x**4 - 1, x).eval(sp.I) == 0
