"""
Object IV -- the graded normal form (Prop 5.1 `prop:normal` and its table).

The three diagnostic factorizations, rebuilt independently from the operators
and factored over Z:
  * phi(x)phi      = (x+1)^2 (x^2-3x+1),                 grow M = phi^2
  * phi^4(x)phi^4  = (x-1)^2 (x^2-47x+1),                grow M = phi^8 = 46.98
  * K(x)K          = (x^2+5)^4 (x^2-5x-5)^2 (x^2+5x-5)^2, grow M = 5.854
Plus: on-circle part divides x^4-1; growth lives in the off-circle real factor;
and the identity a*b = sqrt5 that produces the (x^2+5) imaginary sector.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp
import ge_helpers as H
from ge_helpers import x, phi_seed, K_seed, tensor, sq, mahler, on_circle_roots, is_root_of_unity, PHI

PHI_V = PHI()


def _factor_dict(P):
    """{sympy factor expr: multiplicity} for a monic integer poly."""
    c, facs = sp.factor_list(P.as_expr(), x)
    return {sp.expand(f): m for f, m in facs}


def test_normal_form_phi_tensor_phi():
    """Prop 5.1 table row 1: phi(x)phi = (x+1)^2 (x^2-3x+1),
    on-circle Phi_2^2 (q=2) x grow M=phi^2."""
    P = tensor(phi_seed(), phi_seed())
    assert _factor_dict(P) == {x + 1: 2, x**2 - 3 * x + 1: 1}
    assert abs(mahler(sp.Poly(x**2 - 3 * x + 1, x)) - PHI_V**2) < mp.mpf(10)**(-30)
    assert abs(mahler(P) - PHI_V**2) < mp.mpf(10)**(-30)


def test_normal_form_phi4_tensor_phi4():
    """Prop 5.1 table row 2: phi^4(x)phi^4 = (x-1)^2 (x^2-47x+1),
    on-circle Phi_1^2 (q=0) x grow M=46.98."""
    phi4 = sq(sq(phi_seed()))
    P = tensor(phi4, phi4)
    assert _factor_dict(P) == {x - 1: 2, x**2 - 47 * x + 1: 1}
    grow_M = mahler(sp.Poly(x**2 - 47 * x + 1, x))
    assert abs(grow_M - PHI_V**8) < mp.mpf(10)**(-25)
    assert round(float(grow_M), 2) == 46.98


def test_normal_form_K_tensor_K():
    """Prop 5.1 table row 3: K(x)K = (x^2+5)^4 (x^2-5x-5)^2 (x^2+5x-5)^2,
    imaginary q=1,3 sector x grow M=5.854."""
    P = tensor(K_seed(), K_seed())
    assert _factor_dict(P) == {x**2 + 5: 4, x**2 - 5 * x - 5: 2, x**2 + 5 * x - 5: 2}
    grow_M = mahler(sp.Poly(x**2 - 5 * x - 5, x))
    assert abs(grow_M - (5 + 3 * mp.sqrt(5)) / 2) < mp.mpf(10)**(-30)
    assert round(float(grow_M), 3) == 5.854


def test_ab_equals_sqrt5_imaginary_sector():
    """Prop 5.1 discussion: the K-sector (x^2+5) has roots +-i*sqrt5, i.e. the
    product of K's real half-modulus a and imaginary half-modulus b satisfies
    a*b = sqrt5.  With a^2=(-5+3sqrt5)/2, b^2=(5+3sqrt5)/2 -> a^2 b^2 = 5."""
    r5 = sp.sqrt(5)
    a2 = (-5 + 3 * r5) / 2
    b2 = (5 + 3 * r5) / 2
    assert sp.simplify(a2 * b2 - 5) == 0
    # and x^2+5 indeed has roots +-i sqrt5
    assert sp.simplify(sp.Poly(x**2 + 5, x).eval(sp.I * r5)) == 0


def test_on_circle_part_divides_x4_minus_1_and_growth_off_circle():
    """Prop 5.1: 'in all cases the on-circle part divides x^4-1 and the growth
    lives entirely in the off-circle factor.'  For each object: every modulus-1
    root is a root of x^4-1 (a root of unity), and each factor carrying Mahler
    growth ( M>1 ) has all its roots off the unit circle."""
    phi4 = sq(sq(phi_seed()))
    for P in [tensor(phi_seed(), phi_seed()), tensor(phi4, phi4), tensor(K_seed(), K_seed())]:
        # on-circle roots are 4th roots of unity (roots of x^4-1)
        for r in on_circle_roots(P):
            assert is_root_of_unity(r)
            assert abs(complex(r)**4 - 1) < 1e-18
        # growth factors are off-circle
        for f, m in sp.factor_list(P.as_expr(), x)[1]:
            fp = sp.Poly(f, x)
            if mahler(fp) > 1 + mp.mpf(10)**(-25):
                for rr in H.roots_mp(fp):
                    assert abs(abs(rr) - 1) > mp.mpf(10)**(-6)


def test_grow_factor_measure_equals_object_measure():
    """Prop 5.1: 'M(G) = M(P)' -- the whole Mahler measure sits in the off-circle
    grow factor because the on-circle part has measure 1."""
    phi4 = sq(sq(phi_seed()))
    for P, grow in [(tensor(phi_seed(), phi_seed()), sp.Poly(x**2 - 3 * x + 1, x)),
                    (tensor(phi4, phi4), sp.Poly(x**2 - 47 * x + 1, x))]:
        assert abs(mahler(P) - mahler(grow)) < mp.mpf(10)**(-25)
