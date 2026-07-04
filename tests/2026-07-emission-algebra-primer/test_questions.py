"""Part 6 -- The three questions: self-coupling, the flip, the partner (Sec. 6).

Independent re-derivations of the three self-coupling obstructions, the flip
discriminant/Gram-determinant identity, and the orthogonal partner
K(x)=x^4+5x^2-5 (root split, magnitude, charge completion, parity criterion).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from _eap_helpers import (R, I2, phi, psi, sqrt5, is_zero, mat_eq, comm,
                      charge_one, charge_set, outside_unit, mahler)


# ----------------------------------------------------------------------------
# Prop 6.1 -- Obstruction 1: the generative grade vanishes
# ----------------------------------------------------------------------------
def test_obstruction1_characteristic_operator_annihilates():
    """Prop 6.1: Phi_R(R)=R^2-Tr(R)R+det(R)I = R^2-R-I = 0, and [R,R]=0."""
    PhiR_at_R = R * R - sp.trace(R) * R + R.det() * I2
    assert mat_eq(PhiR_at_R, sp.zeros(2))
    assert mat_eq(comm(R, R), sp.zeros(2))


# ----------------------------------------------------------------------------
# Prop 6.2 -- Obstruction 2: the charge sublattice is closed
# ----------------------------------------------------------------------------
def test_obstruction2_charge_sublattice_closed():
    """Prop 6.2: {0,2} closed under union, sum mod4, doubling; never hits {1,3}."""
    G = {0, 2}
    assert G | G == {0, 2}
    assert {(a + b) % 4 for a in G for b in G} == {0, 2}
    assert {(2 * a) % 4 for a in G} == {0}
    # complement {1,3} is unreachable
    assert ({0, 2} & {1, 3}) == set()


# ----------------------------------------------------------------------------
# Prop 6.3 -- Obstruction 3: the Salem slot is empty
# ----------------------------------------------------------------------------
def test_obstruction3_trace_down_regimes_disjoint():
    """Prop 6.3: trace-down t=theta+1/theta is in [-2,2] on |theta|=1 but
    t>2 for real theta>1, with excess (theta-1)^2/theta; regimes disjoint."""
    th = sp.symbols('theta', positive=True)
    excess = (th + 1 / th) - 2
    assert is_zero(sp.simplify(excess - (th - 1)**2 / th))
    # on the unit circle theta=e^{i a}: t=2cos a in [-2,2]
    a = sp.symbols('a', real=True)
    t_circle = sp.simplify(sp.exp(sp.I * a) + sp.exp(-sp.I * a))
    assert is_zero(t_circle - 2 * sp.cos(a))
    # real growth axis value tau0=beta+1/beta > 2 > phi for any beta>1 (sample)
    for beta in [sp.Rational(11, 10), 2, 5]:
        tau0 = beta + 1 / beta
        assert tau0 > 2
        assert tau0 > phi


def test_obstruction3_all_golden_combos_stay_real():
    """Prop 6.3: unions/products/powers of {phi,psi} are real (never a complex
    unit-circle eigenvalue); sample the closure."""
    vals = [phi, psi, phi**2, phi * psi, psi**2, phi**3, phi**2 * psi]
    for v in vals:
        assert is_zero(sp.im(sp.nsimplify(v)))


# ----------------------------------------------------------------------------
# Prop 6.4 / Ex 6.5 -- the flip
# ----------------------------------------------------------------------------
def test_flip_discriminant_and_gram_determinant():
    """Prop 6.4: g_C(x)=x^2+x-C has discriminant D=1+4C; the trace-form Gram
    determinant of Q(theta) equals the same D=4C+1; flip at C=-1/4."""
    C, x = sp.symbols('C x')
    D = sp.discriminant(x**2 + x - C, x)
    assert is_zero(D - (1 + 4 * C))
    # Gram matrix [[2,-1],[-1,2C+1]] from field traces Tr(1)=2,Tr(theta)=-1,
    # Tr(theta^2)=Tr(C-theta)=2C+1
    gram = sp.Matrix([[2, -1], [-1, 2 * C + 1]])
    assert is_zero(gram.det() - (4 * C + 1))
    # flip (double root) at C=-1/4
    assert sp.solve(sp.Eq(1 + 4 * C, 0), C) == [sp.Rational(-1, 4)]


def test_flip_example_both_sides():
    """Ex 6.5: at C=1, D=5>0, g_1 has real roots {phi^{-1},-phi} (golden gate);
    at C=-1, D=-3<0, g_{-1} has complex pair (primitive cube roots)."""
    x = sp.symbols('x')
    r1 = sp.solve(x**2 + x - 1, x)
    assert set(sp.nsimplify(v) for v in r1) == {sp.nsimplify(1 / phi), sp.nsimplify(-phi)}
    assert sp.discriminant(x**2 + x - 1, x) == 5
    r2 = sp.solve(x**2 + x + 1, x)
    assert set(sp.nsimplify(v) for v in r2) == {sp.nsimplify((-1 + sp.I * sp.sqrt(3)) / 2),
                                                sp.nsimplify((-1 - sp.I * sp.sqrt(3)) / 2)}
    assert sp.discriminant(x**2 + x + 1, x) == -3
    # they are primitive cube roots of unity
    for v in r2:
        assert is_zero(sp.expand(v**3) - 1)


# ----------------------------------------------------------------------------
# Thm 6.6 -- the orthogonal partner K
# ----------------------------------------------------------------------------
def test_K_quadratic_in_y_root_split():
    """Thm 6.6: y=x^2 turns K into y^2+5y-5 with roots
    y_+=(-5+3sqrt5)/2>0 and y_-=(-5-3sqrt5)/2<0 (straddling zero)."""
    y = sp.symbols('y')
    yroots = sp.solve(y**2 + 5 * y - 5, y)
    yp = sp.nsimplify((-5 + 3 * sqrt5) / 2)
    ym = sp.nsimplify((-5 - 3 * sqrt5) / 2)
    assert set(sp.nsimplify(v) for v in yroots) == {yp, ym}
    assert float(yp) > 0 and float(ym) < 0
    # decimal illustrations
    assert abs(float(yp) - 0.854) < 0.001
    assert abs(float(ym) - (-5.854)) < 0.001


def test_K_terrain_root_equals_5qrt_over_phi():
    """Thm 6.6: real roots +/-K with K=sqrt(y_+)=5^{1/4}/phi (~0.924)."""
    yp = (-5 + 3 * sqrt5) / 2
    K = sp.sqrt(yp)
    assert is_zero(sp.simplify(K - 5**sp.Rational(1, 4) / phi))
    assert abs(float(K) - 0.924) < 0.001
    assert float(K) < 1              # real roots are inside the unit circle


def test_K_rotation_root_and_magnitude():
    """Thm 6.6: imaginary roots +/- i*beta with beta=sqrt(|y_-|) (~2.420) and
    magnitude M(K)=beta^2=(5+3sqrt5)/2=phi^2 sqrt5 (~5.854)."""
    ym = (-5 - 3 * sqrt5) / 2
    beta2 = -ym
    assert is_zero(sp.simplify(beta2 - phi**2 * sqrt5))
    assert is_zero(sp.simplify(beta2 - (5 + 3 * sqrt5) / 2))
    assert abs(float(sp.sqrt(beta2)) - 2.420) < 0.001
    assert abs(float(beta2) - 5.854) < 0.001


def test_K_roots_solve_the_quartic():
    """Thm 6.6: the four stated roots (+/-K real, +/- i beta imaginary) are
    exactly the roots of K(x)=x^4+5x^2-5."""
    x = sp.symbols('x')
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2
    K = sp.sqrt(yp)
    beta = sp.sqrt(-ym)
    for root in [K, -K, sp.I * beta, -sp.I * beta]:
        assert is_zero(sp.expand(root**4 + 5 * root**2 - 5))
    # and these are all four roots
    quartic_roots = sp.solve(x**4 + 5 * x**2 - 5, x)
    assert len(quartic_roots) == 4


def test_K_magnitude_via_definition():
    """Thm 6.6: M(K)=beta^2 by the Mahler definition (only the two imaginary
    roots +/- i beta lie outside the unit circle)."""
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2
    K = sp.sqrt(yp)
    beta = sp.sqrt(-ym)
    roots = [K, -K, sp.I * beta, -sp.I * beta]
    out = outside_unit(roots)
    assert len(out) == 2             # the two imaginary roots
    exact = sp.nsimplify(sp.prod([sp.Abs(v) for v in out]))
    assert is_zero(sp.simplify(exact - phi**2 * sqrt5))


def test_K_charge_completes_Z4():
    """Thm 6.6: chi(K)={0,1,2,3}; real roots give {0,2}, imaginary +/- i beta
    (angles +/- pi/2) give the missing {1,3}."""
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2
    K = sp.sqrt(yp)
    beta = sp.sqrt(-ym)
    roots = [K, -K, sp.I * beta, -sp.I * beta]
    assert charge_set(roots) == {0, 1, 2, 3}
    assert charge_set([K, -K]) == {0, 2}
    assert charge_set([sp.I * beta, -sp.I * beta]) == {1, 3}


# ----------------------------------------------------------------------------
# Prop 6.7 -- the parity criterion
# ----------------------------------------------------------------------------
def test_parity_criterion_even_quartic_imaginary_roots():
    """Prop 6.7: an even real quartic (polynomial in x^2) with a negative
    x^2-root places its complex roots on the imaginary axis (charges {1,3}).
    Test several even quartics with negative y-root."""
    x = sp.symbols('x')
    for (b, c) in [(5, -5), (1, -3), (2, -7), (4, -1)]:
        # x^4 + b x^2 + c, need a negative y-root: y^2+by+c=0 has y_-<0
        y = sp.symbols('y')
        yroots = sp.solve(y**2 + b * y + c, y)
        negs = [r for r in yroots if sp.im(r) == 0 and r < 0]
        assert negs, (b, c)
        for yv in negs:
            imag_root = sp.I * sp.sqrt(-yv)
            # purely imaginary -> charge 1 (or 3 for the conjugate)
            assert charge_one(imag_root) in (1, 3)
            assert is_zero(sp.re(sp.nsimplify(imag_root)))


def test_parity_criterion_noneven_has_generic_angle():
    """Prop 6.7 (converse illustration): a NON-even quartic over Q(5^{1/4})
    places a complex root off the imaginary axis (charge-lattice miss).
    Example x^4 - 5^{1/4} x^3 + ... use a concrete non-even quartic with a
    genuinely off-axis complex root."""
    x = sp.symbols('x')
    # x^4 + x + 1 (non-even) has complex roots at generic angles
    roots = sp.nroots(x**4 + x + 1)
    complex_roots = [r for r in roots if abs(sp.im(r)) > 1e-9]
    assert complex_roots
    # at least one complex root is NOT on the real or imaginary axis
    offaxis = [r for r in complex_roots
               if abs(sp.re(r)) > 1e-9 and abs(sp.im(r)) > 1e-9]
    assert offaxis, "expected an off-axis complex root for a non-even quartic"


# ----------------------------------------------------------------------------
# Ex 6.8 -- parity floor counterpart
# ----------------------------------------------------------------------------
def test_exercise_x4_minus_1_charge_completion_at_floor():
    """Ex 6.8: x^4-1=(x-1)(x+1)(x^2+1) has roots {1,-1,i,-i}, charges
    {0,2,1,3}=Z/4Z but M=1 (all on the unit circle) -- charge completion at
    the floor, opposite side from K."""
    x = sp.symbols('x')
    roots = sp.solve(x**4 - 1, x)
    assert set(sp.nsimplify(r) for r in roots) == {1, -1, sp.I, -sp.I}
    assert charge_set(roots) == {0, 1, 2, 3}
    assert mahler(roots) == 1
