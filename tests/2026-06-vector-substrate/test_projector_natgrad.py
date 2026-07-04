"""
Independent verification of the projector/residual calculus and the
natural-gradient reading (Sec. 3.1: Def. 3.1, Prop. 3.2, Ex. 3.4, Ex. 3.5;
Sec. 4.1 Pythagoras; Sec. 7.6: Prop. 7.13, Ex. 7.15).

The projector P=B(B^T G B)^{-1} B^T G is rebuilt and its idempotence,
G-orthogonality, capture criterion, and one-step Newton (natural-gradient)
identity are checked from scratch.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import vsub_nf as nf

x = sp.symbols('x')
G5 = sp.Matrix([[2, 1], [1, 3]])                 # golden trace form


def test_golden_G_from_scratch():
    """The trace form used throughout equals the rebuilt Gram."""
    C = nf.companion_from_poly(x**2 - x - 1)
    assert nf.gram([(1, 0), (0, 1)], C) == G5


def test_projector_phi_line():
    """Ex. 3.4: B=(0,1)^T (line Q.phi), B^T G B=3, P=(1/3)[[0,0],[1,3]]."""
    B = sp.Matrix([0, 1])
    assert (B.T * G5 * B)[0] == 3
    P = nf.projector(B, G5)
    assert P == sp.Rational(1, 3) * sp.Matrix([[0, 0], [1, 3]])


def test_projector_idempotent_and_orthogonal():
    """Prop. 3.2 (i),(ii): P^2=P and B^T G r = 0."""
    B = sp.Matrix([0, 1])
    P = nf.projector(B, G5)
    assert sp.simplify(P * P - P) == sp.zeros(2, 2)
    r = nf.residual((-1, 2), B, G5)
    assert (B.T * G5 * r) == sp.Matrix([0])


def test_capture_2phi_and_novelty_sqrt5():
    """Ex. 3.4: 2phi=(0,2) is captured (r=0); sqrt5=2phi-1=(-1,2) projects to
    (0,5/3) with r=(-1,1/3) and ||r||_G^2 = 5/3."""
    B = sp.Matrix([0, 1])
    P = nf.projector(B, G5)
    assert P * sp.Matrix([0, 2]) == sp.Matrix([0, 2])          # captured
    assert nf.residual((0, 2), B, G5) == sp.Matrix([0, 0])
    Px = P * sp.Matrix([-1, 2])
    assert Px == sp.Matrix([0, sp.Rational(5, 3)])
    r = nf.residual((-1, 2), B, G5)
    assert r == sp.Matrix([-1, sp.Rational(1, 3)])
    assert nf.gnorm2(r, G5) == sp.Rational(5, 3)


def test_capture_criterion_equivalences():
    """Prop. 3.2 (iii),(iv): P x = x  <=>  x in col(B)  <=>  ||r||_G^2 = 0."""
    B = sp.Matrix([0, 1])
    P = nf.projector(B, G5)
    for x_in in [sp.Matrix([0, 5]), sp.Matrix([0, -2])]:      # multiples of phi
        assert P * x_in == x_in
        assert nf.gnorm2(x_in - P * x_in, G5) == 0
    x_out = sp.Matrix([1, 0])                                 # the constant 1: off-line
    assert P * x_out != x_out
    assert nf.gnorm2(x_out - P * x_out, G5) > 0


def test_complex_projector_Qi():
    """Ex. 3.5: Q(i), G_2=2I, B=<1>=(1,0), x=i=(0,1): B^* G_2 B=2, B^* G_2 x=0,
    P x = 0, r=i, ||r||_{G2}^2 = 2 (i is pure off-axis novelty vs <1>)."""
    G2 = sp.Matrix([[2, 0], [0, 2]])
    B = sp.Matrix([1, 0])
    xi = sp.Matrix([0, 1])
    assert (B.T * G2 * B)[0] == 2
    assert (B.T * G2 * xi)[0] == 0
    P = nf.projector(B, G2)
    assert P * xi == sp.Matrix([0, 0])
    r = xi - P * xi
    assert r == sp.Matrix([0, 1])
    assert nf.gnorm2(r, G2) == 2


def test_pythagoras_split():
    """Eq. (17): ||x||_G^2 = ||P x||_G^2 + ||r||_G^2. Check on x=sqrt5 vs line Q.phi."""
    B = sp.Matrix([0, 1])
    P = nf.projector(B, G5)
    xv = sp.Matrix([-1, 2])
    Px = P * xv
    r = xv - Px
    assert nf.gnorm2(xv, G5) == nf.gnorm2(Px, G5) + nf.gnorm2(r, G5)


def test_natural_gradient_one_step():
    """Prop. 7.13 / Ex. 7.15: single Newton step a*=(B^T G B)^{-1} B^T G x lands
    at B a*=P x with residual r; in Q(sqrt5), B=<phi>, x=sqrt5 gives a*=5/3."""
    B = sp.Matrix([0, 1])
    xv = sp.Matrix([-1, 2])
    astar = (B.T * G5 * B).inv() * (B.T * G5 * xv)
    assert astar == sp.Matrix([sp.Rational(5, 3)])
    P = nf.projector(B, G5)
    assert B * astar == P * xv
    assert xv - B * astar == nf.residual((-1, 2), B, G5)
    # gradient/Hessian identities: grad L(0) = -B^T G x, Hess = B^T G B
    grad0 = B.T * G5 * (B * sp.Matrix([0]) - xv)
    hess = B.T * G5 * B
    assert (-hess.inv() * grad0) == astar
    # a Euclidean step would use B^T B = 1 not B^T G B = 3, landing off-projection
    assert (B.T * B)[0] == 1 and (B.T * G5 * B)[0] == 3
