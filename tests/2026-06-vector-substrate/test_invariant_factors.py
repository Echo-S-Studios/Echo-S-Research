"""
Independent verification of the invariant-factor / Smith-normal-form claims
(Sec. 2.3: Cor. 2.9, Thm. 2.10, Ex. 2.11, Rem. 2.12, Ex. 2.13) and the
coordinates-to-minimal-polynomial bridge (Thm. 4.2, Ex. 4.4).

Invariant factors are computed independently via determinantal divisors
(monic gcds of i x i minors of xI-A); non-similarity is cross-checked by
similarity invariants (annihilating polynomials, ranks) that do not use the
invariant-factor machinery.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import vsub_nf as nf

x = sp.symbols('x')


def test_bridge_rho_2sqrt6_matrix():
    """Ex. 4.4: rho(2sqrt6)=rho(theta)^2-5I equals the stated 4x4 matrix,
    where 2sqrt6=theta^2-5 has coordinate vector (-5,0,1,0)."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    R = C**2 - 5*sp.eye(4)
    assert R == sp.Matrix([[-5, 0, -1, 0],
                           [0, -5, 0, -1],
                           [1, 0, 5, 0],
                           [0, 1, 0, 5]])
    # rho of the coordinate vector (-5,0,1,0) gives the same matrix
    assert nf.rho((-5, 0, 1, 0), C) == R


def test_bridge_invariant_factors_and_minpoly():
    """Thm. 4.2 / Ex. 4.4: SNF(xI-rho(2sqrt6))=diag(1,1,x^2-24,x^2-24), so the
    invariant factors are (x^2-24,x^2-24) and the largest (=m_{2sqrt6}) is
    x^2-24; matches (2sqrt6)^2=24 and the doubled factor reflects n/d=2."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    R = C**2 - 5*sp.eye(4)
    ifs = nf.invariant_factors(R)
    assert ifs == [x**2 - 24, x**2 - 24]
    assert sp.expand(nf.matrix_minpoly(R)) == x**2 - 24
    # largest invariant factor == minimal polynomial of the algebraic number
    assert sp.minimal_polynomial(2*sp.sqrt(6), x) == x**2 - 24
    # elementary check (2 sqrt6)^2 = 24
    assert sp.simplify((2*sp.sqrt(6))**2 - 24) == 0
    # product of invariant factors = characteristic polynomial
    assert sp.expand(ifs[0] * ifs[1]) == sp.expand(R.charpoly(x).as_expr())


def test_bridge_monic_integer_admits_seed():
    """Thm. 4.2: alpha in O_K => m_alpha monic over Z (2sqrt6 admitted)."""
    p = sp.Poly(sp.minimal_polynomial(2*sp.sqrt(6), x), x)
    assert p.LC() == 1
    assert all(c.is_integer for c in p.all_coeffs())


def test_charpoly_insufficient_phiphi():
    """Ex. 2.11: rho(phi)+rho(phi) and C((x^2-x-1)^2) share charpoly (x^2-x-1)^2
    and trace 2, but are NOT similar: invariant factors (x^2-x-1,x^2-x-1) vs the
    single (x^2-x-1)^2. Largest invariant factor separates them."""
    Cphi = nf.companion_from_poly(x**2 - x - 1)
    blk = sp.diag(Cphi, Cphi)
    comp = nf.companion_from_poly((x**2 - x - 1)**2)
    # (x^2-x-1)^2 expands as stated
    assert sp.expand((x**2 - x - 1)**2) == x**4 - 2*x**3 - x**2 + 2*x + 1
    # same charpoly and trace
    assert sp.factor(blk.charpoly(x).as_expr()) == (x**2 - x - 1)**2
    assert sp.factor(comp.charpoly(x).as_expr()) == (x**2 - x - 1)**2
    assert blk.trace() == 2 and comp.trace() == 2
    # different invariant factors
    assert nf.invariant_factors(blk) == [x**2 - x - 1, x**2 - x - 1]
    assert nf.invariant_factors(comp) == [x**4 - 2*x**3 - x**2 + 2*x + 1]
    # independent non-similarity witness: block satisfies x^2-x-1, companion does not
    assert blk*blk - blk - sp.eye(4) == sp.zeros(4, 4)
    assert comp*comp - comp - sp.eye(4) != sp.zeros(4, 4)


def test_jordan_full_list_needed():
    """Ex. 2.13: J2(+)J2 and J2(+)J1(+)J1 share charpoly x^4 and minpoly x^2 but
    have invariant factors (x^2,x^2) vs (x,x,x^2): not similar."""
    J2 = sp.Matrix([[0, 1], [0, 0]])
    J1 = sp.Matrix([[0]])
    A1 = sp.diag(J2, J2)
    A2 = sp.diag(J2, J1, J1)
    for A in (A1, A2):
        assert sp.expand(A.charpoly(x).as_expr()) == x**4
        assert sp.expand(nf.matrix_minpoly(A)) == x**2
    assert nf.invariant_factors(A1) == [x**2, x**2]
    assert nf.invariant_factors(A2) == [x, x, x**2]
    # independent non-similarity witness: rank is a similarity invariant
    assert A1.rank() == 2 and A2.rank() == 1


def test_invariant_factor_completeness_property():
    """Thm. 2.10: product of invariant factors = charpoly; largest = minpoly.
    Verified on the derogatory rho(3) inside Q(sqrt5)."""
    C = nf.companion_from_poly(x**2 - x - 1)
    R = nf.rho((3, 0), C)                       # 3 I : charpoly (x-3)^2, minpoly x-3
    ifs = nf.invariant_factors(R)
    assert ifs == [x - 3, x - 3]
    prod = sp.expand(sp.prod(ifs))
    assert prod == sp.expand(R.charpoly(x).as_expr())
    assert sp.expand(ifs[-1]) == sp.expand(nf.matrix_minpoly(R))
