"""
Independent verification of the regular-representation / spectral claims of
"The Vector Substrate" (Sec. 2: Prop. 2.2/regrep, Ex. 2.3, Prop. 2.6/spectral,
Ex. 2.7, Ex. 2.8, Rem. 2.5/mahlerspectral).

We build companion matrices from minimal polynomials in the paper's convention
and check trace/norm/charpoly/minpoly and the Minkowski diagonalisation from
scratch (sympy exact + mpmath high precision).
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import mpmath as mp
import vsub_nf as nf

x = sp.symbols('x')


def test_rho_phi_companion_trace_norm():
    """Ex. 2.3: rho(phi)=C(x^2-x-1)=[[0,1],[1,1]] with Tr rho(phi)=1=Tr(phi)
    and det rho(phi)=-1=N(phi)."""
    C = nf.companion_from_poly(x**2 - x - 1)
    assert C == sp.Matrix([[0, 1], [1, 1]])
    assert C.trace() == 1                     # field trace of phi
    assert C.det() == -1                      # field norm of phi
    # trace = -(coeff of x^{n-1}) ; norm = (-1)^n c_0 : cross-check identity
    assert C.trace() == -(sp.Poly(x**2 - x - 1, x).all_coeffs()[1])
    assert C.det() == (-1)**2 * sp.Poly(x**2 - x - 1, x).all_coeffs()[2]


def test_rho_sqrt5_from_phi():
    """Ex. 2.3: for alpha=sqrt5=2*phi-1, rho(alpha)=2*rho(phi)-I=[[-1,2],[2,1]],
    charpoly=minpoly=x^2-5."""
    C = nf.companion_from_poly(x**2 - x - 1)
    R = 2 * C - sp.eye(2)
    assert R == sp.Matrix([[-1, 2], [2, 1]])
    assert sp.expand(R.charpoly(x).as_expr()) == x**2 - 5
    assert sp.expand(nf.matrix_minpoly(R)) == x**2 - 5
    # sqrt5 = 2 phi - 1 has coordinate vector (-1, 2); trace 0, norm -5
    assert nf.field_trace((-1, 2), C) == 0
    assert nf.field_norm((-1, 2), C) == -5


def test_rho_rational_is_derogatory():
    """Ex. 2.3: rho(3)=3I has charpoly (x-3)^2 but minpoly (x-3):
    the n/d multiplicity of Eq. (2) made visible."""
    C = nf.companion_from_poly(x**2 - x - 1)
    R = nf.rho((3, 0), C)                       # the rational 3 = 3*1 + 0*phi
    assert R == 3 * sp.eye(2)
    assert sp.expand(R.charpoly(x).as_expr()) == sp.expand((x - 3)**2)
    assert sp.expand(nf.matrix_minpoly(R)) == x - 3


def test_charmin_multiplicity_eq2():
    """Prop. 2.2 Eq.(2): charpoly(rho(x))=m_x^{n/d}, minpoly(rho(x))=m_x.
    Take x=sqrt2 inside K=Q(sqrt2+sqrt3) (n=4, d=2): charpoly=(x^2-2)^2, minpoly=x^2-2."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)   # theta=sqrt2+sqrt3
    # sqrt2 = (theta^3 - 9 theta)/2  -> coords (0,-9/2,0,1/2)
    s2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
    R = nf.rho(s2, C)
    assert sp.expand(R.charpoly(x).as_expr()) == sp.expand((x**2 - 2)**2)
    assert sp.expand(nf.matrix_minpoly(R)) == x**2 - 2
    # sanity: this element really squares to 2 (rho(sqrt2)^2 = 2 I)
    assert R * R == 2 * sp.eye(4)


def test_spectral_diagonalisation_golden():
    """Ex. 2.7 / Prop. 2.6: M rho(phi) M^{-1}=diag(phi,phi') with
    M=[[1,phi],[1,phi']], and charpoly(rho(phi))=x^2-x-1."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    C = nf.companion_from_poly(x**2 - x - 1)
    M = sp.Matrix([[1, phi], [1, phic]])
    D = sp.simplify(M * C * M.inv())
    assert sp.simplify(D - sp.diag(phi, phic)) == sp.zeros(2, 2)
    assert sp.expand(C.charpoly(x).as_expr()) == x**2 - x - 1


def test_spectral_diagonalises_every_element():
    """Prop. 2.6: the SAME M diagonalises rho(x) for every x. Check on sqrt5:
    M rho(sqrt5) M^{-1}=diag(sqrt5,-sqrt5)=diag(sigma_1,sigma_2)."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    C = nf.companion_from_poly(x**2 - x - 1)
    M = sp.Matrix([[1, phi], [1, phic]])
    R = nf.rho((-1, 2), C)                      # sqrt5 = 2 phi - 1
    D = sp.simplify(M * R * M.inv())
    assert sp.simplify(D - sp.diag(sp.sqrt(5), -sp.sqrt(5))) == sp.zeros(2, 2)


def test_spectrum_quartic_are_conjugates():
    """Ex. 2.8: rho(sqrt2+sqrt3)=C(x^4-10x^2+1), charpoly=m_theta, eigenvalues
    the four conjugates {+-sqrt2+-sqrt3} ~ {+-3.146,+-0.318}."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    assert sp.expand(C.charpoly(x).as_expr()) == x**4 - 10*x**2 + 1
    mp.mp.dps = 40
    roots = sorted(mp.polyroots([1, 0, -10, 0, 1], extraprec=200), key=lambda r: r.real)
    conj = sorted([mp.sqrt(2) + mp.sqrt(3), mp.sqrt(2) - mp.sqrt(3),
                   -mp.sqrt(2) + mp.sqrt(3), -mp.sqrt(2) - mp.sqrt(3)],
                  key=lambda r: r.real)
    for a, b in zip(roots, conj):
        assert abs(a - b) < mp.mpf(10) ** (-30)
    # magnitudes ~ 3.146 and 0.318, exactly two outside the unit circle
    mags = sorted(abs(r) for r in roots)
    assert abs(mags[0] - mp.mpf("0.3178372451")) < mp.mpf("1e-8")
    assert abs(mags[-1] - mp.mpf("3.1462643699")) < mp.mpf("1e-8")
    assert sum(1 for r in roots if abs(r) > 1) == 2


def test_mahler_spectral_sqrt2_sqrt3():
    """Rem. 2.5 / Ex. 2.8: Mahler(sqrt2+sqrt3)=product of eigenvalues outside
    the unit circle = (sqrt2+sqrt3)^2 = 5+2sqrt6 ~ 9.899."""
    # closed form identity
    assert sp.simplify((sp.sqrt(2) + sp.sqrt(3))**2 - (5 + 2*sp.sqrt(6))) == 0
    # rebuilt from the companion spectrum
    m = nf.mahler_measure_mp(x**4 - 10*x**2 + 1, dps=45)
    assert abs(m - (5 + 2*mp.sqrt(6))) < mp.mpf(10) ** (-30)
    assert abs(m - mp.mpf("9.89897948556635")) < mp.mpf("1e-12")


def test_mahler_plastic_number_is_muS():
    """Rem. 2.5 / Fig. 2: for the plastic number (root of x^3-x-1) only the real
    eigenvalue exceeds 1, so Mahler = mu_S = 1.3247179..."""
    m = nf.mahler_measure_mp(x**3 - x - 1, dps=45)
    mp.mp.dps = 45
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    assert abs(m - muS) < mp.mpf(10) ** (-30)
    assert abs(m - mp.mpf("1.32471795724475")) < mp.mpf("1e-12")
