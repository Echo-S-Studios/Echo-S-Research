"""
Independent verification of Section 2 ("The keystone: the golden generator")
of "The Dissolved Helix and Its Orthogonal Partner".

Covers: Definition 2.1 (Golden data), Proposition 2.2 (The keystone) incl.
eq. (1) [eq:root5], Proposition 2.3 (The self-action) incl. eq. (2) [eq:adspec],
and Definition 2.4 (The two characters: Mahler measure and the Z/4Z phase charge).

Each test rebuilds the quantity from the stated premises (R, the Fibonacci
companion; phi, psi the roots of x^2-x-1) and then compares to the paper's value.
"""
import sympy as sp

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
psi = (1 - sqrt5) / 2
R = sp.Matrix([[0, 1], [1, 1]])          # Def. 2.1: Fibonacci companion
I2 = sp.eye(2)
lam = sp.symbols('lambda')


def chi_of(z):
    """Character II (phase): floor(2*arg(z)/pi mod 4), round-half-up = nearest int.
    Def. 2.4: chi(lambda) = round(2 arg lambda / pi) mod 4."""
    a = sp.arg(z)
    return int(sp.floor(2 * a / sp.pi + sp.Rational(1, 2)) % 4)


# ---------------------------------------------------------------- Prop. 2.2
def test_charpoly_is_fibonacci():
    """Prop. 2.2: charpoly(R) = lambda^2 - lambda - 1."""
    cp = sp.expand(R.charpoly(lam).as_expr())
    assert sp.simplify(cp - (lam**2 - lam - 1)) == 0


def test_trace_and_determinant():
    """Prop. 2.2: Tr R = 1, det R = -1."""
    assert sp.trace(R) == 1
    assert R.det() == -1


def test_R_squared_equals_R_plus_I():
    """Prop. 2.2 (Cayley-Hamilton): R^2 = R + I."""
    assert sp.simplify(R**2 - (R + I2)) == sp.zeros(2, 2)


def test_spectrum_is_phi_psi():
    """Prop. 2.2: spec(R) = {phi, psi}."""
    eig = R.eigenvals()
    got = {sp.nsimplify(sp.radsimp(k)) for k in eig}
    assert got == {sp.nsimplify(phi), sp.nsimplify(psi)}
    assert all(m == 1 for m in eig.values())


def test_vieta_sum_and_product():
    """Prop. 2.2: phi + psi = 1 and phi*psi = -1."""
    assert sp.simplify(phi + psi - 1) == 0
    assert sp.simplify(phi * psi + 1) == 0


def test_root5_identity_eq1():
    """Eq. (1): phi - psi = sqrt5 = phi + phi^{-1}."""
    assert sp.simplify((phi - psi) - sqrt5) == 0
    assert sp.simplify((phi + 1 / phi) - sqrt5) == 0


def test_psi_alternate_forms():
    """Def. 2.1: psi = 1 - phi = -phi^{-1}."""
    assert sp.simplify(psi - (1 - phi)) == 0
    assert sp.simplify(psi - (-1 / phi)) == 0


# ---------------------------------------------------------------- Prop. 2.3
def _ad_R_matrix():
    """Build ad_R = [R, .] as a 4x4 matrix on M_2 ~ R^4 (row-major vec)."""
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]),
             sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    vec = lambda M: [M[0, 0], M[0, 1], M[1, 0], M[1, 1]]
    cols = [sp.Matrix(4, 1, vec(R * E - E * R)) for E in basis]
    return sp.Matrix.hstack(*cols)


def test_adR_spectrum_eq2():
    """Eq. (2): spec(ad_R) = {0, 0, +sqrt5, -sqrt5}."""
    adR = _ad_R_matrix()
    eig = adR.eigenvals()
    eig = {sp.simplify(sp.radsimp(k)): v for k, v in eig.items()}
    assert eig.get(sp.Integer(0)) == 2
    assert eig.get(sqrt5) == 1
    assert eig.get(-sqrt5) == 1
    # characteristic polynomial is lambda^2 (lambda^2 - 5)
    assert sp.simplify(adR.charpoly(lam).as_expr() - lam**2 * (lam**2 - 5)) == 0


def test_coupling_is_sqrt5():
    """Prop. 2.3: the coupling sqrt5 = phi + phi^{-1} = phi - psi is the sole
    irrational in spec(ad_R)."""
    assert sp.simplify((phi + 1 / phi) - (phi - psi)) == 0
    assert sp.nsimplify(sqrt5).is_irrational


# ---------------------------------------------------------------- Def. 2.4
def test_golden_object_mahler():
    """Def. 2.4: Mah(A_phi) = phi, since |psi| = phi^{-1} < 1."""
    assert sp.simplify(sp.Abs(psi) - 1 / phi) == 0
    assert (sp.Abs(psi) < 1) == True         # noqa: E712
    assert (sp.Abs(phi) > 1) == True         # noqa: E712
    mah = sp.Max(1, sp.Abs(phi)) * sp.Max(1, sp.Abs(psi))
    assert sp.simplify(mah - phi) == 0


def test_golden_object_charge():
    """Def. 2.4: chi(A_phi) = {0, 2}  (phi>0 -> 0, psi<0 -> 2)."""
    assert {chi_of(phi), chi_of(psi)} == {0, 2}


def test_operator_laws_on_golden_object():
    """Def. 2.4: Adams power psi^2 squares Mahler and doubles chi mod 4.
    On A_phi: Mah phi -> phi^2 ; chi {0,2} -> {2*0,2*2} mod4 = {0}."""
    # squaring the eigenvalues: {phi^2, psi^2}, both derived from A_phi
    mah_sq = sp.Max(1, sp.Abs(phi**2)) * sp.Max(1, sp.Abs(psi**2))
    assert sp.simplify(mah_sq - phi**2) == 0
    doubled = {(2 * c) % 4 for c in (0, 2)}
    assert doubled == {0}
