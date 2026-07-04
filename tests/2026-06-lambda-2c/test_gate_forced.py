"""
Independent verification of Section 7: "The gate is forced".
Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - prop:ternary : distinct signed channels of ad_R for a degree-d seed = d^2-d+1,
                   which equals 3 iff d=2.
  - lem:mincost  : min Mahler measure exceeding 1 over irreducible integer quadratics
                   is phi = 1.61803..., attained ONLY by x^2-x-1 and x^2+x-1 (disc 5).
  - rem:cubic    : plastic mu_S = 1.32471... (root x^3-x-1) < phi; a cubic seed gives
                   3^2-3+1 = 7 channels (excluded by the ternary lock).
  - prop:firewallimage : x^4+5x^2-5 -> y^2+5y-5 under y=x^2, discriminant 45 = 9*5.
  - lem:keystone : only golden companion has R^2=R+I; tau gives R^2=I-R; radicand
                   seed gives R^2=D*I.
  - thm:gateforced : forced value c=sqrt5/2, lambda=2c=sqrt5=phi-psi.
"""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
x, y = sp.symbols('x y')


# ---------------------------------------------------------------------------
# prop:ternary  --  channel count d^2 - d + 1 = 3 iff d = 2
# ---------------------------------------------------------------------------

def test_ternary_iff_degree_two_algebra():
    """prop:ternary: 'signed channels = d^2 - d + 1 = 3 iff d=2'.
    Solve d^2-d+1 = 3 over positive integers.
    """
    d = sp.symbols('d', integer=True, positive=True)
    sols = sp.solve(sp.Eq(d**2 - d + 1, 3), d)
    assert sols == [2]
    # d=2 -> 3 channels, d=3 -> 7 channels (rem:cubic)
    f = lambda dd: dd**2 - dd + 1
    assert f(2) == 3
    assert f(3) == 7


def test_channel_count_equals_distinct_eigenvalue_differences():
    """prop:ternary: the channels of ad_R are the eigenvalue differences
    {mu_i - mu_j}.  For a degree-d matrix with GENERIC eigenvalues (all pairwise
    differences distinct) there are exactly d^2-d+1 distinct differences
    (d^2-d nonzero + one zero).  We realise the generic case with Golomb rulers
    (perfect difference sets), whose pairwise differences are all distinct.
    """
    golomb = {2: [0, 1], 3: [0, 1, 3], 4: [0, 1, 3, 7]}
    for d in (2, 3, 4):
        mus = [sp.Integer(p) for p in golomb[d]]
        diffs = {sp.simplify(a - b) for a in mus for b in mus}
        assert len(diffs) == d**2 - d + 1


# ---------------------------------------------------------------------------
# lem:mincost / lem:tie  --  min Mahler over integer quadratics
# ---------------------------------------------------------------------------

def _mahler_from_roots(coeffs):
    """Mahler measure of a monic integer polynomial given numpy coeffs."""
    roots = np.roots(coeffs)
    M = 1.0
    for r in roots:
        if abs(r) > 1 + 1e-12:
            M *= abs(r)
    return M


def test_min_mahler_over_integer_quadratics_is_phi():
    """lem:mincost/lem:tie: 'the minimum Mahler measure exceeding 1 over irreducible
    integer quadratics is phi, attained only by x^2-x-1 and x^2+x-1 (disc 5)'.
    Enumerate monic irreducible integer quadratics x^2+bx+c and find the min>1.
    """
    phi = float((1 + mp.sqrt(5)) / 2)
    best = None
    attained = []
    for b in range(-8, 9):
        for c in range(-8, 9):
            poly = sp.Poly(x**2 + b * x + c, x)
            if not poly.is_irreducible:
                continue
            M = _mahler_from_roots([1, b, c])
            if M > 1 + 1e-9:
                if best is None or M < best - 1e-9:
                    best = M
                    attained = [(b, c)]
                elif abs(M - best) < 1e-9:
                    attained.append((b, c))
    assert abs(best - phi) < 1e-9
    assert set(attained) == {(-1, -1), (1, -1)}          # x^2-x-1 and x^2+x-1
    # both have discriminant 5
    for (b, c) in attained:
        assert sp.discriminant(sp.Poly(x**2 + b * x + c, x)) == 5


def test_min_mahler_is_exactly_phi_high_precision():
    """lem:mincost: the attained minimum is exactly phi = (1+sqrt5)/2.
    Mahler(x^2-x-1) = |phi| (only phi lies outside the unit circle).
    """
    phi = (1 + mp.sqrt(5)) / 2
    roots = mp.polyroots([1, -1, -1])
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    assert abs(M - phi) < mp.mpf(10) ** (-30)


# ---------------------------------------------------------------------------
# rem:cubic  --  plastic number floor
# ---------------------------------------------------------------------------

def test_plastic_number_is_root_of_cubic_and_below_phi():
    """rem:cubic: 'the plastic number mu_S = 1.32471... (root of x^3-x-1) < phi'.
    """
    mu = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    assert abs(mu - mp.mpf('1.324717957244746025960908854')) < mp.mpf(10) ** (-25)
    phi = (1 + mp.sqrt(5)) / 2
    assert mu < phi


def test_min_mahler_over_cubics_is_plastic():
    """rem:cubic / appendix drop-degree witness: dropping degree=2 (allow degree 3)
    lowers the minimum Mahler measure to the plastic number mu_S ~ 1.32472.
    Enumerate irreducible integer cubics and confirm the min>1 is the plastic number.
    """
    best = None
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                poly = sp.Poly(x**3 + a * x**2 + b * x + c, x)
                if not poly.is_irreducible:
                    continue
                M = _mahler_from_roots([1, a, b, c])
                if M > 1 + 1e-9 and (best is None or M < best - 1e-12):
                    best = M
    mu = float(mp.findroot(lambda z: z**3 - z - 1, 1.3))
    assert abs(best - mu) < 1e-8


# ---------------------------------------------------------------------------
# prop:firewallimage  --  squaring firewall x^4+5x^2-5 -> y^2+5y-5
# ---------------------------------------------------------------------------

def test_squaring_firewall_kformation():
    """prop:firewallimage: 'x^4+5x^2-5 -> y^2+5y-5 under y=x^2, of discriminant
    45 = 9*5, landing in Q(sqrt5)'.
    """
    f = x**4 + 5 * x**2 - 5
    g_expected = y**2 + 5 * y - 5
    # verify the substitution y=x^2 recovers f exactly
    assert sp.expand(g_expected.subs(y, x**2) - f) == 0
    disc = sp.discriminant(sp.Poly(g_expected, y))
    assert disc == 45
    assert sp.factorint(45) == {3: 2, 5: 1}       # 45 = 9 * 5
    # its roots lie in Q(sqrt5): (-5 +- 3 sqrt5)/2
    roots = sp.solve(g_expected, y)
    for r in roots:
        assert sp.simplify(r).has(sp.sqrt(5))


def _companion(coeffs):
    """Companion matrix of monic x^n + c_{n-1}x^{n-1}+...+c_0 given [c_0,...,c_{n-1}].
    With this convention companion(x^2-x-1) = [[0,1],[1,1]] = the golden keystone.
    """
    n = len(coeffs)
    Cm = sp.zeros(n)
    for i in range(1, n):
        Cm[i, i - 1] = 1
    for i in range(n):
        Cm[i, n - 1] = -coeffs[i]
    return Cm


def test_firewall_fixed_point_and_derogatory_witness():
    """def:firewall / sec:engine table: 'the spectral-lift firewall Phi = companion o
    charpoly is spectrum-preserving and collapses every matrix of a given charpoly
    onto one companion; R(R)=R on companions; phi (+) phi is the derogatory witness'.

      (a) Phi fixes companion matrices: Phi(companion(x^2-x-1)) = companion(x^2-x-1).
      (b) phi (+) phi (block-diag of two golden companions) is derogatory: charpoly
          (x^2-x-1)^2 but minpoly x^2-x-1, so Phi is many-to-one (it and the 4x4
          companion share a charpoly yet differ).
    """
    Rphi = _companion([-1, -1])                 # companion of x^2 - x - 1
    assert Rphi == sp.Matrix([[0, 1], [1, 1]])
    # (a) fixed point of Phi = companion o charpoly
    cp = Rphi.charpoly(x).as_expr()
    Phi_Rphi = _companion([sp.Poly(cp, x).all_coeffs()[-1],
                           sp.Poly(cp, x).all_coeffs()[-2]])
    assert Phi_Rphi == Rphi
    assert sp.expand(cp - (x**2 - x - 1)) == 0
    # (b) phi (+) phi derogatory witness
    S = sp.Matrix(sp.BlockDiagMatrix(Rphi, Rphi))
    I4 = sp.eye(4)
    # S satisfies x^2-x-1 (degree 2), so its minimal polynomial has degree <= 2,
    assert S**2 - S - I4 == sp.zeros(4)
    # while its characteristic polynomial has degree 4 -> minpoly != charpoly:
    assert sp.factor(S.charpoly(x).as_expr()) == (x**2 - x - 1)**2
    assert not S.is_diagonal()                    # non-scalar -> minpoly degree exactly 2
    # hence derogatory (many-to-one collapse under Phi): 2 < 4.


# ---------------------------------------------------------------------------
# lem:keystone  --  R^2 = R + I unique to golden; tau: I-R; radicand: D*I
# ---------------------------------------------------------------------------

def test_keystone_relations_by_seed():
    """lem:keystone: 'only the golden companion satisfies R^2=R+I; the conjugate
    gate tau gives R^2=I-R; each radicand seed gives R^2=D*I'.
    """
    I2 = sp.eye(2)
    # golden keystone companion of x^2-x-1
    Rphi = sp.Matrix([[0, 1], [1, 1]])
    assert Rphi**2 == Rphi + I2
    # conjugate gate tau, companion of x^2+x-1
    Rtau = sp.Matrix([[0, 1], [1, -1]])
    assert Rtau**2 == I2 - Rtau
    # radicand seeds companion of x^2 - D  ->  R^2 = D I
    for D in (2, 3, 5):
        Rrad = sp.Matrix([[0, D], [1, 0]])
        assert Rrad**2 == D * I2


def test_keystone_forces_golden_eigenvalues():
    """lem:keystone: 'R^2=R+I forces eigenvalues {phi, psi}'.
    A matrix satisfying R^2=R+I has minimal/char roots solving lam^2=lam+1.
    """
    lam = sp.symbols('lam')
    roots = sp.solve(sp.Eq(lam**2, lam + 1), lam)
    phi = (1 + sp.sqrt(5)) / 2
    psi = (1 - sp.sqrt(5)) / 2
    assert {sp.simplify(r) for r in roots} == {sp.simplify(phi), sp.simplify(psi)}


def test_gate_forced_value():
    """thm:gateforced / eq:goldenvalue: 'C=1, c=sqrt5/2, lambda=2c=sqrt5=phi-psi'."""
    C = sp.Integer(1)
    c = sp.sqrt(5) / 2
    lam = 2 * c
    phi = (1 + sp.sqrt(5)) / 2
    psi = (1 - sp.sqrt(5)) / 2
    assert sp.simplify(lam - sp.sqrt(5)) == 0
    assert sp.simplify(lam - (phi - psi)) == 0
    assert sp.simplify(2 * c * C - sp.sqrt(1 + 4 * C)) == 0     # gate balance holds
