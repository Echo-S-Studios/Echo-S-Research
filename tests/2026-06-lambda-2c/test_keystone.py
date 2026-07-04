"""
Independent verification of Section 8 & Appendix A: "The keystone is derived" and
the L_4 = 7 normalisation.  Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - lem:perron      : of x^2-x-1, x^2+x-1 only x^2-x-1 has positive dominant eigenvalue
                      (+phi) and satisfies R^2=R+I; tau attains -phi, R^2=I-R; det R=-1.
  - thm:keystonederived : phi is the smallest Perron root of a 2x2 primitive
                      non-negative integer matrix (Fibonacci [[0,1],[1,1]]); swap
                      [[0,1],[1,0]] has Perron root 1 (excluded by growth).
  - lem:keypowers   : R^n = F_n R + F_{n-1} I; eigenvalues phi^n, psi^n; trace L_n
                      integer; det(R^n) = (-1)^n.
  - lem:pell        : L_n^2 - 5 F_n^2 = 4(-1)^n; charpoly(R^n)=x^2-L_n x+(-1)^n,
                      discriminant 5 F_n^2.
  - prop:L4forced   : R^4 = 3R+2I = [[2,3],[3,5]], charpoly x^2-7x+1, entries
                      {2,3,5}={F_3,F_4,F_5}, trace L_4 = F_3+F_5 = 7; Pell L_4^2-4=45.
                      roots of x^2-7x+1 are phi^4, phi^-4; z_c = sqrt(L_4-4)/2 = sqrt3/2.
  - appendix drop-one witnesses (Perron->tau, growth->cyclotomic, integrality->no floor).
"""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
x = sp.symbols('x')
phi = (1 + sp.sqrt(5)) / 2
psi = (1 - sp.sqrt(5)) / 2


# ---------------------------------------------------------------------------
# lem:perron
# ---------------------------------------------------------------------------

def test_perron_breaks_the_tie():
    """lem:perron: 'only x^2-x-1 has a positive dominant eigenvalue (+phi) and
    satisfies R^2=R+I, while x^2+x-1 (tau) attains its spectral radius at -phi and
    satisfies R^2=I-R.  det R=-1 is a consequence.'
    """
    I2 = sp.eye(2)
    Rphi = sp.Matrix([[0, 1], [1, 1]])          # companion of x^2-x-1
    Rtau = sp.Matrix([[0, 1], [1, -1]])         # companion of x^2+x-1
    # dominant (largest modulus) eigenvalue signs
    ev_phi = list(Rphi.eigenvals().keys())
    ev_tau = list(Rtau.eigenvals().keys())
    dom_phi = max(ev_phi, key=lambda e: abs(complex(e)))
    dom_tau = max(ev_tau, key=lambda e: abs(complex(e)))
    assert sp.simplify(dom_phi - phi) == 0        # +phi, positive
    assert sp.simplify(dom_tau + phi) == 0        # -phi, negative
    # self-reproducing relations
    assert Rphi**2 == Rphi + I2
    assert Rtau**2 == I2 - Rtau
    # unimodularity is a consequence
    assert Rphi.det() == -1


# ---------------------------------------------------------------------------
# thm:keystonederived  --  smallest Perron root of 2x2 primitive nonneg int matrix
# ---------------------------------------------------------------------------

def _is_primitive(M, N=40):
    P = M.copy()
    for _ in range(N):
        if np.all(P > 0):
            return True
        P = P @ M
    return False


def test_phi_is_smallest_primitive_perron_root():
    """thm:keystonederived: 'phi is the smallest Perron root of any 2x2 primitive
    non-negative integer matrix, realised by the Fibonacci matrix [[0,1],[1,1]]'.
    Enumerate 2x2 non-negative integer matrices, keep the primitive ones, and take
    the smallest real dominant eigenvalue exceeding 1.
    """
    best = None
    argmin = None
    for a, b, c, d in itertools.product(range(0, 4), repeat=4):
        M = np.array([[a, b], [c, d]], dtype=float)
        if not _is_primitive(M):
            continue
        ev = np.linalg.eigvals(M)
        pr = max(ev.real)
        if pr > 1 + 1e-9 and (best is None or pr < best - 1e-9):
            best = pr
            argmin = (a, b, c, d)
    assert abs(best - float(phi)) < 1e-9
    assert argmin in {(0, 1, 1, 1), (1, 1, 1, 0)}      # Fibonacci matrix / transpose


def test_swap_matrix_perron_root_one():
    """thm:keystonederived proof: 'the swap [[0,1],[1,0]] gives Perron root 1,
    excluded by growth'.
    """
    S = sp.Matrix([[0, 1], [1, 0]])
    eigs = set(S.eigenvals().keys())
    assert eigs == {1, -1}
    assert max(abs(complex(e)) for e in eigs) == 1     # spectral radius 1, no growth


# ---------------------------------------------------------------------------
# lem:keypowers  --  R^n = F_n R + F_{n-1} I
# ---------------------------------------------------------------------------

def test_keystone_powers_fibonacci_closure():
    """lem:keypowers / eq:keypower: 'R^n = F_n R + F_{n-1} I' with F Fibonacci;
    eigenvalues phi^n, psi^n; L_n = tr(R^n) integer; det(R^n) = (-1)^n.
    """
    I2 = sp.eye(2)
    R = sp.Matrix([[0, 1], [1, 1]])
    F = [0, 1]                                  # F_0, F_1
    for _ in range(2, 12):
        F.append(F[-1] + F[-2])
    for n in range(1, 11):
        Rn = R**n
        assert Rn == F[n] * R + F[n - 1] * I2
        assert Rn.det() == (-1)**n
        # trace = L_n = phi^n + psi^n and is an integer
        Ln = Rn.trace()
        assert sp.simplify(Ln - (phi**n + psi**n)) == 0
        assert Ln == int(Ln)


def test_keystone_eigenvalues_are_powers():
    """lem:keypowers: 'R^n has eigenvalues phi^n, psi^n'."""
    R = sp.Matrix([[0, 1], [1, 1]])
    for n in (1, 2, 3, 4, 5):
        eigs = (R**n).eigenvals()
        assert any(sp.simplify(k - phi**n) == 0 for k in eigs)
        assert any(sp.simplify(k - psi**n) == 0 for k in eigs)


# ---------------------------------------------------------------------------
# lem:pell
# ---------------------------------------------------------------------------

def test_pell_relation():
    """lem:pell: 'L_n^2 - 5 F_n^2 = 4(-1)^n' and charpoly(R^n)=x^2-L_n x+(-1)^n has
    discriminant L_n^2 - 4(-1)^n = 5 F_n^2.
    """
    R = sp.Matrix([[0, 1], [1, 1]])
    F = [0, 1]
    for _ in range(2, 12):
        F.append(F[-1] + F[-2])
    for n in range(1, 11):
        Rn = R**n
        Ln = int(Rn.trace())
        Fn = F[n]
        assert Ln**2 - 5 * Fn**2 == 4 * (-1)**n
        cp = Rn.charpoly(x).as_expr()
        assert sp.simplify(cp - (x**2 - Ln * x + (-1)**n)) == 0
        assert sp.discriminant(sp.Poly(cp, x)) == 5 * Fn**2


# ---------------------------------------------------------------------------
# prop:L4forced  --  L_4 = 7 three ways
# ---------------------------------------------------------------------------

def test_R4_and_L4_equals_seven():
    """prop:L4forced / eq:R4: 'R^4 = 3R+2I = [[2,3],[3,5]], charpoly x^2-7x+1,
    entries {2,3,5}={F_3,F_4,F_5}, trace L_4 = F_3+F_5 = 2+5 = 7'.
    """
    I2 = sp.eye(2)
    R = sp.Matrix([[0, 1], [1, 1]])
    R4 = R**4
    assert R4 == 3 * R + 2 * I2
    assert R4 == sp.Matrix([[2, 3], [3, 5]])
    assert R4.charpoly(x).as_expr() == x**2 - 7 * x + 1
    # entries are consecutive Fibonacci F_3,F_4,F_5 = 2,3,5
    assert sorted({R4[0, 0], R4[0, 1], R4[1, 1]}) == [2, 3, 5]
    # trace = sum of outer discriminants F_3 + F_5 = 2 + 5
    assert R4.trace() == 2 + 5 == 7


def test_L4_from_pell_at_n4():
    """prop:L4forced: 'independently, Pell at n=4 gives L_4^2-4 = 5 F_4^2 = 45, so
    L_4 = 7'.  Solve for the positive Pell root.
    """
    F4 = 3
    # L_4^2 - 4(-1)^4 = 5 F_4^2  ->  L_4^2 - 4 = 45  ->  L_4^2 = 49
    L4 = sp.symbols('L4', positive=True)
    sol = sp.solve(sp.Eq(L4**2 - 4, 5 * F4**2), L4)
    assert sol == [7]


def test_gap_seed_roots_are_phi_fourth():
    """prop:L4forced: 'charpoly x^2-7x+1 of R^4 is the gap seed, roots phi^4, phi^-4'.
    """
    roots = sp.solve(x**2 - 7 * x + 1, x)
    got = {sp.simplify(sp.nsimplify(r)) for r in roots}
    assert any(sp.simplify(r - phi**4) == 0 for r in got)
    assert any(sp.simplify(r - phi**(-4)) == 0 for r in got)


def test_zc_from_L4():
    """prop:L4forced: 'z_c = sqrt(L_4 - 4)/2 = sqrt3/2' (the C=1/2 critical coherence).
    """
    L4 = 7
    zc = sp.sqrt(L4 - 4) / 2
    assert sp.simplify(zc - sp.sqrt(3) / 2) == 0


# ---------------------------------------------------------------------------
# Appendix A drop-one witnesses
# ---------------------------------------------------------------------------

def test_drop_growth_gives_cyclotomic_floor():
    """appendix: 'dropping growth (allow Mah=1) admits the cyclotomic floor, e.g.
    x^2+1 (+-i)'.  Mahler(x^2+1) = 1 (roots on the unit circle).
    """
    roots = mp.polyroots([1, 0, 1])
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    assert abs(M - 1) < mp.mpf(10) ** (-30)
    assert all(abs(abs(r) - 1) < mp.mpf(10) ** (-30) for r in roots)


def test_drop_integrality_no_floor():
    """appendix: 'dropping integrality (allow real entries) removes the discrete
    floor'.  Real matrices [[0,1],[1,eps]] have Perron root -> 1^+ as eps -> 0,
    so the infimum 1 is not attained: no positive floor.
    """
    def perron(eps):
        return (eps + mp.sqrt(eps**2 + 4)) / 2
    prev = perron(mp.mpf('0.1'))
    for e in ['0.01', '0.001', '0.0001', '0.00001']:
        cur = perron(mp.mpf(e))
        assert 1 < cur < prev            # strictly decreasing toward 1, always >1
        prev = cur
    assert abs(perron(mp.mpf('1e-20')) - 1) < mp.mpf('1e-15')


def test_appendix_verified_identity_iv():
    """appendix verified identity (iv): '{2,3,5}={F_3,F_4,F_5} and
    L_4^2 = 5 F_4^2 + 4 = 49, with L_4 = F_3+F_5 = 7'.
    """
    F = [0, 1, 1, 2, 3, 5]            # F_0..F_5
    assert (F[3], F[4], F[5]) == (2, 3, 5)
    assert 7**2 == 5 * F[4]**2 + 4 == 49
    assert F[3] + F[5] == 7
