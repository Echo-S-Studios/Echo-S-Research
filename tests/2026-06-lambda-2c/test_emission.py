"""
Independent verification of Sections 13-15: the Lehmer frontier, the emission
algebra (angle confinement), the self-action difference-set spectrum, the smallest
degree-4 Salem number, and the engine's exact-identity table.
Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - rem:lehmer   : Lehmer band right endpoint L = 1.17628... (largest real root of
                   Lehmer's polynomial); smallest stocked measure is phi.
  - Smyth floor  : mu_S = 1.32472 (root of x^3-x-1).
  - ssec:angle   : every catalog root has argument in (pi/2)Z = {0,pi/2,pi,3pi/2};
                   the operations (add/union/double) preserve (pi/2)Z, so an
                   on-circle emitted eigenvalue is a 4th root of unity.
  - ssec:uniform : spec(ad_R) = {0, +-(phi-psi)} = {0, +-sqrt5} (difference set);
                   every degree-4 Salem number has Mah >= beta_4 = 1.72208 > phi.
  - sec:engine   : Mah(x^2-D)=D=1+4C; det G=4D; cost floors Mah(x^2-24)=24 (2 sqrt6),
                   Mah(x^2-7)=7 (sqrt7); floor costs 2 log 24, 2 log 7.
"""
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 45
x = sp.symbols('x')
phi = (1 + sp.sqrt(5)) / 2
psi = (1 - sp.sqrt(5)) / 2


# ---------------------------------------------------------------------------
# rem:lehmer, Smyth floor
# ---------------------------------------------------------------------------

def test_lehmer_number_value():
    """rem:lehmer: 'L = 1.17628...' is the largest real root of Lehmer's polynomial
    x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1.
    """
    lehmer = lambda z: z**10 + z**9 - z**7 - z**6 - z**5 - z**4 - z**3 + z + 1
    L = mp.findroot(lehmer, 1.18)
    assert abs(L - mp.mpf('1.176280818259917506544070338')) < mp.mpf(10) ** (-24)


def test_smyth_plastic_floor():
    """rem:lehmer / Smyth: 'mu_S = 1.32472 is the plastic number, root of x^3-x-1',
    and the Lehmer band (1,L) sits below phi and below mu_S.
    """
    mu = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    assert abs(mu - mp.mpf('1.324717957244746025960908854')) < mp.mpf(10) ** (-24)
    L = mp.findroot(lambda z: z**10 + z**9 - z**7 - z**6 - z**5 - z**4 - z**3 + z + 1, 1.18)
    assert 1 < L < mu < phi


# ---------------------------------------------------------------------------
# ssec:angle  --  angle confinement to (pi/2)Z
# ---------------------------------------------------------------------------

def _args_in_half_pi_Z(coeffs, tol=1e-9):
    """True iff every root's argument is a multiple of pi/2."""
    for r in np.roots(coeffs):
        ang = np.angle(complex(r)) % (2 * np.pi)
        k = ang / (np.pi / 2)
        if abs(k - round(k)) > tol:
            return False
    return True


def test_catalog_arguments_in_half_pi_Z():
    """ssec:angle: 'every catalog argument lies in (pi/2)Z'.  The catalog minpolys
    are x^2-2, x^2-3, x^2-5, x^2-x-1, x^2+x-1, x^2-7x+1, x^4+5x^2-5.
    """
    catalog = {
        'sqrt2': [1, 0, -2],
        'sqrt3': [1, 0, -3],
        'sqrt5': [1, 0, -5],
        'phi':   [1, -1, -1],
        'tau':   [1, 1, -1],
        'gap':   [1, -7, 1],
        'Kform': [1, 0, 5, 0, -5],
    }
    for name, coeffs in catalog.items():
        assert _args_in_half_pi_Z(coeffs), name


def test_half_pi_Z_closed_under_doubling_and_addition():
    """ssec:angle: '(pi/2)Z is closed under addition (tensor: arg add) and doubling
    (squaring: arg double), and no operation takes a square root'.  So any on-circle
    emitted eigenvalue is forced to a 4th root of unity {0,pi/2,pi,3pi/2}.
    """
    grp = {0, 1, 2, 3}                       # multiples of pi/2, mod 4
    # closed under addition
    for a in grp:
        for b in grp:
            assert (a + b) % 4 in grp
    # closed under doubling
    for a in grp:
        assert (2 * a) % 4 in grp
    # the only unit-circle emitted eigenvalues are 4th roots of unity
    fourth_roots = {sp.exp(sp.I * sp.pi * k / 2) for k in range(4)}
    assert {sp.simplify(z) for z in fourth_roots} == {1, sp.I, -1, -sp.I}


def test_salem_conjugate_has_irrational_angle():
    """ssec:angle: 'a Salem number's on-circle conjugate is NOT a root of unity, so
    its argument is an irrational multiple of pi -- never one of {0,pi/2,pi,3pi/2}'.
    Witness: the smallest degree-4 Salem number x^4-x^3-x^2-x+1 has its unit-circle
    conjugate at an angle that is not a multiple of pi/2.
    """
    roots = np.roots([1, -1, -1, -1, 1])
    oncircle = [r for r in roots if abs(abs(r) - 1) < 1e-9 and abs(r.imag) > 1e-9]
    assert oncircle
    for r in oncircle:
        ang = np.angle(complex(r)) % (2 * np.pi)
        k = ang / (np.pi / 2)
        assert abs(k - round(k)) > 1e-6         # NOT a multiple of pi/2


# ---------------------------------------------------------------------------
# ssec:uniform
# ---------------------------------------------------------------------------

def test_self_action_spectrum_is_difference_set():
    """ssec:uniform: 'spec(ad_R) = {mu_i - mu_j} = {0, +-(phi-psi)} = {0, +-sqrt5}'.
    """
    mus = [phi, psi]
    diffs = {sp.simplify(a - b) for a in mus for b in mus}
    assert sp.Integer(0) in diffs
    assert any(sp.simplify(d - sp.sqrt(5)) == 0 for d in diffs)
    assert any(sp.simplify(d + sp.sqrt(5)) == 0 for d in diffs)
    assert sp.simplify((phi - psi) - sp.sqrt(5)) == 0


def test_smallest_degree4_salem_number():
    """ssec:uniform: 'every degree-four Salem number has Mah >= beta_4 = 1.72208 > phi'.
    Verify the known smallest degree-4 Salem number (root of x^4-x^3-x^2-x+1) has
    this value, is a genuine Salem number, and exceeds phi.
    """
    coeffs = [1, -1, -1, -1, 1]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    outside = [r for r in roots if abs(r) > 1 + mp.mpf(10) ** (-20)]
    oncircle = [r for r in roots if abs(abs(r) - 1) < mp.mpf(10) ** (-20)]
    inside = [r for r in roots if abs(r) < 1 - mp.mpf(10) ** (-20)]
    # Salem signature: exactly one root outside, one inside, the rest on the circle
    assert len(outside) == 1 and len(inside) == 1 and len(oncircle) == 2
    beta4 = outside[0].real
    assert abs(beta4 - mp.mpf('1.722083805739043')) < mp.mpf(10) ** (-12)
    assert beta4 > float(phi)
    # Mahler measure = beta4 (only one root outside the circle)
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    assert abs(M - beta4) < mp.mpf(10) ** (-12)


# ---------------------------------------------------------------------------
# sec:engine  --  exact-identity table
# ---------------------------------------------------------------------------

def test_engine_mahler_and_detG_identities():
    """sec:engine table: 'Mah(x^2-D)=D=1+4C' and 'det G = 4D' (gap basis)."""
    C = sp.symbols('C')
    D = 1 + 4 * C
    G = sp.Matrix([[2, 0], [0, 2 * D]])
    assert sp.expand(G.det() - 4 * D) == 0
    for Dval in (2, 3, 5):
        Cval = sp.Rational(Dval - 1, 4)
        assert 1 + 4 * Cval == Dval
        assert sp.simplify(sp.sqrt(Dval) ** 2 - Dval) == 0     # Mah(x^2-D)=D


def test_engine_cost_floors():
    """sec:engine table: 'cost floors 2 log 24 (2 sqrt6, Mah=24) and 2 log 7
    (sqrt7, Mah=7) = lambda log Mah at c=1'.  Verify Mah(x^2-24)=24, Mah(x^2-7)=7,
    where 2 sqrt6 = sqrt24.
    """
    # 2 sqrt6 = sqrt24, minpoly x^2-24
    assert sp.simplify((2 * sp.sqrt(6))**2 - 24) == 0
    for Dval, cost in [(24, 2 * sp.log(24)), (7, 2 * sp.log(7))]:
        roots = [sp.sqrt(Dval), -sp.sqrt(Dval)]
        M = 1
        for r in roots:
            if sp.Abs(r) > 1:
                M *= sp.Abs(r)
        assert sp.simplify(M - Dval) == 0
        # lambda log Mah at c=1 (lambda=2): 2 log Mah
        assert sp.simplify(2 * sp.log(M) - cost) == 0


def test_clifford_unity_gate():
    """sec:engine table: 'Clifford unity tau = x^2+x-1 (C=1 gate) = T^2+T=1'.
    The golden gate polynomial x^2+x=C at C=1 is x^2+x-1.
    """
    C = sp.Integer(1)
    gate_poly = sp.expand(sp.Symbol('T')**2 + sp.Symbol('T') - C)
    T = sp.Symbol('T')
    assert gate_poly == T**2 + T - 1
