"""
Independent verification of Section 6 (the golden limit and its rate):
the golden limit theorem, the proof-sketch identities, the linear rate
(slope phi^-1, curvature magnitude sqrt5-2), and the geometric-rate table.

NOTE: the *sign* of the quadratic term in Prop 6.4's displayed expansion
of (sqrt5 - tau0) does not reproduce from the paper's own derivatives; it
is captured by an xfail test below and recorded in NOTES.md.
"""
import sympy as sp
import mpmath as mp
import pytest

x = sp.symbols('x')
mp.mp.dps = 50
PHI = (1 + mp.sqrt(5)) / 2
SQRT5 = mp.sqrt(5)
phi = (1 + sp.sqrt(5)) / 2


def _beta_n(n):
    """Salem factor beta_n = largest real root of S_n = x^n P - P*,
    P = x^2-x-1, i.e. S_n(x) = x^{n+2} - x^{n+1} - x^n + x^2 + x - 1.
    beta_n -> phi from below."""
    f = lambda z: z**(n + 2) - z**(n + 1) - z**n + z**2 + z - 1
    return mp.findroot(f, PHI - mp.mpf('1e-3'))


def _gap(n):
    b = _beta_n(n)
    return SQRT5 - (b + 1/b)


# ---- proof-sketch identities (Thm 6.5) --------------------------------
def test_proofsketch_reciprocal_polynomial():
    """Thm 6.5 proof: P*(x) = x^2 P(1/x) = 1 - x - x^2 for P = x^2-x-1."""
    P = x**2 - x - 1
    Pstar = sp.expand(x**2 * P.subs(x, 1/x))
    assert sp.simplify(Pstar - (1 - x - x**2)) == 0


def test_proofsketch_Pstar_at_phi():
    """Thm 6.5 proof: P*(phi) = -2 phi."""
    Pstar = 1 - x - x**2
    assert sp.simplify(Pstar.subs(x, phi) - (-2*phi)) == 0


def test_proofsketch_Pprime_at_phi():
    """Thm 6.5 proof: P'(phi) = 2 phi - 1 = sqrt5."""
    P = x**2 - x - 1
    dP = sp.diff(P, x)
    assert sp.simplify(dP.subs(x, phi) - (2*phi - 1)) == 0
    assert sp.simplify(dP.subs(x, phi) - sp.sqrt(5)) == 0


# ---- the golden limit (Thm 6.2) ---------------------------------------
def test_golden_limit_value():
    """Thm 6.2: beta_n -> phi and tau0(beta_n) -> phi + 1/phi = sqrt5."""
    b = _beta_n(30)
    assert abs(b - PHI) < 1e-5
    assert b < PHI                       # approaches from below
    tau0 = b + 1/b
    assert abs(tau0 - SQRT5) < 1e-5


# ---- the linear rate (Prop 6.4) ---------------------------------------
def test_linear_rate_slope():
    """Prop 6.4: trace'(beta)=1-1/beta^2, and trace'(phi)=1-1/phi^2
    = (sqrt5-1)/2 = 1/phi."""
    b = sp.symbols('b', positive=True)
    dtr = sp.diff(b + 1/b, b)
    assert sp.simplify(dtr - (1 - 1/b**2)) == 0
    val = sp.simplify(dtr.subs(b, phi))
    assert sp.simplify(val - (sp.sqrt(5) - 1)/2) == 0
    assert sp.simplify(val - 1/phi) == 0


def test_linear_rate_curvature_magnitude():
    """Prop 6.4: the curvature coefficient (1/2) trace''(phi) = sqrt5-2
    = 1/phi^3 (magnitude of the quadratic coefficient)."""
    b = sp.symbols('b', positive=True)
    d2tr = sp.diff(b + 1/b, b, 2)
    curv = sp.simplify(d2tr.subs(b, phi) / 2)
    assert sp.simplify(curv - (sp.sqrt(5) - 2)) == 0
    assert sp.simplify(curv - 1/phi**3) == 0


def test_linear_rate_expansion_leading_term():
    """Prop 6.4: the leading term of (sqrt5 - tau0) in (phi-beta) is
    +phi^-1 (phi-beta); verified independently by series."""
    b, u = sp.symbols('b u')
    expr = sp.sqrt(5) - (b + 1/b)
    ser = sp.series(expr.subs(b, phi - u), u, 0, 3).removeO()
    assert sp.simplify(ser.coeff(u, 0)) == 0
    assert sp.simplify(ser.coeff(u, 1) - 1/phi) == 0


def test_linear_rate_expansion_quadratic_sign():
    """Prop 6.4 (corrected 2026-07-04): sqrt5 - tau0 = phi^-1 (phi-beta)
    - (sqrt5-2)(phi-beta)^2 + O((phi-beta)^3). The paper originally displayed a
    + on the quadratic term; the independent Taylor coefficient is -(sqrt5-2)
    (= -phi^-3). The linear slope phi^-1 and curvature magnitude sqrt5-2 are
    unchanged."""
    b, u = sp.symbols('b u')
    expr = sp.sqrt(5) - (b + 1/b)
    ser = sp.series(expr.subs(b, phi - u), u, 0, 3).removeO()
    c2 = sp.simplify(ser.coeff(u, 2))
    # corrected coefficient is -(sqrt5 - 2)
    assert sp.simplify(c2 - (-(sp.sqrt(5) - 2))) == 0


# ---- the geometric rate table (Thm 6.6) -------------------------------
def test_geometric_rate_table():
    """Thm 6.6 table: for n in {9,...,27}, gap = sqrt5 - tau0(beta_n),
    gap*phi^n, and the consecutive ratio gap(n-1)/gap(n)."""
    expected = {
        # n : (gap, gap*phi^n, consecutive_ratio = gap(n-1)/gap(n))
        9:  (1.289085e-2, 0.979875, 1.68711),
        12: (2.858434e-3, 0.920407, 1.64037),
        15: (6.613990e-4, 0.902149, 1.62501),
        18: (1.551834e-4, 0.896649, 1.62012),
        21: (3.656841e-5, 0.895048, 1.61863),
        24: (8.628276e-6, 0.894597, 1.61820),
        27: (2.036577e-6, 0.894473, 1.61808),
    }
    for n, (gexp, gphi_exp, ratio_exp) in expected.items():
        g = _gap(n)
        assert abs(g - mp.mpf(gexp)) < abs(mp.mpf(gexp)) * 1e-5
        assert abs(g * PHI**n - mp.mpf(gphi_exp)) < 5e-6
        ratio = _gap(n - 1) / g
        assert abs(ratio - mp.mpf(ratio_exp)) < 5e-5


def test_geometric_rate_asymptotics():
    """Thm 6.6: gap*phi^n -> 2/sqrt5 = 0.894427 and the consecutive ratio
    -> phi. Also (phi - beta_n)*phi^{n-1} -> 2/sqrt5."""
    limit = 2 / SQRT5
    assert abs(limit - mp.mpf('0.894427')) < 5e-7
    g27, g28 = _gap(27), _gap(28)
    assert abs(g27 * PHI**27 - limit) < 1e-4
    assert abs(g27 / g28 - PHI) < 1e-3
    # height-gap asymptotic (phi - beta_n) ~ (2/sqrt5) phi^{1-n}
    b27 = _beta_n(27)
    assert abs((PHI - b27) * PHI**(27 - 1) - limit) < 1e-4
