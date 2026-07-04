"""
Independent verification of Section 3 (redirection map / occupant), the
entropy trade (Sec 5), the golden trace identity and the 'forced' identity
(Sec 6-7) of "The Occupant of the Salem Slot".
"""
import sympy as sp
import mpmath as mp

x, t = sp.symbols('x t')
phi = (1 + sp.sqrt(5)) / 2


def test_amgm_square_identity():
    """Lem 3.1 proof: beta + 1/beta - 2 = (sqrt(beta) - 1/sqrt(beta))^2 >= 0."""
    b = sp.symbols('b', positive=True)
    lhs = b + 1/b - 2
    rhs = (sp.sqrt(b) - 1/sp.sqrt(b))**2
    assert sp.simplify(lhs - rhs) == 0


def test_amgm_bound_and_limit():
    """Lem 3.1: beta + 1/beta > 2 for beta>1, -> 2+ as beta -> 1+."""
    b = sp.symbols('b', positive=True)
    assert sp.limit(b + 1/b, b, 1) == 2
    mp.mp.dps = 40
    for beta in [mp.mpf('1.0001'), mp.mpf('1.1762808'), mp.mpf('1.7220838'), mp.mpf('5')]:
        assert beta + 1/beta > 2


def test_quadratic_redirection_identity():
    """Prop 4.2: tau0 - 2 = (beta-1)^2 / beta  exactly."""
    b = sp.symbols('b', positive=True)
    lhs = (b + 1/b) - 2
    rhs = (b - 1)**2 / b
    assert sp.simplify(lhs - rhs) == 0


def test_quadratic_redirection_near_floor_expansion():
    """Prop 4.2: with beta = 1 + delta, tau0 = 2 + delta^2/(1+delta)
    = 2 + delta^2 + O(delta^3)."""
    d = sp.symbols('d', positive=True)
    tau0 = (1 + d) + 1/(1 + d)
    assert sp.simplify(tau0 - (2 + d**2/(1 + d))) == 0
    ser = sp.series(tau0, d, 0, 3).removeO()
    # 2 + delta^2 + O(delta^3): the delta^0 term is 2, delta^1 term is 0, delta^2 term is 1
    assert ser.coeff(d, 0) == 2
    assert ser.coeff(d, 1) == 0
    assert ser.coeff(d, 2) == 1


def test_occupant_tau0_above_two_above_phi():
    """Thm 3.3: tau0 > 2 > phi for every Salem; chain 2 > phi."""
    assert sp.simplify(2 - phi) > 0                      # 2 > phi
    b = sp.symbols('b', positive=True)
    # tau0 - 2 = (b-1)^2/b > 0 for b>1  =>  tau0 > 2
    assert sp.simplify(((b + 1/b) - 2) - (b - 1)**2/b) == 0


def test_occupant_mahler_ge_tau0():
    """Thm 3.3: Mah(T) = prod max(1,|tau_j|) >= tau0 > 2, illustrated on
    Lehmer's trace-down. Mah is built independently from the roots."""
    mp.mp.dps = 40
    T = t**5 + t**4 - 5*t**3 - 5*t**2 + 4*t + 3
    roots = [complex(r) for r in sp.nroots(T, n=40)]
    mah = 1.0
    for r in roots:
        mah *= max(1.0, abs(r))
    tau0 = max(r.real for r in roots if abs(r.imag) < 1e-20)
    assert mah >= tau0 > 2


def test_entropy_trade():
    """Cor 5.1: log(beta + 1/beta) > log(beta) for beta>1, because
    (beta+1/beta) - beta = 1/beta > 0."""
    b = sp.symbols('b', positive=True)
    assert sp.simplify((b + 1/b) - b - 1/b) == 0     # difference is exactly 1/b>0
    mp.mp.dps = 40
    for beta in [mp.mpf('1.1762808'), mp.mpf('1.5061357'), mp.mpf('1.7220838')]:
        assert mp.log(beta + 1/beta) > mp.log(beta)


def test_golden_trace_identity():
    """Lem 6.1: phi + 1/phi = sqrt5."""
    assert sp.simplify((phi + 1/phi) - sp.sqrt(5)) == 0


def test_forced_selfaction_equals_tracedown():
    """Prop 7.1: seed x^2 - a x - 1 has roots phi_a>1>0>psi_a with
    product -1, so psi_a = -1/phi_a. Then the ad_R grow eigenvalue
    phi_a - psi_a equals the trace-down phi_a + 1/phi_a = sqrt(a^2+4)."""
    a = sp.symbols('a', positive=True)
    roots = sp.solve(x**2 - a*x - 1, x)
    phi_a = [r for r in roots if sp.simplify(sp.limit(r, a, 1) - phi) == 0]
    phi_a = phi_a[0] if phi_a else max(roots, key=lambda r: float(r.subs(a, 1)))
    psi_a = min(roots, key=lambda r: float(r.subs(a, 1)))
    # product of roots is -1  => psi_a = -1/phi_a
    assert sp.simplify(phi_a * psi_a + 1) == 0
    assert sp.simplify(psi_a + 1/phi_a) == 0
    # grow eigenvalue = trace-down = sqrt(a^2+4)
    assert sp.simplify((phi_a - psi_a) - (phi_a + 1/phi_a)) == 0
    assert sp.simplify((phi_a - psi_a) - sp.sqrt(a**2 + 4)) == 0
    # golden seed a=1 -> sqrt5
    assert sp.simplify((phi_a - psi_a).subs(a, 1) - sp.sqrt(5)) == 0
