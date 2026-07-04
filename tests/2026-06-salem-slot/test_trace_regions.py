"""
Independent verification of Section 2 (downstairs / trace-down) and the
lift/interval structure of "The Occupant of the Salem Slot".

Each test re-derives the claimed object from the paper's stated premises
(the trace substitution t = theta + 1/theta) rather than restating the
paper's number.
"""
import sympy as sp
import mpmath as mp

x, t = sp.symbols('x t')
phi = (1 + sp.sqrt(5)) / 2


def test_discriminant_is_flip_t2_minus_4():
    """Def 2.1: the flip D(t)=t^2-4 is the discriminant of x^2 - t x + 1."""
    D = sp.discriminant(x**2 - t*x + 1, x)
    assert sp.simplify(D - (t**2 - 4)) == 0


def test_tracedown_definition_lehmer():
    """Def 2.1 / benchmark: for a reciprocal monic P of degree 2m, the
    trace-down T (degree m) satisfies P(x) = x^m * T(x + 1/x).
    Independently substitute the paper's stated Lehmer trace-down
    T = t^5+t^4-5t^3-5t^2+4t+3 and confirm it reconstructs Lehmer's
    polynomial x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1 exactly."""
    T = t**5 + t**4 - 5*t**3 - 5*t**2 + 4*t + 3
    reconstructed = sp.expand(x**5 * T.subs(t, x + 1/x))
    reconstructed = sp.expand(reconstructed)  # clears x^-k terms since P is reciprocal
    lehmer = (x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1)
    assert sp.simplify(reconstructed - lehmer) == 0


def test_region_on_circle_maps_to_interval():
    """Lem 2.2: theta = e^{i psi} on the unit circle -> t = 2 cos psi in [-2,2]."""
    mp.mp.dps = 40
    for k in range(1, 20):
        psi = mp.pi * k / 21          # psi in (0, pi)
        theta = mp.e**(1j*psi)
        tval = theta + 1/theta
        assert abs(tval.imag) < mp.mpf(10)**(-30)
        assert abs(tval.real - 2*mp.cos(psi)) < mp.mpf(10)**(-30)
        assert -2 <= tval.real <= 2


def test_region_grow_theta_gt_1():
    """Lem 2.2: theta > 1 -> t > 2 (grow), with t -> 2+ as theta -> 1+."""
    mp.mp.dps = 40
    for theta in [mp.mpf('1.0001'), mp.mpf('1.1'), mp.mpf('1.5'),
                  mp.mpf('2'), mp.mpf('10')]:
        assert theta + 1/theta > 2
    # limit at the flip
    assert abs((mp.mpf('1') + mp.mpf('1e-8')) + 1/(mp.mpf('1') + mp.mpf('1e-8')) - 2) < 1e-15


def test_region_decay_theta_lt_minus1():
    """Lem 2.2: theta < -1 -> t < -2 (decay)."""
    mp.mp.dps = 40
    for theta in [mp.mpf('-1.0001'), mp.mpf('-1.5'), mp.mpf('-3'), mp.mpf('-10')]:
        assert theta + 1/theta < -2


def test_no_fourth_channel_disjoint():
    """Cor 2.3: 'on-circle and expanding' = |t|<2 AND t>2 is impossible;
    the intervals (-2,2) and (2, oo) are disjoint."""
    onc = sp.Interval.open(-2, 2)
    grow = sp.Interval.open(2, sp.oo)
    assert onc.intersect(grow) == sp.EmptySet


def test_lehmer_tracedown_totally_real_split():
    """Cor 2.3: a Salem's trace-down T is totally real with one root >2
    and m-1 roots in (-2,2). For Lehmer (m=5): 5 real roots, 1 past 2,
    4 inside (-2,2)."""
    T = t**5 + t**4 - 5*t**3 - 5*t**2 + 4*t + 3
    roots = sp.nroots(T, n=40)
    reals = [complex(r) for r in roots]
    assert all(abs(r.imag) < 1e-30 for r in reals)          # totally real
    reals = sorted(r.real for r in reals)
    past2 = [r for r in reals if r > 2]
    inside = [r for r in reals if -2 < r < 2]
    assert len(past2) == 1
    assert len(inside) == 4


def test_lift_on_lattice_roots_of_unity():
    """Sec 8 / entry: L(T)(x) = x^{deg T} T(x+1/x). On the captured lattice
    {-2,0,2} the lift is a root of unity square:
    L(t-2)=(x-1)^2, L(t)=x^2+1, L(t+2)=(x+1)^2."""
    def L(Tt):
        d = sp.degree(Tt, t)
        return sp.expand(x**d * Tt.subs(t, x + 1/x))
    assert sp.simplify(L(t - 2) - (x - 1)**2) == 0
    assert sp.simplify(L(t) - (x**2 + 1)) == 0
    assert sp.simplify(L(t + 2) - (x + 1)**2) == 0


def test_interval_endpoints_lift():
    """Lem 4.4: t=2 lifts to x^2-2x+1 (root x=1); t=sqrt5 lifts to
    x^2 - sqrt5 x + 1 with roots {phi, 1/phi}."""
    # t = 2 endpoint
    p2 = sp.expand((x + 1/x - 2) * x)
    assert sp.simplify(p2 - (x - 1)**2) == 0
    assert sp.solve(sp.Eq(x**2 - 2*x + 1, 0), x) == [1]
    # t = sqrt5 endpoint -> golden pair
    P = x**2 - sp.sqrt(5)*x + 1
    roots = sp.solve(P, x)
    roots = [sp.nsimplify(sp.radsimp(r)) for r in roots]
    prod = sp.simplify(roots[0] * roots[1])
    ssum = sp.simplify(roots[0] + roots[1])
    assert prod == 1                       # golden pair {phi, 1/phi}
    assert sp.simplify(ssum - sp.sqrt(5)) == 0
    # one of the roots is exactly phi
    assert any(sp.simplify(r - phi) == 0 for r in roots)


def test_trace_strictly_increasing_injective():
    """Redirection is injective because trace is strictly increasing on
    beta>1: d/dbeta (beta+1/beta) = 1 - 1/beta^2 > 0."""
    b = sp.symbols('b', positive=True)
    dtr = sp.diff(b + 1/b, b)
    assert sp.simplify(dtr - (1 - 1/b**2)) == 0
    # strictly positive for b>1
    assert sp.simplify((1 - 1/b**2).subs(b, sp.Rational(3, 2))) > 0
