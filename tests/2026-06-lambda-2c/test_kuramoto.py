"""
Independent verification of Section 11: "Kuramoto division of labour".
Paper: "The Exchange Rate lambda = 2c ...".

Claims (sec:kuramoto, rem:division):
  - critical coupling K_c = 2/(pi g(0)); for a Lorentzian frequency density this
    equals 2 gamma.
  - order parameter r ~ sqrt(K - K_c): the mean-field 1/2 exponent.
  - the critical coherence pinned is z_c = sqrt3/2 (the C=1/2 gate).
"""
import sympy as sp


def test_lorentzian_threshold_is_two_gamma():
    """sec:kuramoto: 'K_c = 2/(pi g(0)) = 2 gamma for a Lorentzian frequency density'.
    Lorentzian g(omega) = gamma / (pi (gamma^2 + omega^2)); g(0) = 1/(pi gamma), so
    K_c = 2/(pi g(0)) = 2 gamma.
    """
    gamma, omega = sp.symbols('gamma omega', positive=True)
    g = gamma / (sp.pi * (gamma**2 + omega**2))
    g0 = g.subs(omega, 0)
    assert sp.simplify(g0 - 1 / (sp.pi * gamma)) == 0
    Kc = 2 / (sp.pi * g0)
    assert sp.simplify(Kc - 2 * gamma) == 0
    # sanity: the density integrates to 1 over the real line
    assert sp.integrate(g, (omega, -sp.oo, sp.oo)) == 1


def test_mean_field_half_exponent():
    """sec:kuramoto: 'order parameter r ~ sqrt(K - K_c)', the mean-field 1/2 exponent.
    Near threshold the self-consistency r = f(r) with f(r) = a r - b r^3 (Landau form,
    a = (K-K_c)*const>0) has the nontrivial branch r ~ sqrt(a/b) ~ sqrt(K-K_c).
    """
    r, a, b = sp.symbols('r a b', positive=True)
    # steady states of dr/dt = a r - b r^3
    sols = sp.solve(sp.Eq(a * r - b * r**3, 0), r)
    nontrivial = [s for s in sols if s != 0]
    # the positive nontrivial branch is sqrt(a/b)
    assert any(sp.simplify(s - sp.sqrt(a / b)) == 0 for s in nontrivial)
    # with a proportional to (K - K_c): r ~ sqrt(K-K_c) -> exponent 1/2
    K, Kc, k0 = sp.symbols('K K_c k0', positive=True)
    r_branch = sp.sqrt((k0 * (K - Kc)) / b)
    # the exponent of (K-K_c) in r_branch is exactly 1/2
    assert sp.simplify(r_branch - sp.sqrt(k0 / b) * (K - Kc)**sp.Rational(1, 2)) == 0


def test_critical_coherence_zc():
    """sec:kuramoto / prop:square: 'critical coherence z_c = sqrt3/2 is the C=1/2
    gate' (z = sqrt D / 2, D = 1+4C; at C=1/2, D=3, z = sqrt3/2).
    """
    Cval = sp.Rational(1, 2)
    D = 1 + 4 * Cval
    z = sp.sqrt(D) / 2
    assert sp.simplify(z - sp.sqrt(3) / 2) == 0
    # perfect-square pinning (eq:square): D = 4 z^2, C = z^2 - 1/4
    zsym = sp.symbols('z', real=True)
    assert sp.simplify((4 * zsym**2) - (2 * zsym)**2) == 0
    assert sp.simplify((zsym**2 - sp.Rational(1, 4)) - ((4 * zsym**2 - 1) / 4)) == 0
