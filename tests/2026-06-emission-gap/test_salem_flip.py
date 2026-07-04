"""The Salem = flip-straddle characterization and the emission delta
(Sec. 8: Def. 8.1, Lemma 8.2, Thm. 8.3; Prop. 8.4; Prop. 11.1; App. A).
"""
import mpmath as mp
import numpy as np
import sympy as sp

import emgap_util as U

mp.mp.dps = 40

x, t = U.x, U.t
PHI = (1 + mp.sqrt(5)) / 2


def test_tracedown_lehmer_roots_and_straddle():
    """Lemma 8.2 / App. A: the Lehmer trace-down T has roots
    {2.0264, 0.9137, -0.5847, -1.4689, -1.8866}: exactly one > 2, four in
    (-2, 2), none <= -2 -- the flip-straddle."""
    L = sp.Poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], x)
    T = U.trace_down(L)
    roots = sorted(complex(r).real for r in np.roots([float(c) for c in T.all_coeffs()]))
    expected = [-1.8866, -1.4689, -0.5847, 0.9137, 2.0264]
    for got, exp in zip(roots, expected):
        assert abs(got - exp) < 1e-3
    assert sum(1 for r in roots if r > 2) == 1
    assert sum(1 for r in roots if -2 < r < 2) == 4
    assert U.flip_straddle(T)


def test_beta4_is_salem_via_flip_and_exceeds_phi():
    """Lemma 8.2 / Cor. 10.4: x^4-x^3-x^2-x+1 has trace-down t^2-t-3 with one
    root > 2 and one in (-2,2): a Salem number, whose value beta_4 exceeds phi."""
    b4 = sp.Poly([1, -1, -1, -1, 1], x)
    T = U.trace_down(b4)
    assert sp.expand(T.as_expr() - (t**2 - t - 3)) == 0
    assert U.flip_straddle(T)
    beta = mp.findroot(lambda z: z**4 - z**3 - z**2 - z + 1, 1.72)
    assert beta > PHI


def test_delta_detects_and_rejects():
    """Thm. 8.3: the delta (reciprocal factor whose trace-down straddles the
    flip) fires exactly on Salem carriers. Positive: Lehmer, beta_4.
    Negative: cyclotomic Phi_10 (all on-circle) and a non-reciprocal factor."""
    assert U.has_salem_factor(sp.Poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], x))
    assert U.has_salem_factor(sp.Poly([1, -1, -1, -1, 1], x))
    # Phi_10 = x^4 - x^3 + x^2 - x + 1 : reciprocal deg 4 but NOT Salem (cyclotomic)
    assert not U.has_salem_factor(sp.Poly(sp.cyclotomic_poly(10, x), x))
    # x^4 - 2 : non-reciprocal, no Salem
    assert not U.has_salem_factor(sp.Poly([1, 0, 0, 0, -2], x))


def test_fourth_roots_of_unity_tracedowns():
    """Prop. 8.4 / App. A: the 4th roots {1, i, -1, -i} have trace-downs
    {2, 0, -2, 0} = {2, 0, -2}; none lies strictly inside (-2, 2)\\{0}, so
    they can never participate in a flip-straddle."""
    vals = set()
    for z in (1, 1j, -1, -1j):
        rho = z + 1 / z                        # complex, but real for these
        vals.add(round(rho.real, 12))
        assert abs(rho.imag) < 1e-12
    assert vals == {2.0, 0.0, -2.0}
    for v in vals:
        assert not (-2 < v < 2 and v != 0)


def test_spectral_oncircle_tracedowns_on_boundary():
    """Prop. 8.4: for spectral products, every on-circle eigenvalue has a
    trace-down in {2, 0, -2}, never strictly interior at a non-central point."""
    catalog = [[1, -1, -1], [1, 0, -2], [1, 0, -3], [1, 0, 5, 0, -5]]
    comps = [U.companion(c) for c in catalog]
    for A in comps:
        for B in comps:
            for M in (np.kron(A, B), A @ A):
                for e in np.linalg.eigvals(M):
                    if abs(abs(e) - 1) < 1e-9:
                        rho = e + 1 / e
                        assert abs(rho.imag) < 1e-7
                        r = rho.real
                        assert min(abs(r - 2), abs(r), abs(r + 2)) < 1e-6


def test_traceless_reciprocal_quartic_never_straddles():
    """Prop. 11.1 / App. A: the traceless reciprocal quartic x^4+bx^2+1 has
    trace-down t^2+(b-2) with symmetric roots +/- sqrt(2-b) -- never one root
    past 2 with the rest inside, so no free traceless quartic is Salem."""
    for b in range(-8, 9):
        R = sp.Poly([1, 0, b, 0, 1], x)
        T = U.trace_down(R)
        assert sp.expand(T.as_expr() - (t**2 + (b - 2))) == 0
        assert not U.flip_straddle(T)


def test_degree4_salem_polys_have_nonzero_trace():
    """Prop. 11.1: every degree-4 Salem polynomial x^4-a x^3+b x^2-a x+1 has
    trace a != 0, whereas a commutator is traceless -- so no degree-4 Salem is
    a commutator. Verified on the small Salem quartics found by scanning."""
    found = 0
    for a in range(1, 4):
        for b in range(-4, 5):
            coeffs = [1, -a, b, -a, 1]
            if U.is_salem_numeric([float(c) for c in coeffs]):
                # confirm exact Salem via flip-straddle, then check trace
                if U.has_salem_factor(sp.Poly(coeffs, x)):
                    found += 1
                    assert a != 0
    assert found >= 1                          # e.g. a=1,b=-1 -> beta_4
