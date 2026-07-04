"""
Independent verification of Section 4 (the flip as a branch point):
the square-root edge series (Lem 4.1) and the monodromy swap (Sec 8/entry).
"""
import sympy as sp
import mpmath as mp

t = sp.symbols('t')


def test_sqrt_edge_series():
    """Lem 4.1: beta(t) = (t + sqrt(t^2-4))/2 expands, with s = sqrt(t-2), as
    beta = 1 + s + s^2/2 + s^3/8 + O(s^4),
    i.e. 1 + sqrt(t-2) + (t-2)/2 + (t-2)^{3/2}/8 + O((t-2)^2)."""
    s = sp.symbols('s', positive=True)
    t_of_s = 2 + s**2                       # so s = sqrt(t-2)
    beta = (t_of_s + sp.sqrt(t_of_s**2 - 4)) / 2
    ser = sp.series(beta, s, 0, 4).removeO()
    assert ser.coeff(s, 0) == 1
    assert ser.coeff(s, 1) == 1
    assert sp.nsimplify(ser.coeff(s, 2)) == sp.Rational(1, 2)
    assert sp.nsimplify(ser.coeff(s, 3)) == sp.Rational(1, 8)


def test_sqrt_edge_leading_asymptotic():
    """Lem 4.1: beta(t) - 1 ~ sqrt(t-2) as t -> 2+  (ratio -> 1)."""
    mp.mp.dps = 40
    def beta(tt):
        return (tt + mp.sqrt(tt**2 - 4)) / 2
    ratios = []
    for eps in [mp.mpf('1e-4'), mp.mpf('1e-6'), mp.mpf('1e-8')]:
        tt = 2 + eps
        ratios.append((beta(tt) - 1) / mp.sqrt(eps))
    # monotonically approaching 1 from above
    assert abs(ratios[-1] - 1) < 1e-3
    assert ratios[0] > ratios[-1] > 1


def test_log_growth_matches_edge():
    """Lem 4.1: log beta ~ beta - 1 ~ sqrt(t-2) near the flip."""
    mp.mp.dps = 40
    def beta(tt):
        return (tt + mp.sqrt(tt**2 - 4)) / 2
    tt = 2 + mp.mpf('1e-8')
    b = beta(tt)
    assert abs(mp.log(b) / (b - 1) - 1) < 1e-3
    assert abs((b - 1) / mp.sqrt(tt - 2) - 1) < 1e-3


def test_monodromy_swaps_beta_and_inverse():
    """Sec entry (monodromy): the two branches of x^2 - t x + 1 = 0 are
    beta and 1/beta (product of roots = constant term = 1). At t=2.3,
    beta=1.71789 and its partner 0.58211 = 1/beta. Encircling the branch
    point t=2 flips the sqrt sign, sending beta -> 1/beta."""
    mp.mp.dps = 40
    tt = mp.mpf('2.3')
    beta_plus = (tt + mp.sqrt(tt**2 - 4)) / 2
    beta_minus = (tt - mp.sqrt(tt**2 - 4)) / 2       # sign flip = monodromy image
    assert abs(beta_plus - mp.mpf('1.71789')) < 5e-6
    assert abs(beta_minus - mp.mpf('0.58211')) < 5e-6
    assert abs(beta_plus * beta_minus - 1) < mp.mpf(10)**(-30)   # product is 1
    assert abs(beta_minus - 1/beta_plus) < mp.mpf(10)**(-30)
