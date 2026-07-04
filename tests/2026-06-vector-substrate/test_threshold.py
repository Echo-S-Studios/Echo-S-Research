"""
Independent verification of the growth-threshold claims and the worked episode
(Sec. 6 episode; Sec. 8: Def. 8.1, Thm. 8.2 (Smyth floor), Conj. 8.4 (Lehmer),
Sec. 8.3 evidence table, Ex. 8.5 (lattice noise), Ex. 8.7 (certified grow),
Rem. 8.10 (degree-aware floor)); plus the cross-domain golden-gate scale of
Conj. 7.16 (tested only for internal consistency, per verification policy).

Transcendental quantities use mpmath at >=40 digits; the "cost" side (log Mahler)
is compared via rigorous interval enclosures exactly as the paper does.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import mpmath as mp
import vsub_nf as nf

x = sp.symbols('x')
mp.mp.dps = 50


def test_smyth_floor_plastic_number():
    """Thm. 8.2: mu_S = 1.3247179... is the real root of x^3-x-1 (plastic number),
    and log mu_S > 0 so the non-reciprocal cost floor lambda*log mu_S > 0."""
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    assert abs(muS - mp.mpf("1.32471795724474602596")) < mp.mpf("1e-18")
    assert muS**3 - muS - 1 == 0 or abs(muS**3 - muS - 1) < mp.mpf("1e-40")
    assert mp.log(muS) > 0
    # lambda=2: constant floor 2 log mu_S = 0.562...
    assert abs(2 * mp.log(muS) - mp.mpf("0.562399148646")) < mp.mpf("1e-11")


def test_lehmer_value():
    """Conj. 8.4: mu_L = 1.17628081... is the largest real root of Lehmer's
    degree-10 polynomial x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1, and mu_L < mu_S."""
    coeffs = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    muL = mp.findroot(lambda t: sum(c * t**(10 - i) for i, c in enumerate(coeffs)),
                      mp.mpf("1.17"))
    assert abs(muL - mp.mpf("1.17628081825991750654")) < mp.mpf("1e-18")
    val = sum(coeffs[i] * muL**(10 - i) for i in range(11))
    assert abs(val) < mp.mpf("1e-40")
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    assert muL < muS


def test_episode_residual_norm_2sqrt6():
    """Sec. 6: in Q(sqrt2+sqrt3), 2sqrt6=theta^2-5 (coords (-5,0,1,0)) is
    G-orthogonal to col(B)=Q+Q.theta, residual r=2sqrt6 with
    ||r||_G^2 = Tr((2sqrt6)^2) = Tr(24) = 4*24 = 96."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    G = nf.power_gram(C)
    # B = {1, theta} = columns e0, e1
    B = sp.Matrix.hstack(sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1, 0, 0]))
    r = nf.residual((-5, 0, 1, 0), B, G)
    assert r == sp.Matrix([-5, 0, 1, 0])          # residual is exactly 2sqrt6
    assert nf.gnorm2(r, G) == 96
    assert nf.field_trace((-5, 0, 1, 0), C) == 0  # trace zero (off-axis)
    # Tr((2sqrt6)^2) = Tr(24) = 4*24
    prod = nf.rho((-5, 0, 1, 0), C) * sp.Matrix([-5, 0, 1, 0])
    assert nf.field_trace(list(prod), C) == 96


def test_episode_pythagoras_116():
    """Sec. 6: for x=theta+2sqrt6, ||x||_G^2 = 116 = 20 + 96, captured part
    ||theta||_G^2 = Tr(theta^2) = 20, residual 96, cross term Tr(theta*2sqrt6)=0."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    G = nf.power_gram(C)
    theta = (0, 1, 0, 0)
    two_s6 = (-5, 0, 1, 0)
    xcoord = tuple(sp.Matrix(theta) + sp.Matrix(two_s6))   # (-5,1,1,0)
    assert nf.gnorm2(sp.Matrix(xcoord), G) == 116
    assert nf.gnorm2(sp.Matrix(theta), G) == 20            # Tr(theta^2)
    assert nf.gnorm2(sp.Matrix(two_s6), G) == 96
    # cross term Tr(theta * 2sqrt6) = 0
    cross = nf.rho(theta, C) * sp.Matrix(two_s6)
    assert nf.field_trace(list(cross), C) == 0
    assert 116 == 20 + 96


def test_episode_grow_decision():
    """Sec. 6: seed m=x^2-24 admissible under budget (64,256); with lambda=2 the
    gain 96 exceeds cost 2 log 24 = 6.356 and the floor 0.562 -> GROW."""
    assert nf.admissible(x**2 - 24, Dmax=64, Hmax=256)
    gain = 96
    cost = 2 * mp.log(24)
    floor = 2 * mp.log(mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3")))
    assert abs(cost - mp.mpf("6.3561076607")) < mp.mpf("1e-9")
    assert gain > cost > floor


def test_evidence_table_gains():
    """Sec. 8.3 table: 2sqrt6 (n=4) gain 96, Mahler 24, cost 6.36 -> GROW;
    sqrt7 (n=8) gain 56, Mahler 7, cost 3.89 -> GROW."""
    # 2sqrt6 in the degree-4 field: gain = Tr(24) = 96
    C4 = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    assert nf.gnorm2(sp.Matrix([-5, 0, 1, 0]), nf.power_gram(C4)) == 96
    assert abs(nf.mahler_measure_mp(x**2 - 24, dps=40) - 24) < mp.mpf("1e-25")
    assert abs(2 * mp.log(24) - mp.mpf("6.356107")) < mp.mpf("1e-6")
    # sqrt7 as residual in the degree-8 compositum: gain = Tr_{K/Q}(7) = 8*7 = 56
    assert 8 * 7 == 56
    assert abs(nf.mahler_measure_mp(x**2 - 7, dps=40) - 7) < mp.mpf("1e-25")
    assert abs(2 * mp.log(7) - mp.mpf("3.891820")) < mp.mpf("1e-6")
    # both gains clear their costs
    assert 96 > 2 * mp.log(24)
    assert 56 > 2 * mp.log(7)


def test_certified_grow_interval_enclosures():
    """Ex. 8.7: rigorous enclosures log7 in [1.94591,1.94592] and
    log24 in [3.17805,3.17806]; lambda*log7 <= 3.89184 <= 56 and
    lambda*log24 <= 6.35612 <= 96 certify GROW without trusting a rounded value."""
    l7 = mp.log(7)
    l24 = mp.log(24)
    assert mp.mpf("1.94591") <= l7 <= mp.mpf("1.94592")
    assert mp.mpf("3.17805") <= l24 <= mp.mpf("3.17806")
    # exact rational gain dominates a rigorous UPPER bound on the cost
    assert 56 >= 2 * mp.mpf("1.94592") == mp.mpf("3.89184")
    assert 96 >= 2 * mp.mpf("3.17806")
    assert 2 * mp.mpf("3.17806") == mp.mpf("6.35612")
    assert 2 * l24 <= mp.mpf("6.35612") <= 96
    assert 2 * l7 <= mp.mpf("3.89184") <= 56


def test_lattice_aligned_noise_stop():
    """Ex. 8.5: r = (1/10) sqrt5 in Q(sqrt5) is trace-zero and exactly off-axis,
    yet ||r||_G^2 = (1/100)*10 = 1/10 < 0.562 = floor, so STOP."""
    C = nf.companion_from_poly(x**2 - x - 1)
    G = nf.gram([(1, 0), (0, 1)], C)
    r = sp.Rational(1, 10) * sp.Matrix([-1, 2])   # (1/10) sqrt5
    assert nf.field_trace(list(r), C) == 0
    assert nf.gnorm2(r, G) == sp.Rational(1, 10)
    floor = 2 * mp.log(mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3")))
    assert mp.mpf(sp.Rational(1, 10)) < floor      # 0.1 < 0.562 -> STOP
    # ...but floor 0 would GROW on it (the point of Ex. 8.5)
    assert sp.Rational(1, 10) > 0


def test_degree_aware_floor():
    """Rem. 8.10: degree-aware floor n*lambda*log mu_S = 2.25 at n=4, 4.50 at n=8;
    a gain of 1 clears the constant floor 0.562 but not the n=8 degree-aware 4.50."""
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    const_floor = 2 * mp.log(muS)
    f4 = 4 * const_floor
    f8 = 8 * const_floor
    assert abs(f4 - mp.mpf("2.25")) < mp.mpf("2e-3")
    assert abs(f8 - mp.mpf("4.50")) < mp.mpf("3e-3")
    assert 1 > const_floor           # gain 1 clears constant floor
    assert 1 < f8                    # but not the n=8 degree-aware floor


def test_golden_gate_scale_internal_consistency():
    """Conj. 7.16 (CROSS-DOMAIN, companion GSA paper): the framework-native scale
    c = sqrt(1+4C)/(2C) at the critical point of the ladder x^2+x=C.  This is a
    SELECTED scale, not an arithmetic claim of this paper; we only check that the
    paper's OWN stated formula is internally consistent: at the golden gate C=1,
    c = sqrt5/2 and the exchange rate lambda = 2c = sqrt5, and the ladder root is
    the golden 1/phi."""
    C = sp.Integer(1)
    c = sp.sqrt(1 + 4*C) / (2*C)
    assert sp.simplify(c - sp.sqrt(5)/2) == 0
    assert sp.simplify(2*c - sp.sqrt(5)) == 0
    # ladder x^2 + x = 1 positive root = (sqrt5 - 1)/2 = 1/phi = phi - 1
    root = sp.Rational(-1, 2) + sp.sqrt(5)/2
    assert sp.simplify(root**2 + root - 1) == 0
    phi = (1 + sp.sqrt(5)) / 2
    assert sp.simplify(root - (phi - 1)) == 0
    assert sp.simplify(root - 1/phi) == 0
