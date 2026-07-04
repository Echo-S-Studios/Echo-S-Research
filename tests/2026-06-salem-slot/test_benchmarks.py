"""
Independent verification of the benchmark table in Section 5 (and the
Lehmer numeric in Sec 4). For each Salem the redirection tau0, the
trace-down T, the Mahler measure Mah(T) and the log(beta)->log(tau0)
entries are rebuilt from scratch and compared to the paper's stated values.
"""
import sympy as sp
import mpmath as mp

t = sp.symbols('t')
mp.mp.dps = 45


def _mahler_of_poly(coeffs):
    """Mahler measure of a monic integer polynomial from its roots:
    prod max(1,|root|)."""
    roots = [complex(r) for r in sp.nroots(sp.Poly(coeffs, t), n=40)]
    m = 1.0
    for r in roots:
        m *= max(1.0, abs(r))
    return m


def _dominant_real_root(poly):
    roots = [complex(r) for r in sp.nroots(poly, n=40)]
    return max(r.real for r in roots if abs(r.imag) < 1e-20)


def test_lehmer_benchmark():
    """Sec 4/5: Lehmer beta=1.1762808, trace-down t^5+t^4-5t^3-5t^2+4t+3,
    tau0=2.026418, Mah(T)=5.615601, logs 0.162->0.706, and the Prop 4.2
    identity tau0-2=(beta-1)^2/beta=0.026418."""
    # beta = largest root of Lehmer's polynomial (built independently)
    lehmer = mp.polyroots([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], maxsteps=200, extraprec=80)
    beta = max((r.real for r in lehmer if abs(r.imag) < 1e-20), default=None)
    beta = mp.mpf(beta)
    assert abs(beta - mp.mpf('1.1762808')) < 5e-7

    tau0 = beta + 1/beta
    assert abs(tau0 - mp.mpf('2.026418')) < 5e-7

    # trace-down must reconstruct Lehmer, and its dominant root must equal tau0
    T = t**5 + t**4 - 5*t**3 - 5*t**2 + 4*t + 3
    assert abs(_dominant_real_root(T) - float(tau0)) < 1e-9

    mah = _mahler_of_poly([1, 1, -5, -5, 4, 3])
    assert abs(mah - 5.615601) < 5e-6

    # Prop 4.2 numeric
    assert abs((beta - 1) - mp.mpf('0.176281')) < 5e-6
    assert abs((beta - 1)**2 / beta - mp.mpf('0.026418')) < 5e-7
    assert abs((tau0 - 2) - (beta - 1)**2 / beta) < mp.mpf(10)**(-30)

    assert abs(mp.log(beta) - mp.mpf('0.162')) < 5e-4
    assert abs(mp.log(tau0) - mp.mpf('0.706')) < 5e-4


def test_beta4_benchmark():
    """Sec 5: deg-4 Salem, trace-down t^2-t-3, tau0=2.302776,
    Mah(T)=3.000000, beta=1.7220838, logs 0.544->0.834. Cross-check that
    the lift of tau0 is a root of x^4-x^3-x^2-x+1."""
    # tau0 = dominant root of the paper's stated trace-down t^2 - t - 3
    tau0 = mp.findroot(lambda z: z**2 - z - 3, mp.mpf('2.3'))
    assert abs(tau0 - mp.mpf('2.302776')) < 5e-7

    beta = (tau0 + mp.sqrt(tau0**2 - 4)) / 2
    assert abs(beta - mp.mpf('1.7220838')) < 5e-7

    mah = _mahler_of_poly([1, -1, -3])
    assert abs(mah - 3.0) < 5e-7

    # cross-consistency: beta solves the degree-4 Salem polynomial
    x = sp.symbols('x')
    salem4 = beta**4 - beta**3 - beta**2 - beta + 1
    assert abs(salem4) < mp.mpf(10)**(-25)

    assert abs(mp.log(beta) - mp.mpf('0.544')) < 5e-4
    assert abs(mp.log(tau0) - mp.mpf('0.834')) < 5e-4


def test_deg6_benchmark():
    """Sec 5: deg-6 Salem, trace-down t^3-t^2-3t+1, tau0=2.170086,
    Mah(T)=3.214320, beta=1.5061357, logs 0.410->0.775. Cross-check that
    the lift of tau0 is a root of x^6-x^5-x^3-x+1 (the S_6 Salem factor)."""
    # tau0 = dominant root of the paper's stated trace-down t^3 - t^2 - 3t + 1
    tau0 = mp.findroot(lambda z: z**3 - z**2 - 3*z + 1, mp.mpf('2.17'))
    assert abs(tau0 - mp.mpf('2.170086')) < 5e-7

    beta = (tau0 + mp.sqrt(tau0**2 - 4)) / 2
    assert abs(beta - mp.mpf('1.5061357')) < 5e-7

    mah = _mahler_of_poly([1, -1, -3, 1])
    assert abs(mah - 3.214320) < 5e-6

    # cross-consistency with the upstairs Salem factor x^6-x^5-x^3-x+1
    salem6 = beta**6 - beta**5 - beta**3 - beta + 1
    assert abs(salem6) < mp.mpf(10)**(-25)

    assert abs(mp.log(beta) - mp.mpf('0.410')) < 5e-4
    assert abs(mp.log(tau0) - mp.mpf('0.775')) < 5e-4
