"""
Independent verification of the floor / bound claims (Secs 2,4,6):
  * Real reciprocal-unit bound  M >= phi^2  (Lem 2.6).
  * Totally-positive (1,2) gap  (Thm 4.6 / ledger G)  -- finite-window replay.
  * Pentagon degree-4 sector    M in {1} U [phi^4, ...)  (Thm 6.5 / ledger I,J).
  * Realification bound 2^{1/5} < mu_S                    (Thm 6.4).
  * Emission-gap sanity on all constructed admissible objects (Prop 2.2).
"""

import sys
import pathlib
import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _cmc_helpers as H

mp.mp.dps = 50
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
TOL = mp.mpf(10) ** (-25)
x = sp.symbols("x")


# --------------------------------------------------------------------------
# Real reciprocal-unit bound  (Lem 2.6)
# --------------------------------------------------------------------------
def test_reciprocal_unit_trace_bound():
    """Lem 2.6: a real reciprocal pair {r,1/r}, r>1, with integer trace has
    trace>=3 (trace 2 forces r=1); hence r>=phi^2, Mahler contribution >=phi^2."""
    # trace 2  ->  x^2-2x+1=(x-1)^2  ->  r=1  (forces r=1, no unit outside circle)
    assert sp.solve(x**2 - 2 * x + 1, x) == [1]
    # trace 3  ->  r = (3+sqrt5)/2 = phi^2  (smallest admissible)
    r = max(sp.solve(x**2 - 3 * x + 1, x))
    assert sp.simplify(r - (1 + sp.sqrt(5)) ** 2 / 4) == 0     # = phi^2
    assert abs(mp.mpf(str(sp.N(r, 45))) - PHI**2) < TOL
    # for every integer trace t>=3 the larger root is >= phi^2
    for t in range(3, 12):
        rr = max(sp.solve(x**2 - t * x + 1, x))
        assert mp.mpf(str(sp.N(rr, 40))) >= PHI**2 - TOL


# --------------------------------------------------------------------------
# Totally-positive (1,2) gap  (Thm 4.6 / ledger G) -- finite-window replay
# --------------------------------------------------------------------------
def test_totally_positive_gap_finite_window():
    """Thm 4.6 / ledger G (computed, deg<=7 in paper; here deg<=4, |c|<=6):
    no totally-positive integer polynomial has Mahler measure in (1,2); the
    least value exceeding 1 is exactly 2, attained by (x-1)(x-2)."""
    K = 6
    measures = []
    # degrees 1..4
    ranges = {
        1: [(a,) for a in range(-K, K + 1)],
        2: [(a, b) for a in range(-K, K + 1) for b in range(-K, K + 1)],
        3: [(a, b, c) for a in range(-K, K + 1) for b in range(-K, K + 1)
            for c in range(-K, K + 1)],
        4: [(a, b, c, d) for a in range(-K, K + 1) for b in range(-K, K + 1)
            for c in range(-K, K + 1) for d in range(-K, K + 1)],
    }
    for deg, tail in ranges.items():
        for t in tail:
            coeffs = [1] + list(t)
            if coeffs[-1] == 0:
                continue                                # p(0)=0 excluded
            if not H.all_real_positive_np(coeffs):
                continue
            m = H.mahler_np(coeffs)
            measures.append(m)
            # the gap: nothing strictly inside (1,2)
            assert not (1 + 1e-6 < m < 2 - 1e-6), (coeffs, m)
    # 2 is attained (by x^2-3x+2 = (x-1)(x-2)) and is the least value > 1
    above_one = [m for m in measures if m > 1 + 1e-6]
    assert min(above_one) == min(above_one)            # non-empty
    assert abs(min(above_one) - 2.0) < 1e-6            # least > 1 is exactly 2


# --------------------------------------------------------------------------
# Pentagon degree-4 sector  (Thm 6.5)
# --------------------------------------------------------------------------
def test_pentagon_construction_expands_to_minimiser():
    """Thm 6.5: with s=phi^2 at +-72 and t=phi^{-2} at +-144, the Galois-coupled
    product (x^2 - s(phi-1)x + s^2)(x^2 + t*phi*x + t^2) = x^4-x^3+6x^2+4x+1."""
    phi = (1 + sp.sqrt(5)) / 2
    s = phi**2
    t = phi**-2
    O = (x**2 - s * (phi - 1) * x + s**2) * (x**2 + t * phi * x + t**2)
    target = x**4 - x**3 + 6 * x**2 + 4 * x + 1
    assert sp.simplify(sp.expand(O) - target) == 0


def test_pentagon_minimiser_root_geometry():
    """Thm 6.5: x^4-x^3+6x^2+4x+1 has a +-72 pair of modulus phi^2 and a
    +-144 pair of modulus phi^{-2}; hence M = (phi^2)^2 = phi^4."""
    rs = H.roots_mp([1, -1, 6, 4, 1])
    two_pi = 2 * mp.pi
    outer = [r for r in rs if abs(r) > 1]
    inner = [r for r in rs if abs(r) < 1]
    assert len(outer) == 2 and len(inner) == 2
    for r in outer:
        assert abs(abs(r) - PHI**2) < TOL
        assert abs(abs(mp.arg(r)) / two_pi - mp.mpf(1) / 5) < 1e-30    # +-72 deg
    for r in inner:
        assert abs(abs(r) - PHI**-2) < TOL
        assert abs(abs(mp.arg(r)) / two_pi - mp.mpf(2) / 5) < 1e-30    # +-144 deg
    assert abs(H.mahler([1, -1, 6, 4, 1]) - PHI**4) < TOL


def test_pentagon_quartic_scan_only_1_and_phi4():
    """Ledger J (computed, |c|<=10 in paper; here |c|<=6): the distinct Mahler
    measures among charge-Z/5 irreducible-or-not quartics are exactly {1, phi^4};
    none lands in (1, phi^4)."""
    K = 6
    phi4 = PHI**4
    found = set()
    for a in range(-K, K + 1):
        for b in range(-K, K + 1):
            for c in range(-K, K + 1):
                for d in range(-K, K + 1):
                    if d == 0:
                        continue
                    coeffs = [1, a, b, c, d]
                    rs = np.roots(coeffs)
                    # fast double-precision screen for charge Z/5 geometry
                    ok = True
                    for r in rs:
                        t5 = 5 * (np.angle(r) / (2 * np.pi))
                        if abs(t5 - round(t5)) > 1e-6:
                            ok = False
                            break
                    if not ok:
                        continue
                    # confirm at high precision
                    if H.charge_group(coeffs, Nmax=60) != 5:
                        continue
                    m = H.mahler(coeffs)
                    if abs(m - 1) < 1e-9:
                        found.add(1)
                    elif abs(m - phi4) < 1e-9:
                        found.add("phi4")
                    else:
                        found.add(("OTHER", coeffs, float(m)))
                    # gap: nothing in (1, phi^4)
                    assert not (1 + 1e-6 < m < float(phi4) - 1e-6), (coeffs, m)
    assert found == {1, "phi4"}                        # exactly {1, phi^4}


# --------------------------------------------------------------------------
# Realification bound 2^{1/5} improved to mu_S  (Thm 6.4)
# --------------------------------------------------------------------------
def test_realification_bound_below_muS():
    """Thm 6.4: the realification bound 2^{1/5}=1.1487 is strictly weaker than
    the Smyth floor mu_S=1.3247 that Thm 6.4 installs."""
    two_fifth = mp.mpf(2) ** (mp.mpf(1) / 5)
    muS = mp.findroot(lambda z: z**3 - z - 1, mp.mpf("1.3"))
    assert abs(two_fifth - mp.mpf("1.1487")) < 5e-5
    assert two_fifth < muS                              # 1.1487 < 1.3247


# --------------------------------------------------------------------------
# Emission-gap sanity on the constructed admissible objects  (Prop 2.2)
# --------------------------------------------------------------------------
def test_emission_gap_holds_on_constructions():
    """Prop 2.2: admissible objects satisfy M in {1} U [phi, inf).  Sanity
    check over every admissible construction used in the paper."""
    objs = [
        [1, 1, 1, 1, 1],                       # Phi_5 : M = 1
        [1, 0, 0, -2],                         # x^3-2 : M = 2
        [1, 0, 1, 0, -1],                      # q_2   : M = phi
        [1, -1, 6, 4, 1],                      # pentagon quartic : phi^4
        [1, 0, 0, 0, 0, -2],                   # x^5-2 : M = 2
        [1, -3, 1],                            # x^2-3x+1 : phi^2
    ]
    for c in objs:
        m = H.mahler(c)
        assert (abs(m - 1) < 1e-20) or (m >= PHI - 1e-20), (c, m)
