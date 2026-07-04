"""
Independent verification of the named constants and closed-form identities in
"The Charge-Measure Coupling on a Spectral Semiring" (v4).

Every value is rebuilt from its DEFINING polynomial / identity (never copied
from the paper's decimal) and only THEN compared to the paper's stated value.
"""

import mpmath as mp
import sympy as sp

mp.mp.dps = 60

s5 = sp.sqrt(5)
PHI = (1 + s5) / 2            # golden ratio, root of x^2 - x - 1
PHI_C = (1 - s5) / 2          # golden conjugate


# --------------------------------------------------------------------------
# Golden ratio and its conjugate  (Sec. 1.3 "golden seed", Sec. 4.3)
# --------------------------------------------------------------------------
def test_phi_is_root_of_x2_x_1():
    """Sec 1.4 / 4.3: phi = (1+sqrt5)/2 is the larger root of x^2 - x - 1."""
    assert sp.simplify(PHI**2 - PHI - 1) == 0


def test_phi_conjugate_equals_minus_inverse():
    """Sec 4.3: phi' = (1-sqrt5)/2 = -1/phi, and phi*phi' = -1."""
    assert sp.simplify(PHI_C - (-1 / PHI)) == 0
    assert sp.simplify(PHI * PHI_C + 1) == 0


def test_phi_conjugate_is_negative_at_argument_pi():
    """Lem 4.5 (pi-ray): phi' < 0, i.e. it lies at argument pi."""
    assert PHI_C.is_real
    assert (PHI_C < 0) == sp.true
    # numeric value ~ -0.618
    assert abs(mp.mpf(str(sp.N(PHI_C, 50))) - mp.mpf("-0.6180339887")) < 1e-9


def test_sqrt5_identity():
    """Sec 4.3: sqrt5 = phi + phi^{-1} = phi - phi'."""
    assert sp.simplify(PHI + 1 / PHI - s5) == 0
    assert sp.simplify(PHI - PHI_C - s5) == 0


def test_phi_squared_closed_form():
    """Lem 2.6 / ledger G: phi^2 = (3+sqrt5)/2 = phi + 1 ~ 2.618."""
    assert sp.simplify(PHI**2 - (3 + s5) / 2) == 0
    assert sp.simplify(PHI**2 - (PHI + 1)) == 0
    assert abs(mp.mpf(str(sp.N(PHI**2, 50))) - mp.mpf("2.618")) < 5e-4


def test_phi_fourth_closed_form():
    """Thm 6.5 / ledger I: phi^4 = (7+3sqrt5)/2 = 3phi+2 ~ 6.854."""
    assert sp.simplify(PHI**4 - (7 + 3 * s5) / 2) == 0
    assert sp.simplify(PHI**4 - (3 * PHI + 2)) == 0
    assert abs(mp.mpf(str(sp.N(PHI**4, 50))) - mp.mpf("6.854")) < 5e-4


# --------------------------------------------------------------------------
# Pentagon cosines  (Prop 6.1)
# --------------------------------------------------------------------------
def test_pentagon_cosines_closed_forms():
    """Prop 6.1: 2cos72 = phi-1 and 2cos144 = -phi."""
    assert sp.simplify(2 * sp.cos(2 * sp.pi / 5) - (PHI - 1)) == 0
    assert sp.simplify(2 * sp.cos(4 * sp.pi / 5) - (-PHI)) == 0


def test_pentagon_cosines_are_roots_of_x2_plus_x_minus_1():
    """Prop 6.1: phi-1 and -phi are the two roots of x^2 + x - 1,
    Galois-conjugate under sqrt5 -> -sqrt5."""
    x = sp.symbols("x")
    r1, r2 = PHI - 1, -PHI
    for r in (r1, r2):
        assert sp.simplify(r**2 + r - 1) == 0
    # sum = -1 (== -coeff), product = -1 (== const): both roots, distinct
    assert sp.simplify((r1 + r2) - (-1)) == 0
    assert sp.simplify(r1 * r2 - (-1)) == 0
    # Galois conjugation sqrt5 -> -sqrt5 swaps them
    swap = {s5: -s5}
    assert sp.simplify(r1.subs(swap) - r2) == 0


# --------------------------------------------------------------------------
# Plastic number mu_S  (Lem 2.5 Smyth)
# --------------------------------------------------------------------------
def test_plastic_number_value_and_defining_cubic():
    """Lem 2.5: mu_S = 1.3247... is the real root of x^3 - x - 1."""
    mu = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    assert abs(mu**3 - mu - 1) < mp.mpf(10) ** (-45)
    assert abs(mu - mp.mpf("1.3247")) < 5e-5           # paper's quoted digits


def test_plastic_number_is_pisot():
    """Lem 2.5: mu_S is the smallest Pisot number => real >1 with both
    conjugates strictly inside the unit circle."""
    rs = mp.polyroots([1, 0, -1, -1], extraprec=200)
    reals = [r for r in rs if abs(mp.im(r)) < 1e-40]
    comps = [r for r in rs if abs(mp.im(r)) >= 1e-40]
    assert len(reals) == 1 and mp.re(reals[0]) > 1
    assert len(comps) == 2
    for c in comps:
        assert abs(c) < 1            # Pisot: other conjugates inside unit disk


# --------------------------------------------------------------------------
# Lehmer's number tau  (Prop 8.2 / ledger L)
# --------------------------------------------------------------------------
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]   # x^10+x^9-x^7-...-x^3+x+1


def test_lehmer_number_value_in_gap():
    """Prop 8.2: Mahler measure of Lehmer's polynomial tau = 1.17628... in (1, phi)."""
    rs = mp.polyroots([mp.mpf(c) for c in LEHMER], extraprec=200)
    tau = mp.mpf(1)
    for r in rs:
        if abs(r) > 1:
            tau *= abs(r)
    assert abs(tau - mp.mpf("1.17628")) < 5e-6         # paper's quoted digits
    assert 1 < tau < mp.mpf(str(sp.N(PHI, 40)))        # strictly below phi


# --------------------------------------------------------------------------
# Minimal degree-4 Salem number beta_4  (ledger F)
# --------------------------------------------------------------------------
def test_beta4_value_in_odd_gap():
    """Ledger F: beta_4 = x^4-x^3-x^2-x+1 has Mahler 1.72208... in (phi, 2)."""
    rs = mp.polyroots([1, -1, -1, -1, 1], extraprec=200)
    m = mp.mpf(1)
    for r in rs:
        if abs(r) > 1:
            m *= abs(r)
    assert abs(m - mp.mpf("1.72208")) < 5e-6
    assert mp.mpf(str(sp.N(PHI, 40))) < m < 2          # inside the (phi,2) gap


def test_phi_numeric_quote():
    """Ledger H: phi = 1.61803..."""
    assert abs(mp.mpf(str(sp.N(PHI, 50))) - mp.mpf("1.61803")) < 5e-6
