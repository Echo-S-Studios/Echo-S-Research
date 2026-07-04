"""
Object III -- the gap / cost floor phi (Prop 4.1 `prop:basegap`,
Cor 4.2 `cor:floor`, and minimality Prop 7.1(3)).

Claims verified by an INDEPENDENT enumeration of monic integer quadratics:
  * no monic integer quadratic x^2+bx+c has Mahler measure in the open band (1,phi);
  * the smallest measure exceeding 1 is exactly phi, attained at x^2-x-1;
  * the sorted realised values begin  phi, 2, 1+sqrt2, phi^2, 1+sqrt3, 3;
  * phi is the smallest Perron number (Prop 7.1(3)), witnessed on the quadratics;
  * {1} U [phi,infty) is closed under multiplication and squaring (Cor 4.2 floor
    propagates).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import ge_helpers as H
from ge_helpers import x, as_poly, mahler, PHI

mp.mp.dps = 50
PHI_V = PHI()
BOUND = 6                       # paper's stated enumeration range |b|,|c| <= 6


def _quadratic_measures():
    """All (measure, b, c) with measure > 1 over |b|,|c| <= 6."""
    out = []
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            m = mahler(as_poly(x**2 + b * x + c))
            if m > 1 + mp.mpf(10)**(-30):
                out.append((m, b, c))
    return out


def test_band_1_to_phi_is_empty():
    """Prop 4.1: 'No monic integer quadratic has Mahler measure in (1,phi).'"""
    for m, b, c in _quadratic_measures():
        assert not (1 < m < PHI_V - mp.mpf(10)**(-30)), \
            f"x^2+{b}x+{c} has measure {m} inside the forbidden band"


def test_smallest_measure_is_phi_at_golden_polynomial():
    """Prop 4.1: 'the smallest measure exceeding 1 is exactly phi, at x^2-x-1.'"""
    ms = _quadratic_measures()
    m_min, b, c = min(ms, key=lambda t: t[0])
    assert abs(m_min - PHI_V) < mp.mpf(10)**(-30)
    assert (b, c) == (-1, -1)                       # x^2 - x - 1
    # independent closed form for phi
    assert abs(m_min - (1 + mp.sqrt(5)) / 2) < mp.mpf(10)**(-40)


def test_first_realized_values_match():
    """Prop 4.1: 'realised values begin phi, 2, 1+sqrt2, phi^2, 1+sqrt3, 3'."""
    vals = sorted({mp.nstr(m, 30) for m, _, _ in _quadratic_measures()},
                  key=lambda s: mp.mpf(s))
    first6 = [mp.mpf(v) for v in vals[:6]]
    expected = [PHI_V, mp.mpf(2), 1 + mp.sqrt(2), PHI_V**2, 1 + mp.sqrt(3), mp.mpf(3)]
    for got, exp in zip(first6, expected):
        assert abs(got - exp) < mp.mpf(10)**(-25)


def test_phi_is_smallest_perron_number_on_quadratics():
    """Prop 7.1(3): 'phi is the smallest Perron number.'  A Perron number is a
    real algebraic integer >1 strictly dominating its conjugate.  Among monic
    integer quadratics, the smallest such is phi (= root of x^2-x-1)."""
    from ge_helpers import roots_mp
    perron_vals = []
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            rts = roots_mp(as_poly(x**2 + b * x + c))
            if len(rts) != 2:
                continue
            reals = [r for r in rts if abs(r.imag) < mp.mpf(10)**(-20)]
            if len(reals) != 2:
                continue
            vals = sorted(abs(r) for r in reals)
            dom = max(reals, key=lambda r: abs(r))
            # Perron: dominant root real, >1, strictly bigger than the other |root|
            if dom.real > 1 and vals[1] > vals[0] + mp.mpf(10)**(-25):
                perron_vals.append(dom.real)
    assert perron_vals, "expected some Perron quadratics"
    assert abs(min(perron_vals) - PHI_V) < mp.mpf(10)**(-25)


def test_floor_set_closed_under_product_and_square():
    """Cor 4.2: the operators map {1} U [phi,infty) into itself (measures multiply
    under (+), square under ()^2), so the floor phi propagates at every size."""
    sample = [mp.mpf(1), PHI_V, PHI_V**2, mp.mpf(2), mp.mpf(5), 46 + mp.mpf('0.98')]

    def in_floor_set(v):
        return abs(v - 1) < mp.mpf(10)**(-25) or v >= PHI_V - mp.mpf(10)**(-25)

    for a in sample:
        assert in_floor_set(a * a)                 # squaring stays in the set
        for b in sample:
            assert in_floor_set(a * b)             # products stay in the set
    # and the band is genuinely a gap: nothing in the set lands in (1,phi)
    for a in sample:
        assert not (1 < a < PHI_V - mp.mpf(10)**(-25))
