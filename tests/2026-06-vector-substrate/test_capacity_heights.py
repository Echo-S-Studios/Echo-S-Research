"""
Independent verification of the height / capacity-gate claims
(Sec. 5: Thm. 5.2 (Landau/Northcott), Prop. 5.5 (admissible-set size),
Ex. 5.6 (Northcott count), Ex. 5.7 (exact heights), Rem. 5.9 (float gate)).
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import mpmath as mp
import vsub_nf as nf
from itertools import product

x = sp.symbols('x')


def test_exact_heights_table():
    """Ex. 5.7: (minpoly, deg, ch, sum c_i^2) for the four worked seeds."""
    rows = {
        "phi":    (x**2 - x - 1, 2, 1, 3),
        "2sqrt6": (x**2 - 24,    2, 24, 577),
        "sqrt7":  (x**2 - 7,     2, 7, 50),
        "phi^4":  (x**2 - 7*x + 1, 2, 7, 51),
    }
    for name, (poly, deg, ch, s2) in rows.items():
        assert sp.Poly(poly, x).degree() == deg, name
        assert nf.coeff_height(poly) == ch, name
        assert nf.two_norm_sq(poly) == s2, name


def test_phi4_minpoly_derivation():
    """Ex. 5.7: gap=phi^4 has minimal polynomial x^2-7x+1 (derived, not restated)."""
    phi = (1 + sp.sqrt(5)) / 2
    assert sp.minimal_polynomial(phi**4, x) == x**2 - 7*x + 1
    # phi^4 = 3 phi + 2 = (7+3 sqrt5)/2 ; trace 7, norm 1
    assert sp.simplify(phi**4 - (3*phi + 2)) == 0


def test_landau_inequality_and_tightness():
    """Thm. 5.2(b): Mahler(p) <= ||p||_2, hence Mahler^2 <= sum c_i^2 (integer);
    Ex. 5.7: for 2sqrt6, Mahler=24, 24^2=576 <= 577 = sum c_i^2."""
    for poly in [x**2 - x - 1, x**2 - 24, x**2 - 7, x**2 - 7*x + 1, x**4 - 10*x**2 + 1]:
        m = nf.mahler_measure_mp(poly, dps=45)
        assert m <= mp.sqrt(nf.two_norm_sq(poly)) + mp.mpf("1e-30")
        assert m**2 <= nf.two_norm_sq(poly) + mp.mpf("1e-25")
    # tightness for 2sqrt6
    m = nf.mahler_measure_mp(x**2 - 24, dps=45)
    assert abs(m - 24) < mp.mpf("1e-25")
    assert 24**2 == 576 and 576 <= 577 == nf.two_norm_sq(x**2 - 24)


def test_mahler_coefficient_bound():
    """Thm. 5.2 proof: |c_j| <= binom(d,j) Mahler(p) (Mahler's coefficient bound)."""
    for poly in [x**2 - 24, x**2 - 7, x**4 - 10*x**2 + 1, x**2 - 7*x + 1]:
        p = sp.Poly(poly, x)
        d = p.degree()
        M = nf.mahler_measure_mp(poly, dps=45)
        coeffs = p.all_coeffs()                  # [c_d=1, c_{d-1}, ..., c_0]
        for j in range(d + 1):
            cj = abs(int(coeffs[d - j]))         # coefficient of x^j
            assert cj <= math.comb(d, j) * M + mp.mpf("1e-25")


def test_northcott_count_D2_H1():
    """Ex. 5.6 / Prop. 5.5: |{monic, deg<=2, ch<=1}| = sum_{d=1}^2 3^d = 12,
    matching the explicit enumeration (3 linear + 9 quadratic)."""
    Dmax, Hmax = 2, 1
    formula = sum((2*Hmax + 1)**d for d in range(1, Dmax + 1))
    assert formula == 12
    # explicit enumeration of monic integer polys, non-leading coeffs in [-Hmax,Hmax]
    count = 0
    coeff_range = range(-Hmax, Hmax + 1)
    for d in range(1, Dmax + 1):
        count += sum(1 for _ in product(coeff_range, repeat=d))
    assert count == 12
    linear = list(product(coeff_range, repeat=1))
    quad = list(product(coeff_range, repeat=2))
    assert len(linear) == 3 and len(quad) == 9


def test_northcott_count_general_formula():
    """Prop. 5.5: size formula holds for another budget too (D=3,H=2 -> 5+25+125)."""
    Dmax, Hmax = 3, 2
    formula = sum((2*Hmax + 1)**d for d in range(1, Dmax + 1))
    assert formula == 5 + 25 + 125 == 155
    count = sum(len(list(product(range(-Hmax, Hmax + 1), repeat=d)))
                for d in range(1, Dmax + 1))
    assert count == formula


def test_northcott_irreducibility_classification():
    """Ex. 5.6: among the nine quadratics x^2+bx+c (b,c in {-1,0,1}):
    x^2-x-1 and x^2+x-1 irreducible; x^2-1, x^2+x reducible; x^2+1 cuts out Q(i)."""
    assert sp.Poly(x**2 - x - 1, x).is_irreducible
    assert sp.Poly(x**2 + x - 1, x).is_irreducible
    assert not sp.Poly(x**2 - 1, x).is_irreducible
    assert not sp.Poly(x**2 + x, x).is_irreducible
    assert sp.factor(x**2 - 1) == (x - 1) * (x + 1)
    assert sp.factor(x**2 + x) == x * (x + 1)
    assert sp.Poly(x**2 + 1, x).is_irreducible           # cuts out Q(i)
    assert sp.minimal_polynomial(sp.I, x) == x**2 + 1


def test_admissibility_gate_is_integer_and_budgeted():
    """Def. 5.4 / Ex. 5.6: 2sqrt6 (x^2-24) needs H_max>=24; admitted at 24,
    rejected at 23; the decision is pure integer (deg,ch)."""
    poly = x**2 - 24
    assert not nf.admissible(poly, Dmax=64, Hmax=23)
    assert nf.admissible(poly, Dmax=64, Hmax=24)
    # episode budget (64,256) admits it comfortably
    assert nf.admissible(poly, Dmax=64, Hmax=256)
    # degree budget bites too: deg 2 needs D_max>=2
    assert not nf.admissible(poly, Dmax=1, Hmax=256)


def test_float_mahler_never_crosses_integer_gate():
    """Rem. 5.9: Mahler(sqrt7)=7 exactly, but a naive float root-finder gives a
    value != 7 (here math.sqrt(7)**2 = 7.000000000000001, strictly above 7), so a
    floating '<=7' test is unreliable; the integer (deg,ch) gate is robust."""
    # exact Mahler is exactly 7
    assert abs(nf.mahler_measure_mp(x**2 - 7, dps=45) - 7) < mp.mpf("1e-30")
    # the specific float artifact quoted in the paper reproduces on this platform
    fv = math.sqrt(7) ** 2
    assert fv != 7.0
    assert fv > 7.0
    assert 0 < fv - 7.0 < 2e-15                           # ~1 ULP above 7
    assert repr(fv) == "7.000000000000001"
    # the exact integer gate admits sqrt7 regardless of the float noise
    assert nf.coeff_height(x**2 - 7) == 7
    assert nf.admissible(x**2 - 7, Dmax=64, Hmax=7)       # ch=7 <= 7, exact
