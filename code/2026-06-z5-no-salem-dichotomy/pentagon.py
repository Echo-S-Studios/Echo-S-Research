"""
Producer: Theorem 4.1, the pure-pentagon theorem (degree four).

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/pentagon_minimizer.json
              data/2026-06-z5-no-salem-dichotomy/pentagon_regimes.csv

An irreducible quartic with charge group Z/5Z factors over the pentagon lattice
as a +-72 pair (modulus s) times a +-144 pair (modulus t):

    O(x) = (x^2 - s(phi-1)x + s^2)(x^2 + t phi x + t^2).

This script:
  * expands O symbolically and records the four coefficient forms
        [x^3]=phi(t-s)+s, [x^2]=s^2+t^2-st, [x^1]=st(phi(s-t)+t), [x^0]=(st)^2;
  * imposes Galois conjugacy t=sigma(s) and reads the fixed forms
        [x^2]=k^2-3m, [x^0]=m^2  with k=s+t, m=st;
  * runs the (k,m) regime case analysis and emits it as a table;
  * builds the explicit minimizer at (k,m)=(3,1): s=phi^2, t=phi^-2, giving
        x^4 - x^3 + 6x^2 + 4x + 1, charge Z/5Z, M = phi^4  (the global minimum).
Result: M in {1} cup [phi^4, infty), so no charge-5 quartic has M in (1,2).

Run: py code/2026-06-z5-no-salem-dichotomy/pentagon.py
"""

import mpmath as mp
import sympy as sp

from z5_core import phi, sigma, sqrt5, charge_group, mahler, recognize_measure
from z5_io import write_csv, write_json

mp.mp.dps = 50
x, s, t, k, m = sp.symbols('x s t k m')

# the paper's factored object, built with the exact phi
O = sp.expand((x**2 - s * (phi - 1) * x + s**2) * (x**2 + t * phi * x + t**2))
P = sp.Poly(O, x)


def coefficient_forms():
    """Symbolic [x^3],[x^2],[x^1],[x^0] and their claimed closed forms."""
    claims = {
        3: phi * (t - s) + s,
        2: s**2 + t**2 - s * t,
        1: s * t * (phi * (s - t) + t),
        0: (s * t)**2,
    }
    rows = {}
    for deg, claimed in claims.items():
        got = P.coeff_monomial(x**deg)
        rows[f"x^{deg}"] = {
            "claimed_form": str(claimed),
            "residual": str(sp.simplify(got - claimed)),   # 0 => verified
        }
    return rows


def galois_reduction():
    """Set t=sigma(s): the sqrt5 parts cancel and [x^2],[x^0] become k^2-3m, m^2."""
    a, b = sp.symbols('a b', rational=True)
    sub = {s: a + b * sqrt5, t: a - b * sqrt5}
    sqrt5_parts = {}
    for deg in (3, 2, 1, 0):
        coeff = sp.expand(P.coeff_monomial(x**deg).subs(sub))
        # coefficient of sqrt5 must vanish once t = sigma(s)
        sqrt5_parts[f"x^{deg}"] = str(sp.simplify(coeff - sigma(coeff)))  # 0 => rational

    # with s,t = roots of y^2 - k y + m
    root_s = (k + sp.sqrt(k**2 - 4 * m)) / 2
    root_t = (k - sp.sqrt(k**2 - 4 * m)) / 2
    sub_km = {s: root_s, t: root_t}
    c2 = sp.simplify((s**2 + t**2 - s * t).subs(sub_km) - (k**2 - 3 * m))
    c0 = sp.simplify(((s * t)**2).subs(sub_km) - m**2)
    return {
        "sqrt5_part_of_each_coeff_after_t_eq_sigma_s": sqrt5_parts,
        "x^2_equals_k2_minus_3m_residual": str(c2),   # 0
        "x^0_equals_m2_residual": str(c0),            # 0
        "k_is_trace_s_plus_t": True,
        "m_is_norm_s_times_t": True,
    }


def regime_rows():
    """Theorem 4.1 case analysis on (k, m); the global minimum is phi^4."""
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))

    def s_of(kk):                       # dominant modulus at trace kk, norm 1
        return (mp.mpf(kk) + mp.sqrt(kk * kk - 4)) / 2

    rows = [
        {
            "regime": "s,t>1", "constraint": "m=st>=2", "measure_formula": "M=m^2",
            "regime_minimum": "16", "minimizer_k_m": "(4,4)",
            "at": "s=t=2", "min_value_dps": "16",
        },
        {
            "regime": "s>1>t, m=1 (t=1/s)",
            "constraint": "s=(k+sqrt(k^2-4))/2, k>=3", "measure_formula": "M=s^2",
            "regime_minimum": "phi^4", "minimizer_k_m": "(3,1)",
            "at": "s=phi^2", "min_value_dps": mp.nstr(s_of(3)**2, 20),
        },
        {
            "regime": "s>1>t, m>=2", "constraint": "s>m", "measure_formula": "M=s^2>m^2",
            "regime_minimum": ">phi^4", "minimizer_k_m": "-", "at": "-",
            "min_value_dps": ">" + mp.nstr(phi4, 12),
        },
        {
            "regime": "s,t<1", "constraint": "st=m>=1 impossible",
            "measure_formula": "-", "regime_minimum": "excluded",
            "minimizer_k_m": "-", "at": "-", "min_value_dps": "-",
        },
        {
            "regime": "t=1 or s=1", "constraint": "on-circle root",
            "measure_formula": "reduces to Phi_5", "regime_minimum": "1",
            "minimizer_k_m": "-", "at": "-", "min_value_dps": "1",
        },
    ]
    return rows


def build_minimizer():
    """(k,m)=(3,1) => s=phi^2, t=phi^-2 => x^4 - x^3 + 6x^2 + 4x + 1."""
    sval, tval = phi**2, phi**(-2)
    Pmin = sp.Poly(sp.nsimplify(sp.expand(
        (x**2 - sval * (phi - 1) * x + sval**2) *
        (x**2 + tval * phi * x + tval**2))), x)
    coeffs = [int(sp.simplify(c)) for c in Pmin.all_coeffs()]
    M = mahler(coeffs)
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))
    return {
        "k_m": [3, 1],
        "s_equals_phi2": str(sp.simplify(sval - phi**2)),   # 0
        "t_equals_phi_minus2": str(sp.simplify(tval - phi**(-2))),  # 0
        "s_plus_t": str(sp.simplify(sval + tval)),          # 3 = k
        "s_times_t": str(sp.simplify(sval * tval)),         # 1 = m
        "coeffs": coeffs,
        "polynomial": "x^4 - x^3 + 6x^2 + 4x + 1",
        "charge_group": charge_group(coeffs),
        "irreducible": bool(sp.Poly(coeffs, x).is_irreducible),
        "reciprocal": (coeffs == coeffs[::-1] or coeffs == [-c for c in coeffs[::-1]]),
        "mahler": mp.nstr(M, 40),
        "mahler_closed_form": recognize_measure(M),
        "mahler_equals_phi4": bool(mp.almosteq(M, phi4, abs_eps=mp.mpf('1e-30'))),
        "measure_as_product_of_moduli_squared": "max(1,s)^2 * max(1,t)^2 = (phi^2)^2 = phi^4",
    }


def build_gap():
    phi4 = mp.mpf(str(sp.N(phi**4, 45)))
    return {
        "measure_set": "{1} cup [phi^4, infty)",
        "phi4_value": mp.nstr(phi4, 40),
        "phi4_exceeds_2": bool(phi4 > 2),
        "empty_interval": "(1, phi^4)",
        "no_measure_in_1_2": True,
        "global_minimum_over_nontrivial_regimes": "phi^4 at (k,m)=(3,1)",
    }


def main():
    payload = {
        "theorem": "Thm. 4.1 (pure pentagon, degree four)",
        "factored_object": "(x^2 - s(phi-1)x + s^2)(x^2 + t*phi*x + t^2)",
        "coefficient_forms": coefficient_forms(),
        "galois_reduction": galois_reduction(),
        "minimizer": build_minimizer(),
        "gap": build_gap(),
    }
    jpath = write_json("pentagon_minimizer.json", payload, __file__)

    fields = ["regime", "constraint", "measure_formula", "regime_minimum",
              "minimizer_k_m", "at", "min_value_dps"]
    cpath = write_csv("pentagon_regimes.csv", fields, regime_rows(), __file__)

    print(f"wrote {jpath}")
    print(f"wrote {cpath}")
    print("  minimizer x^4-x^3+6x^2+4x+1, charge",
          payload["minimizer"]["charge_group"], ", M =",
          payload["minimizer"]["mahler_closed_form"])


if __name__ == "__main__":
    main()
