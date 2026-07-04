"""
PRODUCER: heights, Mahler measure, Landau/Northcott, and the admissibility gate.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex, Sec. 5
  - Thm. 5.2 (Mahler = height^deg; Landau Mah <= ||p||_2; Northcott finiteness)
  - Prop. 5.5 (admissible-set size sum_{d=1}^{Dmax}(2Hmax+1)^d)
  - Ex. 5.6 (Northcott count 12; irreducibility of the nine quadratics)
  - Ex. 5.7 (exact heights table; phi^4 has minpoly x^2-7x+1)
  - Rem. 5.9 (the float Mahler 7.000000000000001 never crosses the integer gate)

Emits:
  data/heights.csv          -- (seed, minpoly, deg, ch, sum c_i^2, Mahler, Landau)
  data/northcott_counts.csv -- admissible-set sizes: formula vs enumeration
  data/admissibility.json   -- irreducibility classification + gate + float gate
"""
import math

import mpmath as mp
import sympy as sp
from itertools import product

import vsub_core as vc

x = vc.x
SCRIPT = "capacity_heights.py"


def heights_rows():
    """Ex. 5.7 + Thm. 5.2(b): exact (deg, ch, sum c_i^2) with Mahler measure and
    the Landau certificate Mah^2 <= sum c_i^2."""
    seeds = [
        ("phi", x**2 - x - 1),
        ("2sqrt6", x**2 - 24),
        ("sqrt7", x**2 - 7),
        ("phi^4", x**2 - 7 * x + 1),
    ]
    rows = []
    for name, poly in seeds:
        p = sp.Poly(poly, x)
        s2 = vc.two_norm_sq(poly)
        mah = vc.mahler_measure_mp(poly, dps=45)
        gap = s2 - mah ** 2                          # Landau slack: sum c_i^2 - Mahler^2 >= 0
        rows.append([
            name,
            sp.sstr(poly),
            p.degree(),
            vc.coeff_height(poly),
            s2,
            mp.nstr(mah, 12),
            mp.nstr(gap, 6),                         # 2sqrt6: 577-576 = 1.0 (Landau near-tight)
        ])
    return rows


def phi4_minpoly():
    """Ex. 5.7: phi^4 has minimal polynomial x^2 - 7x + 1 (= 3 phi + 2)."""
    phi = (1 + sp.sqrt(5)) / 2
    return {
        "phi4_minpoly": vc.sval(sp.minimal_polynomial(phi**4, x)),
        "phi4_equals_3phi_plus_2": bool(sp.simplify(phi**4 - (3 * phi + 2)) == 0),
    }


def mahler_coefficient_bound():
    """Thm. 5.2 proof: |c_j| <= binom(d,j) Mahler(p) for every coefficient."""
    out = {}
    for poly in [x**2 - 24, x**2 - 7, x**4 - 10 * x**2 + 1, x**2 - 7 * x + 1]:
        p = sp.Poly(poly, x)
        d = p.degree()
        M = vc.mahler_measure_mp(poly, dps=45)
        coeffs = p.all_coeffs()                  # [c_d=1, ..., c_0]
        checks = []
        for j in range(d + 1):
            cj = abs(int(coeffs[d - j]))
            checks.append(bool(cj <= math.comb(d, j) * M + mp.mpf("1e-25")))
        out[sp.sstr(poly)] = {"all_bounds_hold": all(checks)}
    return out


def northcott_rows():
    """Prop. 5.5 / Ex. 5.6: |{monic, deg<=Dmax, ch<=Hmax}| = sum (2Hmax+1)^d,
    formula vs explicit enumeration, for two budgets."""
    rows = []
    for Dmax, Hmax in [(2, 1), (3, 2)]:
        formula = sum((2 * Hmax + 1) ** d for d in range(1, Dmax + 1))
        count = sum(len(list(product(range(-Hmax, Hmax + 1), repeat=d)))
                    for d in range(1, Dmax + 1))
        rows.append([Dmax, Hmax, formula, count, "OK" if formula == count else "MISMATCH"])
    return rows


def irreducibility_classification():
    """Ex. 5.6: among the nine quadratics x^2+bx+c (b,c in {-1,0,1})."""
    out = {}
    for b in (-1, 0, 1):
        for c in (-1, 0, 1):
            poly = x**2 + b * x + c
            p = sp.Poly(poly, x)
            out[sp.sstr(poly)] = {
                "irreducible": bool(p.is_irreducible),
                "factored": sp.sstr(sp.factor(poly)),
            }
    return out


def admissibility_gate():
    """Def. 5.4 / Ex. 5.6: 2sqrt6 (x^2-24) admitted at Hmax>=24, rejected at 23;
    the decision is pure integer (deg, ch)."""
    poly = x**2 - 24
    return {
        "seed": "x^2 - 24 (2sqrt6)",
        "reject_at_Hmax_23": not vc.admissible(poly, Dmax=64, Hmax=23),
        "admit_at_Hmax_24": vc.admissible(poly, Dmax=64, Hmax=24),
        "admit_episode_budget_64_256": vc.admissible(poly, Dmax=64, Hmax=256),
        "reject_degree_budget_Dmax_1": not vc.admissible(poly, Dmax=1, Hmax=256),
        "coeff_height": vc.coeff_height(poly),
    }


def float_gate_artifact():
    """Rem. 5.9: Mahler(sqrt7)=7 exactly, but math.sqrt(7)**2 = 7.000000000000001
    (~1 ULP above 7), so a naive float '<=7' test is unreliable; the integer gate
    admits sqrt7 regardless."""
    fv = math.sqrt(7) ** 2
    return {
        "exact_mahler_sqrt7": 7,
        "float_sqrt7_squared_repr": repr(fv),
        "float_strictly_above_7": fv > 7.0,
        "ulp_gap": fv - 7.0,
        "integer_gate_admits_sqrt7": vc.admissible(x**2 - 7, Dmax=64, Hmax=7),
    }


def main():
    p1 = vc.write_csv(
        "heights.csv",
        ["seed", "min_poly", "degree", "coeff_height", "sum_ci_sq",
         "mahler_measure", "landau_gap_sumcisq_minus_M2"],
        heights_rows(),
        SCRIPT,
    )
    p2 = vc.write_csv(
        "northcott_counts.csv",
        ["Dmax", "Hmax", "formula_count", "enumerated_count", "check"],
        northcott_rows(),
        SCRIPT,
    )
    payload = {
        "phi4_minpoly": phi4_minpoly(),
        "mahler_coefficient_bound": mahler_coefficient_bound(),
        "irreducibility_classification": irreducibility_classification(),
        "admissibility_gate": admissibility_gate(),
        "float_gate_artifact": float_gate_artifact(),
    }
    p3 = vc.write_json("admissibility.json", payload, SCRIPT)
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print(f"wrote {p3}")
    print("  2sqrt6 sum c_i^2 = 577 (Landau tight: 24^2=576<=577)")
    print("  float artifact:", payload["float_gate_artifact"]["float_sqrt7_squared_repr"])


if __name__ == "__main__":
    main()
