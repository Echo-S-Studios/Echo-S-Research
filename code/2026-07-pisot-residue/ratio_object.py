"""Producer: the ratio object Rat_p and its structural invariants.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: Proposition 2.1 (ratio object) -- degree n^2, integer coefficients,
(x-1)-multiplicity = sum m_alpha^2; the two-route companion identity
charpoly(C (x) C^{-1}) = Rat_p corroborating the ledger (Section 2); and the
Kronecker-square factorization charpoly(C_{x^4-x+1} (x) C_{x^4-x+1}) =
S_6^2 * (x^4+2x^2-x+1) quoted in Section 2.

Emits:
  data/2026-07-pisot-residue/ratio_object.json
"""
import sympy as sp
from sympy import symbols, Poly, expand, div, factor_list, kronecker_product, Matrix

from pisot_lib import (rat_object, cyclotomic_scan, multiplicity_of_factor,
                       companion, poly_str, scan_json, write_json)

x, y = symbols('x y')

# Canonical monic integer polynomials (squarefree, p(0) != 0) exercised in Sec 2.
CASES = {
    "x^2-x-1":         [1, -1, -1],
    "x^3-x-1":         [1, 0, -1, -1],
    "x^4-2":           [1, 0, 0, 0, -2],
    "x^4-x+1":         [1, 0, 0, -1, 1],
    "x^4-x^3-1":       [1, -1, 0, 0, -1],
    "Z*=x^4-3x^2+1":   [1, 0, -3, 0, 1],
    "x^5-2x^4-2x^3-2x^2-2x-2": [1, -2, -2, -2, -2, -2],
}


def ratio_object_record(coeffs):
    """Build Rat_p and report Proposition 2.1 invariants for one polynomial."""
    R = rat_object(coeffs)
    n = len(coeffs) - 1
    all_int = all(c.is_Integer for c in R.all_coeffs())
    mult1 = multiplicity_of_factor(R.as_expr(), x - 1)
    return {
        "poly": poly_str(coeffs),
        "coeffs_hi_to_lo": [int(c) for c in coeffs],
        "n": n,
        "deg_Rat_p": R.degree(),
        "deg_equals_n_squared": R.degree() == n * n,
        "integer_coefficients": bool(all_int),
        "x_minus_1_multiplicity": mult1,
        "squarefree_input_so_mult_equals_n": mult1 == n,
        "scan": scan_json(cyclotomic_scan(R)),
    }


def repeated_root_witness():
    """Prop 2.1 multiplicity claim with a NON-squarefree input:
    p = (x-1)^2 (x-3) has (x-1)-multiplicity in Rat_p equal to
    sum m_alpha^2 = 2^2 + 1^2 = 5 (not n = 3)."""
    p = expand((x - 1) ** 2 * (x - 3))                 # x^3 - 5x^2 + 7x - 3
    coeffs = [int(k) for k in Poly(p, x).all_coeffs()]
    R = rat_object(coeffs)
    return {
        "poly": "(x-1)^2*(x-3) = " + poly_str(coeffs),
        "coeffs_hi_to_lo": coeffs,
        "n": 3,
        "x_minus_1_multiplicity": multiplicity_of_factor(R.as_expr(), x - 1),
        "expected_sum_m_alpha_squared": 2 ** 2 + 1 ** 2,
    }


def two_route_identity():
    """Ledger corroboration (Sec 2): charpoly(C (x) C^{-1}) = Rat_p, since both
    carry the root multiset {alpha_i / alpha_j}.  The beta_4 mechanism."""
    rows = []
    for name, c in CASES.items():
        C = companion(c)
        K = Matrix(kronecker_product(C, C.inv()))
        cp = Poly(expand(K.charpoly(x).as_expr()), x).primitive()[1]
        R = rat_object(c)
        matches = (cp == R) or (cp == Poly(-R.as_expr(), x))
        rows.append({"poly": poly_str(c), "charpoly_Ckron_Cinv_equals_Rat_p": bool(matches)})
    return rows


def kronecker_square_x4_minus_x_plus_1():
    """Sec 2 quoted identity: charpoly(C (x) C)_{x^4-x+1}
    = S_6^2 * (x^4+2x^2-x+1) with S_6 = x^6-x^4-x^3-x^2+1 the Salem sextic."""
    c = [1, 0, 0, -1, 1]                                # x^4 - x + 1
    C = companion(c)
    K = Matrix(kronecker_product(C, C))
    cp = expand(K.charpoly(x).as_expr())
    quartic = x ** 4 + 2 * x ** 2 - x + 1
    q, r = div(cp, quartic, x)
    facs = factor_list(expand(q))[1]
    fac, mult = facs[0]
    sextic = [int(k) for k in Poly(fac, x).all_coeffs()]
    return {
        "input": "x^4-x+1",
        "quartic_factor_divides": sp.simplify(r) == 0,
        "quartic_factor": poly_str([1, 0, 2, -1, 1]),
        "sextic_factor": poly_str(sextic),
        "sextic_multiplicity": int(mult),
        "sextic_degree": Poly(fac, x).degree(),
        "sextic_irreducible": bool(Poly(fac, x).is_irreducible),
        "sextic_self_reciprocal": sextic == sextic[::-1],
        "sextic_is_S6": sextic == [1, 0, -1, -1, -1, 0, 1],
    }


def main():
    payload = {
        "description": "Proposition 2.1 ratio object Rat_p = prim Res_y(p(y),p(xy)); "
                       "Section 2 two-route and Kronecker-square corroborations.",
        "ratio_objects": [ratio_object_record(c) for c in CASES.values()],
        "repeated_root_multiplicity_witness": repeated_root_witness(),
        "two_route_companion_identity": two_route_identity(),
        "kronecker_square_x4_minus_x_plus_1": kronecker_square_x4_minus_x_plus_1(),
    }
    path = write_json("ratio_object.json", payload, "ratio_object.py")
    print(f"wrote {path}")
    for rec in payload["ratio_objects"]:
        print(f"  {rec['poly']:<28} deg Rat_p={rec['deg_Rat_p']:>3}  "
              f"(x-1)^{rec['x_minus_1_multiplicity']}  scan={rec['scan']['pretty']}")


if __name__ == "__main__":
    main()
