"""Producer: the hypothesis-necessity witness Z* and the x^4-2 sharpness example.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: Proposition 3.4 (the Z* witness) and Remark 3.5 (parent x^4-2 example).

  Z* = f1 f2 = (x^2-x-1)(x^2+x-1) = x^4-3x^2+1 carries the torsion ratio -1
  between phi (root of f1) and -phi (root of f2) at the shared modulus phi;
  its complete scan is {Phi_1^4, Phi_2^4}.  This shows irreducibility cannot be
  dropped from the modulus-pinning theorem, and the disjoint-modulus clause
  cannot be dropped from its mixed form.

  x^4-2 is irreducible with all four roots on one shell 2^{1/4}, ratios = powers
  of i, so its scan is {Phi_1^4, Phi_2^4, Phi_4^4}: the pinning hypothesis
  (uniquely-attained modulus) is necessary for an irreducible input.

Emits:
  data/2026-07-pisot-residue/zstar_sharpness.json
"""
import sympy as sp
from sympy import symbols, expand, simplify, Poly, sqrt, roots, Abs

from pisot_lib import rat_object, cyclotomic_scan, scan_json, poly_str, write_json

x = symbols('x')
PHI = (1 + sqrt(5)) / 2


def zstar_record():
    f1 = x ** 2 - x - 1
    f2 = x ** 2 + x - 1
    prod = expand(f1 * f2)
    r1 = list(roots(f1, x).keys())
    r2 = list(roots(f2, x).keys())
    ratio = simplify(PHI / (-PHI))                 # phi / (-phi) = -1
    R = rat_object([1, 0, -3, 0, 1])
    scan = cyclotomic_scan(R)
    return {
        "Z_star": poly_str([1, 0, -3, 0, 1]),
        "factorization": "(x^2-x-1)(x^2+x-1)",
        "factorization_exact": expand(prod) == x ** 4 - 3 * x ** 2 + 1,
        "f1_irreducible": bool(Poly(f1, x).is_irreducible),
        "f2_irreducible": bool(Poly(f2, x).is_irreducible),
        "f1_roots": [str(r) for r in r1],
        "f2_roots": [str(r) for r in r2],
        "shared_modulus_phi": simplify(Abs(PHI) - Abs(-PHI)) == 0,
        "cross_factor_ratio_phi_over_minus_phi": str(ratio),
        "ratio_is_primitive_2nd_root_of_unity": (ratio == -1) and (simplify(ratio ** 2 - 1) == 0),
        "deg_Rat": R.degree(),
        "scan": scan_json(scan),
        "scan_is_phi1_4_phi2_4": scan == {1: 4, 2: 4},
    }


def x4_minus_2_record():
    R = rat_object([1, 0, 0, 0, -2])
    scan = cyclotomic_scan(R)
    return {
        "polynomial": poly_str([1, 0, 0, 0, -2]),
        "irreducible": bool(Poly(x ** 4 - 2, x).is_irreducible),
        "deg_Rat": R.degree(),
        "scan": scan_json(scan),
        "scan_is_phi1_4_phi2_4_phi4_4": scan == {1: 4, 2: 4, 4: 4},
    }


def main():
    payload = {
        "description": "Proposition 3.4 witness Z* = x^4-3x^2+1 (scan {Phi_1^4,Phi_2^4}) "
                       "and Remark 3.5 parent example x^4-2 (scan {Phi_1^4,Phi_2^4,Phi_4^4}); "
                       "together they isolate every hypothesis of the modulus-pinning theorem.",
        "Z_star_witness": zstar_record(),
        "x4_minus_2_sharpness": x4_minus_2_record(),
    }
    path = write_json("zstar_sharpness.json", payload, "zstar_sharpness.py")
    print(f"wrote {path}")
    print(f"  Z*    scan={payload['Z_star_witness']['scan']['pretty']}")
    print(f"  x^4-2 scan={payload['x4_minus_2_sharpness']['scan']['pretty']}")


if __name__ == "__main__":
    main()
