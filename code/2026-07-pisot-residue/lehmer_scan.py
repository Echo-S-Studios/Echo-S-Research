"""Producer: the canonical Lehmer scan (Section 2 ledger corroboration).

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the Section 2 canonical corroboration row -- Lehmer's degree-10 Salem
polynomial L has deg Rat_L = 100, scan bound 2 n^4 = 20000, and the complete
contact scan returns exactly {Phi_1^10}: Lehmer's number is relationally inert
(no non-trivial root-of-unity ratio among its conjugates).

Emits:
  data/2026-07-pisot-residue/lehmer_scan.json
"""
from sympy import symbols, Poly

from pisot_lib import rat_object, cyclotomic_scan, scan_json, poly_str, write_json

x = symbols('x')

# Lehmer polynomial x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]


def main():
    R = rat_object(LEHMER)
    n = len(LEHMER) - 1
    scan = cyclotomic_scan(R)
    payload = {
        "description": "Section 2 canonical Lehmer scan: Rat_L degree 100, scan {Phi_1^10}, "
                       "Lehmer's number relationally inert.",
        "polynomial": poly_str(LEHMER),
        "coeffs_hi_to_lo": LEHMER,
        "reciprocal_salem": LEHMER == LEHMER[::-1],
        "irreducible": bool(Poly(LEHMER, x).is_irreducible),
        "degree": n,
        "deg_Rat_L": R.degree(),
        "scan_bound_2n4": 2 * n ** 4,
        "scan": scan_json(scan),
        "relationally_inert": scan == {1: n},
    }
    path = write_json("lehmer_scan.json", payload, "lehmer_scan.py")
    print(f"wrote {path}")
    print(f"  deg Rat_L={R.degree()}, bound={2 * n ** 4}, scan={payload['scan']['pretty']}")


if __name__ == "__main__":
    main()
