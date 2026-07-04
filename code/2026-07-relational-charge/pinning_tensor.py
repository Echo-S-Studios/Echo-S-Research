r"""
Producer: modulus pinning, sharpness, and the tensor square (Sections 6-8).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/pinning_instances.csv
              data/2026-07-relational-charge/beta4_tensor_beta4.json

Refactors Theorem 6.15 (thm:pinning: an irreducible p with a uniquely attained
root modulus has no torsion ratio between distinct roots) and its corollaries on
Salem/Pisot inertness, the sharpness witness x^4-2 (all roots on one shell, so
pinning is silent and torsion ratios zeta_4^j appear), the twisted-shell
non-inert witness x^4+x^2+2 = q(x^2) (Example 7.21), and the non-inert tensor
square beta_4 (x) beta_4 (Example 8.5 / ledger S: charpoly F of the 16x16
Kronecker square has (x-1)-multiplicity 4, deg gcd(F,F')=7, exactly 3 distinct
positive real roots -- offset cancellation manufactures a rational block).

Run: py code/2026-07-relational-charge/pinning_tensor.py
"""

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv, write_json

x = C.x

# (name, poly, pinning applies?, expected verdict)
PIN = [
    ("beta4", C.B4, True),
    ("S6", C.S6, True),
    ("S8", C.S8, True),
    ("Lehmer", C.LEHMER, True),
    ("plastic x^3-x-1", C.PLASTIC, True),
    ("x^4-2 (sharpness)", x**4 - 2, False),
    ("x^4+x^2+2 (twisted shell)", C.TWISTSHELL, False),
]


def pin_rows():
    rows = []
    for name, p, pinned in PIN:
        P = sp.Poly(p, x)
        n = P.degree()
        mods = sorted((abs(r) for r in C.roots_mp(p)), reverse=True)
        dominant_gap = mods[0] - mods[1]
        unique_dominant = dominant_gap > mp.mpf(10) ** -12
        all_one_shell = (max(mods) - min(mods)) < mp.mpf(10) ** -18
        Rp = C.ratio_poly(p)
        sig = C.cyclotomic_contacts(Rp)
        rows.append({
            "object": name,
            "degree": n,
            "irreducible": "yes" if P.is_irreducible else "no",
            "dominant_modulus_gap": mp.nstr(dominant_gap, 6),
            "unique_dominant_modulus": "yes" if unique_dominant else "no",
            "all_roots_one_shell": "yes" if all_one_shell else "no",
            "pinning_applies": "yes" if pinned else "no",
            "contact_signature": C.signature_str(sig),
            "inert": "yes" if sig == {1: n} else "no",
        })
    return rows


def tensor_square():
    """Example 8.5 / ledger S: beta_4 (x) beta_4 is NON-inert."""
    Cm = C.companion_matrix(C.B4)
    K = sp.Matrix(sp.kronecker_product(Cm, Cm))
    F = sp.Poly(K.charpoly(x).as_expr(), x)
    Fp = F.diff(x)
    gcd_deg = sp.gcd(F, Fp).degree()
    sqfree = F.quo(sp.gcd(F, Fp))
    distinct_real = sqfree.count_roots(-sp.oo, sp.oo)
    distinct_pos = sqfree.count_roots(0, sp.oo)
    distinct_neg = sqfree.count_roots(-sp.oo, 0)

    # the rational block: positive-real pairwise products {1,1,1,1,tau^2,tau^-2}
    rts = C.roots_mp(C.B4)
    prods = [a * b for a in rts for b in rts]
    positive_reals = [q for q in prods
                      if abs(q.imag) < mp.mpf(10) ** -18 and q.real > mp.mpf(10) ** -18]
    ones = [q for q in positive_reals if abs(q - 1) < mp.mpf(10) ** -15]

    return {
        "object": "beta_4 (x) beta_4",
        "kronecker_charpoly_degree": F.degree(),
        "phi1_multiplicity": C.phi1_multiplicity(F),
        "deg_gcd_F_Fprime": gcd_deg,
        "multiplicity_pattern": "1^4 . (four doubles)",
        "distinct_real_roots": int(distinct_real),
        "distinct_positive_real_roots": int(distinct_pos),
        "distinct_negative_real_roots": int(distinct_neg),
        "rational_block_size": len(positive_reals),
        "diagonal_ones_count": len(ones),
        "rational_block": "{1,1,1,1, tau^2, tau^-2}",
        "inert": False,
        "mechanism": "offset cancellation (Prop 4.9 iii) manufactures a rational block from an inert factor",
        "status": "[forced] structure; single-engine run (ledger S)",
    }


def main():
    rows = pin_rows()
    fields = ["object", "degree", "irreducible", "dominant_modulus_gap",
              "unique_dominant_modulus", "all_roots_one_shell", "pinning_applies",
              "contact_signature", "inert"]
    p1 = write_csv("pinning_instances.csv", fields, rows, __file__)
    print(f"wrote {p1}  ({len(rows)} pinning instances, Thm 6.15 + sharpness)")

    p2 = write_json("beta4_tensor_beta4.json", tensor_square(), __file__)
    print(f"wrote {p2}  (Example 8.5 / ledger S: non-inert tensor square)")


if __name__ == "__main__":
    main()
