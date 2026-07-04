r"""
Producer: the quadratic case of the emission gap (EG), Appendix D.

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/quadratic_eg_sweep.csv
              data/2026-07-relational-charge/quadratic_eg_floors.json

Refactors Appendix D: no charge-admissible object of degree <= 2 has Mahler
measure in the open gap (1, phi).  Emits an exhaustive sweep of the irreducible
quadratics x^2 - t x + n over a small integer box, recording for each its
discriminant type (real / complex pair), Mahler measure, and whether it lands in
(1, phi) -- none do.  The floor JSON records the four case values: the golden
attainment M = phi at (t,n)=(1,-1) [x^2-x-1], the n=+1 floor phi^2 = (3+sqrt5)/2,
and the identities (1+sqrt5)/2 = phi, (3+sqrt5)/2 = phi^2 (ledger R).

Run: py code/2026-07-relational-charge/quadratic_eg.py
"""

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv, write_json

x = C.x
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))


def sweep_rows(box=12):
    rows = []
    for t in range(-box, box + 1):
        for n in range(-box, box + 1):
            if n == 0:
                continue
            p = x**2 - t * x + n
            P = sp.Poly(p, x)
            if not P.is_irreducible:
                continue
            disc = t * t - 4 * n
            M = C.mahler_measure(p)
            in_gap = (M > 1 + mp.mpf(10) ** -20) and (M < PHI - mp.mpf(10) ** -20)
            rows.append({
                "t": t,
                "n": n,
                "discriminant": disc,
                "root_type": "real_pair" if disc > 0 else "complex_pair",
                "mahler": mp.nstr(M, 16),
                "mahler_ge_phi": "yes" if M >= PHI - mp.mpf(10) ** -18 else "no",
                "mahler_eq_1": "yes" if abs(M - 1) < mp.mpf(10) ** -18 else "no",
                "in_open_gap_1_phi": "yes" if in_gap else "no",
            })
    return rows


def floors():
    phi = (1 + sp.sqrt(5)) / 2
    return {
        "identities_ledger_R": {
            "(1+sqrt5)/2 = phi": bool(sp.simplify(phi - (1 + sp.sqrt(5)) / 2) == 0),
            "(3+sqrt5)/2 = phi^2": bool(sp.simplify(phi**2 - (3 + sp.sqrt(5)) / 2) == 0),
            "phi^2 = phi + 1": bool(sp.simplify(phi**2 - phi - 1) == 0),
        },
        "case_real_units_n_minus_1": {
            "M(t) = (|t| + sqrt(t^2+4))/2": True,
            "floor_at_abs_t_1": "phi = (1+sqrt5)/2",
            "attained_by": "x^2-x-1",
            "increasing_in_t": True,
        },
        "case_real_units_n_plus_1": {
            "M(t) = (|t| + sqrt(t^2-4))/2 for |t|>=3": True,
            "floor_at_abs_t_3": "phi^2 = (3+sqrt5)/2",
            "phi2_gt_phi": bool(phi**2 > phi),
        },
        "case_complex_pair": {
            "M = n = |alpha|^2": True,
            "n_1_root_of_unity_M_1": True,
            "n_ge_2_M_eq_n": True,
        },
        "conclusion": "every quadratic Mahler value lies in {1} U [phi, infinity)",
        "phi_attained_by": "x^2-x-1",
        "status": "[forced] (quadratic case of EG); general EG is [computed]/[plausible]",
    }


def main():
    rows = sweep_rows()
    fields = ["t", "n", "discriminant", "root_type", "mahler", "mahler_ge_phi",
              "mahler_eq_1", "in_open_gap_1_phi"]
    p1 = write_csv("quadratic_eg_sweep.csv", fields, rows, __file__)
    gap_hits = sum(1 for r in rows if r["in_open_gap_1_phi"] == "yes")
    print(f"wrote {p1}  ({len(rows)} irreducible quadratics, {gap_hits} in (1,phi))")

    p2 = write_json("quadratic_eg_floors.json", floors(), __file__)
    print(f"wrote {p2}  (Appendix D floor cases + ledger R identities)")


if __name__ == "__main__":
    main()
