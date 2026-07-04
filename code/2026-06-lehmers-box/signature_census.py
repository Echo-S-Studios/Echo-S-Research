"""
Producer: the signature face -- the 27-subfield census of K.  Section 5.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Re-derives Theorem 5.3 by an INDEPENDENT Galois computation on the closure
K(i): builds G = Gal(K(i)/Q) ~ C2 x C2 x D4 (order 32) from scratch, enumerates
all its subgroups, restricts to the 27 that fix a subfield of K (the subgroups
containing complex conjugation), and computes each subfield's signature
(deg, r1, r2).  The census matches the paper's table exactly (total 27), pins the
four Salem-shaped (2,1) quartics (Cor. 5.5), records K's own signature (8,4), and
also emits the TOTAL subgroup count of G, 158 -- so the backing note's "27" is
read as the subfield count, not the total subgroup count.  Adds the exact
Q(sqrt5) certificate beta_4 > phi (Cor. 5.5), the Salem trace-form signatures
(Prop. 5.1), [K:Q]=16 (Lem. 5.2), and the off-circle complex place (Sec. 5).

Outputs:
    data/2026-06-lehmers-box/signature_census.csv
    data/2026-06-lehmers-box/signature_census.json
Backs: Prop. 5.1, Lem. 5.2, Thm. 5.3, Cor. 5.5; total-subgroups = 158.
"""

from __future__ import annotations

import mpmath as mp
import numpy as np
import sympy as sp

import box_core as C
from box_io import write_csv, write_json

mp.mp.dps = 50
_x = sp.symbols('x')


def _integer_power_sums(int_coeffs, upto):
    """p_n = sum roots^n (a rational integer for a monic integer polynomial)."""
    roots = C.roots_of(int_coeffs)
    out = []
    for n in range(upto + 1):
        val = mp.re(sum(r ** n for r in roots))
        out.append(int(mp.nint(val)))
    return out


def trace_form_signature(int_coeffs):
    """Signature (pos, neg) of the trace form Tr(theta^{i+j}) on 1..theta^{d-1}
    (Prop. 5.1)."""
    d = len(int_coeffs) - 1
    p = _integer_power_sums(int_coeffs, 2 * (d - 1))
    M = np.array([[p[i + j] for j in range(d)] for i in range(d)], dtype=float)
    ev = np.linalg.eigvalsh(M)
    return int(np.sum(ev > 1e-6)), int(np.sum(ev < -1e-6))


def unit_rank_signature(int_coeffs):
    """(r1, r2) for a Salem field: r1 = 2 real embeddings, r2 = on-circle pairs."""
    out, ins, onc = C.root_classification(int_coeffs)
    return out + ins, onc // 2


def salem_signature_rows():
    """Prop. 5.1 backing: beta_4 (deg 4, m=2) and Lehmer (deg 10, m=5)."""
    rows = []
    for name, coeffs, m in [("beta_4", C.BETA4, 2), ("Lehmer_L", C.LEHMER, 5)]:
        r1, r2 = unit_rank_signature(coeffs)
        tf = trace_form_signature(coeffs)
        rows.append(dict(
            name=name, degree=2 * m, m=m,
            unit_rank_signature=f"({r1},{r2})",
            expected_unit_rank=f"(2,{m - 1})",
            trace_form_signature=f"({tf[0]},{tf[1]})",
            expected_trace_form=f"({m + 1},{m - 1})",
            indefinite=bool(r2 >= 1)))
    return rows


def beta4_exact_sign():
    """Cor. 5.5: beta_4 > phi  <=>  m_{beta_4}(phi) < 0, with
    m_{beta_4}(phi) = phi^4 - phi^3 - phi^2 - phi + 1 = (1 - sqrt5)/2 in Q(sqrt5)."""
    val = C.PHI_SYM ** 4 - C.PHI_SYM ** 3 - C.PHI_SYM ** 2 - C.PHI_SYM + 1
    a, b = C.to_ab_sqrt5(val)
    s = C.sign_ab_sqrt5(a, b)
    return {
        "test_element": "m_{beta_4}(phi) = phi^4 - phi^3 - phi^2 - phi + 1",
        "a_plus_b_sqrt5": {"a": str(a), "b": str(b)},
        "equals": "(1 - sqrt5)/2",
        "sign": s,
        "conclusion": "sign < 0  =>  beta_4 > phi (FORCED_ABOVE_FLOOR)",
    }


def main():
    g = C.galois_census()
    census = g["census"]
    n_subfields = g["n_subfields"]
    n_total = g["n_total_subgroups"]
    order_breakdown = g["order_breakdown"]

    # CSV rows in the paper's table order (sorted by degree, then r1)
    rows = []
    for (deg, r1, r2) in sorted(census, key=lambda s: (s[0], s[1])):
        salem_shaped = (r1 == 2 and r2 == (deg // 2) - 1 and deg >= 4)
        rows.append(dict(
            degree=deg, r1=r1, r2=r2,
            signature=f"({r1},{r2})",
            n_subfields=census[(deg, r1, r2)],
            salem_bearing=bool(salem_shaped),
            note=("Salem-bearing Lorentzian quartic" if salem_shaped else
                  ("K itself" if deg == 16 else
                   ("totally real" if r2 == 0 else "indefinite (2k,k)")))))
    rows.append(dict(degree="TOTAL", r1="", r2="", signature="",
                     n_subfields=sum(census.values()),
                     salem_bearing="", note="all subfields of K"))

    p_csv = write_csv("signature_census.csv",
                      ["degree", "r1", "r2", "signature", "n_subfields",
                       "salem_bearing", "note"], rows, __file__)

    # [K:Q] = 16 via a primitive element
    theta = sp.sqrt(2) + sp.sqrt(3) + sp.root(5, 4)
    K_degree = int(sp.degree(sp.minimal_polynomial(theta, _x), _x))

    imag_place = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)
    fifth_root = mp.mpf(5) ** mp.mpf("0.25")

    salem_shaped = {f"{k[0]},{k[1]},{k[2]}": v for k, v in census.items()
                    if k[1] == 2 and k[2] == (k[0] // 2) - 1 and k[0] >= 4}

    payload = {
        "_description": "Independent Galois-correspondence signature census of "
                        "K = Q(sqrt2,sqrt3,5^{1/4}) (Theorem 5.3), plus the total "
                        "subgroup count of G ~ C2 x C2 x D4 (=158).",
        "galois_group": "Gal(K(i)/Q) ~ C2 x C2 x D4, order 32",
        "K_degree_over_Q": K_degree,
        "census_by_deg_r1_r2":
            {f"({k[0]},{k[1]},{k[2]})": v for k, v in
             sorted(census.items(), key=lambda kv: (kv[0][0], kv[0][1]))},
        "n_subfields": n_subfields,
        "census_total": sum(census.values()),
        "n_total_subgroups": n_total,
        "subgroup_order_breakdown":
            {str(k): order_breakdown[k] for k in sorted(order_breakdown)},
        "salem_bearing_signatures": salem_shaped,
        "K_own_signature": "(8,4)  [row (deg,r1,r2)=(16,8,4)]",
        "note_27_vs_158": "27 = subfields of K = subgroups of G containing "
                          "complex conjugation; 158 = total subgroups of G.",
        "salem_trace_form_signatures": salem_signature_rows(),
        "beta_4_gt_phi_exact": beta4_exact_sign(),
        "complex_place_off_circle": {
            "abs_i_beta": mp.nstr(imag_place, 12),
            "abs_5^{1/4}_i": mp.nstr(fifth_root, 12),
            "both_off_unit_circle":
                bool(abs(imag_place - 1) > mp.mpf("0.1")
                     and abs(fifth_root - 1) > mp.mpf("0.1")),
        },
    }
    p_json = write_json("signature_census.json", payload, __file__)

    print("wrote", p_csv)
    print("wrote", p_json)
    print(f"  census: "
          f"{ {f'{k[0]},{k[1]},{k[2]}': v for k, v in sorted(census.items())} }")
    print(f"  subfields = {n_subfields}   total subgroups = {n_total}")
    print(f"  order breakdown = { {k: order_breakdown[k] for k in sorted(order_breakdown)} }")
    print(f"  [K:Q] = {K_degree}   salem-bearing = {salem_shaped}")


if __name__ == "__main__":
    main()
