"""Producer: psi^n axis-dilation, the coupled Casimir, and the deviation ladder.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 5.4 - 5.5):
  * Thm 5.7   psi^n dilates the Cartan axis: gap sqrt5 -> F_n sqrt5,
              ad_{R^n} = F_n ad_R, index multiplicative
  * Prop 5.10 Casimir on coupled targets: V1xV2=V1+V3 (3/2, 15/2),
              V2xV2=V0+V2+V4 (0,4,12); tops Cas(V3)=15/2, Cas(V4)=12
  * Thm 5.11  the deviation ladder
        (1/2)Tr_{V_m}(X_n^2) = 5 F_n^2 C(m+2,3) = (5 F_n^2/3) dim(V_m) Cas(V_m)
        with V_3 -> 50 F_n^2, V_4 -> 100 F_n^2

Outputs:
  data/.../deviation_ladder.csv    -- (n,m) grid, both closed forms
  data/.../psi_dilation.json       -- gap dilation, ad_{R^n}=F_n ad_R, coupled Casimir

Unlike tests/test_traceform.py (asserts the two ladder forms agree), this script
evaluates them on the (n,m) grid and emits the table + the dilation record.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, sqrt5, phi, psi, fib, casimir, dim_V, cg_decomp,
                      half_trace_Vm, comm, is_zero, mat_eq, H, S, J)
from eap_io import write_csv, write_json


def ladder_rows():
    rows = []
    for n in range(-4, 5):
        for m in range(0, 7):
            lhs = half_trace_Vm(n, m)
            binom_form = 5 * fib(n) ** 2 * sp.binomial(m + 2, 3)
            casimir_form = sp.Rational(5, 3) * fib(n) ** 2 * dim_V(m) * casimir(m)
            rows.append({
                "n": n, "m": m,
                "half_trace_Vm_X_n^2": int(lhs),
                "5_F_n^2_binom(m+2,3)": int(binom_form),
                "(5F_n^2/3)dim_Cas": int(casimir_form),
                "binom(m+2,3)": int(sp.binomial(m + 2, 3)),
                "forms_agree": bool(is_zero(lhs - binom_form) and is_zero(lhs - casimir_form)),
            })
    return rows


def dilation_payload():
    # (a) spectral gap dilation
    gap_ok = all(is_zero((phi ** n - psi ** n) - fib(n) * sqrt5) for n in range(-8, 9))
    # (b) ad_{R^n} = F_n ad_R on test matrices
    testmats = [sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[1, 2], [3, 4]]), H, S, J]
    adRn_ok = all(mat_eq(comm(R ** n, Y), fib(n) * comm(R, Y))
                  for n in range(-6, 7) for Y in testmats)
    return {
        "a_spectrum": {
            "gap_dilation": "phi - psi = sqrt5  ->  phi^n - psi^n = F_n sqrt5",
            "verified_n_-8_8": bool(gap_ok),
        },
        "b_operator": {
            "ad_Rn_eq_Fn_adR": bool(adRn_ok),
            "hence": "H = X_1  ->  X_n = F_n H; index multiplicative psi^n: X_m -> X_{nm}",
        },
        "c_invariant": {
            "seed": "(1/2)Tr(H^2) = 5 = L_1^2 - 4(-1)^1",
            "lift": "(1/2)Tr(X_n^2) = 5 F_n^2 = L_n^2 - 4(-1)^n",
        },
        "coupled_casimir": {
            "V1xV2": {"summands": [f"V{j}" for j in cg_decomp(1, 2)],
                      "casimirs": [sp.sstr(casimir(j)) for j in cg_decomp(1, 2)]},
            "V2xV2": {"summands": [f"V{j}" for j in cg_decomp(2, 2)],
                      "casimirs": [sp.sstr(casimir(j)) for j in cg_decomp(2, 2)]},
            "Cas_V3": sp.sstr(casimir(3)),      # 15/2
            "Cas_V4": sp.sstr(casimir(4)),      # 12
        },
        "ladder_scalings": {
            "V3_factor_binom(5,3)": int(sp.binomial(5, 3)),   # 10  -> 50 F_n^2
            "V4_factor_binom(6,3)": int(sp.binomial(6, 3)),   # 20  -> 100 F_n^2
        },
    }


def main():
    rows = ladder_rows()
    p1 = write_csv("deviation_ladder.csv", list(rows[0].keys()), rows, __file__)
    p2 = write_json("psi_dilation.json", dilation_payload(), __file__)
    allok = all(r["forms_agree"] for r in rows)
    print(f"wrote {p1}  ({len(rows)} (n,m) cells, both ladder forms agree: {allok})")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
