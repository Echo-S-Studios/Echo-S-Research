"""Producer: the deviation operator and the Trace-Form Duality.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 5.1 - 5.3):
  * Prop 5.4  the deviation operator collapses: X_n = 2R^n - L_n I = F_n H
  * Thm 5.6   the Trace-Form Duality
        (1/2) Tr(X_n^2) = 5 F_n^2 = L_n^2 - 4(-1)^n = (phi^n - psi^n)^2
    verified across n in [-8, 8]
  * Prop 5.2  the Lorentzian metric diag(10,10,-2) (the form the duality uses)

Outputs:
  data/.../traceform_duality.csv   -- one row per n, all four equal columns
  data/.../traceform_metric.json   -- the (2,1) trace-form Gram matrix

Unlike tests/test_traceform.py (asserts the four are equal), this script builds
X_n from the actual matrix power, evaluates all four expressions independently,
and emits them side by side with an equality flag.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, I2, phi, psi, fib, luc, X, H, S, J, trace_form,
                      is_zero, mat_eq)
from eap_io import write_csv, write_json

N_LO, N_HI = -8, 8


def duality_rows():
    rows = []
    for n in range(N_LO, N_HI + 1):
        Xn = X(n)                                        # 2R^n - L_n I
        collapses = mat_eq(Xn, fib(n) * H)               # = F_n H
        half_tr = sp.Rational(1, 2) * sp.trace(Xn * Xn)  # (1/2)Tr(X_n^2)
        five_f2 = 5 * fib(n) ** 2
        lucas_form = luc(n) ** 2 - 4 * (-1) ** (n % 2)
        binet_gap = sp.expand((phi ** n - psi ** n) ** 2)
        all_equal = (is_zero(half_tr - five_f2)
                     and is_zero(five_f2 - lucas_form)
                     and is_zero(five_f2 - binet_gap))
        rows.append({
            "n": n,
            "F_n": int(fib(n)),
            "L_n": int(luc(n)),
            "X_n_eq_F_n_H": bool(collapses),
            "half_trace_X_n^2": int(half_tr),
            "5_F_n^2": int(five_f2),
            "L_n^2-4(-1)^n": int(lucas_form),
            "(phi^n-psi^n)^2": int(sp.nsimplify(binet_gap)),
            "all_four_equal": bool(all_equal),
        })
    return rows


def metric_payload():
    names = ["H", "S", "J"]
    basis = [H, S, J]
    gram = [[int(trace_form(basis[i], basis[j])) for j in range(3)] for i in range(3)]
    invariance = all(
        is_zero(sp.trace((basis[a] * basis[b] - basis[b] * basis[a]) * basis[c])
                + sp.trace(basis[b] * (basis[a] * basis[c] - basis[c] * basis[a])))
        for a in range(3) for b in range(3) for c in range(3)
    )
    return {
        "basis": names,
        "gram_matrix_Tr(XY)": gram,
        "diagonal": [gram[i][i] for i in range(3)],
        "signature": "(2,1)  -- diag(10,10,-2)",
        "space_like": "H, S  (<X,X> = 10 > 0)",
        "time_like": "J  (<J,J> = -2 < 0)",
        "invariance_verified": bool(invariance),
        "seed_length": "(1/2)Tr(H^2) = 5 = L_1^2 - 4(-1)^1",
    }


def main():
    rows = duality_rows()
    p1 = write_csv("traceform_duality.csv", list(rows[0].keys()), rows, __file__)
    p2 = write_json("traceform_metric.json", metric_payload(), __file__)
    allok = all(r["all_four_equal"] for r in rows)
    print(f"wrote {p1}  ({len(rows)} rows, four-way duality holds: {allok})")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
