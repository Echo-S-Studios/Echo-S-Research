"""Producer: the sl2 representation theory (dimensions are the exponents).

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 4):
  * Thm 4.1 / Prop 4.3  irreducibles V_m: dim = m+1, weights {m,..,-m},
    Casimir (1/2)m(m+2)
  * Prop 4.3            the symmetric-power eigenvalue tower
    spec(R|V_n) = {phi^{n-k} psi^k}, and H/sqrt5 integer weights n-2k
  * Thm 4.4 / Ex 4.5    Clebsch-Gordan: V_a x V_b = (+) V_j, dim (a+1)(b+1)=sum(j+1)

Outputs:
  data/.../irreps.csv          -- V_0 .. V_6
  data/.../sympow_tower.csv    -- eigenvalue tower + H-weights, n = 0..5
  data/.../clebsch_gordan.csv  -- decompositions for 0<=a,b<=5

Unlike tests/test_reps.py (asserts the formulae), this script evaluates them
and emits the three tables; each CG row carries the checked dimension identity.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (phi, psi, sqrt5, casimir, dim_V, weights_V, sympow_tower,
                      cg_decomp, is_zero)
from eap_io import write_csv


def irrep_rows():
    rows = []
    for m in range(0, 7):
        rows.append({
            "m": m,
            "dim_V_m": dim_V(m),
            "weights": " ".join(str(w) for w in weights_V(m)),
            "casimir_half_m_m+2": sp.sstr(casimir(m)),
        })
    return rows


def tower_rows():
    rows = []
    for n in range(0, 6):
        tower = sympow_tower(n)
        for k in range(n + 1):
            ev = sp.nsimplify(tower[k])
            hweight = n - 2 * k                       # eigenvalue of H/sqrt5
            hval = (n - k) * sqrt5 + k * (-sqrt5)     # eigenvalue of H itself
            rows.append({
                "n": n, "k": k,
                "eigenvalue_phi^(n-k)psi^k": sp.sstr(ev),
                "H_eigenvalue": sp.sstr(sp.nsimplify(hval)),
                "H_over_sqrt5_weight": hweight,
                "weight_check": bool(is_zero(hval / sqrt5 - hweight)),
            })
    return rows


def cg_rows():
    rows = []
    for a in range(0, 6):
        for b in range(0, 6):
            js = cg_decomp(a, b)
            lhs = (a + 1) * (b + 1)
            rhs = sum(j + 1 for j in js)
            rows.append({
                "a": a, "b": b,
                "summands_V_j": "+".join(f"V{j}" for j in js),
                "dim_(a+1)(b+1)": lhs,
                "sum_(j+1)": rhs,
                "dim_identity_holds": bool(lhs == rhs),
            })
    return rows


def main():
    r1, r2, r3 = irrep_rows(), tower_rows(), cg_rows()
    p1 = write_csv("irreps.csv", list(r1[0].keys()), r1, __file__)
    p2 = write_csv("sympow_tower.csv", list(r2[0].keys()), r2, __file__)
    p3 = write_csv("clebsch_gordan.csv", list(r3[0].keys()), r3, __file__)
    print(f"wrote {p1}  (V_0..V_6)")
    print(f"wrote {p2}  ({len(r2)} tower entries)")
    print(f"wrote {p3}  ({len(r3)} CG pairs; all dim identities: "
          f"{all(r['dim_identity_holds'] for r in r3)})")


if __name__ == "__main__":
    main()
