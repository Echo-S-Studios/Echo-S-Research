"""Producer: the self-action spectrum, root spaces, and the semiring grading.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces:
  * Prop 2.2 / Thm 2.5  spec(ad_R) = {0,0,+sqrt5,-sqrt5}; root vectors N_+, N_-;
    the 0-eigenspace = span{I,R}
  * Prop 3.6            the grading laws: how (M, chi) transform under the three
    operators superposition (+), coupling (x), Adams scaling psi^n
  * Rem 3.7 / Ex 3.8    the tropical coupling subtlety (M(AxB)=phi^4, not phi^3)
  * Prop 3.10           psi^n pulls back to the core (Binet; psi^2 is the keystone)

Outputs:
  data/.../adR_selfaction.json
  data/.../semiring_grading.csv

Unlike tests/test_lie_layer.py + test_semiring.py (asserts), this script
assembles the 4x4 ad_R map, diagonalises it, and emits the spectrum, the root
data, and a machine-readable grading table with worked golden-object values.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, I2, phi, psi, sqrt5, fib, comm, is_zero, mat_eq,
                      mat_to_rows, adR_matrix_on_M2, H, Np, Nm,
                      mahler_exact, charge_set, outside_unit)
from eap_io import write_json, write_csv


def selfaction_payload():
    M = adR_matrix_on_M2()
    eigs = {sp.sstr(sp.simplify(val)): int(mult) for val, mult in M.eigenvals().items()}

    # 0-eigenspace: solve RX = XR, confirm 2-dimensional = span{I, R}
    p, q, r, s = sp.symbols("p q r s")
    Xg = sp.Matrix([[p, q], [r, s]])
    sol = sp.linsolve(list(comm(R, Xg)), [p, q, r, s])
    (solset,) = sol
    free = sorted(sp.Matrix(solset).free_symbols, key=lambda t: t.name)

    return {
        "adR_spectrum_on_M2": eigs,
        "adR_spectrum_note": "{0: 2, sqrt5: 1, -sqrt5: 1}; the irrational "
                             "eigenvalue sqrt5 = phi - psi = phi + 1/phi is the coupling",
        "root_vector_N_plus": mat_to_rows(Np),
        "root_vector_N_minus": mat_to_rows(Nm),
        "eigen_eq_N_plus": "[R, N_+] = +sqrt5 N_+  : " + str(mat_eq(comm(R, Np), sqrt5 * Np)),
        "eigen_eq_N_minus": "[R, N_-] = -sqrt5 N_- : " + str(mat_eq(comm(R, Nm), -sqrt5 * Nm)),
        "zero_space_dimension": len(free),
        "zero_space_is_span_I_R": bool(len(free) == 2),
        "decomposition": "M_2 = span{I,R} (0-space, dim 2) (+) C N_+ (+sqrt5) (+) C N_- (-sqrt5)",
    }


def grading_rows():
    """Prop 3.6 grading table + worked golden-object instances."""
    rows = []
    # law rows (symbolic statements of the transform)
    rows.append({"operator": "superposition (+)", "magnitude_law": "M(A)*M(B)",
                 "charge_law": "chi(A) U chi(B)",
                 "worked_example": "M({phi^2,psi^2}+{phi,psi}) = phi^2*phi = phi^3"})
    rows.append({"operator": "coupling (x)", "magnitude_law": "prod max(1,|a_i b_j|) [tropical]",
                 "charge_law": "chi(A)+chi(B) mod 4",
                 "worked_example": "M({phi^2,psi^2} x {phi,psi}) = phi^4  (NOT naive phi^3)"})
    rows.append({"operator": "Adams psi^n", "magnitude_law": "M(A)^n",
                 "charge_law": "n*chi(A) mod 4",
                 "worked_example": "M(psi^3{phi,psi}) = phi^3; chi = 3*{0,2} = {0,2}"})
    return rows


def semiring_facts():
    """Concrete numbers that the grading table refers to (for provenance JSON)."""
    A2, A1 = [phi ** 2, psi ** 2], [phi, psi]
    # tropical coupling
    coupled = [a * b for a in A2 for b in A1]
    out = outside_unit(coupled)
    trop = sp.nsimplify(sp.prod([sp.Abs(v) for v in out]))
    naive = phi ** 2 * phi
    return {
        "golden_object_A_phi": {
            "magnitude_M": sp.sstr(mahler_exact([phi, psi])),   # phi
            "charge_chi": charge_set([phi, psi]),               # {0,2}
        },
        "superposition_factors": {
            "M_union": sp.sstr(mahler_exact(A2 + A1)),           # phi^3
            "equals_product_M_A_M_B": bool(is_zero(mahler_exact(A2 + A1) - phi ** 3)),
        },
        "coupling_tropical": {
            "coupled_spectrum": ["phi^3", "-phi", "phi^-1", "-phi^-3"],
            "num_outside_unit_circle": len(out),
            "M_tropical": sp.sstr(trop),                        # phi^4
            "M_naive_product": sp.sstr(sp.nsimplify(naive)),    # phi^3
            "tropical_over_naive": sp.sstr(sp.nsimplify(trop / naive)),  # phi
            "note": "the cross term psi^2*phi = phi^-1 falls THROUGH the unit "
                    "circle and is floored to 1; naive M(A)M(B) is wrong by phi",
        },
        "adams_pullback": {
            "psi2_is_keystone": bool(is_zero(phi ** 2 - (phi + 1))),
            "binet_n_0_8": bool(all(is_zero(phi ** n - (fib(n) * phi + fib(n - 1)))
                                    for n in range(0, 9))),
            "multiplicative_not_additive": "psi^2(phi+1) - [psi^2(phi)+psi^2(1)] = 2phi",
            "non_additivity_gap": sp.sstr(sp.expand((phi + 1) ** 2 - (phi ** 2 + 1))),  # 2*phi
        },
    }


def main():
    p1 = write_json("adR_selfaction.json", selfaction_payload(), __file__)
    rows = grading_rows()
    p2 = write_csv("semiring_grading.csv",
                   ["operator", "magnitude_law", "charge_law", "worked_example"],
                   rows, __file__)
    p3 = write_json("semiring_grading_facts.json", semiring_facts(), __file__)
    print(f"wrote {p1}")
    print(f"wrote {p2}  ({len(rows)} operators)")
    print(f"wrote {p3}")


if __name__ == "__main__":
    main()
