"""Producer: the sl2 bracket table, the (2,1) metric, and the null-frame transition.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 2.3 - 2.6):
  * Thm 2.6 / Prop 2.9  the sl2 bracket tables in both bases:
        rational  {H,S,J}:  [H,S]=10J, [H,J]=2S, [S,J]=-2H
        split     {H,N+,N-}: [H,N+]=2sqrt5 N+, [H,N-]=-2sqrt5 N-, [N+,N-]=sqrt5 H
  * Prop 2.7 / 5.2      Lorentzian signature (2,1): H^2=S^2=5I, J^2=-I, N+^2=N-^2=0,
        and the trace-form Gram matrix diag(10,10,-2)
  * Thm 2.10            null frame N+/- = (S +/- sqrt5 J)/2; transition matrices
        M_{R->E}, M_{E->R} (det = -2/sqrt5); conjugator V diagonalising H

Outputs:
  data/.../sl2_bracket_table.csv   -- every bracket in both bases + the metric
  data/.../transition.json         -- rational constants, transition matrices, V

Unlike tests/test_lie_layer.py (asserts each bracket), this script computes the
full antisymmetric bracket table and the Gram matrix and emits them.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, I2, phi, psi, sqrt5, comm, is_zero, mat_eq, mat_to_rows,
                      trace_form, H, S, J, Np, Nm, V)
from eap_io import write_csv, write_json


def bracket_rows():
    """Full [X,Y] table over both bases, as canonical linear combinations."""
    rows = []

    def linear_form(M, basis, names):
        """Express traceless M as an integer/algebraic combo of basis matrices."""
        coeffs = sp.symbols(f"c0:{len(basis)}")
        expr = sum((c * B for c, B in zip(coeffs, basis)), sp.zeros(2, 2))
        sol = sp.solve([sp.Eq(sp.expand(expr[i, j]), M[i, j])
                        for i in range(2) for j in range(2)], coeffs, dict=True)
        if not sol:
            return sp.sstr(sp.nsimplify(M))
        sol = sol[0]
        terms = []
        for c, nm in zip(coeffs, names):
            v = sp.nsimplify(sol.get(c, 0))
            if v != 0:
                terms.append(f"({sp.sstr(v)})*{nm}")
        return " + ".join(terms) if terms else "0"

    rat = [H, S, J]
    rat_names = ["H", "S", "J"]
    for i in range(3):
        for j in range(3):
            br = comm(rat[i], rat[j])
            rows.append({"basis": "rational {H,S,J}",
                         "X": rat_names[i], "Y": rat_names[j],
                         "[X,Y]": linear_form(br, rat, rat_names)})

    spl = [H, Np, Nm]
    spl_names = ["H", "N+", "N-"]
    for i in range(3):
        for j in range(3):
            br = comm(spl[i], spl[j])
            rows.append({"basis": "split {H,N+,N-}",
                         "X": spl_names[i], "Y": spl_names[j],
                         "[X,Y]": linear_form(br, spl, spl_names)})
    return rows


def metric_rows():
    """Trace-form Gram matrix on {H,S,J}: diag(10,10,-2), signature (2,1)."""
    names = ["H", "S", "J"]
    basis = [H, S, J]
    rows = []
    for i in range(3):
        for j in range(3):
            rows.append({"basis": "trace-form <X,Y>=Tr(XY)",
                         "X": names[i], "Y": names[j],
                         "[X,Y]": sp.sstr(sp.nsimplify(trace_form(basis[i], basis[j])))})
    return rows


def transition_payload():
    MRE = sp.Matrix([[1, 0, 0], [0, 1, 1 / sqrt5], [0, 1, -1 / sqrt5]])
    MER = sp.Matrix([[1, 0, 0],
                     [0, sp.Rational(1, 2), sp.Rational(1, 2)],
                     [0, sqrt5 / 2, -sqrt5 / 2]])
    Vi = V.inv()
    return {
        "rational_structure_constants": {
            "[H,S]": "10 J", "[H,J]": "2 S", "[S,J]": "-2 H",
            "HS": mat_to_rows(H * S), "SH": mat_to_rows(S * H),
        },
        "splitting_obstruction": {
            "adH_on_span_S_J": [[0, 2], [10, 0]],
            "char_poly": "lambda^2 - 20",
            "roots": ["2*sqrt5", "-2*sqrt5"],
            "note": "not diagonalisable over Q; minimal field is Q(sqrt20)=Q(sqrt5)",
        },
        "signature_2_1": {
            "H^2": mat_to_rows(H * H), "S^2": mat_to_rows(S * S),
            "J^2": mat_to_rows(J * J),
            "N+^2_is_zero": bool(mat_eq(Np * Np, sp.zeros(2))),
            "N-^2_is_zero": bool(mat_eq(Nm * Nm, sp.zeros(2))),
            "trace_form_gram": "diag(10, 10, -2)",
            "signature": "(2,1): two space-like boosts H,S; one time-like rotation J",
        },
        "null_frame": {
            "N+": "(S + sqrt5 J)/2", "N-": "(S - sqrt5 J)/2",
            "N+_check": bool(mat_eq(Np, (S + sqrt5 * J) / 2)),
            "N-_check": bool(mat_eq(Nm, (S - sqrt5 * J) / 2)),
            "M_R_to_E": mat_to_rows(MRE),
            "M_E_to_R": mat_to_rows(MER),
            "mutually_inverse": bool(mat_eq(MRE * MER, sp.eye(3))),
            "det_M_R_to_E": sp.sstr(sp.nsimplify(MRE.det())),  # -2/sqrt5
            "det_is_irrational": "det = -2/sqrt5 certifies the frame change leaves Q",
        },
        "conjugator_V": {
            "V": mat_to_rows(V),
            "V^-1 H V": mat_to_rows(Vi * H * V),          # diag(sqrt5,-sqrt5)
            "V^-1 N+ V": mat_to_rows(Vi * Np * V),        # [[0,(5-sqrt5)/2],[0,0]]
            "V^-1 N- V": mat_to_rows(Vi * Nm * V),        # [[0,0],[(5+sqrt5)/2,0]]
        },
    }


def main():
    rows = bracket_rows() + metric_rows()
    p1 = write_csv("sl2_bracket_table.csv", ["basis", "X", "Y", "[X,Y]"], rows, __file__)
    p2 = write_json("transition.json", transition_payload(), __file__)
    print(f"wrote {p1}  ({len(rows)} bracket/metric rows)")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
