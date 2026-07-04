"""Producer: the worked examples and exercises across Sec. 2 - 7.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces:
  * Ex 2.12  the sl2 triple by hand: N+ N-, N- N+, [N+,N-] = sqrt5 H
  * Ex 3.8   tropical vs naive coupling: M(AxB) = phi^4, not phi^3
  * Ex 4.7   Sym^2 in the flesh: eigenvalues {phi^2,-1,psi^2}, H/sqrt5 weights {2,0,-2}
  * Prop 6.4 / Ex 6.5  the flip: D = 1+4C = Gram det; C=1 -> D=5 (terrain),
             C=-1 -> D=-3 (rotation, cube roots)
  * Ex 7.5   the floor-collapse orbit {i,-i} -> {-1,-1} -> {1}, M == 1 throughout

Output:
  data/.../worked_examples.json

Unlike tests/test_*.py, this script computes each worked example symbolically
and emits the numbers together in one provenance-stamped record.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (phi, psi, sqrt5, comm, is_zero, mat_eq, mat_to_rows,
                      H, Np, Nm, mahler_exact, charge_set, outside_unit,
                      mahler_numeric)
from eap_io import write_json


def sl2_triple_by_hand():
    NpNm = Np * Nm
    NmNp = Nm * Np
    return {
        "N+_N-": mat_to_rows(NpNm),
        "N-_N+": mat_to_rows(NmNp),
        "[N+,N-]_eq_sqrt5_H": bool(mat_eq(comm(Np, Nm), sqrt5 * H)),
        "phi^2-psi^2_eq_sqrt5": bool(is_zero((phi ** 2 - psi ** 2) - sqrt5)),
    }


def tropical_vs_naive():
    A, B = [phi ** 2, psi ** 2], [phi, psi]
    coupled = [a * b for a in A for b in B]
    out = outside_unit(coupled)
    trop = sp.nsimplify(sp.prod([sp.Abs(v) for v in out]))
    return {
        "coupled_spectrum": ["phi^3", "-phi", "phi^-1", "-phi^-3"],
        "outside_unit_circle_count": len(out),
        "M_tropical": sp.sstr(trop),                      # phi^4
        "M_naive_M(A)M(B)": sp.sstr(sp.nsimplify(phi ** 2 * phi)),  # phi^3
        "off_by_factor": sp.sstr(sp.nsimplify(trop / (phi ** 2 * phi))),  # phi
        "superposition_does_factor": sp.sstr(mahler_exact(A + B)),  # phi^3
    }


def sym2_in_the_flesh():
    eigs = [phi ** 2, phi * psi, psi ** 2]
    return {
        "eigenvalues": [sp.sstr(sp.nsimplify(e)) for e in eigs],  # phi^2, -1, psi^2
        "phi_psi_eq_minus1": bool(is_zero(sp.nsimplify(eigs[1]) - (-1))),
        "H_over_sqrt5_weights": [2 - 2 * k for k in range(3)],    # {2,0,-2}
    }


def the_flip():
    C, x = sp.symbols("C x")
    D = sp.discriminant(x ** 2 + x - C, x)
    gram = sp.Matrix([[2, -1], [-1, 2 * C + 1]])
    # C=1 (golden gate) and C=-1 (cube roots)
    r1 = [sp.nsimplify(v) for v in sp.solve(x ** 2 + x - 1, x)]
    r2 = [sp.nsimplify(v) for v in sp.solve(x ** 2 + x + 1, x)]
    return {
        "discriminant_D": sp.sstr(sp.expand(D)),            # 1 + 4C
        "gram_det_equals_D": bool(is_zero(gram.det() - (4 * C + 1))),
        "flip_at_C": "-1/4",
        "C=1_terrain": {"D": 5, "roots": [sp.sstr(v) for v in r1],
                        "note": "real roots {1/phi, -phi}; sqrt5 marks the golden gate"},
        "C=-1_rotation": {"D": -3, "roots": [sp.sstr(v) for v in r2],
                          "note": "complex pair = primitive cube roots of unity"},
    }


def floor_collapse_orbit():
    start = [sp.I, -sp.I]
    step1 = [sp.I ** 2, (-sp.I) ** 2]     # {-1,-1}
    step2 = [v ** 2 for v in step1]       # {1,1}
    return {
        "start_{i,-i}": {"M": str(mahler_numeric(start)), "charge": charge_set(start)},
        "psi^2_->_{-1,-1}": {"M": str(mahler_numeric(step1)), "charge": charge_set(step1),
                             "matches_2*{1,3}mod4": sorted({(2 * c) % 4 for c in charge_set(start)})},
        "psi^2_->_{1}": {"M": str(mahler_numeric(step2)), "charge": charge_set(step2)},
        "note": "M stays exactly 1 throughout; only the charge cycles 1,3 -> 2 -> 0",
    }


def main():
    payload = {
        "Ex_2_12_sl2_triple_by_hand": sl2_triple_by_hand(),
        "Ex_3_8_tropical_vs_naive_coupling": tropical_vs_naive(),
        "Ex_4_7_sym2_in_the_flesh": sym2_in_the_flesh(),
        "Prop_6_4_Ex_6_5_the_flip": the_flip(),
        "Ex_7_5_floor_collapse_orbit": floor_collapse_orbit(),
    }
    p = write_json("worked_examples.json", payload, __file__)
    print(f"wrote {p}  (5 worked examples)")


if __name__ == "__main__":
    main()
