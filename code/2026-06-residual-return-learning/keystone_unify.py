r"""
Producer: Section 6 "The unifying principle and the phi keystone", Proposition 6.2
of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits the phi-keystone identity: the same object seen three independent ways --
  (1) the exact integer companion of x^2-x-1 (the substrate's first field Q(phi)),
  (2) the symmetric matrix R of the language return operator (R^2 = R + I), and
  (3) the Clifford element cl(1/2, 1, -1/2, 0),
all the void law x^2 = x + 1, with Mahler measure = phi.

  -> data/keystone_identity.json

Also emits the residual-return correspondence table (the five-part structure that
is byte-identical across the learning and language faces, Table 5).

  -> data/correspondence_table.csv

Run:  py code/2026-06-residual-return-learning/keystone_unify.py
"""
import sympy as sp
from sympy import Matrix, Rational as Q, eye

import rrl_core as core

SCRIPT = "keystone_unify.py"
x = core.X


def produce_keystone():
    """Prop 6.2: the three routes coincide as one 2x2 matrix / one minimal poly."""
    route1 = core.companion([1, -1, -1])                 # C(x^2-x-1)
    route2 = core.KEYSTONE_R                             # mat(1/2,1,-1/2,0), R^2=R+I
    route3 = core.mat(Q(1, 2), 1, Q(-1, 2), 0)          # cl(1/2,1,-1/2,0)
    mah = core.mahler_measure([1, -1, -1])
    phi_mp = core.phi_mpf()
    phi = (1 + sp.sqrt(5)) / 2
    payload = {
        "proposition": "the phi-keystone is one object (Prop 6.2)",
        "route1_companion_of_golden_law": core.mat_to_strlist(route1),   # [[0,1],[1,1]]
        "route2_language_keystone_R": core.mat_to_strlist(route2),
        "route3_clifford_cl_half_1_mhalf_0": core.mat_to_strlist(route3),
        "all_three_equal": bool(route1 == route2 == route3),
        "minimal_polynomial": core.poly_str(sp.expand(route1.charpoly(x).as_expr())),
        "minimal_polynomial_coeffs": core.poly_coeffs(route1.charpoly(x).as_expr()),
        "golden_law_R2_eq_R_plus_I": bool(route2 * route2 == route2 + eye(2)),
        "mahler_measure": str(mah),
        "mahler_equals_phi": bool(abs(mah - phi_mp) < core.mp.mpf("1e-30")),
        "phi_value": str(sp.nsimplify(phi)),
        "TYPE_R_gen": bool(route2 * route2 == route2 + eye(2)),
        "det_R": int(route2.det()),                      # -1
        "scalar_part_tau_R": str(core.cl_coords(route2)[0]),  # 1/2
    }
    core.write_json(SCRIPT, "keystone_identity.json", payload)
    return payload


def produce_correspondence():
    """Table 5: the residual-return correspondence, one mechanism in two carriers."""
    rows = [
        ["loss/residual", "r = x - P x  (trace form G)", "L(X) = R X + X R - X"],
        ["captured iff", "||r||_G^2 = 0 (exact over Q)", "L(X) = 0 (exact over Q)"],
        ["landing space", "col(B) subset K", "ker L = span{e1+2e2, i}"],
        ["model state", "forced basis B (growing)", "lexicon subset ker L (growing)"],
        ["sole mutator", "confirm (adjoin seed)", "commit / lexicon add"],
        ["generalisation", "same residual direction captured once",
         "same exact residue => one entry"],
        ["witness", "SHA-256 chain (eq chain)", "the SAME SHA-256 chain"],
        ["exactness gate", "Q/Z decision; float refused", "Q decision; float refused"],
    ]
    core.write_csv(SCRIPT, "correspondence_table.csv",
                   ["role", "learning_number_field_K", "language_Cl_2_0"], rows)
    return rows


def main():
    ks = produce_keystone()
    produce_correspondence()
    print(f"[{SCRIPT}] three routes equal: {ks['all_three_equal']}, "
          f"minpoly {ks['minimal_polynomial']}")
    print(f"[{SCRIPT}] Mahler measure = {ks['mahler_measure']} "
          f"(= phi: {ks['mahler_equals_phi']})")
    print(f"[{SCRIPT}] wrote 2 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
