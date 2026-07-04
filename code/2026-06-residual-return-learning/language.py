r"""
Producer: Section 5 "A residual-valued language" over Cl(2,0) = M_2(R) of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits:
  * Sec 5.1-5.2 (the exact carrier + keystone): the matrix isomorphism relations,
      the keystone R (golden law R^2=R+I), Cayley-Hamilton, and the gate P0
      -> data/carrier_readings.json
  * Prop 5.8 + Thm 5.4 (the return operator and its kernel): the 4x4 matrix of
      L(X)=RX+XR-X, L(R) (R not in kernel), ker L = span{e1+2e2, i} (dim 2),
      and the Sylvester eigenvalues {sqrt5,0,0,-sqrt5}
      -> data/return_operator_kernel.json
  * Rem 5.5 (the disproof): L(H(0,-1,1,0)) = H(-3,0,0,0) (NOT in kernel), while
      cl([[0,-1],[1,0]]) = i IS in the kernel
      -> included in return_operator_kernel.json
  * eq (proj) / Thm 5.6 (the commit projector): the exact idempotent orthogonal
      projector onto ker L
      -> data/commit_projector.csv
  * Thm 5.6 / Ex 5.7 (generalisation + growing lexicon): five tokens -> four exact
      residues (E1 and E1+1 merge; i, 2i, -i distinct); sign-aware words
      -> data/lexicon.csv
  * Table 6 (the jurisdiction firewall counts): 20 + 5 = 25 wired, + 2 = 27 bank
      -> data/firewall.csv

Run:  py code/2026-06-residual-return-learning/language.py

All computation is through the exact Clifford matrix isomorphism (rrl_core.mat /
cl_coords) with sympy rationals; the producer emits values, not assertions.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, eye, zeros

import rrl_core as core

SCRIPT = "language.py"
x = core.X
R = core.KEYSTONE_R

# kernel basis columns  e1+2e2 = (0,1,2,0),  i = (0,0,0,1)  in coords (1,e1,e2,i)
KERNEL_K = Matrix([[0, 0], [1, 0], [2, 0], [0, 1]])
COMMIT = KERNEL_K * (KERNEL_K.T * KERNEL_K).inv() * KERNEL_K.T
KERNEL_LABELS = ["e1+2e2", "i"]


def word(coord4):
    """Sign-aware sparse code of commit(X) over the kernel basis, exact L1-relative
    threshold tau=1/4 (Thm 5.6 / eq proj)."""
    committed = COMMIT * Matrix(coord4)
    c = (KERNEL_K.T * KERNEL_K).inv() * KERNEL_K.T * committed   # kernel-basis coords
    l1 = sum(abs(v) for v in c)
    tokens = []
    for j, cj in enumerate(c):
        if l1 != 0 and abs(cj) >= Q(1, 4) * l1:
            tokens.append(("+" if cj > 0 else "-") + KERNEL_LABELS[j])
    return "(" + "".join(tokens) + ")"


def produce_carrier():
    """Sec 5.1-5.2: the exact carrier readings and the keystone golden law."""
    e1, e2, i, one = core.mat(0, 1, 0, 0), core.mat(0, 0, 1, 0), core.mat(0, 0, 0, 1), eye(2)
    P0 = core.mat(Q(1, 2), 0, Q(1, 2), 0)          # gate P0 = 1/2(1+e2)
    generic = core.mat(Q(3, 7), -2, Q(5, 3), 4)
    ch = lambda M: bool(sp.expand(M * M - M.trace() * M + M.det() * eye(2)) == zeros(2))
    payload = {
        "matrix_iso": "mat(a+b e1+c e2+d i) = [[a+c, b-d],[b+d, a-c]]",
        "generator_relations": {
            "e1^2=1": bool(e1 * e1 == one), "e2^2=1": bool(e2 * e2 == one),
            "i^2=-1": bool(i * i == -one), "i=e1e2": bool(e1 * e2 == i),
        },
        "keystone_R": {
            "mat": core.mat_to_strlist(R),                 # [[0,1],[1,1]]
            "symmetric": bool(R == R.T),
            "golden_law_R2_eq_R_plus_I": bool(R * R == R + eye(2)),
            "charpoly": core.poly_str(sp.expand(R.charpoly(x).as_expr())),  # x^2-x-1
            "scalar_part_tau": str(core.cl_coords(R)[0]),  # 1/2
            "det": int(R.det()),                           # -1
            "TYPE": ["gen"],                               # satisfies golden law
        },
        "gate_P0": {
            "mat": core.mat_to_strlist(P0),                # [[1,0],[0,0]]
            "idempotent": bool(P0 * P0 == P0),
            "symmetric": bool(P0 == P0.T),
            "TYPE": ["idem", "rest"],
        },
        "cayley_hamilton_X2_minus_trX_plus_detI": {
            "on_keystone": ch(R), "on_generic_holding": ch(generic),
        },
    }
    core.write_json(SCRIPT, "carrier_readings.json", payload)
    return payload


def produce_return_operator():
    """Prop 5.8 / Thm 5.4 / Rem 5.5: L, its kernel, Sylvester eigenvalues, disproof."""
    L = core.return_operator_matrix()
    LR = core.return_operator(R)                           # L(R) = 5/2 + e1 - 1/2 e2
    eigs = L.eigenvals()
    ns = L.nullspace()
    # disproof (Rem 5.5)
    disproof = core.return_operator(core.mat(0, -1, 1, 0))     # -> H(-3,0,0,0)
    N = Matrix([[0, -1], [1, 0]])                              # cl(N) = i in kernel
    payload = {
        "return_operator": "L(X) = R X + X R - X",
        "L_matrix_basis_1_e1_e2_i": core.mat_to_strlist(L),
        "L_of_keystone_coords": [str(v) for v in core.cl_coords(LR)],  # (5/2,1,-1/2,0)
        "keystone_in_kernel": bool(LR == zeros(2)),                    # False
        "kernel_dimension": len(ns),                                  # 2
        "kernel_basis": ["e1+2e2 = H(0,1,2,0)", "i = H(0,0,0,1)"],
        "kernel_check_L_e1p2e2_zero": bool(L * Matrix([0, 1, 2, 0]) == zeros(4, 1)),
        "kernel_check_L_i_zero": bool(L * Matrix([0, 0, 0, 1]) == zeros(4, 1)),
        "sylvester_eigenvalues": {str(k): int(v) for k, v in eigs.items()},
        "zero_eigenvalue_multiplicity": int(eigs.get(sp.Integer(0), 0)),
        "disproof": {
            "claim": "N=H(0,-1,1,0) in ker L is FALSE",
            "L_of_H_0_m1_1_0_coords": [str(v) for v in core.cl_coords(disproof)],  # (-3,0,0,0)
            "correct_image_cl_of_[[0,-1],[1,0]]": [str(v) for v in core.cl_coords(N)],  # (0,0,0,1)=i
            "correct_image_in_kernel": bool(core.return_operator(N) == zeros(2)),
        },
    }
    core.write_json(SCRIPT, "return_operator_kernel.json", payload)
    return payload


def produce_commit_projector():
    """eq (proj) / Thm 5.6: the commit is the exact orthogonal projector onto ker L."""
    rows = [[i] + [str(COMMIT[i, j]) for j in range(4)] for i in range(4)]
    core.write_csv(SCRIPT, "commit_projector.csv",
                   ["row_coord", "1", "e1", "e2", "i"], rows)
    return COMMIT


def produce_lexicon():
    """Thm 5.6 / Ex 5.7: five tokens -> four distinct exact residues, with words."""
    tokens = [
        ("E1",   [0, 1, 0, 0]),
        ("E1+1", [1, 1, 0, 0]),
        ("i",    [0, 0, 0, 1]),
        ("2i",   [0, 0, 0, 2]),
        ("-i",   [0, 0, 0, -1]),
    ]
    seen, rows = {}, []
    for name, coord in tokens:
        committed = tuple(COMMIT * Matrix(coord))
        val_str = "H(" + ",".join(str(v) for v in committed) + ")"
        w = word(coord)
        if committed in seen:
            entry = seen[committed]
            status = f"merged into #{entry}"
        else:
            entry = len(seen) + 1
            seen[committed] = entry
            status = f"new entry #{entry}"
        rows.append([name, val_str, w, status])
    core.write_csv(SCRIPT, "lexicon.csv",
                   ["token", "committed_value", "word", "lexicon"], rows)
    return rows, len(seen)


def produce_firewall():
    """Table 6: the jurisdiction firewall counts."""
    theorem, computed, interpretive, false_as_stated = 20, 5, 2, 0
    rows = [
        ["theorem", "yes", theorem, "proved/exact-verified fact"],
        ["computed", "yes", computed, "derived/learned (lexicon entry, reading)"],
        ["interpretive", "no", interpretive, "glossed reading, on request, tagged"],
        ["false_as_stated", "no", false_as_stated, "defined for completeness; unused"],
        ["WIRED_TOTAL", "-", theorem + computed, "crosses by default (20+5)"],
        ["BANK_TOTAL", "-", theorem + computed + interpretive, "all statements (20+5+2)"],
    ]
    core.write_csv(SCRIPT, "firewall.csv",
                   ["jurisdiction", "wired", "count", "meaning"], rows)
    return rows


def main():
    car = produce_carrier()
    ro = produce_return_operator()
    produce_commit_projector()
    lex, n_entries = produce_lexicon()
    produce_firewall()
    print(f"[{SCRIPT}] keystone golden law: "
          f"{car['keystone_R']['golden_law_R2_eq_R_plus_I']}, "
          f"charpoly {car['keystone_R']['charpoly']}")
    print(f"[{SCRIPT}] ker L dim = {ro['kernel_dimension']}, "
          f"Sylvester eigs = {ro['sylvester_eigenvalues']}")
    print(f"[{SCRIPT}] disproof L(H(0,-1,1,0)) = "
          f"H({','.join(ro['disproof']['L_of_H_0_m1_1_0_coords'])})")
    print(f"[{SCRIPT}] lexicon: 5 tokens -> {n_entries} entries "
          f"({[(r[0], r[2]) for r in lex]})")
    print(f"[{SCRIPT}] wrote 5 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
