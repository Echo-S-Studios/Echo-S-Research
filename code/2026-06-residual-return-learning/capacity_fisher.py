r"""
Producer: Section 4 "Capacity and the applied growth gate" of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits:
  * Ex 4.2 (the gate decides on exact integers): for the seeds x^2-24 and x^2-7,
      degree, coefficient height, sum c_i^2 (Landau L2), exact Mahler measure, and
      the GROW/REJECT/STOP verdicts under several integer budgets
      -> data/northcott_gate.csv  and  data/northcott_gate.json
  * Ex 4.4 (exact Fisher matrices): Fisher_exp = (1/n)(G - (1/n) t t^T) for
      Q(sqrt5), Q(sqrt2,sqrt3), Q(sqrt2,sqrt3,sqrt7), and the trace-zero identity
      ||sqrt5||_G^2 = 10 = n*Fisher, with G = n*Fisher on every trace-zero direction
      -> data/fisher_matrices.json

Run:  py code/2026-06-residual-return-learning/capacity_fisher.py

Grams are rebuilt from field traces / the regular representation, and the gate is
an integer-only decision function -- the producer computes verdicts, not asserts.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt

import rrl_core as core

SCRIPT = "capacity_fisher.py"
x = core.X


def capacity_decision(deg, height, gain, Dmax, Hmax, floor=0):
    """The shipped exact gate (integers only): STOP if gain<=floor, else REJECT if
    inadmissible (deg>Dmax or height>Hmax), else GROW."""
    if gain <= floor:
        return "STOP"
    if deg > Dmax or height > Hmax:
        return "REJECT"
    return "GROW"


def produce_gate():
    """Ex 4.2: Northcott/Landau admissibility certificate + gate verdicts."""
    seeds = [
        ("2*sqrt6", [1, 0, -24], 96),
        ("sqrt7",   [1, 0, -7],  56),
    ]
    rows, records = [], []
    for name, coeffs, gain in seeds:
        deg = len(coeffs) - 1
        height = max(abs(c) for c in coeffs)
        sumsq = sum(c * c for c in coeffs)
        mah = core.mahler_measure(coeffs)
        rec = {
            "seed": name,
            "minpoly_coeffs": coeffs,
            "degree": deg,
            "coeff_height": height,
            "sum_of_squares_L2": sumsq,
            "landau_bound_Mah_le": str(sp.sqrt(sumsq)),
            "exact_mahler_measure": str(mah),
            "gain_residual_norm": gain,
        }
        records.append(rec)
        rows.append([name, str(coeffs).replace(" ", ""), deg, height, sumsq,
                     str(round(float(mah), 6)), gain])
    # gate verdict table for the x^2-24 seed under three budgets (Ex 4.2)
    verdicts = {
        "Budget(Dmax=64,Hmax=256), gain=96": capacity_decision(2, 24, 96, 64, 256),
        "Budget(Dmax=4,Hmax=10),  gain=96":  capacity_decision(2, 24, 96, 4, 10),
        "Budget(Dmax=64,Hmax=256), gain=0":  capacity_decision(2, 24, 0, 64, 256),
    }
    core.write_csv(SCRIPT, "northcott_gate.csv",
                   ["seed", "minpoly_coeffs", "degree", "coeff_height",
                    "sum_of_squares", "mahler_measure", "gain"], rows)
    core.write_json(SCRIPT, "northcott_gate.json",
                    {"seeds": records, "gate_verdicts_x2_minus_24": verdicts})
    return records, verdicts


def fisher_exp(G, tvec):
    """Fisher_exp = (1/n)(G - (1/n) t t^T), t the trace vector."""
    n = G.shape[0]
    t = Matrix(tvec)
    return (G - (t * t.T) / n) / n


def produce_fisher():
    """Ex 4.4: exact Fisher matrices and the trace-zero identity G = n*Fisher."""
    # Q(sqrt5): G rebuilt from companion(x^2-x-1) power traces; t=(Tr1,Tr phi)=(2,1)
    C = core.companion([1, -1, -1])
    G5 = Matrix([[(C**0).trace(), (C**1).trace()],
                 [(C**1).trace(), (C**2).trace()]])
    F5 = fisher_exp(G5, [2, 1])
    # residual sqrt5 = 2 phi - 1 = (-1,2), trace-zero, ||.||_G^2 = 10 = 2*Fisher
    v = Matrix([-1, 2])
    t5 = Matrix([2, 1])
    norm5 = int((v.T * G5 * v)[0])
    n_fisher5 = int(2 * (v.T * F5 * v)[0])

    # Q(sqrt2,sqrt3): G=diag(4,8,12,24), t=(4,0,0,0)
    G4 = sp.diag(4, 8, 12, 24)
    F4 = fisher_exp(G4, [4, 0, 0, 0])
    # Q(sqrt2,sqrt3,sqrt7): Kronecker order G=diag(8,56,16,112,24,168,48,336)
    G8 = Matrix(sp.kronecker_product(G4, sp.diag(2, 14)))
    F8 = fisher_exp(G8, [8, 0, 0, 0, 0, 0, 0, 0])

    payload = {
        "definition": "Fisher_exp = (1/n)(G - (1/n) t t^T), t = trace vector",
        "Q_sqrt5": {
            "n": 2, "G": core.mat_to_strlist(G5),
            "trace_vector_t": [2, 1],
            "Fisher_exp": core.mat_to_strlist(F5),
            "Fisher_diag_or_full": "full [[0,0],[0,5/4]]",
        },
        "Q_sqrt2_sqrt3": {
            "n": 4, "G_diag": [int(G4[i, i]) for i in range(4)],
            "trace_vector_t": [4, 0, 0, 0],
            "Fisher_exp_diag": [str(F4[i, i]) for i in range(4)],  # 0,2,3,6
        },
        "Q_sqrt2_sqrt3_sqrt7": {
            "n": 8, "G_diag": [int(G8[i, i]) for i in range(8)],
            "trace_vector_t": [8, 0, 0, 0, 0, 0, 0, 0],
            "Fisher_exp_diag": [str(F8[i, i]) for i in range(8)],  # 0,7,2,14,3,21,6,42
        },
        "trace_zero_identity": {
            "residual": "sqrt5 = 2*phi - 1 = (-1,2)",
            "is_trace_zero": bool((t5.T * v)[0] == 0),
            "residual_norm_G": norm5,             # 10
            "n_times_Fisher": n_fisher5,          # 10
            "identity_holds": bool(norm5 == n_fisher5 == 10),
            "note": "G = n*Fisher_exp on the trace-zero subspace only "
                    "(they differ at the constant); requires 1 in col(B)",
        },
    }
    core.write_json(SCRIPT, "fisher_matrices.json", payload)
    return payload


def main():
    recs, verdicts = produce_gate()
    fish = produce_fisher()
    print(f"[{SCRIPT}] gate seeds: "
          f"{[(r['seed'], r['degree'], r['coeff_height'], r['sum_of_squares_L2']) for r in recs]}")
    print(f"[{SCRIPT}] x^2-24 verdicts: {verdicts}")
    print(f"[{SCRIPT}] Q(sqrt2,sqrt3) Fisher diag = "
          f"{fish['Q_sqrt2_sqrt3']['Fisher_exp_diag']}")
    print(f"[{SCRIPT}] trace-zero identity holds: "
          f"{fish['trace_zero_identity']['identity_holds']}")
    print(f"[{SCRIPT}] wrote 3 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
