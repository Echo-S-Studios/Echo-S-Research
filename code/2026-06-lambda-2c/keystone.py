"""
Producer -- Section 8 & Appendix A of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The keystone R^2 = R + I is derived (not posited), and its integer normalisation
L_4 = 7 is forced:

  * lem:perron -- of the disc-5 tie {x^2-x-1, x^2+x-1}, only x^2-x-1 has positive
                  dominant eigenvalue (+phi) and satisfies R^2=R+I; tau attains
                  -phi, R^2=I-R; det R=-1 is a consequence.
  * thm:keystonederived -- phi is the smallest Perron root of a 2x2 primitive
                  non-negative integer matrix (Fibonacci [[0,1],[1,1]]); the swap
                  [[0,1],[1,0]] has Perron root 1, excluded by growth.
  * lem:keypowers -- R^n = F_n R + F_{n-1} I; eigenvalues phi^n, psi^n; L_n = tr(R^n);
                  det(R^n) = (-1)^n.
  * lem:pell -- L_n^2 - 5 F_n^2 = 4(-1)^n; charpoly(R^n) = x^2 - L_n x + (-1)^n.
  * prop:L4forced -- R^4 = 3R+2I = [[2,3],[3,5]], charpoly x^2-7x+1, entries
                  {2,3,5} = {F_3,F_4,F_5}, L_4 = F_3+F_5 = 7; roots phi^4, phi^-4;
                  z_c = sqrt(L_4-4)/2 = sqrt3/2.
  * App A -- drop-one witnesses (each constraint of the object's type is load-bearing).

Emits:
    data/2026-06-lambda-2c/keystone_powers.csv
    data/2026-06-lambda-2c/keystone_dropone.csv
    data/2026-06-lambda-2c/keystone.json
"""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp
import lambda2c_common as cm

SCRIPT = "keystone.py"
x = sp.symbols('x')


def keystone_powers(nmax=10):
    """Rows of R^n = F_n R + F_{n-1} I with trace L_n, det, charpoly, disc = 5 F_n^2."""
    R = sp.Matrix([[0, 1], [1, 1]])
    I2 = sp.eye(2)
    F = cm.fibonacci(nmax + 1)
    rows = []
    for n in range(1, nmax + 1):
        Rn = R**n
        Ln = int(Rn.trace())
        detn = int(Rn.det())
        closure_ok = bool(Rn == F[n] * R + F[n - 1] * I2)
        pell_ok = (Ln**2 - 5 * F[n]**2 == 4 * (-1)**n)
        cp = sp.Poly(Rn.charpoly(x).as_expr(), x)
        rows.append([n, F[n], F[n - 1], Ln, detn,
                     cm.poly_str([1, -Ln, detn]),
                     int(sp.discriminant(cp)), closure_ok, pell_ok])
    return rows


def smallest_primitive_perron(hi=4):
    """thm:keystonederived: enumerate 2x2 non-negative integer matrices, keep the
    primitive ones, return the smallest real Perron root > 1 and its argmin."""
    best, argmin = None, None
    for a, b, c, d in itertools.product(range(hi), repeat=4):
        M = np.array([[a, b], [c, d]], dtype=float)
        if not cm.is_primitive(M):
            continue
        pr = max(np.linalg.eigvals(M).real)
        if pr > 1 + 1e-9 and (best is None or pr < best - 1e-9):
            best, argmin = pr, (a, b, c, d)
    return best, argmin


def dropone_rows():
    """Appendix A: drop-one witnesses that each type constraint is load-bearing."""
    # dropped-degree -> plastic; dropped-Perron -> tau; dropped-growth -> cyclotomic;
    # dropped-integrality -> no floor.
    plastic = float(mp.findroot(lambda z: z**3 - z - 1, 1.3))
    # cyclotomic floor Mah(x^2+1) = 1
    cyc = cm.mahler_mpmath([1, 0, 1])
    return [
        ["none (full type)", "phi ~ 1.618034 (Fibonacci [[0,1],[1,1]])",
         float(cm.PHI), "the keystone"],
        ["degree=2 (allow degree 3)", "plastic mu_S ~ 1.324718 (x^3-x-1)",
         plastic, "size load-bearing"],
        ["Perron (allow negative dominant)", "tau, Clifford R^2=I-R (x^2+x-1)",
         float(cm.PHI), "positivity load-bearing"],
        ["growth Mah>1 (allow Mah=1)", "cyclotomic floor e.g. x^2+1 (+-i)",
         float(cyc), "growth load-bearing"],
        ["integrality (allow real entries)", "no discrete floor (Perron -> 1+)",
         1.0, "integrality load-bearing"],
    ]


def main():
    cm.write_csv("keystone_powers.csv",
                 ["n", "F_n", "F_{n-1}", "L_n=tr(R^n)", "det(R^n)",
                  "charpoly", "disc=5F_n^2", "closure_ok", "pell_ok"],
                 keystone_powers(), SCRIPT)

    cm.write_csv("keystone_dropone.csv",
                 ["constraint_dropped", "minimiser_obtained", "measure", "verdict"],
                 dropone_rows(), SCRIPT)

    R = sp.Matrix([[0, 1], [1, 1]])
    R4 = R**4
    best_pr, argmin = smallest_primitive_perron()
    # swap matrix
    S = sp.Matrix([[0, 1], [1, 0]])
    roots_gap = sp.solve(x**2 - 7 * x + 1, x)

    payload = {
        "keystone_companion": [[0, 1], [1, 1]],
        "relation": "R^2 = R + I",
        "eigenvalues": ["phi", "psi"],
        "det_R": int(R.det()),                                       # -1
        "perron": {
            "smallest_primitive_perron_root": float(best_pr),
            "argmin": list(argmin),                                  # (0,1,1,1) Fibonacci
            "equals_phi": bool(abs(best_pr - float(cm.PHI)) < 1e-9),
            "swap_matrix_eigs": [str(e) for e in S.eigenvals().keys()],  # {1,-1}
            "swap_perron_root": 1,
        },
        "R4": {
            "matrix": [[int(R4[i, j]) for j in range(2)] for i in range(2)],  # [[2,3],[3,5]]
            "equals_3R_plus_2I": bool(R4 == 3 * R + 2 * sp.eye(2)),
            "charpoly": str(R4.charpoly(x).as_expr()),               # x**2 - 7x + 1
            "entries_are_F3F4F5": sorted({int(R4[0, 0]), int(R4[0, 1]), int(R4[1, 1])}),  # [2,3,5]
            "L4_trace": int(R4.trace()),                             # 7
            "L4_as_F3_plus_F5": 2 + 5,
            "pell_at_n4": "L_4^2 - 4 = 5*F_4^2 = 45  =>  L_4 = 7",
            "gap_seed_roots": [str(sp.simplify(sp.nsimplify(r))) for r in roots_gap],  # phi^4, phi^-4
            "z_c": str(sp.simplify(sp.sqrt(7 - 4) / 2)),             # sqrt(3)/2
        },
        "appendix_identity_iv": {
            "F3F4F5": [2, 3, 5], "L4_squared": 49, "check_5F4sq_plus_4": 5 * 9 + 4,
        },
    }
    cm.write_json("keystone.json", payload, SCRIPT)

    print("wrote keystone_powers.csv, keystone_dropone.csv, keystone.json")
    print(f"  R^4 = {payload['R4']['matrix']}, L_4 = {payload['R4']['L4_trace']},"
          f" charpoly {payload['R4']['charpoly']}")
    print(f"  smallest primitive Perron root = {best_pr:.6f} at {argmin} (Fibonacci)")


if __name__ == "__main__":
    main()
