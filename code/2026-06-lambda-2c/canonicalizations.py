"""
Producer -- Sections 4 & 6 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The three canonicalizations of the conformal scale c and the frame-shift value:

  * canon table -- c in {1 (Jeffreys volume-matching), n (degree-invariant
                   significance), sqrt(1+4C)/(2C) (frame-shift / eigenframe gap)};
                   cost floors 2 log muS, 2n log muS, and the spectral gap itself.
  * def:frameshift / eq:frameshift -- solving the gate balance 2 c C = sqrt(1+4C)
                   gives c = sqrt(1+4C)/(2C); at the golden gate C=1,
                   c = sqrt5/2 and lambda = 2c = sqrt5 = phi - psi (the spectral gap).

Emits:
    data/2026-06-lambda-2c/canonicalizations.csv
    data/2026-06-lambda-2c/frameshift.json
"""
import sympy as sp
import lambda2c_common as cm

SCRIPT = "canonicalizations.py"


def frameshift_solution():
    """Solve 2 c C = sqrt(1+4C) for c (eq:frameshift)."""
    Cpos = sp.symbols('C', positive=True)
    c = sp.symbols('c', positive=True)
    sol = sp.solve(sp.Eq(2 * c * Cpos, sp.sqrt(1 + 4 * Cpos)), c)
    return Cpos, sol


def main():
    Cpos, sol = frameshift_solution()
    c_frame = sp.sqrt(1 + 4 * Cpos) / (2 * Cpos)
    solved_matches = len(sol) == 1 and sp.simplify(sol[0] - c_frame) == 0

    c_at_1 = sp.simplify(c_frame.subs(Cpos, 1))
    lam_at_1 = sp.simplify(2 * c_at_1)

    # canon table (CSV)
    rows = [
        ["Jeffreys volume-matching", "1", "2 log muS"],
        ["degree-invariant significance", "n", "2n log muS"],
        ["frame-shift (eigenframe gap)", "sqrt(1+4C)/(2C)", "spectral gap; lambda=sqrt5 at C=1"],
    ]
    cm.write_csv("canonicalizations.csv",
                 ["principle", "selects_c", "cost_floor"], rows, SCRIPT)

    payload = {
        "frameshift_equation": "2 c C = sqrt(1+4C)",
        "frameshift_solution": str(sol[0]) if sol else None,      # sqrt(4C+1)/(2C)
        "frameshift_matches_closed_form": bool(solved_matches),
        "golden_gate": {
            "C": 1,
            "c": str(c_at_1),                                     # sqrt(5)/2
            "c_float": cm.approx(c_at_1),
            "lambda": str(lam_at_1),                              # sqrt(5)
            "lambda_float": cm.approx(lam_at_1),
            "lambda_equals_spectral_gap": bool(
                sp.simplify(lam_at_1 - (cm.PHI - cm.PSI)) == 0),  # phi - psi = sqrt5
            "phi_minus_psi": str(sp.simplify(cm.PHI - cm.PSI)),
        },
        "canonicalizations": {
            "jeffreys": {"c": "1", "status": "declared", "cost_floor": "2 log muS"},
            "degree": {"c": "n", "status": "declared", "cost_floor": "2n log muS"},
            "frameshift": {"c": "sqrt(1+4C)/(2C)", "status": "forced (gate at C=1)",
                           "cost_floor": "spectral gap sqrt(1+4C)"},
        },
    }
    cm.write_json("frameshift.json", payload, SCRIPT)

    print("wrote canonicalizations.csv, frameshift.json")
    print("  frame-shift c =", payload["frameshift_solution"],
          "-> at C=1:  c =", payload["golden_gate"]["c"],
          ", lambda =", payload["golden_gate"]["lambda"])


if __name__ == "__main__":
    main()
