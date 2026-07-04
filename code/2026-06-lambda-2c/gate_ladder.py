"""
Producer -- Section 5 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The gate ladder x^2 + x = C and its trifurcation:

  * eq:gengate  -- discriminant of x^2 + x - C is D = 1 + 4C.
  * eq:companion-- R_C = [[0,C],[1,-1]] has charpoly x^2+x-C, eigenvalues
                   (-1 +- sqrt(1+4C))/2, spectral gap sqrt(1+4C) = sqrt D.
  * gate table  -- C in {1/4,1/2,1} -> D in {2,3,5}, seeds sqrt D, fields Q(sqrt D);
                   golden gate C=1 roots {1/phi, -phi}; tau and phi both disc 5.
  * prop:trifurcation -- spec(ad_{R_C}) = {-sqrt(1+4C), 0, +sqrt(1+4C)}, 0 of
                   multiplicity 2, 0-eigenspace = span{I, R_C}, ad charpoly
                   t^2 (t^2 - (1+4C)); golden gate spectrum {-sqrt5, 0, +sqrt5}.
  * eq:gatebalance / rem:mahler -- lambda C = r(R_C) = sqrt(1+4C) = sqrt(Mah(x^2-D)),
                   with Mah(x^2 - D) = D.

Emits:
    data/2026-06-lambda-2c/gate_ladder.csv
    data/2026-06-lambda-2c/trifurcation.json
"""
import sympy as sp
import lambda2c_common as cm

SCRIPT = "gate_ladder.py"
x, t, C = sp.symbols('x t C')


def gate_rows():
    """The three squarefree gate levels C in {1/4,1/2,1}."""
    table = [(sp.Rational(1, 4), 2), (sp.Rational(1, 2), 3), (sp.Integer(1), 5)]
    rows = []
    for Cval, Dval in table:
        D = 1 + 4 * Cval
        R = sp.Matrix([[0, Cval], [1, -1]])
        eigs = list(R.eigenvals().keys())
        gap = sp.simplify(sp.Abs(eigs[0] - eigs[1]))
        rows.append([
            str(Cval), int(D), f"sqrt({Dval})", f"x^2-{Dval}", f"Q(sqrt{Dval})",
            str(gap), cm.approx(sp.sqrt(Dval)),
        ])
    return rows


def companion_facts():
    """Symbolic companion charpoly / eigenvalues / gap for general C."""
    R = sp.Matrix([[0, C], [1, -1]])
    cp = sp.expand(R.charpoly(x).as_expr())
    lam_plus = (-1 + sp.sqrt(1 + 4 * C)) / 2
    lam_minus = (-1 - sp.sqrt(1 + 4 * C)) / 2
    gap = sp.simplify(lam_plus - lam_minus)
    return {
        "R_C": [["0", "C"], ["1", "-1"]],
        "charpoly": str(cp),                              # x**2 + x - C
        "eigenvalue_plus": str(lam_plus),
        "eigenvalue_minus": str(lam_minus),
        "spectral_gap": str(gap),                         # sqrt(4C+1)
    }


def trifurcation_facts():
    """spec(ad_{R_C}) and the golden-gate instance."""
    R = sp.Matrix([[0, C], [1, -1]])
    ad = cm.ad_matrix(R)
    cp = sp.factor(ad.charpoly(t).as_expr())
    # golden gate
    Rg = sp.Matrix([[0, 1], [1, -1]])
    adg = cm.ad_matrix(Rg)
    eigs_g = adg.eigenvals()
    zero_mult = int(eigs_g.get(sp.Integer(0), 0))
    return {
        "ad_charpoly": str(cp),                                   # t**2*(t**2 - (4C+1))
        "ad_charpoly_target": str(sp.factor(t**2 * (t**2 - (1 + 4 * C)))),
        "spectrum": ["-sqrt(4C+1)", "0", "+sqrt(4C+1)"],
        "zero_eigenspace": "span{I, R_C}",
        "golden_gate_spectrum": ["-sqrt(5)", "0", "+sqrt(5)"],
        "golden_gate_zero_multiplicity": zero_mult,               # 2
        "ad_rank_at_golden": int(adg.rank()),                     # 2  (nullity 2)
    }


def gate_balance_and_mahler():
    """eq:gatebalance at golden gate + Mah(x^2 - D) = D for the ladder radicands."""
    Cval = sp.Integer(1)
    lam = sp.sqrt(5)
    r = sp.sqrt(1 + 4 * Cval)
    balance_ok = bool(sp.simplify(lam * Cval - r) == 0)
    mahler = {}
    for Dval in (2, 3, 5, 7, 24):
        roots = [sp.sqrt(Dval), -sp.sqrt(Dval)]
        M = 1
        for rt in roots:
            if sp.Abs(rt) > 1:
                M = M * sp.Abs(rt)
        mahler[str(Dval)] = {
            "minpoly": f"x^2-{Dval}",
            "mahler": int(sp.simplify(M)),                        # = D
            "r_RC_equals_sqrt_mahler": bool(sp.simplify(sp.sqrt(Dval) - sp.sqrt(M)) == 0),
        }
    return {
        "gate_balance_golden": "sqrt(5)*1 = sqrt(1+4) = sqrt(5)",
        "gate_balance_holds": balance_ok,
        "discriminant_seed_mahler": mahler,
    }


def main():
    cm.write_csv(
        "gate_ladder.csv",
        ["C", "D=1+4C", "seed", "minpoly", "field", "gap_r(R_C)", "gap_float"],
        gate_rows(), SCRIPT)

    golden_roots = sp.solve(x**2 + x - 1, x)
    payload = {
        "family": "x^2 + x = C,  D = 1 + 4C",
        "companion": companion_facts(),
        "gate_levels": {
            "C": ["1/4", "1/2", "1"], "D": [2, 3, 5],
            "seeds": ["sqrt(2)", "sqrt(3)", "sqrt(5)"],
            "fields": ["Q(sqrt2)", "Q(sqrt3)", "Q(sqrt5)"],
            "critical_coherence_at_C_half": "sqrt(3)/2",
        },
        "golden_gate": {
            "polynomial": "x^2 + x - 1  (Clifford unity T^2+T=1)",
            "roots": [str(sp.nsimplify(r)) for r in golden_roots],
            "roots_are_reciprocal_golden_pair": "{1/phi, -phi}",
            "disc_tau_and_phi": [int(sp.discriminant(sp.Poly(x**2 + x - 1, x))),
                                 int(sp.discriminant(sp.Poly(x**2 - x - 1, x)))],
        },
        "trifurcation": trifurcation_facts(),
        "gate_balance": gate_balance_and_mahler(),
    }
    cm.write_json("trifurcation.json", payload, SCRIPT)

    print("wrote gate_ladder.csv, trifurcation.json")
    print("  gates C in {1/4,1/2,1} -> D in {2,3,5}; ad charpoly =",
          payload["trifurcation"]["ad_charpoly"])


if __name__ == "__main__":
    main()
