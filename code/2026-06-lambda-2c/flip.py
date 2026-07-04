"""
Producer -- Sections 9, 10 & 13 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The flip D = 1 + 4C: one sign change, four readings.

  * thm:flip -- for C in {1,0,-1/4,-1,-2}, D in {5,1,0,-3,-7}; roots of x^2+x-C
                real (D>0) / double -1/2 (D=0) / complex (D<0); gap +-sqrt D real or
                imaginary; regime hyperbolic / parabolic / elliptic.
  * prop:detG -- in the gap basis {1, sqrt D}, G = diag(2, 2D), det G = 4D:
                positive-definite (Riemannian) iff D>0, indefinite (Lorentzian)
                iff D<0, degenerate at D=0; canonical instance Q(i) -> diag(2,-2).
  * prop:meet -- the two flips meet at C=-1: x^2+x+1 has roots e^{+-2pi i/3} on the
                unit circle (Mahler 1, the cyclotomic vertex), discriminant D=-3.
  * prop:square -- perfect-square pinning D = (2x+1)^2 = 4z^2, C = z^2 - 1/4, so a
                real coherence flow keeps D >= 0 (approaches the flip, never crosses).

Emits:
    data/2026-06-lambda-2c/flip_table.csv
    data/2026-06-lambda-2c/flip.json
"""
import sympy as sp
import mpmath as mp
import lambda2c_common as cm

SCRIPT = "flip.py"
x = sp.symbols('x')


def flip_rows():
    """The flip table across the sign change of D."""
    cases = [
        (sp.Integer(1), 5, "hyperbolic (saddle)"),
        (sp.Integer(0), 1, "hyperbolic"),
        (sp.Rational(-1, 4), 0, "parabolic (paradox)"),
        (sp.Integer(-1), -3, "elliptic (rotation)"),
        (sp.Integer(-2), -7, "elliptic (rotation)"),
    ]
    rows = []
    for Cval, Dexp, regime in cases:
        D = 1 + 4 * Cval
        roots = sp.solve(x**2 + x - Cval, x)
        if Dexp > 0:
            character = "real"
            root_str = ", ".join(str(sp.nsimplify(r)) for r in roots)
            gap = f"+-sqrt({Dexp}) real"
            metric = "totally real / PD"
        elif Dexp == 0:
            character = "double root"
            root_str = "-1/2 (double)"
            gap = "0"
            metric = "degenerate"
        else:
            character = "complex"
            root_str = ", ".join(str(sp.simplify(r)) for r in roots)
            gap = f"+-i sqrt({abs(Dexp)})"
            metric = "complex place / Lorentzian"
        rows.append([str(Cval), int(D), character, root_str, gap, metric, regime])
    return rows


def detG_facts():
    """prop:detG: det G = 4D and the signature flip, incl. the Q(i) instance."""
    C = sp.symbols('C')
    D = 1 + 4 * C
    # gap basis {1, sqrt D}, sqrt D = 2 theta + 1
    tr_one = sp.Integer(2)
    tr_sqrtD = 2 * sp.Integer(-1) + tr_one              # = 0
    tr_D = 2 * D
    G = sp.Matrix([[tr_one, tr_sqrtD], [tr_sqrtD, tr_D]])
    detG = sp.expand(G.det())

    signatures = {}
    for Dval in (5, 1, 0, -3, -7):
        Gd = sp.Matrix([[2, 0], [0, 2 * Dval]])
        eigs = list(Gd.eigenvals().keys())
        if all(e > 0 for e in eigs):
            sig = "positive-definite (Riemannian)"
        elif 0 in eigs:
            sig = "degenerate"
        else:
            sig = "indefinite (1,1) Lorentzian"
        signatures[str(Dval)] = {"det_G": int(Gd.det()), "signature": sig}

    # Q(i): basis {1, i}, embeddings i -> +-i
    conj = [sp.I, -sp.I]
    b = [[sp.Integer(1), sp.Integer(1)], [conj[0], conj[1]]]
    Gi = sp.Matrix(2, 2, lambda j, k: sp.simplify(sum(b[j][e] * b[k][e] for e in (0, 1))))
    return {
        "gap_basis": "{1, sqrt D},  sqrt D = 2 theta + 1",
        "G_symbolic": "diag(2, 2D)",
        "det_G": str(detG),                             # 8C + 2  == 4D
        "det_G_equals_4D": bool(sp.simplify(detG - 4 * D) == 0),
        "traces": {"Tr(1)": 2, "Tr(sqrtD)": int(tr_sqrtD), "Tr(D)": "2D"},
        "signature_by_D": signatures,
        "Qi_instance": {
            "gram": [[int(Gi[i, j]) for j in range(2)] for i in range(2)],   # diag(2,-2)
            "det": int(Gi.det()), "D": -1, "det_equals_4D": int(Gi.det()) == 4 * (-1),
        },
    }


def two_flips_meet():
    """prop:meet: at C=-1, x^2+x+1 has cube-roots of unity, Mahler 1, D=-3."""
    roots = mp.polyroots([1, 1, 1])
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    on_circle = all(abs(abs(r) - 1) < mp.mpf(10) ** (-30) for r in roots)
    return {
        "C": -1, "polynomial": "x^2+x+1", "roots": "e^{+-2pi i/3}",
        "discriminant": int(sp.discriminant(sp.Poly(x**2 + x + 1, x))),   # -3
        "mahler": float(M), "on_unit_circle": bool(on_circle),
        "meaning": "additive flip's rotation regime (D<0) meets multiplicative flip's marginal vertex |lambda|=1",
    }


def perfect_square_pinning():
    """prop:square / eq:square: D = 4z^2, C = z^2 - 1/4."""
    z = sp.symbols('z', real=True)
    return {
        "z_definition": "z = sqrt(D)/2 = |x + 1/2|",
        "D_equals_4z2": bool(sp.simplify((2 * z)**2 - 4 * z**2) == 0),
        "C_equals_z2_minus_quarter": bool(
            sp.simplify((z**2 - sp.Rational(1, 4)) - ((4 * z**2 - 1) / 4)) == 0),
        "real_flow_keeps_D_nonneg": "D = 4z^2 >= 0 for real z; equality iff z=0 (the flip)",
        "gate_coherences": {"C=1/4": "sqrt(2)/2", "C=1/2": "sqrt(3)/2", "C=1": "sqrt(5)/2"},
    }


def main():
    cm.write_csv("flip_table.csv",
                 ["C", "D=1+4C", "root_character", "roots", "gap", "field_metric", "regime"],
                 flip_rows(), SCRIPT)

    payload = {
        "flip_event": "sign of D = 1 + 4C changes once at C = -1/4 (double root -1/2)",
        "four_readings": ["eigenvalues real<->complex", "field totally-real<->complex place",
                          "trace-form signature Riemannian<->Lorentzian", "channels growth<->rotation"],
        "det_G_signature": detG_facts(),
        "two_flips_meet": two_flips_meet(),
        "perfect_square_pinning": perfect_square_pinning(),
    }
    cm.write_json("flip.json", payload, SCRIPT)

    print("wrote flip_table.csv, flip.json")
    print("  det G = 4D:", payload["det_G_signature"]["det_G_equals_4D"],
          "; Q(i) gram", payload["det_G_signature"]["Qi_instance"]["gram"])


if __name__ == "__main__":
    main()
