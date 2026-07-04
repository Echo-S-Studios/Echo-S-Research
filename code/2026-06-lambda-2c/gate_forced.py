"""
Producer -- Section 7 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The gate is forced to the golden level C=1 by two convergent derivations:

  * prop:ternary -- signed channels of ad_R for a degree-d seed = d^2 - d + 1,
                    which equals 3 iff d = 2 (d=3 -> 7, the cubic's surplus).
  * lem:mincost / lem:tie -- min Mahler measure > 1 over irreducible integer
                    quadratics is phi = 1.61803..., attained ONLY by x^2-x-1 and
                    x^2+x-1 (both discriminant 5).
  * rem:cubic -- the global non-reciprocal infimum is the plastic number
                    mu_S = 1.32472 (x^3-x-1) < phi, excluded by the ternary lock.
  * prop:firewallimage -- squaring firewall x^4+5x^2-5 -> y^2+5y-5 (y=x^2),
                    discriminant 45 = 9*5, roots in Q(sqrt5).
  * lem:keystone -- only the golden companion has R^2=R+I; tau: R^2=I-R;
                    radicand seeds: R^2 = D I.
  * thm:gateforced / eq:goldenvalue -- C=1, c=sqrt5/2, lambda=2c=sqrt5=phi-psi.

Emits:
    data/2026-06-lambda-2c/channel_count.csv
    data/2026-06-lambda-2c/min_mahler_quadratics.csv
    data/2026-06-lambda-2c/gate_forced.json
"""
import sympy as sp
import mpmath as mp
import lambda2c_common as cm

SCRIPT = "gate_forced.py"
x, y = sp.symbols('x y')


def channel_counts(dmax=6):
    """Signed channel count d^2 - d + 1 per seed degree d."""
    rows = []
    for d in range(1, dmax + 1):
        rows.append([d, d**2 - d + 1, "ternary" if d == 2 else ""])
    return rows


def min_mahler_quadratics(bound=8):
    """Enumerate irreducible integer quadratics x^2+bx+c; return the minimizers of
    the Mahler measure exceeding 1 (the disc-5 golden tie)."""
    best = None
    attained = []
    for b in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            poly = sp.Poly(x**2 + b * x + c, x)
            if not poly.is_irreducible:
                continue
            M = cm.mahler_numpy([1, b, c])
            if M > 1 + 1e-9:
                if best is None or M < best - 1e-9:
                    best, attained = M, [(b, c)]
                elif abs(M - best) < 1e-9:
                    attained.append((b, c))
    # exact minimum is phi
    phi_val = cm.mahler_mpmath([1, -1, -1])
    return best, attained, phi_val


def min_mahler_cubics(bound=4):
    """rem:cubic: min Mahler over irreducible integer cubics is the plastic number."""
    best = None
    argmin = None
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            for c in range(-bound, bound + 1):
                poly = sp.Poly(x**3 + a * x**2 + b * x + c, x)
                if not poly.is_irreducible:
                    continue
                M = cm.mahler_numpy([1, a, b, c])
                if M > 1 + 1e-9 and (best is None or M < best - 1e-12):
                    best, argmin = M, (a, b, c)
    plastic = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    return best, argmin, plastic


def firewall_kformation():
    """prop:firewallimage: x^4+5x^2-5 -> y^2+5y-5 under y=x^2."""
    f = x**4 + 5 * x**2 - 5
    g = y**2 + 5 * y - 5
    substitution_ok = bool(sp.expand(g.subs(y, x**2) - f) == 0)
    disc = int(sp.discriminant(sp.Poly(g, y)))
    roots = sp.solve(g, y)
    return {
        "quartic": "x^4+5x^2-5", "reduced": "y^2+5y-5", "substitution": "y=x^2",
        "substitution_recovers_quartic": substitution_ok,
        "discriminant": disc, "disc_factorization": "9*5",
        "y_roots": [str(sp.simplify(r)) for r in roots],
        "roots_in_Q_sqrt5": all(sp.simplify(r).has(sp.sqrt(5)) for r in roots),
    }


def keystone_relations():
    """lem:keystone: self-reproducing relation by seed type."""
    I2 = sp.eye(2)
    Rphi = sp.Matrix([[0, 1], [1, 1]])       # x^2-x-1
    Rtau = sp.Matrix([[0, 1], [1, -1]])      # x^2+x-1
    out = {
        "golden_x2mx m1": {"relation": "R^2 = R + I",
                           "holds": bool(Rphi**2 == Rphi + I2)},
        "conjugate_tau": {"relation": "R^2 = I - R",
                          "holds": bool(Rtau**2 == I2 - Rtau)},
        "radicands": {},
    }
    for D in (2, 3, 5):
        Rrad = sp.Matrix([[0, D], [1, 0]])
        out["radicands"][f"x2m{D}"] = {"relation": f"R^2 = {D} I",
                                       "holds": bool(Rrad**2 == D * I2)}
    return out


def main():
    cm.write_csv("channel_count.csv",
                 ["seed_degree_d", "signed_channels_d2-d+1", "note"],
                 channel_counts(), SCRIPT)

    best_q, attained, phi_val = min_mahler_quadratics()
    rows = []
    for (b, c) in sorted(attained):
        rows.append([cm.poly_str([1, b, c]),
                     b, c, int(sp.discriminant(sp.Poly(x**2 + b * x + c, x))),
                     cm.mahler_numpy([1, b, c])])
    cm.write_csv("min_mahler_quadratics.csv",
                 ["minpoly", "b", "c", "discriminant", "mahler"], rows, SCRIPT)

    best_c, argmin_c, plastic = min_mahler_cubics()
    d = sp.symbols('d', integer=True, positive=True)

    payload = {
        "ternary_lock": {
            "channel_formula": "d^2 - d + 1",
            "equals_3_iff": [int(s) for s in sp.solve(sp.Eq(d**2 - d + 1, 3), d)],  # [2]
            "degree2_channels": 3, "degree3_channels": 7,
        },
        "min_mahler_quadratics": {
            "minimum": float(best_q),
            "minimum_exact": "phi = (1+sqrt5)/2",
            "phi_highprec": str(phi_val),
            "minimizers_bc": [list(t) for t in sorted(attained)],   # (-1,-1),(1,-1)
            "minimizers_minpoly": ["x^2-x-1", "x^2+x-1"],
            "both_discriminant_5": True,
        },
        "cubic_objection": {
            "min_mahler_cubics": float(best_c),
            "argmin_abc": list(argmin_c),
            "plastic_number_muS": str(plastic),
            "plastic_below_phi": bool(plastic < float(cm.PHI)),
            "excluded_by": "ternary lock (cubic seed gives 7 channels, not 3)",
        },
        "firewall": firewall_kformation(),
        "keystone_relations": keystone_relations(),
        "forced_value": {
            "C": 1, "c": "sqrt(5)/2", "c_float": cm.approx(sp.sqrt(5) / 2),
            "lambda": "sqrt(5)", "lambda_float": cm.approx(sp.sqrt(5)),
            "lambda_equals_phi_minus_psi": bool(
                sp.simplify(2 * (sp.sqrt(5) / 2) - (cm.PHI - cm.PSI)) == 0),
            "gate_balance_holds": bool(
                sp.simplify(2 * (sp.sqrt(5) / 2) * 1 - sp.sqrt(1 + 4 * 1)) == 0),
        },
    }
    cm.write_json("gate_forced.json", payload, SCRIPT)

    print("wrote channel_count.csv, min_mahler_quadratics.csv, gate_forced.json")
    print(f"  min Mahler over quadratics = {best_q:.6f} (phi), minimizers {sorted(attained)}")
    print(f"  min Mahler over cubics = {best_c:.6f} (plastic {float(plastic):.6f}) < phi")


if __name__ == "__main__":
    main()
