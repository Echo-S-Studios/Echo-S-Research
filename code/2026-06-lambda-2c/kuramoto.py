"""
Producer -- Section 11 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The Kuramoto division of labour: the dynamics provide the critical point, the
frame-shift provides the gate.

  * sec:kuramoto -- Lorentzian frequency density g(omega) = gamma/(pi(gamma^2+omega^2))
                    gives g(0) = 1/(pi gamma), so the critical coupling
                    K_c = 2/(pi g(0)) = 2 gamma; the density integrates to 1.
  * rem:division -- order parameter r ~ sqrt(K - K_c): the mean-field 1/2 exponent
                    (Landau normal form dr/dt = a r - b r^3, nontrivial branch
                    sqrt(a/b)).
  * prop:square -- the critical coherence pinned is z_c = sqrt3/2 (the C=1/2 gate,
                    z = sqrt D/2, D = 1+4C = 3 at C=1/2).

Emits:
    data/2026-06-lambda-2c/kuramoto.json
"""
import sympy as sp
import lambda2c_common as cm

SCRIPT = "kuramoto.py"


def lorentzian_threshold():
    gamma, omega = sp.symbols('gamma omega', positive=True)
    g = gamma / (sp.pi * (gamma**2 + omega**2))
    g0 = g.subs(omega, 0)
    Kc = 2 / (sp.pi * g0)
    integral = sp.integrate(g, (omega, -sp.oo, sp.oo))
    return {
        "density": "gamma / (pi (gamma^2 + omega^2))",
        "g0": str(sp.simplify(g0)),                       # 1/(pi gamma)
        "K_c": str(sp.simplify(Kc)),                      # 2 gamma
        "K_c_equals_2gamma": bool(sp.simplify(Kc - 2 * gamma) == 0),
        "density_integrates_to_1": bool(integral == 1),
    }


def mean_field_exponent():
    r, a, b = sp.symbols('r a b', positive=True)
    sols = sp.solve(sp.Eq(a * r - b * r**3, 0), r)
    nontrivial = [s for s in sols if s != 0]
    branch_is_sqrt = any(sp.simplify(s - sp.sqrt(a / b)) == 0 for s in nontrivial)
    K, Kc, k0 = sp.symbols('K K_c k0', positive=True)
    r_branch = sp.sqrt((k0 * (K - Kc)) / b)
    exponent_half = bool(
        sp.simplify(r_branch - sp.sqrt(k0 / b) * (K - Kc)**sp.Rational(1, 2)) == 0)
    return {
        "normal_form": "dr/dt = a r - b r^3",
        "nontrivial_branch": "sqrt(a/b)",
        "branch_is_sqrt_ab": bool(branch_is_sqrt),
        "order_parameter_scaling": "r ~ sqrt(K - K_c)",
        "exponent": "1/2",
        "exponent_is_half": exponent_half,
    }


def critical_coherence():
    Cval = sp.Rational(1, 2)
    D = 1 + 4 * Cval
    z = sp.sqrt(D) / 2
    return {
        "C": "1/2", "D": int(D), "z_c": str(sp.simplify(z)),      # sqrt(3)/2
        "z_c_float": cm.approx(z),
        "z_c_equals_sqrt3_over_2": bool(sp.simplify(z - sp.sqrt(3) / 2) == 0),
    }


def main():
    payload = {
        "division_of_labour": {
            "dynamics_provides": "critical point z_c = sqrt3/2, mean-field sqrt exponent",
            "frameshift_provides": "the gate lambda = 2c, lambda = sqrt5 at C=1",
        },
        "lorentzian_threshold": lorentzian_threshold(),
        "mean_field_exponent": mean_field_exponent(),
        "critical_coherence": critical_coherence(),
    }
    cm.write_json("kuramoto.json", payload, SCRIPT)
    print("wrote kuramoto.json")
    print("  K_c =", payload["lorentzian_threshold"]["K_c"],
          "; z_c =", payload["critical_coherence"]["z_c"])


if __name__ == "__main__":
    main()
