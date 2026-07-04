"""
Producer: the golden limit and its convergence rate (Section 6).

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * data/prop64_taylor.json  (Prop 6.4 "Linear rate, slope phi^-1", Thm 6.2/6.6)
      - golden limit: beta_n -> phi, tau0(beta_n) -> phi+1/phi = sqrt5   (Thm 6.2)
      - trace'(phi) = 1 - phi^-2 = phi^-1 (linear slope)                 (Prop 6.4)
      - curvature magnitude (1/2) trace''(phi) = sqrt5-2 = phi^-3        (Prop 6.4)
      - the CORRECTED Taylor expansion of (sqrt5 - tau0) in u=(phi-beta):
            sqrt5 - tau0 = phi^-1 * u  -  (sqrt5-2) * u^2 + O(u^3)
        i.e. the QUADRATIC coefficient is  -(sqrt5-2) = -phi^-3  (negative;
        the errata sign correction of 2026-07-04)
      - proof-sketch identities: P* = 1-x-x^2, P*(phi) = -2 phi,
        P'(phi) = 2 phi - 1 = sqrt5                                      (Thm 6.6)
  * data/geometric_rate.csv  (Thm 6.6 table, n = 9,12,...,27 + limit row)
        n, gap = sqrt5 - tau0(beta_n), gap*phi^n, consecutive_ratio gap(n-1)/gap(n)

Run:  py code/2026-06-salem-slot/make_golden_rate.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

x = sc.x
mp = sc.mp
phi = sc.phi_sym


def golden_limit():
    """Thm 6.2: beta_n -> phi from below, tau0 -> sqrt5."""
    mp.mp.dps = 50
    b = sc.beta_n(30)
    tau0 = b + 1 / b
    return {
        "beta_30": mp.nstr(b, 16),
        "phi": mp.nstr(sc.PHI(), 16),
        "approaches_from_below": bool(b < sc.PHI()),
        "tau0_30": mp.nstr(tau0, 16),
        "sqrt5": mp.nstr(sc.SQRT5(), 16),
        "tau0_to_sqrt5_gap": mp.nstr(sc.SQRT5() - tau0, 6),
    }


def linear_rate_and_taylor():
    """Prop 6.4: slope phi^-1, curvature magnitude sqrt5-2, and the corrected
    quadratic sign in the expansion of (sqrt5 - tau0)."""
    b, u = sp.symbols('b u')
    # slope and curvature of the trace map itself
    dtr = sp.diff(b + 1 / b, b)
    d2tr = sp.diff(b + 1 / b, b, 2)
    slope = sp.simplify(dtr.subs(b, phi))
    curv = sp.simplify(d2tr.subs(b, phi) / 2)
    # Taylor of (sqrt5 - tau0) at beta = phi in u = (phi - beta)
    expr = sp.sqrt(5) - (b + 1 / b)
    ser = sp.series(expr.subs(b, phi - u), u, 0, 4).removeO()
    c1 = sp.simplify(ser.coeff(u, 1))
    c2 = sp.simplify(ser.coeff(u, 2))
    c3 = sp.simplify(ser.coeff(u, 3))
    return {
        "trace_prime": "trace'(beta) = 1 - 1/beta^2",
        "trace_prime_at_phi": str(sp.nsimplify(slope)),
        "slope_equals_phi_inverse": sp.simplify(slope - 1 / phi) == 0,
        "slope_equals_sqrt5_minus_1_over_2": sp.simplify(slope - (sp.sqrt(5) - 1) / 2) == 0,
        "curvature_half_trace_second_at_phi": str(sp.nsimplify(curv)),
        "curvature_equals_sqrt5_minus_2": sp.simplify(curv - (sp.sqrt(5) - 2)) == 0,
        "curvature_equals_phi_inverse_cubed": sp.simplify(curv - 1 / phi ** 3) == 0,
        "taylor_of_sqrt5_minus_tau0": {
            "variable": "u = phi - beta",
            "expansion": "sqrt5 - tau0 = phi^-1 * u - (sqrt5-2) * u^2 + O(u^3)",
            "coeff_u1_symbolic": str(sp.nsimplify(c1)),
            "coeff_u1_equals_phi_inverse": sp.simplify(c1 - 1 / phi) == 0,
            "coeff_u2_symbolic": str(sp.nsimplify(c2)),
            "coeff_u2_equals_minus_sqrt5_minus_2": sp.simplify(c2 - (-(sp.sqrt(5) - 2))) == 0,
            "coeff_u2_equals_minus_phi_inverse_cubed": sp.simplify(c2 - (-1 / phi ** 3)) == 0,
            "coeff_u3_symbolic": str(sp.nsimplify(c3)),
            "note": "CORRECTED sign (errata 2026-07-04): the quadratic term is "
                    "NEGATIVE, -(sqrt5-2) = -phi^-3, not +(sqrt5-2). The linear "
                    "slope phi^-1 and the curvature MAGNITUDE sqrt5-2 are unchanged.",
        },
    }


def proof_sketch():
    """Thm 6.6 proof: P* = 1-x-x^2, P*(phi) = -2 phi, P'(phi) = 2phi-1 = sqrt5."""
    P = x ** 2 - x - 1
    Pstar = sp.expand(x ** 2 * P.subs(x, 1 / x))
    dP = sp.diff(P, x)
    return {
        "P": "x^2 - x - 1",
        "Pstar": str(Pstar),
        "Pstar_equals_1_minus_x_minus_x2": sp.simplify(Pstar - (1 - x - x ** 2)) == 0,
        "Pstar_at_phi": str(sp.simplify(Pstar.subs(x, phi))),
        "Pstar_at_phi_equals_minus_2phi": sp.simplify(Pstar.subs(x, phi) - (-2 * phi)) == 0,
        "Pprime_at_phi": str(sp.simplify(dP.subs(x, phi))),
        "Pprime_at_phi_equals_sqrt5": sp.simplify(dP.subs(x, phi) - sp.sqrt(5)) == 0,
    }


def geometric_rate_rows():
    """Thm 6.6 table for n = 9,12,...,27: gap, gap*phi^n, consecutive ratio."""
    mp.mp.dps = 60
    PHI = sc.PHI()
    rows = []
    for n in (9, 12, 15, 18, 21, 24, 27):
        g = sc.golden_gap(n)
        g_prev = sc.golden_gap(n - 1)
        rows.append({
            "n": n,
            "gap_sqrt5_minus_tau0": mp.nstr(g, 7),
            "gap_times_phi_n": mp.nstr(g * PHI ** n, 7),
            "consecutive_ratio": mp.nstr(g_prev / g, 7),
        })
    rows.append({
        "n": "limit",
        "gap_sqrt5_minus_tau0": "0",
        "gap_times_phi_n": mp.nstr(2 / sc.SQRT5(), 7),   # 2/sqrt5 = 0.894427
        "consecutive_ratio": mp.nstr(PHI, 7),            # phi = 1.618034
    })
    return rows


def main():
    payload = {
        "_description": "Golden limit and its rate (Section 6). Includes the CORRECTED "
                        "Prop 6.4 quadratic Taylor coefficient -(sqrt5-2) = -phi^-3.",
        "golden_limit_Thm6_2": golden_limit(),
        "linear_rate_and_taylor_Prop6_4": linear_rate_and_taylor(),
        "proof_sketch_Thm6_6": proof_sketch(),
        "geometric_rate_asymptotics_Thm6_6": {
            "gap_times_phi_n_limit": mp.nstr(2 / sc.SQRT5(), 12) + " = 2/sqrt5",
            "consecutive_ratio_limit": mp.nstr(sc.PHI(), 12) + " = phi",
        },
    }
    p_json = io.write_json("prop64_taylor.json", payload, __file__)

    rows = geometric_rate_rows()
    fields = ["n", "gap_sqrt5_minus_tau0", "gap_times_phi_n", "consecutive_ratio"]
    p_csv = io.write_csv("geometric_rate.csv", fields, rows, __file__)
    print("wrote", p_json)
    print("wrote", p_csv)
    for r in rows:
        print(f"  n={str(r['n']):5s} gap={r['gap_sqrt5_minus_tau0']:>12s} "
              f"gap*phi^n={r['gap_times_phi_n']:>10s} ratio={r['consecutive_ratio']}")


if __name__ == "__main__":
    main()
