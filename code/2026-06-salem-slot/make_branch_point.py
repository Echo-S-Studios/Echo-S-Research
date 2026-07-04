"""
Producer: the branch-point structure at the flip (Section 4).

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * data/branch_point.json
      - square-root edge series of beta(t)=(t+sqrt(t^2-4))/2 in s=sqrt(t-2):
        beta = 1 + s + s^2/2 + s^3/8 + O(s^4)                          (Lem 4.1)
        and the leading asymptotic beta-1 ~ sqrt(t-2), log beta ~ beta-1
      - Prop 4.2 exact identity  tau0 - 2 = (beta-1)^2/beta, with the
        near-floor expansion tau0 = 2 + delta^2 + O(delta^3) and the
        Lehmer numeric (beta-1=0.176281 -> tau0-2=0.026418)            (Prop 4.2)
      - monodromy round t=2 swaps beta <-> 1/beta (product of the two
        branches of x^2 - t x + 1 is 1); numeric at t=2.3               (Sec entry)

Run:  py code/2026-06-salem-slot/make_branch_point.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

mp = sc.mp


def sqrt_edge_series():
    """Lem 4.1: expand beta(t) in s = sqrt(t-2)."""
    s = sp.symbols('s', positive=True)
    t_of_s = 2 + s ** 2
    beta = (t_of_s + sp.sqrt(t_of_s ** 2 - 4)) / 2
    ser = sp.series(beta, s, 0, 4).removeO()
    coeffs = {k: sp.nsimplify(ser.coeff(s, k)) for k in range(4)}
    return {
        "variable": "s = sqrt(t-2)",
        "series": "beta = 1 + s + s^2/2 + s^3/8 + O(s^4)",
        "coeff_s0": str(coeffs[0]), "coeff_s1": str(coeffs[1]),
        "coeff_s2": str(coeffs[2]), "coeff_s3": str(coeffs[3]),
        "matches_paper": (coeffs[0] == 1 and coeffs[1] == 1
                          and coeffs[2] == sp.Rational(1, 2)
                          and coeffs[3] == sp.Rational(1, 8)),
    }


def sqrt_edge_asymptotic():
    """Lem 4.1: beta-1 ~ sqrt(t-2) and log beta ~ beta-1 as t -> 2+."""
    mp.mp.dps = 40
    beta = lambda tt: (tt + mp.sqrt(tt ** 2 - 4)) / 2
    rows = []
    for eps in (mp.mpf('1e-4'), mp.mpf('1e-6'), mp.mpf('1e-8')):
        tt = 2 + eps
        b = beta(tt)
        rows.append({
            "t_minus_2": mp.nstr(eps, 2),
            "beta_minus_1_over_sqrt": mp.nstr((b - 1) / mp.sqrt(eps), 12),
            "logbeta_over_beta_minus_1": mp.nstr(mp.log(b) / (b - 1), 12),
        })
    return {"claim": "(beta-1)/sqrt(t-2) -> 1 and log(beta)/(beta-1) -> 1 as t->2+",
            "samples": rows}


def quadratic_redirection():
    """Prop 4.2: tau0 - 2 = (beta-1)^2/beta exactly; near floor 2+delta^2."""
    b = sp.symbols('b', positive=True)
    d = sp.symbols('d', positive=True)
    exact = sp.simplify(((b + 1 / b) - 2) - (b - 1) ** 2 / b) == 0
    tau0_delta = (1 + d) + 1 / (1 + d)
    ser = sp.series(tau0_delta, d, 0, 3).removeO()
    # Lehmer numeric
    mp.mp.dps = 45
    beta_lehmer = sc.dominant_real_root([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
    tau0 = beta_lehmer + 1 / beta_lehmer
    return {
        "exact_identity": "tau0 - 2 = (beta-1)^2 / beta",
        "exact_holds": exact,
        "near_floor_expansion": "tau0 = 2 + delta^2/(1+delta) = 2 + delta^2 + O(delta^3)",
        "series_coeff_delta0": str(ser.coeff(d, 0)),
        "series_coeff_delta1": str(ser.coeff(d, 1)),
        "series_coeff_delta2": str(ser.coeff(d, 2)),
        "lehmer_beta_minus_1": mp.nstr(beta_lehmer - 1, 8),
        "lehmer_tau0_minus_2": mp.nstr(tau0 - 2, 8),
        "lehmer_quad_check": mp.nstr((beta_lehmer - 1) ** 2 / beta_lehmer, 8),
    }


def monodromy():
    """Sec entry: encircling the branch point t=2 swaps beta <-> 1/beta."""
    mp.mp.dps = 40
    tt = mp.mpf('2.3')
    beta_plus = (tt + mp.sqrt(tt ** 2 - 4)) / 2
    beta_minus = (tt - mp.sqrt(tt ** 2 - 4)) / 2
    return {
        "claim": "analytic continuation of beta(t) once around t=2 sends beta -> 1/beta; "
                 "the two branches are the two roots of x^2 - t x + 1 (product = 1)",
        "t": "2.3",
        "beta_plus": mp.nstr(beta_plus, 10),
        "beta_minus": mp.nstr(beta_minus, 10),
        "product": mp.nstr(beta_plus * beta_minus, 6),
        "beta_minus_is_inverse": abs(beta_minus - 1 / beta_plus) < mp.mpf(10) ** (-30),
    }


def main():
    payload = {
        "_description": "The flip t=2 as a square-root branch point of the trace cover "
                        "(Section 4): sqrt-edge series (Lem 4.1), quadratic redirection "
                        "(Prop 4.2), and the monodromy swap (Section 'entry').",
        "sqrt_edge_series_Lem4_1": sqrt_edge_series(),
        "sqrt_edge_asymptotic_Lem4_1": sqrt_edge_asymptotic(),
        "quadratic_redirection_Prop4_2": quadratic_redirection(),
        "monodromy_swap": monodromy(),
    }
    path = io.write_json("branch_point.json", payload, __file__)
    print("wrote", path)


if __name__ == "__main__":
    main()
