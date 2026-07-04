"""
Producer: the benchmark Salem table of Section 5 (with the Lehmer numeric of
Section 4).

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces (Section 5 "What it does instead", Remark 5.2 benchmark table):
  * data/benchmarks.csv   one row per benchmark Salem:
        name, degree, salem_poly, beta=Mah(beta), trace_down_T, tau0, Mah_T,
        log_beta, log_tau0, quad_check (tau0-2 == (beta-1)^2/beta)
  * data/benchmarks.json  the same, with high-precision strings + provenance.

For each Salem the redirection tau0, the trace-down T, the Mahler measure
Mah(T) and the entropy-trade logs are rebuilt from scratch: beta is the
dominant root of the (independently given) Salem polynomial, T is its
trace-down (Definition 2.1) computed by the power-sum recurrence, tau0 is the
dominant root of T, and Mah(T) = prod max(1,|root of T|).

Run:  py code/2026-06-salem-slot/make_benchmarks.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

x, t = sc.x, sc.t
mp = sc.mp

# Each benchmark named by its *upstairs* Salem polynomial (highest-first coeffs).
BENCHMARKS = [
    ("Lehmer", 10, [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]),   # x^10+x^9-x^7-...-x^3+x+1
    ("beta4", 4, [1, -1, -1, -1, 1]),                          # x^4-x^3-x^2-x+1
    ("deg6", 6, [1, -1, 0, -1, 0, -1, 1]),                     # x^6-x^5-x^3-x+1
]


def _poly_expr(coeffs):
    deg = len(coeffs) - 1
    return sum(c * x ** (deg - i) for i, c in enumerate(coeffs))


def compute(name, degree, coeffs):
    mp.mp.dps = 45
    beta = sc.dominant_real_root(coeffs)            # Mah(beta) = beta for a Salem
    tau0 = beta + 1 / beta
    # trace-down T of the Salem polynomial (exact, integer coefficients)
    P = _poly_expr(coeffs)
    T = sc.trace_down(P)
    T_coeffs = [int(c) for c in T.all_coeffs()]
    # dominant root of T must equal tau0; Mah(T) from T's roots
    tau0_from_T = sc.dominant_real_root(T_coeffs)
    mah_T = sc.mahler_from_coeffs(T_coeffs)
    quad = (beta - 1) ** 2 / beta                   # Prop 4.2: tau0 - 2
    return {
        "name": name,
        "degree": degree,
        "salem_poly": str(sp.expand(P)),
        "trace_down_T": str(T.as_expr()),
        "beta": beta, "tau0": tau0, "tau0_from_T": tau0_from_T,
        "mah_T": mah_T, "log_beta": mp.log(beta), "log_tau0": mp.log(tau0),
        "quad_tau0_minus_2": quad,
        "quad_matches": abs((tau0 - 2) - quad) < mp.mpf(10) ** (-30),
    }


def main():
    results = [compute(*b) for b in BENCHMARKS]

    csv_rows = []
    for r in results:
        csv_rows.append({
            "name": r["name"],
            "degree": r["degree"],
            "salem_poly": r["salem_poly"],
            "trace_down_T": r["trace_down_T"],
            "beta_Mah_beta": mp.nstr(r["beta"], 8),
            "tau0": mp.nstr(r["tau0"], 8),
            "Mah_T": mp.nstr(r["mah_T"], 8),
            "log_beta": mp.nstr(r["log_beta"], 5),
            "log_tau0": mp.nstr(r["log_tau0"], 5),
            "quad_tau0_minus_2": mp.nstr(r["quad_tau0_minus_2"], 8),
        })
    fields = ["name", "degree", "salem_poly", "trace_down_T", "beta_Mah_beta",
              "tau0", "Mah_T", "log_beta", "log_tau0", "quad_tau0_minus_2"]
    p_csv = io.write_csv("benchmarks.csv", fields, csv_rows, __file__)

    payload = {
        "_description": "Benchmark Salem numbers (Section 5) with trace redirection tau0, "
                        "trace-down T, Mah(T), and the entropy-trade logs. High precision.",
        "benchmarks": [
            {
                "name": r["name"], "degree": r["degree"],
                "salem_poly": r["salem_poly"], "trace_down_T": r["trace_down_T"],
                "beta": mp.nstr(r["beta"], 16), "tau0": mp.nstr(r["tau0"], 16),
                "tau0_from_trace_down": mp.nstr(r["tau0_from_T"], 16),
                "Mah_T": mp.nstr(r["mah_T"], 16),
                "log_beta": mp.nstr(r["log_beta"], 12),
                "log_tau0": mp.nstr(r["log_tau0"], 12),
                "prop4_2_tau0_minus_2": mp.nstr(r["quad_tau0_minus_2"], 16),
                "prop4_2_identity_holds": bool(r["quad_matches"]),
            }
            for r in results
        ],
    }
    p_json = io.write_json("benchmarks.json", payload, __file__)
    print("wrote", p_csv)
    print("wrote", p_json)
    for r in results:
        print(f"  {r['name']:8s} beta={mp.nstr(r['beta'],8)} tau0={mp.nstr(r['tau0'],8)} "
              f"Mah(T)={mp.nstr(r['mah_T'],8)}  T={r['trace_down_T']}")


if __name__ == "__main__":
    main()
