"""
Producer: the Pisot-trace accumulation (Section 8) and the upstairs
cyclotomic x grow factorization of S_n (Section 9).

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * data/pisot_traces.csv   (Prop 8.1) accumulation points = Pisot traces:
        name, defining_poly, number, redirection tau0 = number + 1/number
        (the golden phi -> sqrt5, the least Pisot-limit accumulation; the
         smallest Pisot / plastic mu_P -> 2.079596)
  * data/sn_factorizations.csv   (Section 9 table) S_n = x^n P - P*, P=x^2-x-1:
        n, S_n, factorization (cyclotomic x Salem), salem_factor,
        salem_mahler (dominant root, climbing 1.5061 -> 1.6054 -> 1.6134 < phi)

Run:  py code/2026-06-salem-slot/make_pisot_factorization.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

x = sc.x
mp = sc.mp


def pisot_trace_rows():
    """Prop 8.1: golden phi -> sqrt5 (least Pisot-limit trace) and the smallest
    Pisot number, the plastic mu_P (root of x^3-x-1), -> 2.079596."""
    mp.mp.dps = 40
    rows = []
    # golden ratio: least limit point of the Pisot set; trace = sqrt5
    phi = sc.PHI()
    rows.append({
        "name": "golden ratio phi (least Pisot limit point)",
        "defining_poly": "x^2 - x - 1",
        "number": mp.nstr(phi, 10),
        "tau0_redirection": mp.nstr(phi + 1 / phi, 10),
        "note": "tau0 = sqrt5 = 2.236068 (least Pisot-limit accumulation)",
    })
    # plastic number: smallest Pisot number
    mu = mp.findroot(lambda z: z ** 3 - z - 1, mp.mpf('1.3'))
    rows.append({
        "name": "plastic number mu_P (smallest Pisot)",
        "defining_poly": "x^3 - x - 1",
        "number": mp.nstr(mu, 10),
        "tau0_redirection": mp.nstr(mu + 1 / mu, 10),
        "note": "tau0 = 2.079596 > 2 (a redirection into grow)",
    })
    return rows


def _Sn(n):
    return sc.Sn_poly(n)


def sn_factor_rows():
    """Section 9 table: S_n factors as cyclotomic x Salem for n=6,10,12."""
    mp.mp.dps = 40
    phi = sc.PHI()
    # the Salem factor of each S_n (the non-cyclotomic reciprocal factor)
    salem_factors = {
        6: x ** 6 - x ** 5 - x ** 3 - x + 1,
        10: x ** 8 - 2 * x ** 7 + x ** 6 - x ** 4 + x ** 2 - 2 * x + 1,
        12: x ** 12 - x ** 11 - x ** 9 - x ** 7 - x ** 5 - x ** 3 - x + 1,
    }
    rows = []
    for n in (6, 10, 12):
        Sn = _Sn(n)
        fac = sp.factor(Sn)
        salem = salem_factors[n]
        roots = [complex(r) for r in sp.nroots(sp.Poly(salem, x), n=40)]
        beta = max(r.real for r in roots if abs(r.imag) < 1e-18)
        rows.append({
            "n": n,
            "S_n": str(sp.expand(Sn)),
            "factorization": str(fac),
            "salem_factor": str(sp.expand(salem)),
            "salem_mahler_dominant_root": mp.nstr(mp.mpf(beta), 8),
            "below_phi": bool(beta < float(phi)),
        })
    return rows


def main():
    p1 = io.write_csv(
        "pisot_traces.csv",
        ["name", "defining_poly", "number", "tau0_redirection", "note"],
        pisot_trace_rows(), __file__,
    )
    p2 = io.write_csv(
        "sn_factorizations.csv",
        ["n", "S_n", "factorization", "salem_factor",
         "salem_mahler_dominant_root", "below_phi"],
        sn_factor_rows(), __file__,
    )
    print("wrote", p1)
    print("wrote", p2)


if __name__ == "__main__":
    main()
