"""
Producer -- Object IV: the graded normal form.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Proposition 5.1 (prop:normal) and its three-row table.

Rebuilds the three diagnostic objects from the operators, factors each over Z,
classifies every factor by sector (on-circle | x^4-1, imaginary axis, off-circle
grow), and EMITS to data/2026-06-generative-emptiness/object4_normalform.json:
  * phi(x)phi     = (x+1)^2 (x^2-3x+1),                  grow M = phi^2   = 2.618;
  * phi^4(x)phi^4 = (x-1)^2 (x^2-47x+1),                 grow M = phi^8   = 46.98;
  * K(x)K = (x^2+5)^4 (x^2-5x-5)^2 (x^2+5x-5)^2,          grow M = 5.854;
  * for each: on-circle part divides x^4-1, growth lives off-circle, and M(G)=M(P);
  * the identity a*b = sqrt5 that produces the imaginary sector (x^2+5).

Run:
    py code/2026-06-generative-emptiness/object4_normalform.py
"""
import json
import os

import sympy as sp
import mpmath as mp

import ge_core as G
from ge_core import (x, phi_seed, K_seed, tensor, sq, mahler, roots_mp, pretty,
                     PHI, mpf_str)

OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness", "object4_normalform.json",
)

PHI_V = PHI()
X4M1 = sp.Poly(x**4 - 1, x)


def _sector(fp):
    """Classify an irreducible factor by sector."""
    rts = roots_mp(fp)
    on_circle = all(abs(abs(r) - 1) < mp.mpf(10)**(-15) for r in rts)
    imaginary = all(abs(r.real) < mp.mpf(10)**(-15) for r in rts)
    divides_x4m1 = sp.rem(X4M1.as_expr(), fp.as_expr(), x) == 0
    if on_circle and divides_x4m1:
        return "on-circle (divides x^4-1)"
    if imaginary:
        return "imaginary axis (charges 1,3)"
    return "off-circle grow"


def factorization(P):
    """Graded factorization of P: per-factor poly, multiplicity, measure, sector."""
    _, facs = sp.factor_list(P.as_expr(), x)
    out = []
    for f, m in facs:
        fp = sp.Poly(sp.expand(f), x)
        out.append({
            "factor": pretty(fp),
            "multiplicity": m,
            "mahler_measure": mpf_str(mahler(fp), 20),
            "sector": _sector(fp),
        })
    return out


def diagnostic(name, P, grow_expr, grow_closed_name, grow_closed_val):
    facs = factorization(P)
    grow = sp.Poly(grow_expr, x)
    mP, mG = mahler(P), mahler(grow)

    # on-circle part C(x): factors dividing x^4-1 (must have Mahler measure 1)
    on_circle_ok = all(
        (sp.rem(X4M1.as_expr(), sp.Poly(f["factor"], x).as_expr(), x) == 0)
        for f in facs if f["sector"].startswith("on-circle")
    )
    # off-circle part G(x): ALL factors not dividing x^4-1 (Prop 5.1 M(G)=M(P)).
    # For the two phi objects G is a single real grow factor; for K(x)K the
    # measure is distributed over the imaginary sector and two real grow factors.
    off_total = mp.mpf(1)
    for f in facs:
        if not f["sector"].startswith("on-circle"):
            off_total *= mahler(sp.Poly(f["factor"], x)) ** f["multiplicity"]

    return {
        "object": name,
        "poly": pretty(P),
        "factorization": facs,
        "grow_factor": pretty(grow),
        "grow_measure": mpf_str(mG, 20),
        "grow_measure_closed_form": grow_closed_name,
        "grow_measure_rounded": round(float(mG), 3),
        "grow_matches_closed_form":
            bool(abs(mG - grow_closed_val) < mp.mpf(10)**(-20)),
        "object_measure": mpf_str(mP, 20),
        "off_circle_total_measure": mpf_str(off_total, 20),
        "off_circle_measure_equals_object":
            bool(abs(mP - off_total) < mp.mpf(10)**(-15)),
        "on_circle_divides_x4_minus_1": on_circle_ok,
    }


def ab_equals_sqrt5():
    """The K-sector (x^2+5) has roots +-i*sqrt5: a^2 b^2 = 5 with a,b K's moduli."""
    r5 = sp.sqrt(5)
    a2 = (-5 + 3 * r5) / 2      # real half-modulus squared
    b2 = (5 + 3 * r5) / 2       # imaginary half-modulus squared
    return {
        "a_squared": "(-5 + 3*sqrt5)/2",
        "b_squared": "(5 + 3*sqrt5)/2",
        "a2_times_b2_minus_5": int(sp.simplify(a2 * b2 - 5)),
        "identity": "a * b = sqrt5  ->  imaginary sector x^2 + 5 (roots +-i*sqrt5)",
        "x2_plus_5_has_root_i_sqrt5":
            sp.simplify(sp.Poly(x**2 + 5, x).eval(sp.I * r5)) == 0,
    }


def main():
    phi = phi_seed()
    phi4 = sq(sq(phi))
    rows = [
        diagnostic("phi (x) phi", tensor(phi, phi),
                   x**2 - 3 * x + 1, "phi^2", PHI_V**2),
        diagnostic("phi^4 (x) phi^4", tensor(phi4, phi4),
                   x**2 - 47 * x + 1, "phi^8", PHI_V**8),
        diagnostic("K (x) K", tensor(K_seed(), K_seed()),
                   x**2 - 5 * x - 5, "(5+3*sqrt5)/2", (5 + 3 * mp.sqrt(5)) / 2),
    ]
    data = G.provenance("object4_normalform.py")
    data.update({
        "object": "IV -- the graded normal form",
        "paper_result": "Proposition 5.1 (prop:normal)",
        "content_polynomial": "x^4 - 1 = (x-1)(x+1)(x^2+1)",
        "diagnostics": rows,
        "imaginary_sector_identity": ab_equals_sqrt5(),
    })
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"wrote {OUT}")
    for r in rows:
        parts = " * ".join(
            f"({fx['factor']})^{fx['multiplicity']}" if fx["multiplicity"] > 1
            else f"({fx['factor']})" for fx in r["factorization"])
        print(f"  {r['object']:16s} = {parts}")
        print(f"      grow M = {r['grow_measure_rounded']} "
              f"({r['grow_measure_closed_form']}), off-circle M(G)=M(P): "
              f"{r['off_circle_measure_equals_object']}")


if __name__ == "__main__":
    main()
