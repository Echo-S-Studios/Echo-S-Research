"""
Producer -- Object II: the content polynomial x^4 - 1.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Proposition 3.1 (prop:content), Remark 3.2 (rem:algebra-not-circle).

Computes and EMITS, to data/2026-06-generative-emptiness/object2_content.json:
  * the cyclotomic factorization  x^4-1 = Phi_1 Phi_2 Phi_4 = (x-1)(x+1)(x^2+1);
  * the four roots of x^4-1 with their Z/4Z charges (angle lattice 0,pi/2,pi,3pi/2);
  * the realised on-circle roots -- phi(x)phi contributes -1, phi^4(x)phi^4 gives +1;
  * K's imaginary place +-i*beta with modulus beta = 2.4195 (off the unit circle);
  * the sweep showing the full Z/4Z sits in the charge while only Z/2Z={+-1} is on
    the circle (why x^4-1, not x^2-1, is the content polynomial).

Run:
    py code/2026-06-generative-emptiness/object2_content.py
"""
import json
import os

import sympy as sp
import mpmath as mp

import ge_core as G
from ge_core import (x, phi_seed, K_seed, tensor, sq, roots_mp, on_circle_roots,
                     charge_of_root, charges, pretty, mpf_str)

OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness", "object2_content.json",
)


def cyclotomic_factorization():
    """x^4-1 = Phi_1 Phi_2 Phi_4 = (x-1)(x+1)(x^2+1)."""
    Phi1 = sp.cyclotomic_poly(1, x)
    Phi2 = sp.cyclotomic_poly(2, x)
    Phi4 = sp.cyclotomic_poly(4, x)
    return {
        "content_polynomial": "x^4 - 1",
        "factors": {"Phi_1": str(Phi1), "Phi_2": str(Phi2), "Phi_4": str(Phi4)},
        "product_equals_x4_minus_1":
            sp.expand((x**4 - 1) - Phi1 * Phi2 * Phi4) == 0,
    }


def content_roots():
    """The four roots of x^4-1 sit at the Z/4Z angle lattice {0,pi/2,pi,3pi/2}."""
    out = []
    for r in roots_mp(sp.Poly(x**4 - 1, x)):
        ang = float((mp.arg(r) % (2 * mp.pi)))
        out.append({
            "root": [mpf_str(r.real, 12), mpf_str(r.imag, 12)],
            "arg_over_halfpi": round(float((mp.arg(r) % (2 * mp.pi)) / (mp.pi / 2))),
            "arg_radians": round(ang, 12),
            "charge": charge_of_root(r),
        })
    return sorted(out, key=lambda d: d["charge"])


def realised_on_circle():
    """phi(x)phi contributes -1 (charge 2); phi^4(x)phi^4 contributes +1 (charge 0)."""
    pp = tensor(phi_seed(), phi_seed())
    phi4 = sq(sq(phi_seed()))
    p44 = tensor(phi4, phi4)
    return [
        {"object": "phi (x) phi", "poly": pretty(pp),
         "root_minus_1": int(pp.eval(-1)) == 0,
         "root_plus_1": int(pp.eval(1)) == 0, "contributes": "-1 (charge 2)"},
        {"object": "phi^4 (x) phi^4", "poly": pretty(p44),
         "root_minus_1": int(p44.eval(-1)) == 0,
         "root_plus_1": int(p44.eval(1)) == 0, "contributes": "+1 (charge 0)"},
    ]


def K_imaginary_place():
    """K's charge-{1,3} place is +-i*beta with beta = sqrt((5+3sqrt5)/2) = 2.4195."""
    imag = [r for r in roots_mp(K_seed())
            if abs(r.real) < mp.mpf(10)**(-20) and abs(r.imag) > 0]
    beta = abs(imag[0])
    beta_exact = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)
    return {
        "n_imaginary_roots": len(imag),
        "beta": mpf_str(beta, 20),
        "beta_closed_form": "sqrt((5 + 3*sqrt5)/2)",
        "beta_rounded_4dp": round(float(beta), 4),
        "matches_closed_form": bool(abs(beta - beta_exact) < mp.mpf(10)**(-30)),
        "off_unit_circle": bool(abs(beta - 1) > mp.mpf("0.1")),
    }


def charge_vs_circle():
    """Full Z/4Z in the charge, only Z/2Z={+-1} realised on the circle."""
    phi4 = sq(sq(phi_seed()))
    objects = [phi_seed(), K_seed(), tensor(phi_seed(), phi_seed()),
               tensor(phi4, phi4), tensor(K_seed(), K_seed())]
    seen = set()
    for P in objects:
        for r in on_circle_roots(P):
            seen.add((round(float(r.real)), round(float(r.imag))))
    return {
        "on_circle_roots_seen": sorted([list(t) for t in seen]),
        "on_circle_is_Z2": seen <= {(1, 0), (-1, 0)},
        "imaginary_units_on_circle": (0, 1) in seen or (0, -1) in seen,
        "charge_group_via_K": sorted(set(charges(K_seed()))),
        "x2_minus_1_has_root_i": sp.Poly(x**2 - 1, x).eval(sp.I) == 0,
        "x4_minus_1_has_root_i": sp.Poly(x**4 - 1, x).eval(sp.I) == 0,
    }


def main():
    data = G.provenance("object2_content.py")
    data.update({
        "object": "II -- the content polynomial x^4 - 1",
        "paper_result": "Proposition 3.1 (prop:content), Remark 3.2",
        "cyclotomic_factorization": cyclotomic_factorization(),
        "content_roots": content_roots(),
        "realised_on_circle": realised_on_circle(),
        "K_imaginary_place": K_imaginary_place(),
        "charge_vs_circle": charge_vs_circle(),
    })
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"wrote {OUT}")
    print("  x^4-1 =", " * ".join(
        f"({v})" for v in data["cyclotomic_factorization"]["factors"].values()))
    print("  charges of x^4-1 roots:", [d["charge"] for d in data["content_roots"]])
    print("  beta =", data["K_imaginary_place"]["beta_rounded_4dp"],
          "off-circle:", data["K_imaginary_place"]["off_unit_circle"])
    print("  on-circle roots seen:", data["charge_vs_circle"]["on_circle_roots_seen"])


if __name__ == "__main__":
    main()
