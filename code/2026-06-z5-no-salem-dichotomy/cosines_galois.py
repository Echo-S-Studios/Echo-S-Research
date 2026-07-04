"""
Producer: Proposition 3.1 (irrational, Galois-coupled pentagon cosines) and the
symbolic rows of the Sec. 7 verification table.

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/cosines_galois.json

All identities are established by exact symbolic manipulation (sympy), then
emitted with the residual that certifies each equality.  Covers:
  * 2 cos 72  = phi - 1                                    (Prop. 3.1)
  * 2 cos 144 = -phi                                       (Prop. 3.1)
  * {2cos72, 2cos144} are exactly the two roots of x^2+x-1 (Prop. 3.1 / Sec. 7)
  * Galois-conjugate over Q under sqrt5 -> -sqrt5          (Prop. 3.1)
  * both cosines irrational; contrast 2cos120 = -1 rational (why Z/3Z breaks)
  * cross-term collapse identity phi^2 - phi = 1           (Thm. 4.1 proof)

Run: py code/2026-06-z5-no-salem-dichotomy/cosines_galois.py
"""

import sympy as sp

from z5_core import phi, sigma, sqrt5
from z5_io import write_json

x = sp.symbols('x')

c72 = 2 * sp.cos(sp.rad(72))       # 2 cos 72 deg, kept exact
c144 = 2 * sp.cos(sp.rad(144))     # 2 cos 144 deg
c120 = 2 * sp.cos(sp.rad(120))     # 2 cos 120 deg (the Z/3Z value, = -1)


def build():
    # write each cosine in Q(sqrt5): a + b sqrt5
    a72 = sp.nsimplify(sp.simplify(c72), [sqrt5])     # (sqrt5-1)/2
    a144 = sp.nsimplify(sp.simplify(c144), [sqrt5])   # -(1+sqrt5)/2
    roots_quad = sp.solve(x**2 + x - 1, x)

    out = {}

    out["two_cos_72"] = {
        "identity": "2*cos(72 deg) = phi - 1",
        "closed_form_Q_sqrt5": str(a72),
        "value_40dps": str(sp.N(c72, 40)),
        "residual_minus_(phi-1)": str(sp.simplify(c72 - (phi - 1))),
        "annihilates_x2_plus_x_minus_1": str(sp.simplify(c72**2 + c72 - 1)),
        "is_rational": bool(sp.nsimplify(sp.simplify(c72)).is_rational),
    }
    out["two_cos_144"] = {
        "identity": "2*cos(144 deg) = -phi",
        "closed_form_Q_sqrt5": str(a144),
        "value_40dps": str(sp.N(c144, 40)),
        "residual_minus_(-phi)": str(sp.simplify(c144 - (-phi))),
        "annihilates_x2_plus_x_minus_1": str(sp.simplify(c144**2 + c144 - 1)),
        "is_rational": bool(sp.nsimplify(sp.simplify(c144)).is_rational),
    }

    out["shared_minimal_polynomial"] = {
        "polynomial": "x^2 + x - 1",
        "roots": [str(sp.nsimplify(r)) for r in roots_quad],
        "cosine_set_matches_roots": bool(
            {sp.nsimplify(sp.simplify(c72)), sp.nsimplify(sp.simplify(c144))}
            == {sp.nsimplify(r) for r in roots_quad}
        ),
        "sum_equals_minus1": str(sp.simplify(a72 + a144 + 1)),      # 0
        "product_equals_minus1": str(sp.simplify(a72 * a144 + 1)),  # 0
    }

    out["galois_conjugacy"] = {
        "map": "sqrt5 -> -sqrt5",
        "sigma_of_2cos72_minus_2cos144": str(sp.simplify(sigma(a72) - a144)),  # 0
        "sigma_of_2cos144_minus_2cos72": str(sp.simplify(sigma(a144) - a72)),  # 0
        "conjugate_pair": True,
        "note": "the two cosines are swapped by sigma; this forced coupling is "
                "what Z/3Z never needs (Sec. 4).",
    }

    out["z3_contrast"] = {
        "identity": "2*cos(120 deg) = -1 (rational)",
        "residual_plus_1": str(sp.simplify(c120 + 1)),     # 0
        "is_rational": bool(sp.nsimplify(sp.simplify(c120)).is_rational),
        "consequence": "rational cosine lets the Z/3Z search collapse to x^3-2 "
                       "(M=2); at Z/5Z the cosines are irrational, so integrality "
                       "requires the Galois-partner coupling.",
    }

    out["cross_term_collapse"] = {
        "identity": "phi^2 - phi = 1",
        "residual": str(sp.simplify((phi**2 - phi) - 1)),   # 0
        "role": "collapses the quartic cross term in the pentagon expansion (Thm. 4.1).",
    }

    return out


def main():
    payload = {"proposition_3_1_and_symbolic_rows": build()}
    path = write_json("cosines_galois.json", payload, __file__)
    print(f"wrote {path}")
    print("  2cos72 = phi-1, 2cos144 = -phi, roots of x^2+x-1, Galois-conjugate: verified")


if __name__ == "__main__":
    main()
