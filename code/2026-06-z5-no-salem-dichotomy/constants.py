"""
Producer: named constants and closed-form identities of the Z/5Z paper.

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/constants.json

Recomputes every named constant from its DEFINING polynomial / identity (never
transcribed from a printed decimal), emitting the exact closed form and a
40-digit value.  Covers:
  * mu_S (plastic number) = real root of x^3 - x - 1, Pisot        (Lem. 2.4, Table 1)
  * realification value 2^(1/5) = 1.1487                            (Rem. 6.2(1))
  * phi = (1+sqrt5)/2, phi^2 - phi - 1 = 0                          (Prop. 3.1)
  * phi^2 = (3+sqrt5)/2 = phi+1 ~ 2.618                             (Lem. 2.5, Thm. 3.1b)
  * phi^4 = (7+3sqrt5)/2 = 3phi+2 ~ 6.854102                        (Thm. 4.1, Sec. 7)
  * 2cos72 = phi-1, 2cos144 = -phi (roots of x^2+x-1)               (Prop. 3.1)
  * 2+sqrt3 (reciprocal Z/5 window value)                          (Thm. 3.1b, Sec. 7)
  * the forced-floor ordering 2^(1/5) < mu_S < 2                    (Table 1, Thm. 3.1)

Run: py code/2026-06-z5-no-salem-dichotomy/constants.py
"""

import mpmath as mp
import sympy as sp

from z5_core import phi, sqrt5, sigma
from z5_io import write_json

mp.mp.dps = 60
DIGITS = 40


def _dec(expr, digits: int = DIGITS) -> str:
    """Exact sympy expression -> plain decimal string at `digits` digits."""
    return mp.nstr(mp.mpf(str(sp.N(expr, digits + 10))), digits)


def build_constants():
    out = {}

    # --- golden ratio (Prop. 3.1) ---
    out["phi"] = {
        "defining_poly": "x^2 - x - 1",
        "closed_form": "(1+sqrt5)/2",
        "value": _dec(phi),
        "poly_check_zero": str(sp.simplify(phi**2 - phi - 1)),
        "paper_ref": "Prop. 3.1",
    }
    out["phi_squared"] = {
        "closed_form": "(3+sqrt5)/2 = phi+1",
        "value": _dec(phi**2),
        "check_zero": str(sp.simplify(phi**2 - (3 + sqrt5) / 2)),
        "role": "forced Mahler floor of a real reciprocal unit (Lem. 2.5)",
        "paper_ref": "Lem. 2.5 / Thm. 3.1(b)",
    }
    out["phi_fourth"] = {
        "closed_form": "(7+3*sqrt5)/2 = 3*phi+2",
        "value": _dec(phi**4),
        "check_zero": str(sp.simplify(phi**4 - (7 + 3 * sqrt5) / 2)),
        "check_3phi_plus_2": str(sp.simplify(phi**4 - (3 * phi + 2))),
        "role": "degree-4 pure-pentagon Mahler floor (Thm. 4.1); minimizer measure",
        "paper_ref": "Thm. 4.1 / Sec. 7",
    }

    # --- plastic number mu_S, the Smyth floor (Lem. 2.4) ---
    mu = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    ps_roots = mp.polyroots([1, 0, -1, -1], extraprec=200)
    conj_moduli = sorted(
        float(abs(r)) for r in ps_roots if abs(mp.im(r)) >= mp.mpf(10) ** (-30)
    )
    out["mu_S_plastic_number"] = {
        "defining_poly": "x^3 - x - 1",
        "value": mp.nstr(mu, DIGITS),
        "residual_at_root": mp.nstr(mu**3 - mu - 1, 3),
        "is_pisot_conjugate_moduli_below_1": bool(all(m < 1 for m in conj_moduli)),
        "conjugate_moduli": [mp.nstr(mp.mpf(m), 12) for m in conj_moduli],
        "description": "smallest Pisot number; forced non-reciprocal Mahler floor (Smyth)",
        "paper_ref": "Lem. 2.4 (Smyth)",
    }

    # --- realification value 2^(1/5) (Rem. 6.2(1)) ---
    two15 = mp.root(2, 5)
    out["realification_bound_2_pow_1_5"] = {
        "closed_form": "2^(1/5)",
        "value": mp.nstr(two15, DIGITS),
        "fifth_power": mp.nstr(two15**5, 5),
        "role": "prior forced floor via psi^5, improved to mu_S in this paper",
        "paper_ref": "Rem. 6.2(1)",
    }

    # --- pentagon cosines (Prop. 3.1) ---
    c72 = 2 * sp.cos(sp.rad(72))
    c144 = 2 * sp.cos(sp.rad(144))
    out["two_cos_72"] = {
        "closed_form": "phi - 1 = (sqrt5-1)/2",
        "value": _dec(c72),
        "check_zero": str(sp.simplify(c72 - (phi - 1))),
        "root_of": "x^2 + x - 1",
        "paper_ref": "Prop. 3.1",
    }
    out["two_cos_144"] = {
        "closed_form": "-phi = -(1+sqrt5)/2",
        "value": _dec(c144),
        "check_zero": str(sp.simplify(c144 - (-phi))),
        "root_of": "x^2 + x - 1",
        "galois_conjugate_of": "two_cos_72 under sqrt5 -> -sqrt5",
        "paper_ref": "Prop. 3.1",
    }

    # --- 2+sqrt3, extra reciprocal Z/5 window value (Sec. 7 corrected set) ---
    tps3 = 2 + sp.sqrt(3)
    out["two_plus_sqrt3"] = {
        "closed_form": "2 + sqrt3",
        "value": _dec(tps3),
        "arises_as": "M of Phi_5*(x^2-4x+1); trace-4 real reciprocal unit",
        "paper_ref": "Thm. 3.1(b) / Sec. 7 window {1, phi^2, 2+sqrt3}",
    }

    return out


def build_orderings():
    """The forced-floor improvement ordering 2^(1/5) < mu_S < 2 (Table 1)."""
    two15 = mp.root(2, 5)
    mu = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    return {
        "sequence": "2^(1/5) < mu_S < 2",
        "two_pow_1_5": mp.nstr(two15, 20),
        "mu_S": mp.nstr(mu, 20),
        "realized_floor": 2,
        "ordering_holds": bool(two15 < mu < 2),
        "improvement_mu_S_minus_2_pow_1_5": mp.nstr(mu - two15, 12),
        "paper_ref": "Table 1 / Thm. 3.1 (floor improved 2^(1/5) -> mu_S)",
    }


def main():
    payload = {
        "constants": build_constants(),
        "forced_floor_ordering": build_orderings(),
    }
    path = write_json("constants.json", payload, __file__)
    print(f"wrote {path}")
    print(f"  {len(payload['constants'])} constants recomputed from defining data")


if __name__ == "__main__":
    main()
