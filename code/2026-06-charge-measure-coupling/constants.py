"""
Producer: named constants and closed-form identities.

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/constants.json

Recomputes every named constant of the paper from its DEFINING polynomial or
identity (never transcribed from the paper's printed decimals), emitting the
exact closed form and a 40-digit value.  Covers:
  * golden ratio phi and conjugate phi'          (Sec 1.3, Sec 4.3, Lem 4.5)
  * sqrt5 = phi + 1/phi = phi - phi'             (Sec 4.3)
  * phi^2 = (3+sqrt5)/2, phi^4 = (7+3sqrt5)/2    (Lem 2.6, Thm 6.5, ledger G/I)
  * pentagon cosines 2cos72=phi-1, 2cos144=-phi  (Prop 6.1)
  * plastic number mu_S = root of x^3-x-1        (Lem 2.5, Smyth floor)
  * Lehmer's number tau                          (Prop 8.2, ledger L)
  * minimal degree-4 Salem beta_4                (ledger F)
  * 2+sqrt3 (reciprocal Z/5 window value)        (Thm 6.4 window, erratum)

Run: py code/2026-06-charge-measure-coupling/constants.py
"""

import mpmath as mp
import sympy as sp

import cmc_core as core
from cmc_io import write_json

mp.mp.dps = 60
DIGITS = 40

s5 = sp.sqrt(5)
PHI = (1 + s5) / 2
PHI_C = (1 - s5) / 2


def _dec(expr, digits: int = DIGITS) -> str:
    """Exact sympy expression -> plain decimal string at `digits` digits."""
    return mp.nstr(mp.mpf(str(sp.N(expr, digits + 10))), digits)


def build_constants():
    out = {}

    # --- golden ratio and conjugate (Sec 1.3, 4.3, Lem 4.5) ---
    out["phi"] = {
        "defining_poly": "x^2 - x - 1",
        "closed_form": "(1+sqrt5)/2",
        "value": _dec(PHI),
        "poly_check_zero": str(sp.simplify(PHI**2 - PHI - 1)),
        "paper_ref": "Sec 1.4 / 4.3",
    }
    out["phi_conjugate"] = {
        "closed_form": "(1-sqrt5)/2 = -1/phi",
        "value": _dec(PHI_C),
        "is_negative_at_argument_pi": bool(PHI_C < 0),
        "phi_times_conj": str(sp.simplify(PHI * PHI_C)),  # = -1
        "paper_ref": "Sec 4.3 / Lem 4.5",
    }
    out["sqrt5"] = {
        "identity": "sqrt5 = phi + 1/phi = phi - phi'",
        "value": _dec(s5),
        "check_phi_plus_inv": str(sp.simplify(PHI + 1 / PHI - s5)),
        "check_phi_minus_conj": str(sp.simplify(PHI - PHI_C - s5)),
        "paper_ref": "Sec 4.3",
    }

    # --- powers of phi (Lem 2.6, Thm 6.5, ledgers G, H, I) ---
    out["phi_squared"] = {
        "closed_form": "(3+sqrt5)/2 = phi + 1",
        "value": _dec(PHI**2),
        "check": str(sp.simplify(PHI**2 - (3 + s5) / 2)),
        "paper_ref": "Lem 2.6 / ledger G",
    }
    out["phi_fourth"] = {
        "closed_form": "(7+3*sqrt5)/2 = 3*phi + 2",
        "value": _dec(PHI**4),
        "check": str(sp.simplify(PHI**4 - (7 + 3 * s5) / 2)),
        "paper_ref": "Thm 6.5 / ledger I",
    }

    # --- pentagon cosines (Prop 6.1) ---
    c72 = 2 * sp.cos(2 * sp.pi / 5)
    c144 = 2 * sp.cos(4 * sp.pi / 5)
    out["two_cos_72"] = {
        "closed_form": "phi - 1",
        "value": _dec(c72),
        "check": str(sp.simplify(c72 - (PHI - 1))),
        "root_of": "x^2 + x - 1",
        "paper_ref": "Prop 6.1",
    }
    out["two_cos_144"] = {
        "closed_form": "-phi",
        "value": _dec(c144),
        "check": str(sp.simplify(c144 - (-PHI))),
        "root_of": "x^2 + x - 1",
        "galois_conjugate_of": "two_cos_72 under sqrt5 -> -sqrt5",
        "paper_ref": "Prop 6.1",
    }

    # --- plastic number mu_S, the Smyth floor (Lem 2.5) ---
    mu = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    ps_roots = mp.polyroots([1, 0, -1, -1], extraprec=200)
    complex_conj_moduli = sorted(
        float(abs(r)) for r in ps_roots if abs(mp.im(r)) >= mp.mpf(10) ** (-30)
    )
    out["mu_S_plastic_number"] = {
        "defining_poly": "x^3 - x - 1",
        "value": mp.nstr(mu, DIGITS),
        "residual_at_root": mp.nstr(mu**3 - mu - 1, 3),
        "is_pisot_conjugate_moduli_below_1": all(m < 1 for m in complex_conj_moduli),
        "description": "smallest Pisot number; forced non-reciprocal Mahler floor",
        "paper_ref": "Lem 2.5 (Smyth)",
    }

    # --- Lehmer's number tau (Prop 8.2, ledger L) ---
    lehmer = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    tau = core.mahler(lehmer)
    out["lehmer_tau"] = {
        "defining_poly": "x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1",
        "value": mp.nstr(tau, DIGITS),
        "in_gap_1_phi": bool(1 < tau < mp.mpf(str(sp.N(PHI, 45)))),
        "description": "smallest known Salem number; below the emission floor phi",
        "paper_ref": "Prop 8.2 / ledger L",
    }

    # --- minimal degree-4 Salem number beta_4 (ledger F) ---
    beta4 = core.mahler([1, -1, -1, -1, 1])
    out["beta_4"] = {
        "defining_poly": "x^4 - x^3 - x^2 - x + 1",
        "value": mp.nstr(beta4, DIGITS),
        "in_gap_phi_2": bool(mp.mpf(str(sp.N(PHI, 45))) < beta4 < 2),
        "description": "minimal degree-4 Salem number; sits in the odd-charge gap (phi,2)",
        "paper_ref": "ledger F",
    }

    # --- 2+sqrt3, the extra reciprocal Z/5 window value (erratum) ---
    two_plus_sqrt3 = 2 + sp.sqrt(3)
    out["two_plus_sqrt3"] = {
        "closed_form": "2 + sqrt3",
        "value": _dec(two_plus_sqrt3),
        "arises_as": "Mahler measure of Phi_5 * (x^2 - 4x + 1); reciprocal Z/5 window",
        "paper_ref": "Thm 6.4 window {1, phi^2, 2+sqrt3}",
    }

    return out


def main():
    payload = {"constants": build_constants()}
    path = write_json("constants.json", payload, __file__)
    print(f"wrote {path}")
    print(f"  {len(payload['constants'])} constants recomputed from defining data")


if __name__ == "__main__":
    main()
