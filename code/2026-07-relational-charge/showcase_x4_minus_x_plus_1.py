r"""
Producer: the fully-rigid quartic x^4 - x + 1 (Example 6.20, ledger X).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/x4_minus_x_plus_1_showcase.json

Refactors Example 6.20 (ex:x4quartic): p = x^4-x+1 is irreducible,
non-reciprocal, with no real roots, and no uniquely attained modulus (so the
pinning theorem is silent).  Three exact routes settle the full type:
  * Structure: charpoly(C_p (x) C_p) = S_6(x)^2 * (x^4+2x^2-x+1), with S_6 the
    certified Salem sextic and the second factor the irreducible psi^2 image;
  * Shell scan: Rat_p contacts {Phi_1^4} (no within-shell torsion);
  * The nu-criterion verbatim: Rat_{Rat_p} has degree 256 and complete
    signature {Phi_1^28} (multiplicity 16+4+8 = 28), certifying full inertness.
Two by-products: Gal(p) = S_4, and M(x^4-x+1) = tau_{S_6} (smallest deg-6 Salem).

Run: py code/2026-07-relational-charge/showcase_x4_minus_x_plus_1.py
"""

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_json

x = C.x


def main():
    p = C.X4MX1
    P = sp.Poly(p, x)

    # basic structure
    irreducible = P.is_irreducible
    reciprocal = C.is_reciprocal(p)
    real_roots = P.count_roots(-sp.oo, sp.oo)

    # Kronecker square charpoly factorization
    Cm = C.companion_matrix(p)
    K = sp.Matrix(sp.kronecker_product(Cm, Cm))
    F = sp.factor(K.charpoly(x).as_expr())
    expected = (x**6 - x**4 - x**3 - x**2 + 1) ** 2 * (x**4 + 2 * x**2 - x + 1)
    kron_ok = sp.expand(F - expected) == 0
    psi2_irreducible = sp.Poly(x**4 + 2 * x**2 - x + 1, x).is_irreducible

    # shell scan Rat_p
    Rp = C.ratio_poly(p)
    sig_shell = C.cyclotomic_contacts(Rp)

    # nested nu-scan Rat_{Rat_p}
    RRp = C.ratio_poly(Rp.as_expr())
    sig_nested = C.cyclotomic_contacts(RRp)
    phi1_nested = C.phi1_multiplicity(RRp)

    # Galois group
    g, _ = P.galois_group()
    gal_order = g.order()
    gal_transitive = g.is_transitive()

    # Mahler = tau_{S6}
    Mq = C.mahler_measure(p)
    tau_S6 = max(r.real for r in C.roots_mp(C.S6) if abs(r.imag) < mp.mpf(10) ** -20)
    theta0 = max(r.real for r in C.roots_mp(C.PLASTIC) if abs(r.imag) < mp.mpf(10) ** -20)
    phi = (1 + mp.sqrt(5)) / 2

    payload = {
        "object": "x^4 - x + 1",
        "basic_structure": {
            "irreducible": bool(irreducible),
            "reciprocal": bool(reciprocal),
            "real_roots": int(real_roots),
            "pinning_silent": "no uniquely attained modulus (two conjugate-pair shells r != 1/r)",
        },
        "kronecker_square_factorization": {
            "charpoly_C_x_C": "S_6(x)^2 * (x^4+2x^2-x+1)",
            "matches_expected": bool(kron_ok),
            "psi2_factor_irreducible": bool(psi2_irreducible),
        },
        "shell_scan_Rat_p": {
            "rat_degree": Rp.degree(),
            "contact_signature": C.signature_str(sig_shell),
            "inert_shell": sig_shell == {1: 4},
        },
        "nested_nu_scan_Rat_Rat_p": {
            "nested_degree": RRp.degree(),
            "scan_bound_2d2": 2 * RRp.degree() ** 2,
            "contact_signature": C.signature_str(sig_nested),
            "phi1_multiplicity": phi1_nested,
            "predicted_multiplicity_16_4_8": 16 + 4 + 8,
            "matches_prediction": sig_nested == {1: 28},
        },
        "galois_group": {
            "order": int(gal_order),
            "transitive": bool(gal_transitive),
            "is_S4": int(gal_order) == 24 and bool(gal_transitive),
        },
        "mahler_measure": {
            "value": mp.nstr(Mq, 20),
            "equals_tau_S6": bool(abs(Mq - tau_S6) < mp.mpf(10) ** -20),
            "tau_S6": mp.nstr(tau_S6, 20),
            "strictly_between_theta0_and_phi": bool(theta0 < Mq < phi),
            "note": "smallest degree-6 Salem number, in (theta0, phi)",
        },
        "verdict": "fully relationally inert (rigidity imported across (x) via S_6)",
        "status": "[forced] per instance (ledger X)",
    }
    path = write_json("x4_minus_x_plus_1_showcase.json", payload, __file__)
    print(f"wrote {path}")
    print(f"  x^4-x+1: Kronecker={kron_ok}, shell {C.signature_str(sig_shell)}, "
          f"nested deg {RRp.degree()} {C.signature_str(sig_nested)}, "
          f"Gal order {gal_order}")


if __name__ == "__main__":
    main()
