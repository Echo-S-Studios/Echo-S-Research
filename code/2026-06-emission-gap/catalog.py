"""Producer: the emission catalog and its per-seed invariants.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces: data/2026-06-emission-gap/catalog.csv
Paper results: Def. 2.1 (the seven catalog seeds and their minimal polynomials),
Lemma 4.1 (eigenvalue arguments in (pi/2)Z), Cor. 6.1 / App. A (catalog Mahler
measures {phi, phi, 2, 3, 5, phi^4, beta^2}, minimum phi), and Prop. 9.1 (field
and trace-form signatures; the single Lorentzian generator K).
"""
import mpmath as mp
import sympy as sp

import emgap_core as C

SCRIPT = "catalog.py"

# closed forms for the catalog Mahler measures (Cor. 6.1)
_PHI = (1 + mp.sqrt(5)) / 2
CLOSED_FORM = {
    "phi": "phi", "tau": "phi", "sqrt2": "2", "sqrt3": "3", "sqrt5": "5",
    "gap": "phi^4", "K": "beta^2 = (5+3 sqrt5)/2",
}


def poly_str(coeffs):
    return str(sp.Poly(coeffs, C.x).as_expr())


def compute_rows():
    rows = []
    for name, coeffs in C.CATALOG.items():
        roots = C.mp_roots(coeffs)
        args = sorted({round(C.arg_degrees(r)) % 360 for r in roots if abs(r) > mp.mpf(10) ** (-30)})
        # argument confinement: every argument is a multiple of 90 degrees
        conf = all(min(a % 90, 90 - (a % 90)) < 1 for a in args)
        measure = C.mahler(coeffs)
        sig = C.signature_from_minpoly(coeffs)
        tf = C.trace_form_signature(coeffs)
        # forbidden Salem band is (1, mu_S); mu_S = 1.3247179572...
        mu_S = float(mp.findroot(lambda z: z**3 - z - 1, mp.mpf("1.324")))
        in_band = bool(1 + 1e-9 < float(measure) < mu_S - 1e-9)
        rows.append({
            "seed": name,
            "minpoly": poly_str(coeffs),
            "degree": len(coeffs) - 1,
            "eigenvalue_args_deg": ";".join(str(a) for a in args),
            "args_in_half_pi_Z": conf,
            "mahler_measure": C.s(measure, 20),
            "mahler_closed_form": CLOSED_FORM[name],
            "field_signature_r1_r2": f"({sig[0]},{sig[1]})",
            "trace_form_signature": f"({tf[0]},{tf[1]})",
            "kind": "totally_real" if sig[1] == 0 else "lorentzian",
            "mahler_in_salem_band": in_band,
        })
    return rows


def main():
    rows = compute_rows()
    fields = ["seed", "minpoly", "degree", "eigenvalue_args_deg", "args_in_half_pi_Z",
              "mahler_measure", "mahler_closed_form", "field_signature_r1_r2",
              "trace_form_signature", "kind", "mahler_in_salem_band"]
    path = C.write_csv("catalog.csv", fields, rows, SCRIPT)
    print(f"wrote {path}")
    measures = [float(sp.Float(r["mahler_measure"])) for r in rows]
    print(f"  min catalog Mahler measure = {min(measures)} (phi = {float(_PHI)})")
    print(f"  seeds in Salem band (1, mu_S): {sum(r['mahler_in_salem_band'] for r in rows)}")


if __name__ == "__main__":
    main()
