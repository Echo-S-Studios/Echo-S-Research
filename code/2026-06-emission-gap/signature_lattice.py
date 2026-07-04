"""Producer: field signatures, the trace form, and the signature lattice of K.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces: data/2026-06-emission-gap/field_signatures.csv
Paper results: a Salem field of degree 2m has signature (2, m-1) and trace-form
signature (m+1, m-1) (Prop. 9.1); the catalog K-formation is the single
Lorentzian generator, (2,1) / (3,1), with its complex place off the circle; the
signature lattice of K = Q(sqrt2, sqrt3, 5^{1/4}) has every non-totally-real
subfield of shape (2k, k), so the only Salem-shaped subfield is Q(5^{1/4}) at
degree 4 (Thm. 10.3); cyclotomic fields Q(zeta_n) are totally complex (0, phi(n)/2),
disjoint from the Salem signature (Prop. 12.3).
"""
import mpmath as mp
import sympy as sp

import emgap_core as C

SCRIPT = "signature_lattice.py"
x = C.x


def classify(r1, r2):
    if r2 == 0:
        return "totally_real"
    if r1 == 0:
        return "totally_complex_CM"
    if r1 == 2 * r2:
        return "signature_(2k,k)"
    if r1 == 2 and r2 >= 1:
        return "salem_shape_(2,m-1)"
    return "mixed"


def minpoly_rows():
    """Signatures from explicit minimal polynomials (Salem fields + catalog seeds)."""
    entries = [
        ("beta_4 (min deg-4 Salem)", [1, -1, -1, -1, 1], "Prop 9.1 / Cor 10.4"),
        ("Lehmer (deg-10 Salem)", C.LEHMER, "Prop 9.1"),
        ("K-formation x^4+5x^2-5", C.CATALOG["K"], "Prop 9.1 (Lorentzian generator)"),
        ("2K difference x^4+20x^2-80", [1, 0, 20, 0, -80], "Lemma 10.2 (real, not Salem)"),
        ("sqrt2 seed x^2-2", C.CATALOG["sqrt2"], "Prop 9.1 (definite)"),
        ("sqrt5 seed x^2-5", C.CATALOG["sqrt5"], "Prop 9.1 (definite)"),
        ("phi seed x^2-x-1", C.CATALOG["phi"], "Prop 9.1 (definite)"),
    ]
    rows = []
    for label, coeffs, loc in entries:
        r1, r2 = C.signature_from_minpoly(coeffs)
        tf = C.trace_form_signature(coeffs)
        m = (len(coeffs) - 1) // 2
        rows.append({
            "label": label, "source": "minpoly",
            "generator_or_minpoly": str(sp.Poly(coeffs, x).as_expr()),
            "degree": len(coeffs) - 1,
            "signature_r1_r2": f"({r1},{r2})",
            "trace_form_signature": f"({tf[0]},{tf[1]})",
            "trace_form_predicted_m+1_m-1": f"({m+1},{m-1})" if r2 else "n/a",
            "kind": classify(r1, r2),
            "paper_location": loc,
        })
    return rows


def subfield_rows():
    """Signature lattice of K via primitive elements (Thm. 10.3)."""
    gens = [
        ("Q(5^{1/4})", sp.root(5, 4), "Thm 10.3 (only Salem-shaped subfield)"),
        ("Q(sqrt2,5^{1/4})", sp.sqrt(2) + sp.root(5, 4), "Thm 10.3"),
        ("Q(sqrt2,sqrt3,sqrt5)", sp.sqrt(2) + sp.sqrt(3) + sp.sqrt(5), "Thm 10.3 (totally real)"),
        ("K=Q(sqrt2,sqrt3,5^{1/4})", sp.sqrt(2) + sp.sqrt(3) + sp.root(5, 4), "Thm 10.3 (full field)"),
    ]
    rows = []
    for label, gen, loc in gens:
        (r1, r2), deg = C.field_signature(gen)
        rows.append({
            "label": label, "source": "primitive_element",
            "generator_or_minpoly": str(gen),
            "degree": deg,
            "signature_r1_r2": f"({r1},{r2})",
            "trace_form_signature": "",
            "trace_form_predicted_m+1_m-1": "",
            "kind": classify(r1, r2),
            "paper_location": loc,
        })
    return rows


def cyclotomic_rows():
    """Q(zeta_n) is totally complex (0, phi(n)/2), disjoint from Salem (Prop. 12.3)."""
    rows = []
    for n in (5, 7, 8, 12):
        cp = sp.Poly(sp.cyclotomic_poly(n, x), x)
        coeffs = [int(c) for c in cp.all_coeffs()]
        r1, r2 = C.signature_from_minpoly(coeffs)
        rows.append({
            "label": f"Q(zeta_{n})", "source": "cyclotomic",
            "generator_or_minpoly": f"Phi_{n}(x)",
            "degree": cp.degree(),
            "signature_r1_r2": f"({r1},{r2})",
            "trace_form_signature": "",
            "trace_form_predicted_m+1_m-1": "",
            "kind": classify(r1, r2),
            "paper_location": "Prop 12.3 (circulant field, CM)",
        })
    return rows


def main():
    rows = minpoly_rows() + subfield_rows() + cyclotomic_rows()
    fields = ["label", "source", "generator_or_minpoly", "degree", "signature_r1_r2",
              "trace_form_signature", "trace_form_predicted_m+1_m-1", "kind", "paper_location"]
    C.write_csv("field_signatures.csv", fields, rows, SCRIPT)
    print("field_signatures.csv:")
    for r in rows:
        print(f"  {r['label']:26s} deg {r['degree']:2d}  sig {r['signature_r1_r2']:7s}  {r['kind']}")


if __name__ == "__main__":
    main()
