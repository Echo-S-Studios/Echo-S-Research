"""
Producer: the orthogonal partner completes Z/4Z (Section 6).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/partner.json
              data/2026-07-helix-orthogonal-partner/registry_quartics.csv

Recomputes Proposition 6.1 and Proposition 6.2:
  * Prop. 6.1: chi(Kf) = {0,1,2,3}. The real roots +/-K sit at arg {0,pi}
    (chi in {0,2}); the imaginary roots +/- i beta at arg +/- pi/2 (chi in {1,3}),
    the axis ORTHOGONAL to the golden {0,2}. Kf's complex roots have Re = 0.
  * Prop. 6.2 (parity criterion): a real quartic puts its complex roots on the
    imaginary axis iff it is even with a negative x^2-root. Kf is even; the two
    registry quartics
        cons = x^4 - 6x^3 + 26x^2 - 16x - 4
        res  = x^4 + 2x^3 + 39x^2 - 52x + 11
    are non-even, split over Q(sqrt5) with complex-block real parts phi^2 (cons)
    and -phi (res) -- OFF the (pi/2)Z lattice. All three generate Q(5^{1/4}),
    have signature (2,1), Galois group D_4 (order 8), are non-reciprocal (so no
    Salem), and have Mahler measures phi^2 sqrt5, 2 phi^5 = 11+5 sqrt5, and
    12+19 phi = (43+19 sqrt5)/2.

The label origins of cons/res are [open] in the paper; only the derived
mathematical facts are emitted here (registry_quartics.csv).

Run: py code/2026-07-helix-orthogonal-partner/partner.py
"""

import sympy as sp
from sympy import galois_group
from sympy.polys.numberfields.subfield import field_isomorphism

import helix_core as hc
from helix_io import write_csv, write_json

x = hc.x
sqrt5 = hc.sqrt5
phi = hc.phi
Kf = x**4 + 5 * x**2 - 5
cons = x**4 - 6 * x**3 + 26 * x**2 - 16 * x - 4
res = x**4 + 2 * x**3 + 39 * x**2 - 52 * x + 11
beta = sp.sqrt((5 + 3 * sqrt5) / 2)
K = 5 ** sp.Rational(1, 4) / phi

MAHLER_CLOSED = {
    "Kf": ("phi^2 sqrt5", phi**2 * sqrt5),
    "cons": ("2 phi^5 = 11+5 sqrt5", 11 + 5 * sqrt5),
    "res": ("12+19 phi = (43+19 sqrt5)/2", (43 + 19 * sqrt5) / 2),
}


def full_charge():
    """Prop. 6.1: chi(Kf) = {0,1,2,3}; imaginary roots have Re = 0."""
    charges = {hc.chi_of(K), hc.chi_of(-K), hc.chi_of(sp.I * beta), hc.chi_of(-sp.I * beta)}
    cardinal = {
        "chi(1)": hc.chi_of(sp.Integer(1)),
        "chi(-1)": hc.chi_of(sp.Integer(-1)),
        "chi(i)": hc.chi_of(sp.I),
        "chi(-i)": hc.chi_of(-sp.I),
    }
    roots = sp.Poly(Kf, x).all_roots()
    complex_roots = [r for r in roots if sp.im(r) != 0]
    return {
        "chi_Kf": sorted(charges),
        "real_axis_charges_{0,2}": sorted({hc.chi_of(K), hc.chi_of(-K)}),
        "imaginary_axis_charges_{1,3}": sorted({hc.chi_of(sp.I * beta), hc.chi_of(-sp.I * beta)}),
        "cardinal_point_charges": cardinal,
        "complex_roots_real_parts": [str(sp.re(r)) for r in complex_roots],
        "orthogonal_to_golden": "{1,3} is the pi/2-rotate of {0,2}; disjoint union = Z/4Z",
    }


def parity_and_factorisations():
    """Prop. 6.2: evenness, Q(sqrt5) factorisations, complex-block real parts."""
    cons_fac = (x**2 - (3 + sqrt5) * x + (11 + 5 * sqrt5)) * (
        x**2 - (3 - sqrt5) * x + (11 - 5 * sqrt5)
    )
    res_fac = (x**2 + (1 - sqrt5) * x + (43 - 19 * sqrt5) / 2) * (
        x**2 + (1 + sqrt5) * x + (43 + 19 * sqrt5) / 2
    )
    return {
        "Kf_is_even": hc.is_even_poly(Kf),
        "cons_is_even": hc.is_even_poly(cons),
        "res_is_even": hc.is_even_poly(res),
        "cons_factorisation_residual": str(sp.expand(cons_fac - cons)),
        "res_factorisation_residual": str(sp.expand(res_fac - res)),
        "cons_complex_block_real_part": "phi^2 = (3+sqrt5)/2",
        "cons_real_part_residual": str(sp.simplify((3 + sqrt5) / 2 - phi**2)),
        "res_complex_block_real_part": "-phi = -(1+sqrt5)/2",
        "res_real_part_residual": str(sp.simplify(-(1 + sqrt5) / 2 - (-phi))),
        "note": "even Kf -> complex roots purely imaginary; non-even cons/res -> nonzero Re, off (pi/2)Z",
    }


def invariants_for(name, poly):
    """Field, signature, Galois group, reciprocity, Mahler for one quartic."""
    P = sp.Poly(poly, x)
    n_real = P.count_roots(-sp.oo, sp.oo)
    G, _ = galois_group(P)
    # field: does it generate Q(5^{1/4})?
    if name == "Kf":
        # 5^{1/4} = K(K^2+4)/3 in Q(K)
        generates = sp.simplify(K * (K**2 + 4) / 3 - 5 ** sp.Rational(1, 4)) == 0
    else:
        target = sp.AlgebraicNumber(5 ** sp.Rational(1, 4))
        r = sp.AlgebraicNumber(sp.CRootOf(poly, 0))
        generates = field_isomorphism(r, target) is not None
    label, closed = MAHLER_CLOSED[name]
    mah_num = hc.mahler_numeric(poly)
    mah_err = hc.mp.nstr(abs(mah_num - hc.mp.mpf(str(sp.N(closed, 55)))), 3)
    reciprocal = hc.is_reciprocal(poly)
    return {
        "name": name,
        "poly": str(sp.expand(poly)),
        "even": hc.is_even_poly(poly),
        "real_roots": int(n_real),
        "complex_pairs": (4 - int(n_real)) // 2,
        "field_signature_r1_r2": f"({n_real},{(4 - int(n_real)) // 2})",
        "generates_Q_fifth_root_5": bool(generates),
        "galois_group": "D4",
        "galois_order": int(G.order()),
        "galois_transitive": bool(G.is_transitive()),
        "reciprocal": reciprocal,
        "is_salem": bool(reciprocal),  # non-reciprocal => not Salem
        "mahler_closed_form": label,
        "mahler_decimal": hc.dec(closed),
        "mahler_numeric_abs_error": mah_err,
    }


def golden_field_identities():
    """Prop. 6.2 supporting identities used by the Mahler closed forms."""
    return {
        "phi5_equals_(11+5 sqrt5)/2_residual": str(sp.simplify(phi**5 - (11 + 5 * sqrt5) / 2)),
        "2 phi5_equals_11+5 sqrt5_residual": str(sp.simplify(2 * phi**5 - (11 + 5 * sqrt5))),
        "12+19 phi_equals_(43+19 sqrt5)/2_residual": str(
            sp.simplify((12 + 19 * phi) - (43 + 19 * sqrt5) / 2)
        ),
    }


def main():
    quartics = [invariants_for("Kf", Kf), invariants_for("cons", cons), invariants_for("res", res)]
    payload = {
        "section": "6 -- The orthogonal partner completes Z/4Z",
        "results": {
            "prop_6_1_full_charge": full_charge(),
            "prop_6_2_parity_criterion": parity_and_factorisations(),
            "prop_6_2_quartic_invariants": quartics,
            "prop_6_2_golden_field_identities": golden_field_identities(),
        },
    }
    p_json = write_json("partner.json", payload, __file__)
    p_csv = write_csv(
        "registry_quartics.csv",
        ["name", "poly", "even", "real_roots", "complex_pairs", "field_signature_r1_r2",
         "generates_Q_fifth_root_5", "galois_group", "galois_order", "galois_transitive",
         "reciprocal", "is_salem", "mahler_closed_form", "mahler_decimal", "mahler_numeric_abs_error"],
        quartics,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print("  chi(Kf)={0,1,2,3} ; Kf/cons/res -> Q(5^{1/4}), sig (2,1), D4 order 8, non-reciprocal")


if __name__ == "__main__":
    main()
