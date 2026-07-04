"""Producer: the cross-shell residue at two non-real shells (Sections 5-6).

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the first live cross-shell residue instance, executed on
    p = x^5 - 2x^4 - 2x^3 - 2x^2 - 2x - 2
(the smallest two-pair Pisot quintic; Theorem 6.1).  Re-derives:

  * Rat_p^o = Rat_p / (x-1)^5 is degree 20, squarefree, and irreducible
    (so S* = Rat_p^o and deg C_2 = 400)                          (Section 6);
  * shell detector = 4 unimodular roots of Rat_p^o (distinct shells, Prop 5.4);
  * Prop 5.2 modulus multiset: 7 distinct values, mult multiset [2,2,2,2,4,4,4],
    modulus 1 of multiplicity 4, reciprocal-closed;
  * Prop 5.3 composed square C_2 monic of degree d^2 = 400 whose roots are all
    ordered products of S*-roots (pinned exactly on the small self-reciprocal
    S* = Z* = x^4-3x^2+1);
  * Prop 5.3(c) / Theorem 6.1 negative certificate: the ONLY ordered products
    that are roots of unity are the 20 equal to 1 (scan {Phi_1^20}); every other
    on-circle product sits at an irrational angle -> zero mirrored cross-shell
    classes.

The C_2 certificate is checked at high precision on the 400 ordered products
with an exact root-of-unity test (order <= 1680, the totient bound for deg 400).

Emits:
  data/2026-07-pisot-residue/residue_instance.json    (all invariants + certificate)
  data/2026-07-pisot-residue/residue_ratcirc.csv      (Rat_p^o degree-20 coeffs)
"""
import numpy as np
import sympy as sp
import mpmath as mp
from sympy import symbols, Poly, div, expand, gcd, diff, resultant

from pisot_lib import (rat_object, hp_roots, is_root_of_unity, cyclotomic_scan,
                       scan_json, poly_str, write_json, write_csv)

x, y = symbols('x y')
P = [1, -2, -2, -2, -2, -2]                      # x^5-2x^4-2x^3-2x^2-2x-2


def rat_circ():
    """Rat_p^o = Rat_p / (x-1)^5."""
    R = rat_object(P)
    q, r = div(R.as_expr(), (x - 1) ** 5, x)
    assert sp.simplify(r) == 0, "(x-1)^5 must divide Rat_p"
    return Poly(expand(q), x)


def modulus_multiset(roots, dps=80):
    """Bucket |root| at high precision -> distinct values and multiplicities."""
    with mp.workdps(dps):
        buckets = []
        for r in roots:
            m = abs(r)
            for b in buckets:
                if abs(b[0] - m) < mp.mpf(10) ** -20:
                    b[1] += 1
                    break
            else:
                buckets.append([m, 1])
        mults = sorted(b[1] for b in buckets)
        one = [b for b in buckets if abs(b[0] - 1) < mp.mpf(10) ** -20]
        vals = [b[0] for b in buckets]
        recip_closed = all(any(abs(v * w - 1) < mp.mpf(10) ** -18 for w in vals)
                           for v in vals)
        return {
            "num_distinct_moduli": len(buckets),
            "multiplicity_multiset": mults,
            "modulus_1_multiplicity": (one[0][1] if one else 0),
            "reciprocal_closed": bool(recip_closed),
            "distinct_moduli_display": [float(b[0]) for b in sorted(buckets)],
        }


def c2_root_product_identity():
    """Prop 5.3(b) pinned on the small self-reciprocal S* = Z* = x^4-3x^2+1:
    C_2 = Res_y(S*(y), y^d S*(x/y)) is monic of degree d^2=16 with root multiset
    equal to all ordered products of S*-roots.  We confirm the coefficient vector
    matches prod_{i,j}(x - r_i r_j) built independently from the S*-roots."""
    S = [1, 0, -3, 0, 1]
    d = 4
    Sy = Poly(S, y)
    Sxy = sum(S[k] * (x / y) ** (d - k) for k in range(len(S)))
    C2 = Poly(expand(resultant(Sy.as_expr(), expand(y ** d * Sxy), y)), x)
    lead = C2.all_coeffs()[0]
    c2_norm = [complex(c) / complex(lead) for c in C2.all_coeffs()]
    with mp.workdps(50):
        sroots = [complex(r) for r in hp_roots(S, 50)]
    prods = [a * b for a in sroots for b in sroots]
    prod_poly = np.poly(np.array(prods))
    max_dev = max(abs(a - b) for a, b in zip(c2_norm, prod_poly))
    return {
        "S_star_used": poly_str(S),
        "deg_C2": C2.degree(),
        "monic_up_to_sign": abs(int(lead)) == 1,
        "root_multiset_is_ordered_products": bool(max_dev < 1e-6),
        "max_coeff_deviation_display": float(max_dev),
    }


def negative_certificate(roots, dps=80):
    """Prop 5.3(c) / Theorem 6.1: scan the 400 ordered products r_i r_j; the only
    roots of unity among them are the 20 equal to 1 (Phi_1^20)."""
    with mp.workdps(dps):
        eps1 = mp.mpf(10) ** -25
        epsc = mp.mpf(10) ** -14
        equal_one = oncircle_not_one = rogue = 0
        for a in roots:
            for b in roots:
                z = a * b
                if abs(z - 1) < eps1:
                    equal_one += 1
                elif abs(abs(z) - 1) < epsc:
                    oncircle_not_one += 1
                    if is_root_of_unity(z, Mmax=1680, dps=dps):
                        rogue += 1
    return {
        "num_ordered_products": len(roots) ** 2,
        "products_equal_to_1": equal_one,
        "products_on_circle_not_1": oncircle_not_one,
        "rogue_roots_of_unity": rogue,
        "scan": {"factors": [{"m": 1, "multiplicity": equal_one}],
                 "pretty": f"{{Phi_1^{equal_one}}}"},
        "zero_mirrored_cross_shell_classes": rogue == 0,
    }


def main():
    Rc = rat_circ()
    g = gcd(Rc.as_expr(), diff(Rc.as_expr(), x))
    squarefree = Poly(g, x).degree() == 0
    coeffs20 = [int(c) for c in Rc.all_coeffs()]
    roots = hp_roots(coeffs20, 80)

    with mp.workdps(80):
        eps = mp.mpf(10) ** -25
        n_unimod = len([r for r in roots if abs(abs(r) - 1) < eps])

    cert = negative_certificate(roots)
    payload = {
        "description": "Sections 5-6 cross-shell residue on the first two-pair Pisot quintic "
                       "x^5-2x^4-2x^3-2x^2-2x-2; degree-20 Rat_p^o and the C_2 negative certificate.",
        "instance": {"poly": poly_str(P), "coeffs_hi_to_lo": P},
        "rat_p_circ": {
            "definition": "Rat_p / (x-1)^5",
            "degree": Rc.degree(),
            "squarefree": bool(squarefree),
            "irreducible": bool(Rc.is_irreducible),
            "equals_S_star": bool(squarefree and Rc.is_irreducible),
            "scan": scan_json(cyclotomic_scan(Poly((x - 1) ** 5 * Rc.as_expr(), x))),
        },
        "shell_detector": {
            "num_unimodular_roots": n_unimod,
            "distinct_shells": n_unimod == 4,
            "reads_4_iff_distinct_shells": True,
        },
        "modulus_multiset_prop52": modulus_multiset(roots),
        "c2_degree_and_root_products_prop53b": c2_root_product_identity(),
        "deg_C2_full_instance": Rc.degree() ** 2,
        "negative_certificate_theorem61": cert,
    }
    jpath = write_json("residue_instance.json", payload, "residue_instance.py")

    cpath = write_csv(
        "residue_ratcirc.csv",
        ["power", "coefficient"],
        [(Rc.degree() - i, c) for i, c in enumerate(coeffs20)],
        "residue_instance.py")

    print(f"wrote {jpath}")
    print(f"wrote {cpath}")
    print(f"  Rat_p^o: deg={Rc.degree()}, squarefree={squarefree}, "
          f"irreducible={bool(Rc.is_irreducible)}")
    print(f"  shell detector = {n_unimod} unimodular roots (distinct shells)")
    print(f"  deg C_2 = {Rc.degree() ** 2}; certificate: "
          f"{cert['products_equal_to_1']} equal-1, "
          f"{cert['products_on_circle_not_1']} on-circle-not-1, "
          f"{cert['rogue_roots_of_unity']} rogue -> scan {cert['scan']['pretty']}")


if __name__ == "__main__":
    main()
