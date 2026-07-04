"""
Producer: the superselection sector / angle-charge structure (Section 'entry').

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * data/angle_charges.csv   angle charge A of catalog generators and the
        Salem/off-lattice witnesses (Prop 'charge'):
        object, poly, root_args_over_pi, all_on_(pi/2)Z, A
  * data/tensor_products.json   the reducibility witnesses (Section 'entry'):
        phi (x) phi = (x+1)^2 (x^2-3x+1), on-circle root -1 (order 2)
        phi^4 (x) phi^4 = (x-1)^2 (x^2-47x+1), on-circle root +1 (order 1)
        + the (pi/2)Z lattice closure under doubling and addition (Thm superselect)
  * data/entry_routes.csv   the six entry routes and their effect on the charge
        (the "reaches slot?" table of Section 'entry' / Thm 'rejects')

Run:  py code/2026-06-salem-slot/make_superselection.py
"""

from __future__ import annotations

import cmath

import sympy as sp

import salem_core as sc
import salem_io as io

x = sc.x
mp = sc.mp


def _args_over_pi(poly):
    return [cmath.phase(complex(r)) / cmath.pi for r in sp.nroots(sp.Poly(poly, x), n=30)]


def charge_rows():
    """Prop 'charge': A of the catalog generators + Salem/off-lattice witnesses."""
    rows = []
    generators = [
        ("real seed phi", x ** 2 - x - 1),
        ("K (imaginary pair)", x ** 4 + 5 * x ** 2 - 5),
    ]
    for name, poly in generators:
        args = _args_over_pi(poly)
        rows.append({
            "object": name,
            "poly": str(sp.expand(poly)),
            "root_args_over_pi": ";".join(f"{a:.4f}" for a in sorted(args)),
            "all_on_half_pi_Z": all(sc.on_half_pi_lattice(a) for a in args),
            "A": sc.angle_charge(poly),
        })
    # A Salem: Lehmer has on-circle conjugates at irrational angle => A=0
    lehmer = x ** 10 + x ** 9 - x ** 7 - x ** 6 - x ** 5 - x ** 4 - x ** 3 + x + 1
    args = _args_over_pi(lehmer)
    rows.append({
        "object": "Lehmer Salem (deg 10)",
        "poly": "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1",
        "root_args_over_pi": "irrational multiples of pi present",
        "all_on_half_pi_Z": all(sc.on_half_pi_lattice(a) for a in args),
        "A": sc.angle_charge(lehmer),
    })
    return rows


def beta4_offlattice():
    """Sec entry: beta_4 off-lattice trace root (1-sqrt13)/2 lifts to a Salem
    conjugate at arg 0.726 pi (not in (pi/2)Z)."""
    mp.mp.dps = 40
    tr = (1 - mp.sqrt(13)) / 2
    root = (tr + mp.sqrt(tr ** 2 - 4)) / 2   # complex lift, |root|=1
    arg_over_pi = mp.arg(root) / mp.pi
    return {
        "trace_root": mp.nstr(tr, 10),
        "on_unit_circle": abs(abs(root) - 1) < mp.mpf(10) ** (-30),
        "arg_over_pi": mp.nstr(arg_over_pi, 8),
        "off_half_pi_lattice": abs(2 * arg_over_pi - round(float(2 * arg_over_pi))) > 1e-3,
        "note": "off-lattice captured position lifts to an irrational-angle Salem conjugate",
    }


def tensor_products():
    """Sec entry: phi(x)phi and phi^4(x)phi^4 are reducible (cyclotomic x grow)."""
    seed = x ** 2 - x - 1
    roots = list(sp.roots(seed).keys())
    prods = [a * b for a in roots for b in roots]
    poly_pp = sp.expand(sp.prod([x - p for p in prods]))
    # phi^4 minimal polynomial x^2 - 7x + 1
    seed4 = sp.minimal_polynomial(sc.phi_sym ** 4, x)
    roots4 = list(sp.roots(seed4).keys())
    prods4 = [sp.simplify(a * b) for a in roots4 for b in roots4]
    poly_44 = sp.expand(sp.prod([x - p for p in prods4]))
    return {
        "phi_tensor_phi": {
            "factored": str(sp.factor(poly_pp)),
            "equals_(x+1)^2(x^2-3x+1)": sp.factor(poly_pp) == (x + 1) ** 2 * (x ** 2 - 3 * x + 1),
            "on_circle_root": "-1 (2nd root of unity, order 2)",
            "grow_block_roots": "phi^2, phi^-2 (roots of x^2-3x+1)",
        },
        "phi4_tensor_phi4": {
            "phi4_min_poly": str(sp.expand(seed4)),
            "factored": str(sp.factor(poly_44)),
            "equals_(x-1)^2(x^2-47x+1)": sp.factor(poly_44) == (x - 1) ** 2 * (x ** 2 - 47 * x + 1),
            "on_circle_root": "+1 (1st root of unity, order 1)",
        },
    }


def lattice_closure():
    """Thm superselect: (pi/2)Z closed under angle-doubling and angle-addition.
    In units of pi/2 mod 4 (fourth-roots-of-unity lattice)."""
    lattice = {0, 1, 2, 3}
    doubled = {(2 * k) % 4 for k in lattice}
    summed = {(j + k) % 4 for j in lattice for k in lattice}
    return {
        "lattice_units_of_half_pi_mod_4": sorted(lattice),
        "doubling_2k_mod4": sorted(doubled),
        "doubling_closed": doubled <= lattice,
        "addition_jk_mod4": sorted(summed),
        "addition_closed": summed <= lattice,
        "claim": "squaring doubles arguments (2*(pi/2)Z subset (pi/2)Z); "
                 "tensor adds arguments ((pi/2)Z+(pi/2)Z subset (pi/2)Z)",
    }


def entry_routes():
    """The six entry routes of Section 'entry' and their effect on charge A."""
    return [
        {"route": "trace lift L (linear)",
         "mechanism": "x^2-tx+1; lattice t in {-2,0,2} -> roots of unity",
         "effect_on_A": "preserves A=1", "reaches_slot": "no (reducible)"},
        {"route": "direct sum (+) (superposition)",
         "mechanism": "union of spectra",
         "effect_on_A": "preserves A=1", "reaches_slot": "no"},
        {"route": "block-diagonal (suspended superposition)",
         "mechanism": "captured (+) grow",
         "effect_on_A": "preserves A=1, reducible", "reaches_slot": "no"},
        {"route": "monodromy / curvature (non-linear travel)",
         "mechanism": "permutes {beta,1/beta} around branch point t=2",
         "effect_on_A": "preserves A (same root set)", "reaches_slot": "no"},
        {"route": "free commutator [A,B]",
         "mechanism": "off-diagonal coupling",
         "effect_on_A": "breaks the field", "reaches_slot": "not an operation"},
        {"route": "limit (Salem's theorem)",
         "mechanism": "Pisot -> Salem accumulation",
         "effect_on_A": "A=1 -> A=0 in the limit", "reaches_slot": "approached, never attained"},
    ]


def main():
    # angle_charges.csv
    rows = charge_rows()
    fields = ["object", "poly", "root_args_over_pi", "all_on_half_pi_Z", "A"]
    p_csv = io.write_csv("angle_charges.csv", fields, rows, __file__)

    # tensor_products.json (+ lattice closure + beta4 off-lattice witness)
    payload = {
        "_description": "Superselection sector: reducibility of the tensor products "
                        "(cyclotomic x grow), the (pi/2)Z lattice closure of the "
                        "operators, and the beta_4 off-lattice lift (Section 'entry').",
        "tensor_products": tensor_products(),
        "lattice_closure_Thm_superselect": lattice_closure(),
        "beta4_offlattice_lift": beta4_offlattice(),
    }
    p_json = io.write_json("tensor_products.json", payload, __file__)

    # entry_routes.csv
    er_fields = ["route", "mechanism", "effect_on_A", "reaches_slot"]
    p_er = io.write_csv("entry_routes.csv", er_fields, entry_routes(), __file__)

    print("wrote", p_csv)
    print("wrote", p_json)
    print("wrote", p_er)
    for r in rows:
        print(f"  {r['object']:24s} A={r['A']}")


if __name__ == "__main__":
    main()
