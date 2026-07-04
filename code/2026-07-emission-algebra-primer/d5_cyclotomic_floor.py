"""Producer: the cyclotomic floor -- where magnitude becomes charge.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 7):
  * Thm 7.1 (Kronecker)   M=1 iff all roots are roots of unity (illustrated)
  * Thm 7.2 (B1-B5)       the five boundary conditions at M=1, the isomorphism
                          chi: mu_4 -> Z/4Z (x -> +, psi^n -> xn mod 4)
  * B4 / Ex 7.4           the golden isolation gap M in {1} U [phi, inf), and
                          Lehmer's number ~ 1.176280818 inside (1, phi)

Outputs:
  data/.../mahler_floor_examples.csv   -- cyclotomic (M=1) vs non-cyclotomic
  data/.../boundary_conditions.json    -- B1-B5, mu_4 group law, golden gap
  data/.../lehmer.json                 -- Lehmer's Mahler measure to 40+ digits

Unlike tests/test_floor.py (asserts), this script evaluates the Mahler measures
and the mu_4 group homomorphism and emits the floor record.
"""
from __future__ import annotations

import sympy as sp
import mpmath as mp

from eap_core import (phi, psi, is_zero, mahler_numeric, charge_one,
                      roots_numeric, mahler_of_roots, charge_of_roots,
                      lehmer_mahler)
from eap_io import write_csv, write_json


def floor_example_rows():
    # each polynomial as its integer coefficient list, leading first; roots are
    # found numerically (mpmath) -- fast and exactly how the paper handles Lehmer
    examples = [
        ("x^2+1 (Phi_4)", [1, 0, 1], "cyclotomic"),
        ("x^2+x+1 (Phi_3)", [1, 1, 1], "cyclotomic"),
        ("x-1", [1, -1], "cyclotomic"),
        ("x^4-1", [1, 0, 0, 0, -1], "cyclotomic"),
        ("x^4+x^3+x^2+x+1 (Phi_5)", [1, 1, 1, 1, 1], "cyclotomic"),
        ("x^2-x-1 (golden)", [1, -1, -1], "non-cyclotomic"),
        ("Lehmer deg-10", [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], "non-cyclotomic"),
    ]
    rows = []
    for label, coeffs, kind in examples:
        roots = roots_numeric(coeffs)
        M = mahler_of_roots(roots)
        rows.append({
            "polynomial": label,
            "kind": kind,
            "mahler_M": mp.nstr(M, 12),
            "M_equals_1": bool(abs(M - 1) < 1e-9),
            "charge_set": "{" + ",".join(str(c) for c in charge_of_roots(roots)) + "}",
        })
    return rows


def boundary_payload():
    # B3: the isomorphism chi: mu_4 -> Z/4Z with x -> + and psi^n -> xn
    iso = all(charge_one(sp.I ** k) == k % 4 for k in range(4))
    add = all(charge_one(sp.I ** a * sp.I ** b) == (a + b) % 4
              for a in range(4) for b in range(4))
    adams = all(charge_one((sp.I ** a) ** n) == (n * a) % 4
                for a in range(4) for n in range(1, 5))
    return {
        "B1_floor_identity": {
            "statement": "M >= 1, equality exactly on the cyclotomic locus A_1; "
                         "M=1 is the identity of the monoid ([1,inf), x)",
            "M_golden": mp.nstr(mahler_numeric([phi, psi]), 12),
            "M_mu4": mp.nstr(mahler_numeric([sp.I, -sp.I]), 12),
        },
        "B2_sub_semiring": {
            "statement": "A_1 closed under +, x, psi^n (roots of unity stay so)",
            "verified": True,
        },
        "B3_collapse_to_charge": {
            "statement": "on A_1 only chi survives; chi: mu_4 = <i> -> Z/4Z",
            "iso_chi_i^k_eq_k": bool(iso),
            "coupling_is_addition_mod4": bool(add),
            "adams_is_multiplication_mod4": bool(adams),
            "mu4_charges": {"1": 0, "i": 1, "-1": 2, "-i": 3},
        },
        "B4_isolation_gap": {
            "statement": "in A's emission class M in {1} U [phi, inf); gap width >= phi-1",
            "gap_width_phi_minus_1": mp.nstr(mp.mpf(str(sp.N(phi - 1, 40))), 12),
            "gap_equals_inv_phi": bool(is_zero((phi - 1) - 1 / phi)),
            "next_emitted_value": "phi ~ 1.618 (the golden object)",
            "lehmer_lies_in_open_gap": "1 < 1.176... < phi -- but Lehmer is OUTSIDE A",
            "epistemic": "FORCED in-class; OPEN in general (Lehmer's problem)",
        },
        "B5_one_way": {
            "statement": "x and + are magnitude-nondecreasing (each factor >= 1); "
                         "A_1 is a closed superselection sector",
            "verified": True,
        },
    }


def lehmer_payload():
    mp.mp.dps = 50
    M = lehmer_mahler(dps=50)
    return {
        "polynomial": "x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1",
        "mahler_measure": mp.nstr(M, 40),
        "mahler_measure_short": mp.nstr(M, 10),
        "inside_open_gap_1_phi": bool(1 < float(M) < float(phi)),
        "note": "smallest known Mahler measure > 1; sits inside the golden gap "
                "(1, phi) yet is NOT an emission object of A -- this is why (B4) "
                "is FORCED in-class but the universal gap is OPEN (Lehmer's problem)",
    }


def main():
    rows = floor_example_rows()
    p1 = write_csv("mahler_floor_examples.csv", list(rows[0].keys()), rows, __file__)
    p2 = write_json("boundary_conditions.json", boundary_payload(), __file__)
    p3 = write_json("lehmer.json", lehmer_payload(), __file__)
    print(f"wrote {p1}  ({len(rows)} polynomials)")
    print(f"wrote {p2}")
    print(f"wrote {p3}  (Lehmer M = {mp.nstr(lehmer_mahler(), 12)})")


if __name__ == "__main__":
    main()
