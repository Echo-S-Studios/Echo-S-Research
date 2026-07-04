"""Producer: the core is the maximal order of the golden field.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Thm 1.7 and the Intuition after it):
  * disc(x^2-x-1)=5 is fundamental, field disc(Q(sqrt5))=5 (since 5 = 1 mod 4),
    so [O_K : Z[phi]] = 1 and Z[R] = O_{Q(sqrt5)} is the MAXIMAL order;
  * the clumsy seed 2*phi has minimal polynomial x^2-2x-4, disc 20 = 2^2*5,
    giving a defective suborder of index 2.

Output:
  data/.../maximal_order.json

Unlike tests/test_core.py, this script computes the discriminants and index
from sympy's discriminant/minimal_polynomial and emits the arithmetic record.
"""
from __future__ import annotations

import sympy as sp

from eap_io import write_json


def field_disc(d: int) -> int:
    """Discriminant of Q(sqrt d), d squarefree: d if d = 1 (mod 4) else 4d."""
    return d if d % 4 == 1 else 4 * d


def main():
    x = sp.symbols("x")

    # the golden seed R: defining polynomial x^2 - x - 1
    disc_golden = int(sp.discriminant(x ** 2 - x - 1, x))
    fdisc = field_disc(5)
    index_sq_golden = sp.Rational(disc_golden, fdisc)
    index_golden = sp.sqrt(index_sq_golden)

    # the clumsy seed 2*phi = 1 + sqrt5
    twophi = sp.expand(2 * ((1 + sp.sqrt(5)) / 2))
    minpoly = sp.expand(sp.minimal_polynomial(twophi, x))
    disc_clumsy = int(sp.discriminant(minpoly, x))
    index_sq_clumsy = sp.Rational(disc_clumsy, fdisc)
    index_clumsy = sp.sqrt(index_sq_clumsy)

    payload = {
        "golden_seed": {
            "defining_polynomial": "x^2 - x - 1",
            "disc_of_polynomial": disc_golden,
            "field_disc_Q_sqrt5": fdisc,
            "field_disc_rule": "5 = 1 (mod 4) so disc(Q(sqrt5)) = 5 (fundamental)",
            "index_squared": int(index_sq_golden),
            "index_OK_over_ZR": int(index_golden),
            "is_maximal_order": bool(index_golden == 1),
            "conclusion": "Z[R] = Z[phi] = O_{Q(sqrt5)} (nothing missing)",
        },
        "clumsy_seed_2phi": {
            "element": "2*phi = 1 + sqrt5",
            "minimal_polynomial": sp.sstr(minpoly),
            "disc_of_polynomial": disc_clumsy,
            "disc_factored": "20 = 2^2 * 5",
            "index_squared": int(index_sq_clumsy),
            "index_OK_over_suborder": int(index_clumsy),
            "is_maximal_order": bool(index_clumsy == 1),
            "conclusion": "proper suborder of index 2 -- arithmetically defective",
        },
    }
    p = write_json("maximal_order.json", payload, __file__)
    print(f"wrote {p}")
    print(f"  golden index = {int(index_golden)} (maximal), "
          f"clumsy index = {int(index_clumsy)}")


if __name__ == "__main__":
    main()
