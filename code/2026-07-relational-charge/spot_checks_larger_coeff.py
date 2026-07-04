r"""
Producer: larger-coefficient degree-12 Salem spot checks (ledger V).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/spot_checks_larger_coeff.csv

Refactors Remark 6.14 / ledger V: the first five degree-12 Salem twist-classes
with a coefficient outside {-1,0,1} (lexicographic in {-2,...,2}^6) are all
relationally inert.  For each such class this emits the coefficient vector, the
polynomial, the ratio-object degree 144, and the complete contact signature
{Phi_1^12}.

Run: py code/2026-07-relational-charge/spot_checks_larger_coeff.py
"""

import itertools

import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x


def build(c):
    c1, c2, c3, c4, c5, c6 = c
    return sp.Poly([1, c1, c2, c3, c4, c5, c6, c5, c4, c3, c2, c1, 1], x)


def twist(c):
    c1, c2, c3, c4, c5, c6 = c
    return (-c1, c2, -c3, c4, -c5, c6)


def salem_pattern(p):
    a, b, mid, at2, atm2 = C.trace_sturm_pattern(p.as_expr())
    return a + b == 1 and mid == 5 and at2 == 0 and atm2 == 0


def find_first_five():
    """Lexicographic scan over {-2,...,2}^6 for the first five twist-classes
    that (a) use a coefficient outside {-1,0,1} and (b) are certified Salem."""
    seen, found = set(), []
    for c in itertools.product([-2, -1, 0, 1, 2], repeat=6):
        if c in seen:
            continue
        seen.add(c)
        seen.add(twist(c))
        if all(abs(v) <= 1 for v in c):  # must reach outside {-1,0,1}
            continue
        p = build(c)
        if p.eval(1) == 0 or p.eval(-1) == 0:
            continue
        if not salem_pattern(p):
            continue
        if not p.is_irreducible:
            continue
        found.append(c)
        if len(found) == 5:
            break
    return found


def main():
    reps = find_first_five()
    rows = []
    for c in reps:
        p = build(c)
        Rp = C.ratio_poly(p.as_expr())
        sig = C.cyclotomic_contacts(Rp)
        rows.append({
            "coeffs": " ".join(str(v) for v in c),
            "poly": str(p.as_expr()),
            "is_salem": "yes" if C.is_salem_polynomial(p.as_expr()) else "no",
            "rat_degree": Rp.degree(),
            "contact_signature": C.signature_str(sig),
            "inert": "yes" if sig == {1: 12} else "no",
        })
    fields = ["coeffs", "poly", "is_salem", "rat_degree", "contact_signature",
              "inert"]
    path = write_csv("spot_checks_larger_coeff.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} larger-coefficient Salem spot checks (ledger V); "
          f"all inert = {all(r['inert'] == 'yes' for r in rows)}")


if __name__ == "__main__":
    main()
