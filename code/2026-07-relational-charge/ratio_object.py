r"""
Producer: the ratio object Rat_p and its diagonal normalization.

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/ratio_objects.csv

Refactors Definition 6.3 (def:rat) and Lemma 6.6 (lem:diag).  For each named
monic p with p(0) != 0 it emits, in one row: the degree of Rat_p (which is
exactly n^2, "no drop"), the leading coefficient (+-p(0))^n of the raw
resultant, the Phi_1 = (x-1) multiplicity (= n for squarefree p, the diagonal
ratios alpha_i/alpha_i), whether p is squarefree, and -- for |p(0)| = 1 -- the
two-route agreement (ledger G) between the resultant construction and the
companion Kronecker matrix C_p (x) C_p^{-1}.

Run: py code/2026-07-relational-charge/ratio_object.py
"""

import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x

# name, polynomial, degree n
CASES = [
    ("x^3-2", x**3 - 2, 3),
    ("x^3+2", x**3 + 2, 3),
    ("x^4-2", x**4 - 2, 4),
    ("q2=x^4+x^2-1", C.Q2, 4),
    ("x^4+5x^2+5", C.GROUPDROP, 4),
    ("K=x^4+5x^2-5", C.KSEED, 4),
    ("beta4", C.B4, 4),
    ("S6", C.S6, 6),
    ("S8", C.S8, 8),
    ("Lehmer", C.LEHMER, 10),
    ("plastic=x^3-x-1", C.PLASTIC, 3),
    ("x^2+x+2", x**2 + x + 2, 2),
    ("x^4-x+1", C.X4MX1, 4),
    ("x^4+x^2+2", C.TWISTSHELL, 4),
]


def build_rows():
    rows = []
    for name, p, n in CASES:
        P = sp.Poly(p, x)
        Rp = C.ratio_poly(p)
        # raw (un-normalised) resultant, to exhibit the leading coeff (+-p(0))^n
        Rraw = sp.resultant(P.as_expr().subs(x, C.y),
                            sp.expand(P.as_expr().subs(x, x * C.y)), C.y)
        raw_deg = sp.Poly(sp.expand(Rraw), x).degree()
        squarefree = sp.gcd(P, P.diff(x)).degree() == 0
        p0 = int(P.eval(0))
        # leading coefficient of the raw resultant is (+-p(0))^n (no degree drop)
        raw_lead = int(sp.Poly(sp.expand(Rraw), x).LC())
        two_route = ""
        if abs(p0) == 1:
            s_res = C.cyclotomic_contacts(Rp)
            s_kron = C.cyclotomic_contacts(C.ratio_poly_via_kronecker(p))
            two_route = "agree" if s_res == s_kron else "DISAGREE"
        rows.append({
            "object": name,
            "degree_n": n,
            "rat_degree": Rp.degree(),
            "rat_degree_equals_n2": "yes" if Rp.degree() == n * n else "no",
            "raw_resultant_degree": raw_deg,
            "raw_leading_coeff": raw_lead,
            "abs_raw_lead_equals_p0_to_n": "yes" if abs(raw_lead) == abs(p0) ** n else "no",
            "p0": p0,
            "squarefree": "yes" if squarefree else "no",
            "phi1_multiplicity": C.phi1_multiplicity(Rp),
            "phi1_equals_n": "yes" if C.phi1_multiplicity(Rp) == n else "no",
            "two_route_kronecker": two_route,
        })
    return rows


def main():
    rows = build_rows()
    fields = ["object", "degree_n", "rat_degree", "rat_degree_equals_n2",
              "raw_resultant_degree", "raw_leading_coeff",
              "abs_raw_lead_equals_p0_to_n", "p0", "squarefree",
              "phi1_multiplicity", "phi1_equals_n", "two_route_kronecker"]
    path = write_csv("ratio_objects.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} ratio objects (Definition 6.3 / Lemma 6.6 / ledger G)")


if __name__ == "__main__":
    main()
