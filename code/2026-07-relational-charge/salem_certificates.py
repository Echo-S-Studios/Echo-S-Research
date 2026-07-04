r"""
Producer: Salem / Pisot certificates and inertness (Section 7).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/salem_certificates.csv

Refactors Theorem 7.9 (beta_4, Lehmer inert), Theorem 7.10 / ledger O
(beta_4, S_6, S_8, L certified Salem by reciprocity + irreducibility +
trace-polynomial Sturm pattern (1, 0, deg/2-1)), and Example 7.19 / ledger W
(the plastic number x^3-x-1, smallest Pisot, inert).

For each object it emits the full exact Salem/Pisot certificate: reciprocal,
irreducible, the trace-Sturm pattern, the root-modulus profile (one dominant
> 1, one < 1, the rest exactly on the unit circle for a Salem number), the
Mahler measure, the ratio-object degree n^2, the complete contact signature
(the inert diagonal {Phi_1^n} only), and the relational verdict (inert).

Run: py code/2026-07-relational-charge/salem_certificates.py
"""

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x

OBJECTS = [
    ("beta4", C.B4, "salem"),
    ("S6", C.S6, "salem"),
    ("S8", C.S8, "salem"),
    ("Lehmer", C.LEHMER, "salem"),
    ("plastic x^3-x-1", C.PLASTIC, "pisot"),
]


def modulus_profile(p):
    """(#|.|>1, #|.|<1, #|.|=1) at high precision."""
    mods = [abs(r) for r in C.roots_mp(p)]
    tol = mp.mpf(10) ** -18
    outside = sum(1 for m in mods if m > 1 + tol)
    inside = sum(1 for m in mods if m < 1 - tol)
    oncirc = sum(1 for m in mods if abs(m - 1) <= tol)
    return outside, inside, oncirc


def build_rows():
    rows = []
    for name, p, kind in OBJECTS:
        P = sp.Poly(p, x)
        n = P.degree()
        recip = C.is_reciprocal(p)
        irred = P.is_irreducible
        M = C.mahler_measure(p)
        Rp = C.ratio_poly(p)
        sig = C.cyclotomic_contacts(Rp)
        outside, inside, oncirc = modulus_profile(p)
        row = {
            "object": name,
            "kind": kind,
            "degree": n,
            "reciprocal": "yes" if recip else "no",
            "irreducible": "yes" if irred else "no",
            "mahler": mp.nstr(M, 18),
            "modulus_gt1": outside,
            "modulus_lt1": inside,
            "modulus_eq1": oncirc,
            "rat_degree": Rp.degree(),
            "contact_signature": C.signature_str(sig),
            "inert": "yes" if sig == {1: n} else "no",
        }
        if kind == "salem" and n % 2 == 0:
            a, b, mid, at2, atm2 = C.trace_sturm_pattern(p)
            d = n // 2
            row["trace_sturm_pattern"] = f"({a},{b},{mid})"
            row["trace_pattern_ok"] = (
                "yes" if (a == 1 and b == 0 and mid == d - 1 and at2 == 0 and atm2 == 0)
                else "no")
            row["is_salem_certificate"] = "yes" if C.is_salem_polynomial(p) else "no"
        else:
            row["trace_sturm_pattern"] = "n/a (odd degree)"
            row["trace_pattern_ok"] = "n/a"
            row["is_salem_certificate"] = "n/a (Pisot, not Salem)"
        rows.append(row)
    return rows


def main():
    rows = build_rows()
    fields = ["object", "kind", "degree", "reciprocal", "irreducible",
              "trace_sturm_pattern", "trace_pattern_ok", "is_salem_certificate",
              "mahler", "modulus_gt1", "modulus_lt1", "modulus_eq1",
              "rat_degree", "contact_signature", "inert"]
    path = write_csv("salem_certificates.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} certificates (Thm 7.9, Thm 7.10 / ledger O, P, W)")
    for r in rows:
        print(f"    {r['object']:<16} {r['kind']:<6} M={r['mahler'][:9]:<9}"
              f" sig {r['contact_signature']:<14} inert={r['inert']}")


if __name__ == "__main__":
    main()
