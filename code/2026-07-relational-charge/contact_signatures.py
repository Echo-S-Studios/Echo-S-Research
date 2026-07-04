r"""
Producer: complete cyclotomic-contact signatures (Appendix A ledger).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/contact_signatures.csv

Runs the paper's own probe (Appendix B: bounded cyclotomic trial-division to
the complete totient bound 2(deg P)^2 of Lemma 4.5) on the ratio object of every
named instance in Appendix A, and on the mixed ratio objects of Section 8.
Emits, per object: deg Rat_p (or Rat_{p,q}), the complete scan bound 2 deg^2,
and the complete contact signature {Phi_M^mult}.

Reproduces ledger entries:
  A x^3-2 {Phi_1^3,Phi_3^3}      B x^3+2 {Phi_1^3,Phi_3^3} (gauge-blind, = A)
  C x^4-2 {Phi_1^4,Phi_2^4,Phi_4^4}   D q2 {Phi_1^4,Phi_2^4}
  E x^4+5x^2+5 {Phi_1^4,Phi_2^4}      F K=x^4+5x^2-5 {Phi_1^4,Phi_2^4} (=E; type differs)
  G beta4 {Phi_1^4}              H Lehmer {Phi_1^10}
  P S6 {Phi_1^6}, S8 {Phi_1^8}   U x^2+x+2 {Phi_1^2}, x^4+x^2+2 {Phi_1^4,Phi_2^4}
  W plastic {Phi_1^3}           X x^4-x+1 {Phi_1^4}
  I/Q the six mixed Salem pairs -> empty signature (no circle-locking)

Run: py code/2026-07-relational-charge/contact_signatures.py
"""

import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x

# ledger, object-name, polynomial (single-object ratio scans)
SINGLE = [
    ("A", "x^3-2", x**3 - 2),
    ("B", "x^3+2", x**3 + 2),
    ("C", "x^4-2", x**4 - 2),
    ("D", "q2=x^4+x^2-1", C.Q2),
    ("E", "x^4+5x^2+5", C.GROUPDROP),
    ("F", "K=x^4+5x^2-5", C.KSEED),
    ("G", "beta4", C.B4),
    ("H", "Lehmer", C.LEHMER),
    ("P", "S6", C.S6),
    ("P", "S8", C.S8),
    ("U", "q=x^2+x+2", x**2 + x + 2),
    ("U", "p=x^4+x^2+2", C.TWISTSHELL),
    ("W", "plastic=x^3-x-1", C.PLASTIC),
    ("X", "x^4-x+1", C.X4MX1),
]

# ledger, pair-name, (p, q) mixed ratio scans (Section 8: circle-locking)
MIXED = [
    ("I", "L x beta4", C.LEHMER, C.B4),
    ("Q", "beta4 x S6", C.B4, C.S6),
    ("Q", "beta4 x S8", C.B4, C.S8),
    ("Q", "S6 x S8", C.S6, C.S8),
    ("Q", "L x S6", C.LEHMER, C.S6),
    ("Q", "L x S8", C.LEHMER, C.S8),
]


def build_rows():
    rows = []
    for ledger, name, p in SINGLE:
        Rp = C.ratio_poly(p)
        d = Rp.degree()
        sig = C.cyclotomic_contacts(Rp)
        rows.append({
            "ledger": ledger,
            "object": name,
            "kind": "ratio",
            "rat_degree": d,
            "scan_bound_2d2": 2 * d * d,
            "contact_signature": C.signature_str(sig),
            "phi1_multiplicity": C.phi1_multiplicity(Rp),
        })
    for ledger, name, p, q in MIXED:
        coprime = sp.gcd(sp.Poly(p, x), sp.Poly(q, x)).degree() == 0
        Rpq = C.mixed_ratio_poly(p, q)
        d = Rpq.degree()
        sig = C.cyclotomic_contacts(Rpq)
        rows.append({
            "ledger": ledger,
            "object": name + (" (gcd=1)" if coprime else " (gcd>1)"),
            "kind": "mixed",
            "rat_degree": d,
            "scan_bound_2d2": 2 * d * d,
            "contact_signature": C.signature_str(sig),
            "phi1_multiplicity": C.phi1_multiplicity(Rpq),
        })
    return rows


def main():
    rows = build_rows()
    fields = ["ledger", "object", "kind", "rat_degree", "scan_bound_2d2",
              "contact_signature", "phi1_multiplicity"]
    path = write_csv("contact_signatures.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} contact signatures (Appendix A: A-H, P, U, W, X, I, Q)")
    for r in rows:
        print(f"    {r['ledger']:>2} {r['object']:<22} deg {r['rat_degree']:>3}"
              f"  {r['contact_signature']}")


if __name__ == "__main__":
    main()
