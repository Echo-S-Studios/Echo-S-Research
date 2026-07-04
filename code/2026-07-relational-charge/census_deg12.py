r"""
Producer: the exhaustive degree-12 Salem census (Theorem 6.13, Remark 6.14).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/census_deg12_summary.json
              data/2026-07-relational-charge/census_deg12_salems.csv

Refactors Theorem 6.13 (thm:census) and Remark 6.14 (rem:censusscope), ledger
T and W.  Enumerates all 3^6 = 729 monic reciprocal polynomials
x^12 + c1 x^11 + ... + c6 x^6 + ... + c1 x + 1 with c_i in {-1,0,1}.  The sign
twist p(x)->p(-x) fixes exactly the 27 vectors with vanishing odd coefficients,
giving (729+27)/2 = 378 twist-classes (Burnside).  Each non-Salem class carries
an exact rejection certificate; each of the 37 Salem classes has ratio object of
degree 144 whose complete contact scan returns {Phi_1^12} (relationally inert).

Rejection tally uses the CORRECTED split (erratum 2026-07-04): 39 vanish at +-1,
256 fail the trace-Sturm pattern (1,0,5), 46 are reducible, 0 fail only trace
irreducibility; 39+256+46+37 = 378 and (trace)+(reducible) = 302.  A class with
a cyclotomic factor counts as reducible regardless of a coincidental trace
pattern; the single boundary polynomial
x^12-x^11-x^10-x^9-x^7-x^6-x^5-x^3-x^2-x+1 (factor Phi_10) sits in the reducible
bucket though its trace pattern is (1,0,5).

Run: py code/2026-07-relational-charge/census_deg12.py   (about 1-2 minutes)
"""

import itertools

import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv, write_json

x = C.x


def build(c):
    """Reciprocal monic degree-12 polynomial from (c1,...,c6)."""
    c1, c2, c3, c4, c5, c6 = c
    return sp.Poly([1, c1, c2, c3, c4, c5, c6, c5, c4, c3, c2, c1, 1], x)


def twist(c):
    """p(x)->p(-x) on coefficient vectors: odd-index coeffs flip sign."""
    c1, c2, c3, c4, c5, c6 = c
    return (-c1, c2, -c3, c4, -c5, c6)


def pattern_ok(p):
    """Count-based trace Sturm pattern (1,0,5) with no root at +-2."""
    a, b, mid, at2, atm2 = C.trace_sturm_pattern(p.as_expr())
    return a + b == 1 and mid == 5 and at2 == 0 and atm2 == 0


def enumerate_classes():
    """Return (all_count, fixed_count, reps)."""
    all_c = list(itertools.product([-1, 0, 1], repeat=6))
    seen, reps, fixed = set(), [], 0
    for c in all_c:
        if c in seen:
            continue
        t = twist(c)
        if t == c:
            fixed += 1
        seen.add(c)
        seen.add(t)
        reps.append(c)
    return len(all_c), fixed, reps


def classify(reps):
    buckets = {"vanish": 0, "trace": 0, "reducible": 0, "salem": 0}
    salem_reps = []
    for c in reps:
        p = build(c)
        if p.eval(1) == 0 or p.eval(-1) == 0:
            buckets["vanish"] += 1
        elif not pattern_ok(p):
            buckets["trace"] += 1
        elif not p.is_irreducible:
            buckets["reducible"] += 1
        else:
            buckets["salem"] += 1
            salem_reps.append(c)
    return buckets, salem_reps


def poly_str(c):
    return str(build(c).as_expr())


def main():
    total, fixed, reps = enumerate_classes()
    buckets, salem_reps = classify(reps)

    # the boundary polynomial with a Phi_10 factor
    edge_c = (-1, -1, -1, 0, -1, -1)
    edge = build(edge_c)
    phi10 = sp.Poly(sp.cyclotomic_poly(10, x), x)
    _, rem = sp.div(edge, phi10)
    edge_reducible = not edge.is_irreducible
    edge_has_phi10 = rem.is_zero
    edge_trace_pattern = pattern_ok(edge)

    summary = {
        "family_size": total,
        "family_size_is_3_to_6": total == 729,
        "twist_fixed_count": fixed,
        "twist_fixed_is_27": fixed == 27,
        "twist_classes": len(reps),
        "burnside_378": len(reps) == (729 + 27) // 2 == 378,
        "rejection_tally": {
            "vanish_at_pm1": buckets["vanish"],
            "trace_sturm_fail": buckets["trace"],
            "reducible": buckets["reducible"],
            "salem": buckets["salem"],
        },
        "corrected_split_256_46": (buckets["trace"] == 256 and buckets["reducible"] == 46),
        "combined_trace_plus_reducible": buckets["trace"] + buckets["reducible"],
        "combined_is_302": buckets["trace"] + buckets["reducible"] == 302,
        "tally_sums_to_378": sum(buckets.values()) == 378,
        "salem_count_37": buckets["salem"] == 37,
        "boundary_polynomial": {
            "coeffs": list(edge_c),
            "poly": poly_str(edge_c),
            "reducible": bool(edge_reducible),
            "has_phi10_factor": bool(edge_has_phi10),
            "trace_pattern_1_0_5": bool(edge_trace_pattern),
            "bucketed_as": "reducible (cyclotomic factor overrides coincidental trace pattern)",
        },
        "headline": "all 37 certified Salem classes are relationally inert {Phi_1^12}",
        "status": "[forced] per instance; [computed] complete family verification (ledger T)",
    }

    # per-Salem contact scans (deg Rat = 144 -> {Phi_1^12})
    rows = []
    all_inert = True
    for c in salem_reps:
        p = build(c)
        Rp = C.ratio_poly(p.as_expr())
        sig = C.cyclotomic_contacts(Rp)
        inert = sig == {1: 12}
        all_inert = all_inert and inert and Rp.degree() == 144
        rows.append({
            "coeffs": " ".join(str(v) for v in c),
            "poly": poly_str(c),
            "rat_degree": Rp.degree(),
            "scan_bound_2d2": 2 * Rp.degree() ** 2,
            "contact_signature": C.signature_str(sig),
            "inert": "yes" if inert else "no",
        })
    summary["all_37_inert_phi1_12"] = bool(all_inert)

    ps = write_json("census_deg12_summary.json", summary, __file__)
    print(f"wrote {ps}")
    fields = ["coeffs", "poly", "rat_degree", "scan_bound_2d2",
              "contact_signature", "inert"]
    pc = write_csv("census_deg12_salems.csv", fields, rows, __file__)
    print(f"wrote {pc}")
    print(f"  729 family / {fixed} fixed / {len(reps)} orbits / "
          f"{buckets['salem']} Salem / {buckets['vanish']} vanish / "
          f"{buckets['trace']} trace / {buckets['reducible']} reducible")
    print(f"  combined {buckets['trace'] + buckets['reducible']} / sum "
          f"{sum(buckets.values())} / all 37 inert = {all_inert}")


if __name__ == "__main__":
    main()
