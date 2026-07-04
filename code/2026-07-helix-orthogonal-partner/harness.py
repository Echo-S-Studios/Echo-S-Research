"""
Producer: the harness -- constraint-import as verification (Section 8).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/harness.json
              data/2026-07-helix-orthogonal-partner/harness_ledger.csv

Section 8 is largely metalogical; its load-bearing COMPUTATIONAL kernels are:
  * Prop. 8.4 (filter asymmetry): a single exact counterexample refutes a
    universal (delta = dashv), but no finite exact range promotes one (at most
    [computed]). Independent witness: Euler's n^2+n+41 is prime for n=0..39 and
    composite at n=40 = 41^2.
  * Prop. 8.6 (a generator is inert): on the closed {0,2} sub-semigroup the claim
    "no reachable charge lies in {1,3}" has no exact counterexample, so the
    casualty set C_dashv is empty and no partition forms.
  * Rem. 8.9 (the transferable kernel): sqrt5 = phi + phi^{-1} = phi - psi, the
    sqrt-branch of the 2-to-1 fold.

It also reproduces the ledger of Table 2 (survivors / casualties / residue) as a
machine-readable partition -- the paper's own certificate trail.

Run: py code/2026-07-helix-orthogonal-partner/harness.py
"""

import sympy as sp

import helix_core as hc
from helix_io import write_csv, write_json

# Table 2 of the paper: imported claim | verdict | certificate | class.
LEDGER = [
    # survivors -- become internal relations
    ("screw generator = R, R^2 = R+I", "forced", "Prop. 2.2", "survivor"),
    ("winding = chi in Z/4Z ; pitch = M >= phi", "forced", "Def. 2.4, Lem. 3.2", "survivor"),
    ("coupling rate = spec(ad_R) = {0, +/-sqrt5}", "forced", "Prop. 2.3", "survivor"),
    ("terrain/rotation = two faces of D=1+4C", "forced", "Thm. 4.2, Prop. 5.1", "survivor"),
    ("closure via the orthogonal partner Kf", "forced", "Prop. 6.1, Thm. 7.2", "survivor"),
    # casualties -- discarded (not-gamma forced)
    ("a helix couples to a copy of itself", "refuted", "Cor. 3.5: [R,R]=0, empty slot", "casualty"),
    ("phi^-4 ~ 0.146 cutoff ('last power > 0.1')", "refuted", "superseded; lambda=2c never uses it", "casualty"),
    ("the algebra is founded on the helix", "refuted", "order keystone->char->flip", "casualty"),
    # residue -- stays posited
    ("(x) = 'intertwine', (+) = 'stack'", "posited", "law forced; reading not adjudicable", "residue"),
    ("the helix as a reading of the flip", "posited", "survives as interpretation, not object", "residue"),
]


def filter_asymmetry():
    """Prop. 8.4: one counterexample refutes; a finite range never promotes."""
    # dashv side: universal 'every a in {0,2} is 0 mod 4', witness a=2
    domain = {0, 2}
    witnesses = [a for a in domain if (a % 4) != 0]
    # non-promotion side: Euler polynomial
    euler = lambda n: n**2 + n + 41
    finite_range_all_prime = all(sp.isprime(euler(n)) for n in range(40))
    return {
        "refutation_from_one_witness": {
            "universal": "every a in {0,2} satisfies a == 0 (mod 4)",
            "counterexample": witnesses,
            "verdict": "dashv (refuted) -- decisive from a single exact instance",
        },
        "finite_range_never_promotes": {
            "universal": "n^2+n+41 is prime for all n>=0",
            "passes_range_0_to_39": bool(finite_range_all_prime),
            "fails_at_n_40": bool(not sp.isprime(euler(40))),
            "demotion_certificate": f"f(40) = {euler(40)} = 41*41",
            "sound_tag_on_passed_range": "computed (never forced)",
        },
    }


def generator_inert():
    """Prop. 8.6: closed {0,2} has empty casualty set -> a pure generator is inert."""
    reachable = {0, 2}
    for _ in range(8):
        reachable |= {(a + b) % 4 for a in reachable for b in reachable}
        reachable |= {(2 * a) % 4 for a in reachable}
    casualties = reachable & {1, 3}
    # a substantive import DOES have a decidable casualty: [R,R]=0 refutes self-coupling
    self_couple_refuted = (hc.R * hc.R - hc.R * hc.R) == sp.zeros(2, 2)
    return {
        "reachable_charges": sorted(reachable),
        "casualty_set_C_dashv": sorted(casualties),
        "generator_is_inert": casualties == set(),
        "substantive_import_has_a_falsifiable_member": bool(self_couple_refuted),
        "example_refutation": "[R,R]=0 refutes 'a helix couples to a copy of itself' (Cor. 3.5)",
    }


def transferable_kernel():
    """Rem. 8.9: sqrt5 = phi + phi^{-1} = phi - psi (the fold's sqrt-branch)."""
    return {
        "kernel_identity": "sqrt5 = phi + phi^{-1} = phi - psi",
        "phi_plus_inv_residual": str(sp.simplify((hc.phi + 1 / hc.phi) - hc.sqrt5)),
        "phi_minus_psi_residual": str(sp.simplify((hc.phi - hc.psi) - hc.sqrt5)),
        "role": "the single differential fact crossing between constructions (2-to-1 fold boundary)",
    }


def ledger_partition():
    """Table 2: the survivors / casualties / residue partition."""
    rows = [
        {"claim": c, "verdict": v, "certificate": cert, "class": cls}
        for (c, v, cert, cls) in LEDGER
    ]
    counts = {}
    for r in rows:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    return rows, counts


def main():
    rows, counts = ledger_partition()
    payload = {
        "section": "8 -- The harness: constraint-import as verification",
        "results": {
            "prop_8_4_filter_asymmetry": filter_asymmetry(),
            "prop_8_6_generator_inert": generator_inert(),
            "rem_8_9_transferable_kernel": transferable_kernel(),
            "table_2_ledger_counts": counts,
        },
    }
    p_json = write_json("harness.json", payload, __file__)
    p_csv = write_csv(
        "harness_ledger.csv",
        ["claim", "verdict", "certificate", "class"],
        rows,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print(f"  ledger partition: {counts} ; asymmetry C_dashv != empty makes the import substantive")


if __name__ == "__main__":
    main()
