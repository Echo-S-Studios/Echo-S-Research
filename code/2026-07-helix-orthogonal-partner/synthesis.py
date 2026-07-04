"""
Producer: the helix dictionary and the coupling theorem (Section 7).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/helix_dictionary.csv
              data/2026-07-helix-orthogonal-partner/synthesis.json

Reproduces Table 1 (the [posited] helix reading over the [forced] substrate) as
a machine-readable row-per-anatomy record, and records Theorem 7.2 (terrain and
rotation are the two faces of one flip) as the assembled statement plus the exact
substrate fact that closes it: the disjoint union {0,2} u {1,3} = Z/4Z.

The dictionary text is transcribed from the paper's Table 1; each row also
carries the substrate charge it names so the [forced] content is checkable
against the other data artifacts (keystone.json, kformation.json, partner.json).

Run: py code/2026-07-helix-orthogonal-partner/synthesis.py
"""

import helix_core as hc
from helix_io import write_csv, write_json

# Table 1 of the paper: helix anatomy | dissolved helix (golden strand R) |
# orthogonal partner (Kf) | tag.
DICTIONARY = [
    {
        "helix_anatomy": "advance/pitch (terrain)",
        "golden_strand_R": "phi-growth, M>=phi; ad_R grow eig sqrt5",
        "orthogonal_partner_Kf": "real roots +/-K = +/- 5^{1/4}/phi",
        "tag": "forced",
    },
    {
        "helix_anatomy": "rotation/winding",
        "golden_strand_R": "real-axis flip {0,pi} => chi in {0,2}",
        "orthogonal_partner_Kf": "imaginary +/- i beta at +/- pi/2 => chi in {1,3}",
        "tag": "forced",
    },
    {
        "helix_anatomy": "face of the flip D=1+4C",
        "golden_strand_R": "D>0: real growth",
        "orthogonal_partner_Kf": "D<0: rotation (Kf holds both)",
        "tag": "forced",
    },
    {
        "helix_anatomy": "charge reached",
        "golden_strand_R": "{0,2} (real axis)",
        "orthogonal_partner_Kf": "{1,3} (imag.) => completes Z/4Z",
        "tag": "forced",
    },
    {
        "helix_anatomy": "trace-form metric",
        "golden_strand_R": "Riemannian (definite), Q(sqrt5)",
        "orthogonal_partner_Kf": "Lorentzian (2,1), Q(5^{1/4})",
        "tag": "forced",
    },
    {
        "helix_anatomy": "self-coupling",
        "golden_strand_R": "impossible: [R,R]=0, empty Salem slot",
        "orthogonal_partner_Kf": "the partner, not a copy",
        "tag": "forced",
    },
]


def coupling_theorem():
    """Thm. 7.2: assemble the exact substrate fact {0,2} u {1,3} = Z/4Z."""
    golden = {0, 2}
    partner = {1, 3}
    return {
        "statement": "R realises the D>0 (terrain) face; Kf adjoins the D<0 (rotation) face; "
        "their union realises Z/4Z, and the golden mode's only closure is the orthogonal partner.",
        "golden_axis_charges": sorted(golden),
        "partner_axis_charges": sorted(partner),
        "disjoint": golden.isdisjoint(partner),
        "union_is_full_Z4Z": sorted(golden | partner) == [0, 1, 2, 3],
        "helix_reading_tag": "posited (grounded in Section 5)",
        "algebraic_content_tag": "forced",
    }


def main():
    payload = {
        "section": "7 -- Synthesis: the helix dictionary and the coupling theorem",
        "results": {
            "table_1_helix_dictionary_rows": DICTIONARY,
            "thm_7_2_coupling": coupling_theorem(),
        },
    }
    p_json = write_json("synthesis.json", payload, __file__)
    p_csv = write_csv(
        "helix_dictionary.csv",
        ["helix_anatomy", "golden_strand_R", "orthogonal_partner_Kf", "tag"],
        DICTIONARY,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print("  Table 1 reproduced (6 rows); {0,2} u {1,3} = Z/4Z closes the coupling theorem")


if __name__ == "__main__":
    main()
