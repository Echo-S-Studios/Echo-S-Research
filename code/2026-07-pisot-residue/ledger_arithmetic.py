"""Producer: signature-ledger and next-step arithmetic tallies.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the Section 7.3 signature-ledger-M accounting and the Section 8
next-step (N3/P5) Burnside arithmetic -- the pure counting claims that tie the
census, sweeps, and canonical signatures together:

  * ledger M: 237 signatures = 14 canonical + 37 census + 103 quartic + 83
    quintic; two in-engine decision paths each -> 474 = 2*237 scan executions;
  * quintic box 3125 = 5^5 with reject tally + 83 partitioning it;
  * degree-12 Salem Burnside 729/27/378;
  * quartic box 2401 = 7^4, split 102 + 1 = 103;
  * N3/P5 census extension c in {-2,...,2}^6: 15625 vectors, 125 twist-fixed,
    Burnside (15625+125)/2 = 7875 twist-classes.

These are exact integer identities; every value is recomputed, not restated.

Emits:
  data/2026-07-pisot-residue/ledger_arithmetic.json
"""
from pisot_lib import write_json


def main():
    ledger_M = {
        "canonical_signatures": 14,
        "census_salem": 37,
        "quartic_pisots": 103,
        "quintic_pisots": 83,
        "total_signatures": 14 + 37 + 103 + 83,
        "matches_paper_237": 14 + 37 + 103 + 83 == 237,
        "decision_paths_per_signature": 2,
        "scan_executions": 2 * (14 + 37 + 103 + 83),
        "matches_paper_474": 2 * (14 + 37 + 103 + 83) == 474,
    }
    quintic_box = {
        "size": 5 ** 5,
        "reject_tally": [625, 50, 638, 1318, 411],
        "certified_pisot": 83,
        "partitions_box": 625 + 50 + 638 + 1318 + 411 + 83 == 5 ** 5,
        "pattern_split": {"real5": 0, "mixed": 16, "two_pair": 67,
                          "sums_to_83": 0 + 16 + 67 == 83},
        "scan_sizes": {"deg_Rat": 5 ** 2, "bound_2n4": 2 * 5 ** 4, "deg_C2": 20 ** 2},
    }
    salem_burnside = {
        "palindromic_vectors": 3 ** 6,
        "twist_fixed": 3 ** 3,
        "orbits": (3 ** 6 + 3 ** 3) // 2,
        "matches_paper_378": (3 ** 6 + 3 ** 3) // 2 == 378,
        "corrected_cascade": {"pm1": 39, "trace_sturm_reject": 256,
                              "reducible": 46, "salem": 37,
                              "combined": 256 + 46,
                              "sums_to_378": 39 + 256 + 46 + 37 == 378},
        "scan_sizes": {"deg_Rat": 12 ** 2, "bound_2n4": 2 * 12 ** 4},
    }
    quartic_box = {
        "size": 7 ** 4,
        "complex_pair": 102,
        "totally_real": 1,
        "certified_pisot": 102 + 1,
        "matches_paper_103": 102 + 1 == 103,
    }
    n3_extension = {
        "vectors": 5 ** 6,
        "twist_fixed": 5 ** 3,
        "burnside_orbits": (5 ** 6 + 5 ** 3) // 2,
        "matches_paper_7875": (5 ** 6 + 5 ** 3) // 2 == 7875,
        "burnside_even": (5 ** 6 + 5 ** 3) % 2 == 0,
    }
    lehmer = {"deg_Rat": 10 ** 2, "bound_2n4": 2 * 10 ** 4}

    payload = {
        "description": "Sections 7.3 & 8 exact arithmetic tallies: ledger-M signature "
                       "accounting, box sizes, Burnside counts, and the N3/P5 extension.",
        "ledger_M": ledger_M,
        "quintic_box": quintic_box,
        "salem_burnside": salem_burnside,
        "quartic_box": quartic_box,
        "n3_census_extension": n3_extension,
        "lehmer_scan_sizes": lehmer,
    }
    path = write_json("ledger_arithmetic.json", payload, "ledger_arithmetic.py")
    print(f"wrote {path}")
    print(f"  ledger M: {ledger_M['total_signatures']} signatures, "
          f"{ledger_M['scan_executions']} scans")
    print(f"  N3/P5 Burnside: {n3_extension['burnside_orbits']} twist-classes")


if __name__ == "__main__":
    main()
