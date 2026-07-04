"""
Producer: the Section 7 verification window + Proposition 6.1 (the residual
location).

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/window_objects.csv
              data/2026-06-z5-no-salem-dichotomy/window_summary.json

Brute-force enumeration of every monic integer polynomial in the paper's exact
window -- quartics |c|<=10, quintics |c|<=4, sextics |c|<=3 -- keeping the
charge-EXACTLY-5 squarefree objects (Def. 2.1), each with a high-precision
Mahler measure (Def. 2.3), its recognized closed form, reciprocity, and its
Z-factorization.

Reproduces the Sec. 7 table (corrected 2026-07-04) and Prop. 6.1:
  * charge-5 quartics: distinct M = {1, phi^4}; none in (1,2).
  * 13 non-reciprocal charge-5 objects; realized min M = 2 (x^5-2 and kin);
    all obey M >= mu_S; NONE in [mu_S, 2)  (the empty residual window).
  * 4 reciprocal charge-5 objects; distinct measures {1, phi^2, 2+sqrt3};
    none in (1,2); each = Phi_5 x (real reciprocal unit).
  * OVERALL: no charge-5 object has M in (1,2); the realized floor is 2.

Run: py code/2026-06-z5-no-salem-dichotomy/enumerate_window.py
   (a few-hundred-thousand-poly numpy screen; ~1-2 min.)
"""

import mpmath as mp
import sympy as sp

from z5_core import phi, window_objects
from z5_io import write_csv, write_json

mp.mp.dps = 50

MU_S = mp.findroot(lambda tt: tt**3 - tt - 1, mp.mpf('1.3'))
PHI2 = mp.mpf(str(sp.N(phi**2, 45)))
PHI4 = mp.mpf(str(sp.N(phi**4, 45)))


def _in_open_1_2(M) -> bool:
    return bool(mp.mpf('1') + mp.mpf('1e-9') < M < mp.mpf('2') - mp.mpf('1e-9'))


def summarize(objs):
    recip = [o for o in objs if o["reciprocal"]]
    nonrecip = [o for o in objs if not o["reciprocal"]]
    quartics = [o for o in objs if o["degree"] == 4]

    def distinct_closed_forms(items):
        # stable, sorted by numeric value
        seen = {}
        for o in items:
            seen.setdefault(o["measure_closed_form"], o["mahler"])
        return [k for k, _ in sorted(seen.items(), key=lambda kv: float(kv[1]))]

    Mgt1 = [o["mahler"] for o in objs if o["mahler"] > 1 + mp.mpf('1e-9')]
    realized_floor = min(Mgt1) if Mgt1 else None

    summary = {
        "window": "quartics |c|<=10, quintics |c|<=4, sextics |c|<=3",
        "total_charge5_objects": len(objs),
        "quartics": {
            "count": len(quartics),
            "distinct_measures": distinct_closed_forms(quartics),
            "none_in_open_1_2": not any(_in_open_1_2(o["mahler"]) for o in quartics),
        },
        "non_reciprocal": {
            "count": len(nonrecip),
            "realized_minimum": mp.nstr(min(o["mahler"] for o in nonrecip), 20)
                                 if nonrecip else None,
            "realized_minimum_closed_form":
                min(nonrecip, key=lambda o: o["mahler"])["measure_closed_form"]
                if nonrecip else None,
            "x5_minus_2_present":
                any(tuple(o["coeffs"]) == (1, 0, 0, 0, 0, -2) for o in nonrecip),
            "all_obey_smyth_floor_mu_S":
                all(o["mahler"] >= MU_S - mp.mpf('1e-9') for o in nonrecip),
            "count_in_muS_to_2":
                sum(1 for o in nonrecip
                    if MU_S - mp.mpf('1e-9') <= o["mahler"] < 2 - mp.mpf('1e-9')),
            "residual_window_empty":
                not any(MU_S - mp.mpf('1e-9') <= o["mahler"] < 2 - mp.mpf('1e-9')
                        for o in nonrecip),
        },
        "reciprocal": {
            "count": len(recip),
            "distinct_measures": distinct_closed_forms(recip),
            "none_in_open_1_2": not any(_in_open_1_2(o["mahler"]) for o in recip),
            "objects": [
                {"poly": o["poly"], "factorization": o["factorization"],
                 "mahler_closed_form": o["measure_closed_form"]}
                for o in sorted(recip, key=lambda o: (float(o["mahler"]), o["degree"]))
            ],
        },
        "overall": {
            "count_in_open_1_2": sum(1 for o in objs if _in_open_1_2(o["mahler"])),
            "no_charge5_object_in_open_1_2":
                not any(_in_open_1_2(o["mahler"]) for o in objs),
            "realized_floor": mp.nstr(realized_floor, 20) if realized_floor else None,
        },
        "constants_used": {
            "mu_S": mp.nstr(MU_S, 20),
            "phi^2": mp.nstr(PHI2, 20),
            "phi^4": mp.nstr(PHI4, 20),
        },
    }
    return summary


def main():
    print("enumerating the Sec. 7 window (this takes ~1-2 min)...")
    objs = window_objects()
    print(f"  found {len(objs)} charge-exactly-5 squarefree objects")

    fields = ["degree", "coeffs", "poly", "reciprocal", "mahler_str",
              "measure_closed_form", "factorization", "window_bound"]
    rows = []
    for o in sorted(objs, key=lambda o: (o["degree"], float(o["mahler"]),
                                         tuple(o["coeffs"]))):
        rows.append({
            "degree": o["degree"],
            "coeffs": " ".join(str(c) for c in o["coeffs"]),
            "poly": o["poly"],
            "reciprocal": o["reciprocal"],
            "mahler_str": o["mahler_str"],
            "measure_closed_form": o["measure_closed_form"],
            "factorization": o["factorization"],
            "window_bound": o["window_bound"],
        })
    cpath = write_csv("window_objects.csv", fields, rows, __file__)

    payload = {"section_7_window_summary": summarize(objs)}
    jpath = write_json("window_summary.json", payload, __file__)

    s = payload["section_7_window_summary"]
    print(f"wrote {cpath}")
    print(f"wrote {jpath}")
    print(f"  non-reciprocal: {s['non_reciprocal']['count']} (expect 13); "
          f"reciprocal: {s['reciprocal']['count']} (expect 4)")
    print(f"  overall in (1,2): {s['overall']['count_in_open_1_2']} (expect 0); "
          f"floor {s['overall']['realized_floor']}")


if __name__ == "__main__":
    main()
