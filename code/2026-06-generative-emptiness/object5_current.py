"""
Producer -- Object V: the conserved current.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Proposition 6.1 (prop:current), Proposition 6.2 (prop:clean).

Builds a two-generation orbit from the seeds {phi, K} and EMITS:
  * data/2026-06-generative-emptiness/object5_orbit.csv -- every orbit object with
    its Mahler measure, whether it stays on the Z/4Z lattice (block-diagonal /
    conservation), whether its measure avoids the forbidden band (1,phi), and
    whether all on-circle roots are roots of unity (clean radial growth);
  * data/2026-06-generative-emptiness/object5_measures.json -- the seven listed
    orbit measures phi, phi^2, phi^4, 46.98, 76.63, 122.99, 8049.92, each rebuilt
    from an explicit operator word with its closed form.

Run:
    py code/2026-06-generative-emptiness/object5_current.py
"""
import csv
import json
import os

import mpmath as mp

import ge_core as G
from ge_core import (phi_seed, K_seed, tensor, sq, dsum, mahler, charges,
                     on_circle_roots, is_root_of_unity, pretty, PHI, mpf_str)

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness",
)
CSV_OUT = os.path.join(DATA_DIR, "object5_orbit.csv")
JSON_OUT = os.path.join(DATA_DIR, "object5_measures.json")

PHI_V = PHI()


def orbit():
    """Two-generation orbit generated from {phi, K} under (x), ( )^2, (+)."""
    seeds = [phi_seed(), K_seed()]
    objs = list(seeds)
    gen1 = []
    for P in seeds:
        gen1.append(sq(P))
        for Q in seeds:
            gen1.append(tensor(P, Q))
            gen1.append(dsum(P, Q))
    objs += gen1
    for P in gen1:
        objs.append(sq(P))
        objs.append(tensor(P, phi_seed()))
    return objs


def in_floor_set(m):
    return abs(m - 1) < mp.mpf(10)**(-20) or m >= PHI_V - mp.mpf(10)**(-20)


def orbit_records():
    seen = {}
    for P in orbit():
        key = pretty(P)
        if key in seen:
            continue
        m = mahler(P)
        ch = charges(P)
        on_lattice = "OFF" not in ch
        oc = on_circle_roots(P)
        roots_of_unity = all(is_root_of_unity(r) for r in oc)
        seen[key] = {
            "poly": key,
            "mahler_measure": mpf_str(m, 20),
            "mahler_rounded": round(float(m), 4),
            "on_Z4_lattice": on_lattice,
            "avoids_forbidden_band": bool(in_floor_set(m)),
            "n_on_circle_roots": len(oc),
            "on_circle_roots_of_unity": bool(roots_of_unity),
        }
    return sorted(seen.values(), key=lambda d: mp.mpf(d["mahler_measure"]))


def listed_measures():
    """The seven printed orbit measures, each from an explicit operator word."""
    phi = phi_seed()
    phi4 = sq(sq(phi))
    words = [
        ("phi", phi, "phi", PHI_V),
        ("phi^2", sq(phi), "phi^2", PHI_V**2),
        ("phi^4", sq(sq(phi)), "phi^4", PHI_V**4),
        ("phi^4 (x) phi^4", tensor(phi4, phi4), "phi^8 = 46.98", PHI_V**8),
        ("phi (x) K", tensor(phi, K_seed()), "76.63 (M of phi (x) K)", None),
        ("sq(phi) (+) [phi^4 (x) phi^4]",
         dsum(sq(phi), tensor(phi4, phi4)), "phi^10 = 122.99", PHI_V**10),
        ("phi (x) phi (x) K",
         tensor(tensor(phi, phi), K_seed()), "8049.92 (M of phi (x) phi (x) K)", None),
    ]
    out = []
    for word, P, closed, val in words:
        m = mahler(P)
        rec = {
            "operator_word": word,
            "poly": pretty(P),
            "mahler_measure": mpf_str(m, 20),
            "mahler_rounded_2dp": round(float(m), 2),
            "closed_form": closed,
        }
        if val is not None:
            rec["matches_closed_form"] = bool(abs(m - val) < mp.mpf(10)**(-18))
        out.append(rec)
    return out


def write_csv(records):
    with open(CSV_OUT, "w", newline="", encoding="utf-8") as f:
        f.write(G.csv_header("object5_current.py") + "\n")
        f.write("# two-generation orbit from seeds {phi, K}; conservation + "
                "clean radial growth (Prop 6.1, 6.2)\n")
        cols = ["poly", "mahler_measure", "mahler_rounded", "on_Z4_lattice",
                "avoids_forbidden_band", "n_on_circle_roots",
                "on_circle_roots_of_unity"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in records:
            w.writerow(r)


def main():
    records = orbit_records()
    write_csv(records)

    measures = listed_measures()
    data = G.provenance("object5_current.py")
    data.update({
        "object": "V -- the conserved current",
        "paper_result": "Proposition 6.1 (block-diagonal), Proposition 6.2 (clean growth)",
        "orbit_n_objects": len(records),
        "block_diagonal_all_on_lattice": all(r["on_Z4_lattice"] for r in records),
        "clean_growth_all_avoid_band": all(r["avoids_forbidden_band"] for r in records),
        "on_circle_all_roots_of_unity":
            all(r["on_circle_roots_of_unity"] for r in records),
        "listed_measures": measures,
        "orbit_csv": "object5_orbit.csv",
    })
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"wrote {CSV_OUT}")
    print(f"wrote {JSON_OUT}")
    print(f"  orbit objects: {len(records)}; all on-lattice: "
          f"{data['block_diagonal_all_on_lattice']}; all avoid band: "
          f"{data['clean_growth_all_avoid_band']}")
    print("  listed measures:", [m["mahler_rounded_2dp"] for m in measures])


if __name__ == "__main__":
    main()
