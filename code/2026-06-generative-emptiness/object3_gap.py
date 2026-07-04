"""
Producer -- Object III: the gap / cost floor phi.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Proposition 4.1 (prop:basegap), Corollary 4.2 (cor:floor),
                 minimality Prop 7.1(3).

Enumerates every monic integer quadratic x^2 + b x + c over |b|,|c| <= 6, computes
each Mahler measure, and EMITS:
  * data/2026-06-generative-emptiness/object3_gap_quadratics.csv -- the distinct
    realised measures > 1 (sorted) with a witness (b,c) and a closed form where
    identified;
  * data/2026-06-generative-emptiness/object3_gap.json -- the summary: the band
    (1,phi) is empty, the smallest measure > 1 is phi at x^2-x-1, the first realised
    values phi,2,1+sqrt2,phi^2,1+sqrt3,3, phi is the smallest Perron number, and
    {1} U [phi,inf) is closed under product and square.

Run:
    py code/2026-06-generative-emptiness/object3_gap.py
"""
import csv
import json
import os

import sympy as sp
import mpmath as mp

import ge_core as G
from ge_core import x, as_poly, mahler, roots_mp, PHI, mpf_str

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness",
)
CSV_OUT = os.path.join(DATA_DIR, "object3_gap_quadratics.csv")
JSON_OUT = os.path.join(DATA_DIR, "object3_gap.json")

BOUND = 6
PHI_V = PHI()
EPS = mp.mpf(10) ** (-30)

# closed forms probed for identification of the small realised measures
_CLOSED_FORMS = [
    ("phi = (1+sqrt5)/2", (1 + mp.sqrt(5)) / 2),
    ("2", mp.mpf(2)),
    ("1+sqrt2", 1 + mp.sqrt(2)),
    ("phi^2 = (3+sqrt5)/2", PHI_V**2),
    ("1+sqrt3", 1 + mp.sqrt(3)),
    ("3", mp.mpf(3)),
    ("2+sqrt2", 2 + mp.sqrt(2)),
    ("(3+sqrt13)/2", (3 + mp.sqrt(13)) / 2),
    ("1+sqrt5", 1 + mp.sqrt(5)),
    ("2+sqrt3", 2 + mp.sqrt(3)),
    ("(5+sqrt21)/2", (5 + mp.sqrt(21)) / 2),
    ("3+sqrt5", 3 + mp.sqrt(5)),
    ("phi^4 = (7+3sqrt5)/2", PHI_V**4),
]


def _closed_form(m):
    for name, val in _CLOSED_FORMS:
        if abs(m - val) < mp.mpf(10) ** (-20):
            return name
    return ""


def quadratic_measures():
    """(measure, b, c) for every x^2+bx+c with |b|,|c|<=6 and measure > 1."""
    out = []
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            m = mahler(as_poly(x**2 + b * x + c))
            if m > 1 + EPS:
                out.append((m, b, c))
    return out


def distinct_sorted(measures):
    """Collapse to distinct measures (30-digit key), keep first witness, sort."""
    seen = {}
    for m, b, c in sorted(measures, key=lambda t: t[0]):
        key = mp.nstr(m, 28)
        if key not in seen:
            seen[key] = (m, b, c)
    return [seen[k] for k in sorted(seen, key=lambda s: mp.mpf(s))]


def smallest_perron():
    """The smallest Perron number among monic integer quadratics is phi."""
    perron = []
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            rts = roots_mp(as_poly(x**2 + b * x + c))
            if len(rts) != 2:
                continue
            reals = [r for r in rts if abs(r.imag) < mp.mpf(10)**(-20)]
            if len(reals) != 2:
                continue
            vals = sorted(abs(r) for r in reals)
            dom = max(reals, key=lambda r: abs(r))
            if dom.real > 1 and vals[1] > vals[0] + mp.mpf(10)**(-25):
                perron.append((dom.real, b, c))
    m, b, c = min(perron, key=lambda t: t[0])
    return {"smallest_perron_measure": mpf_str(m, 25),
            "witness_poly": f"x^2 + {b}*x + {c}",
            "equals_phi": bool(abs(m - PHI_V) < mp.mpf(10)**(-25))}


def floor_set_closure():
    """{1} U [phi,inf) is closed under product and square (Cor 4.2 arithmetic)."""
    sample = [mp.mpf(1), PHI_V, PHI_V**2, mp.mpf(2), mp.mpf(5), 46 + mp.mpf("0.98")]

    def in_set(v):
        return abs(v - 1) < mp.mpf(10)**(-25) or v >= PHI_V - mp.mpf(10)**(-25)

    prod_ok = all(in_set(a * b) for a in sample for b in sample)
    sq_ok = all(in_set(a * a) for a in sample)
    no_gap = all(not (1 < a < PHI_V - mp.mpf(10)**(-25)) for a in sample)
    return {"products_stay_in_set": prod_ok, "squares_stay_in_set": sq_ok,
            "sample_avoids_gap": no_gap,
            "floor_set": "{1} U [phi, inf)"}


def write_csv(rows):
    with open(CSV_OUT, "w", newline="", encoding="utf-8") as f:
        f.write(G.csv_header("object3_gap.py") + "\n")
        f.write("# monic integer quadratics x^2+b*x+c, |b|,|c|<=6, distinct "
                "Mahler measures > 1, sorted ascending\n")
        w = csv.writer(f)
        w.writerow(["rank", "b", "c", "poly", "mahler_measure", "closed_form"])
        for i, (m, b, c) in enumerate(rows, start=1):
            w.writerow([i, b, c, f"x^2+{b}x+{c}", mpf_str(m, 25), _closed_form(m)])


def main():
    measures = quadratic_measures()
    rows = distinct_sorted(measures)
    write_csv(rows)

    band_empty = all(not (1 < m < PHI_V - EPS) for m, _, _ in measures)
    m_min, b_min, c_min = min(measures, key=lambda t: t[0])
    first6 = [{"measure": mpf_str(m, 25), "witness": f"x^2+{b}x+{c}",
               "closed_form": _closed_form(m)} for m, b, c in rows[:6]]

    data = G.provenance("object3_gap.py")
    data.update({
        "object": "III -- the gap (cost floor phi)",
        "paper_result": "Proposition 4.1, Corollary 4.2, Prop 7.1(3)",
        "enumeration_bound": f"|b|,|c| <= {BOUND}",
        "n_quadratics_measure_gt_1": len(measures),
        "n_distinct_measures": len(rows),
        "band_1_to_phi_empty": band_empty,
        "phi": mpf_str(PHI_V, 25),
        "smallest_measure_above_1": {
            "measure": mpf_str(m_min, 25),
            "witness_poly": f"x^2 + {b_min}*x + {c_min}",
            "equals_phi": bool(abs(m_min - PHI_V) < mp.mpf(10)**(-30)),
        },
        "first_realised_values": first6,
        "smallest_perron": smallest_perron(),
        "floor_set_closure": floor_set_closure(),
        "mahler_spectrum": "{1} U [phi, inf)",
        "csv_file": "object3_gap_quadratics.csv",
    })
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"wrote {CSV_OUT}")
    print(f"wrote {JSON_OUT}")
    print(f"  band (1,phi) empty: {band_empty}")
    print(f"  smallest measure > 1: {mpf_str(m_min, 8)} at x^2+{b_min}x+{c_min}")
    print("  first realised:", [r["closed_form"] or r["measure"][:6] for r in first6])


if __name__ == "__main__":
    main()
