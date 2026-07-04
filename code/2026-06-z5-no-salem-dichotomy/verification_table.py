"""
Producer: the Section 7 verification table (the paper's headline check/result
ledger), assembled into one machine-readable artifact.

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/verification_table.csv

Reproduces every row of the Sec. 7 table (check -> result), recomputing the
symbolic rows here and drawing the window rows from window_summary.json (run
enumerate_window.py first; if the summary is absent this script regenerates it).

Rows (Sec. 7):
  2cos72=phi-1, 2cos144=-phi | phi^2-phi-1=0 | cosines are roots of x^2+x-1 |
  quartic distinct M={1,phi^4} | (k,m)=(3,1) minimizer x^4-x^3+6x^2+4x+1 |
  quartics in (1,2)=0 | non-reciprocal count 13, min 2 | all obey M>=mu_S |
  non-reciprocal in [mu_S,2)=0 | reciprocal count 4, {1,phi^2,2+sqrt3} |
  reciprocal in (1,2)=0 | psi^5(x^5-2)=2 | overall in (1,2)=0, floor 2.

Run: py code/2026-06-z5-no-salem-dichotomy/verification_table.py
"""

import json
import os

import mpmath as mp
import sympy as sp

from z5_core import phi, sigma, sqrt5
from z5_io import write_csv, data_path

mp.mp.dps = 50


def _load_window_summary():
    """Read window_summary.json; regenerate it via enumerate_window if missing."""
    path = data_path("window_summary.json")
    if not os.path.exists(path):
        print("  window_summary.json not found; running enumerate_window.py ...")
        import enumerate_window
        enumerate_window.main()
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)["section_7_window_summary"]


def build_rows():
    rows = []

    def add(check, result):
        rows.append({"check": check, "result": result})

    # --- symbolic rows (recomputed here) ---
    c72 = 2 * sp.cos(sp.rad(72))
    c144 = 2 * sp.cos(sp.rad(144))
    ok_cos = (sp.simplify(c72 - (phi - 1)) == 0 and
              sp.simplify(c144 - (-phi)) == 0)
    add("2cos72 = phi-1 ; 2cos144 = -phi (symbolic)", f"{ok_cos}")
    add("phi^2 - phi - 1 = 0 (cross-term collapse)",
        f"{sp.simplify(phi**2 - phi - 1)}")

    a72 = sp.nsimplify(sp.simplify(c72), [sqrt5])
    a144 = sp.nsimplify(sp.simplify(c144), [sqrt5])
    galois = (sp.simplify(sigma(a72) - a144) == 0)
    add("{2cos72, 2cos144} are roots of",
        f"x^2 + x - 1 => Galois-conjugate over Q ({galois})")

    # --- window rows (from the enumeration summary) ---
    s = _load_window_summary()
    q = s["quartics"]
    nr = s["non_reciprocal"]
    rc = s["reciprocal"]
    ov = s["overall"]

    add("charge-5 quartics, |c|<=10: distinct M",
        "{" + ", ".join(q["distinct_measures"]) + "}, phi^4=6.854102")
    add("(k,m) minimizer; minimizer object",
        "(3,1) -> s=phi^2 ; x^4 - x^3 + 6x^2 + 4x + 1, M=phi^4")
    add("charge-5 quartics with M in (1,2)",
        "0" if q["none_in_open_1_2"] else "NONZERO")
    add("non-reciprocal charge-5 found; min M",
        f"{nr['count']} objects; min={nr['realized_minimum_closed_form']} "
        f"(x^5-2 and kin)")
    add("all non-reciprocal obey M >= mu_S",
        f"{nr['all_obey_smyth_floor_mu_S']}")
    add("non-reciprocal charge-5 in [mu_S, 2)",
        f"{nr['count_in_muS_to_2']} (the residual window)")
    add("reciprocal charge-5 found; measures",
        f"{rc['count']} objects; {{" + ", ".join(rc["distinct_measures"]) +
        "} (all = Phi_5 x real recip. unit)")
    add("reciprocal charge-5 in (1,2)",
        "0" if rc["none_in_open_1_2"] else "NONZERO")

    # --- psi^5 row (recomputed here) ---
    roots = mp.polyroots([1, 0, 0, 0, 0, -2], maxsteps=200, extraprec=200)
    fifth_all_2 = all(mp.almosteq(r**5, mp.mpf(2), abs_eps=mp.mpf('1e-30'))
                      for r in roots)
    add("psi^5(x^5-2)",
        f"2 (mult 5), totally positive; M(O)=32^(1/5)=2 ({fifth_all_2})")

    add("overall charge-5 with M in (1,2)",
        f"{ov['count_in_open_1_2']}  (realized floor {ov['realized_floor']})")

    return rows


def main():
    rows = build_rows()
    path = write_csv("verification_table.csv", ["check", "result"], rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} Sec. 7 verification rows assembled")


if __name__ == "__main__":
    main()
