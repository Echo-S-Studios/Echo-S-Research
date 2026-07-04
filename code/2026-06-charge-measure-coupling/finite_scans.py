"""
Producer: the finite-window "[computed]" scans (Sec 4.6, ledger J, Sec 6 window).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    :
  * data/2026-06-charge-measure-coupling/totally_positive_measures.csv
        Thm 4.6 / ledger G: distinct Mahler measures of totally-positive integer
        polynomials; the (1,2) gap (least value >1 is exactly 2).
  * data/2026-06-charge-measure-coupling/z5_quartic_measures.csv
        ledger J: distinct Mahler measures of charge-Z/5 quartics; only {1, phi^4}.
  * data/2026-06-charge-measure-coupling/z5_reciprocal_window.csv
        Thm 6.4 window: reciprocal charge-Z/5 measures {1, phi^2, 2+sqrt3}.
  * data/2026-06-charge-measure-coupling/scan_summary.json  (gap facts + windows)

A fast double-precision numpy screen is followed by exact 50-digit
re-classification (charge group + Mahler) on the survivors, exactly the
two-stage discipline of Section 7.  Windows are those of the verification suite
(|c|<=6); the paper's stated windows are wider (deg<=7; |c|<=10) and are noted.

Run: py code/2026-06-charge-measure-coupling/finite_scans.py
"""

import mpmath as mp
import numpy as np
import sympy as sp

import cmc_core as core
from cmc_io import write_csv, write_json

mp.mp.dps = 50
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
PHI4 = float(PHI**4)


# --------------------------------------------------------------------------
# fast double-precision screens
# --------------------------------------------------------------------------
def _all_real_positive(coeffs, tol=1e-9):
    for r in np.roots(coeffs):
        if abs(r.imag) > tol * max(1.0, abs(r)) or r.real <= tol:
            return False
    return True


def _mahler_np(coeffs):
    m = abs(coeffs[0])
    for r in np.roots(coeffs):
        if abs(r) > 1.0:
            m *= abs(r)
    return m


def _looks_charge5(coeffs, tol=1e-6):
    for r in np.roots(coeffs):
        t5 = 5 * (np.angle(r) / (2 * np.pi))
        if abs(t5 - round(t5)) > tol:
            return False
    return True


# --------------------------------------------------------------------------
# Thm 4.6 / ledger G : totally-positive (1,2) gap
# --------------------------------------------------------------------------
def totally_positive_scan(K=6, degmax=4):
    from itertools import product
    seen = {}          # rounded measure -> (count, example)
    min_above_one = None
    gap_violation = None
    for deg in range(1, degmax + 1):
        for tail in product(range(-K, K + 1), repeat=deg):
            coeffs = [1] + list(tail)
            if coeffs[-1] == 0 or not _all_real_positive(coeffs):
                continue
            m = _mahler_np(coeffs)
            key = round(m, 6)
            if key not in seen:
                seen[key] = [0, coeffs]
            seen[key][0] += 1
            if m > 1 + 1e-6:
                if min_above_one is None or m < min_above_one:
                    min_above_one = m
                if 1 + 1e-6 < m < 2 - 1e-6:
                    gap_violation = (coeffs, m)
    rows = [{"mahler_measure": f"{k:.6f}", "count": v[0],
             "example_poly_coeffs": " ".join(map(str, v[1]))}
            for k, v in sorted(seen.items())]
    summary = {
        "window": f"degrees 1..{degmax}, |c|<={K} (paper: deg<=7)",
        "distinct_measures": len(rows),
        "least_measure_above_1": round(min_above_one, 6) if min_above_one else None,
        "gap_1_2_empty": gap_violation is None,
        "paper_ref": "Thm 4.6 / ledger G",
    }
    return rows, summary


# --------------------------------------------------------------------------
# ledger J : charge-Z/5 quartics -> {1, phi^4}
# --------------------------------------------------------------------------
def z5_quartic_scan(K=6):
    seen = {}
    gap_violation = None
    for a in range(-K, K + 1):
        for b in range(-K, K + 1):
            for c in range(-K, K + 1):
                for d in range(-K, K + 1):
                    if d == 0:
                        continue
                    coeffs = [1, a, b, c, d]
                    if not _looks_charge5(coeffs):
                        continue
                    if core.charge_group(coeffs, nmax=60) != 5:
                        continue
                    m = core.mahler(coeffs)
                    if abs(m - 1) < mp.mpf(10) ** (-9):
                        label = "1"
                    elif abs(m - PHI**4) < mp.mpf(10) ** (-9):
                        label = "phi^4"
                    else:
                        label = f"OTHER={mp.nstr(m, 12)}"
                    seen.setdefault(label, [0, coeffs])
                    seen[label][0] += 1
                    if 1 + 1e-6 < float(m) < PHI4 - 1e-6:
                        gap_violation = (coeffs, float(m))
    rows = [{"mahler_closed_form": k, "count": v[0],
             "example_poly_coeffs": " ".join(map(str, v[1]))}
            for k, v in sorted(seen.items())]
    summary = {
        "window": f"quartics |c|<={K} (paper: |c|<=10)",
        "distinct_measures": sorted(seen.keys()),
        "only_1_and_phi4": set(seen.keys()) == {"1", "phi^4"},
        "gap_1_phi4_empty": gap_violation is None,
        "paper_ref": "ledger J / Thm 6.5",
    }
    return rows, summary


# --------------------------------------------------------------------------
# Thm 6.4 window : reciprocal charge-Z/5 measures {1, phi^2, 2+sqrt3}
# --------------------------------------------------------------------------
def z5_reciprocal_window(tmax=5):
    _x = sp.symbols("x")
    phi5 = 1 + _x + _x**2 + _x**3 + _x**4
    rows = []
    for t in range(2, tmax + 1):
        poly = sp.Poly(sp.expand(phi5 * (_x**2 - t * _x + 1)), _x)
        coeffs = [int(a) for a in poly.all_coeffs()]
        m = core.mahler(coeffs)
        cg = core.charge_group(coeffs)
        recip = core.is_reciprocal(coeffs)
        if abs(m - 1) < mp.mpf(10) ** (-12):
            cf = "1"
        elif abs(m - PHI**2) < mp.mpf(10) ** (-12):
            cf = "phi^2"
        elif abs(m - (mp.mpf(2) + mp.sqrt(3))) < mp.mpf(10) ** (-12):
            cf = "2+sqrt3"
        else:
            cf = ""
        rows.append({
            "factor_trace_t": t,
            "object": f"Phi_5*(x^2-{t}x+1)",
            "charge_group": f"Z/{cg}" if cg else "none",
            "reciprocal": "yes" if recip else "no",
            "mahler_value": mp.nstr(m, 18),
            "mahler_closed_form": cf,
        })
    return rows


def main():
    tp_rows, tp_summary = totally_positive_scan()
    p1 = write_csv("totally_positive_measures.csv",
                   ["mahler_measure", "count", "example_poly_coeffs"],
                   tp_rows, __file__)

    z5q_rows, z5q_summary = z5_quartic_scan()
    p2 = write_csv("z5_quartic_measures.csv",
                   ["mahler_closed_form", "count", "example_poly_coeffs"],
                   z5q_rows, __file__)

    z5r_rows = z5_reciprocal_window()
    p3 = write_csv("z5_reciprocal_window.csv",
                   ["factor_trace_t", "object", "charge_group", "reciprocal",
                    "mahler_value", "mahler_closed_form"],
                   z5r_rows, __file__)

    p4 = write_json("scan_summary.json", {
        "totally_positive_gap": tp_summary,
        "z5_quartic_sector": z5q_summary,
        "z5_reciprocal_window_measures": sorted(
            {r["mahler_closed_form"] for r in z5r_rows if r["mahler_closed_form"]}),
    }, __file__)

    for p in (p1, p2, p3, p4):
        print(f"wrote {p}")
    print(f"  totally-positive: least measure >1 is "
          f"{tp_summary['least_measure_above_1']}, (1,2) empty="
          f"{tp_summary['gap_1_2_empty']}")
    print(f"  Z/5 quartics: measures {z5q_summary['distinct_measures']}")


if __name__ == "__main__":
    main()
