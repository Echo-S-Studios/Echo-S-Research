"""
Producer: the floor wall -- the empty Mahler strip (1, phi).  Section 3.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Re-runs the exhaustive degree-two scan of Lemma 3.1 (P1-GATE-02/03): every monic
integer quadratic x^2 + b x + c with |b|,|c| <= 12 (625 of them), its exact
Mahler measure, and whether it falls in the forbidden open strip (1, phi).  The
scan confirms the strip is empty and the unique minimiser above 1 is phi at
discriminant 5 (x^2 +/- x - 1).  Also emits the operator-preservation arithmetic
of Lemma 3.2 / Remark 3.3: the multiplicative semigroup {phi,2,3,5,phi^4,beta^2}
stays at or above phi^2 under product and squaring (no accumulation at 1).

Outputs:
    data/2026-06-lehmers-box/floor_strip_scan.csv     (all 625 quadratics)
    data/2026-06-lehmers-box/floor_wall_summary.json
Backs: Lem. 3.1 (empty strip, minimiser), Lem. 3.2 / Rem. 3.3 (operator floor).
"""

from __future__ import annotations

import mpmath as mp

from box_io import write_csv, write_json

mp.mp.dps = 50
PHI = (1 + mp.sqrt(5)) / 2
EPS = mp.mpf(10) ** -30
BOUND = 12


def mahler_quadratic(b, c):
    """Exact-enough Mahler measure of monic x^2 + b x + c (integers b,c)."""
    disc = b * b - 4 * c
    if disc >= 0:
        sq = mp.sqrt(disc)
        r1 = (-b + sq) / 2
        r2 = (-b - sq) / 2
        return max(mp.mpf(1), abs(r1)) * max(mp.mpf(1), abs(r2))
    # complex-conjugate pair, both of modulus sqrt(c)  (c > 0 when disc < 0)
    return max(mp.mpf(1), mp.sqrt(c)) ** 2


def scan():
    """Full |b|,|c| <= BOUND scan.  Returns (rows, minimiser_value, argmins)."""
    best = None
    argmins = []
    # first pass: find the minimum measure strictly above 1
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            m = mahler_quadratic(b, c)
            if m > 1 + EPS:
                if best is None or m < best - EPS:
                    best, argmins = m, [(b, c)]
                elif abs(m - best) < EPS:
                    argmins.append((b, c))

    rows = []
    for b in range(-BOUND, BOUND + 1):
        for c in range(-BOUND, BOUND + 1):
            m = mahler_quadratic(b, c)
            in_strip = bool(1 + EPS < m < PHI - EPS)
            is_min = bool(abs(m - best) < EPS) if best is not None else False
            rows.append(dict(
                b=b, c=c, discriminant=b * b - 4 * c,
                mahler_measure=mp.nstr(m, 18),
                in_open_strip_1_phi=in_strip,
                is_minimiser_above_1=is_min))
    return rows, best, argmins


def operator_semigroup_facts():
    """Lem. 3.2 / Rem. 3.3: generators of the emitted Mahler semigroup and the
    closure arithmetic that keeps them at/above phi^2 (no accumulation at 1)."""
    beta2 = (5 + 3 * mp.sqrt(5)) / 2                  # imag-root^2 of K = 5.854...
    gens = {"phi": PHI, "2": mp.mpf(2), "3": mp.mpf(3), "5": mp.mpf(5),
            "phi^4": PHI ** 4, "beta^2": beta2}
    phi2 = PHI ** 2
    all_ge_floor = all(g >= PHI - EPS for g in gens.values())
    products_ge_phi2 = all(a * b >= phi2 - EPS
                           for a in gens.values() for b in gens.values())
    squares_ge_phi2 = all(g ** 2 >= phi2 - EPS for g in gens.values())
    return {
        "generators": {k: mp.nstr(v, 18) for k, v in gens.items()},
        "phi_squared": mp.nstr(phi2, 18),
        "all_generators_at_or_above_floor": all_ge_floor,
        "all_pairwise_products_ge_phi_squared": products_ge_phi2,
        "all_squares_ge_phi_squared": squares_ge_phi2,
        "note": "product (direct sum) and squaring keep the image in "
                "{1} U [phi, inf); no realised value falls in (1, phi).",
    }


def main():
    rows, best, argmins = scan()
    cols = ["b", "c", "discriminant", "mahler_measure",
            "in_open_strip_1_phi", "is_minimiser_above_1"]
    p_csv = write_csv("floor_strip_scan.csv", cols, rows, __file__)

    hits = [(r["b"], r["c"], r["mahler_measure"])
            for r in rows if r["in_open_strip_1_phi"]]
    argmin_discs = sorted(set(b * b - 4 * c for (b, c) in argmins))
    payload = {
        "_description": "Exhaustive degree-two Mahler scan (Lemma 3.1) plus the "
                        "operator floor arithmetic (Lemma 3.2 / Remark 3.3).",
        "scan_bound_abs_b_c": BOUND,
        "n_quadratics_scanned": len(rows),
        "open_strip": "(1, phi)",
        "phi": mp.nstr(PHI, 30),
        "n_measures_in_open_strip_1_phi": len(hits),
        "measures_in_open_strip_1_phi": hits,          # empty: strip is empty
        "minimiser_measure_above_1": mp.nstr(best, 30),
        "minimiser_equals_phi": bool(abs(best - PHI) < mp.mpf(10) ** -25),
        "argmin_polynomials_b_c": [list(bc) for bc in sorted(argmins)],
        "argmin_discriminants": argmin_discs,
        "operator_floor": operator_semigroup_facts(),
    }
    p_json = write_json("floor_wall_summary.json", payload, __file__)

    print("wrote", p_csv, f"({len(rows)} rows)")
    print("wrote", p_json)
    print(f"  quadratics scanned: {len(rows)}")
    print(f"  measures in open strip (1, phi): {len(hits)}")
    print(f"  minimiser above 1: {mp.nstr(best, 20)} (== phi: "
          f"{abs(best - PHI) < mp.mpf(10) ** -25})")
    print(f"  argmins: {sorted(argmins)}  discriminants: {argmin_discs}")


if __name__ == "__main__":
    main()
