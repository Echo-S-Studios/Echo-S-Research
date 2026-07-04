"""Producer: the degree-2 Mahler gap and its preservation under spectral ops.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces:
  data/2026-06-emission-gap/degree2_quadratic_scan.csv   (Lemma 6.1)
  data/2026-06-emission-gap/mahler_gap_summary.json      (Lemma 6.1, 6.2, App. A)
Paper results: a scan of all 625 integer quadratics x^2+bx+c with |b|,|c|<=12
finds NO Mahler measure in (1, phi), and the smallest measure above 1 is exactly
phi (Lemma 6.1); the reciprocal jump x^2-bx+1 (M=1 for |b|<=2, >=phi^2 for
|b|>=3); direct-sum multiplicativity M(p+q)=M(p)M(q) and squaring M -> M^2
(Lemma 6.2); a Kronecker sample whose minimum nontrivial measure is phi^2.
"""
import mpmath as mp

import emgap_core as C

SCRIPT = "mahler_gap.py"
PHI = (1 + mp.sqrt(5)) / 2
MU_S = mp.findroot(lambda z: z**3 - z - 1, mp.mpf("1.324"))


def degree2_scan():
    """All 625 quadratics x^2+bx+c, |b|,|c|<=12, with Mahler measure & gap flag."""
    rows, min_above_one = [], None
    for b in range(-12, 13):
        for c in range(-12, 13):
            m = C.mahler([1, b, c])
            in_gap = bool(1 + mp.mpf("1e-9") < m < PHI - mp.mpf("1e-9"))
            if m > 1 + mp.mpf("1e-9") and (min_above_one is None or m < min_above_one):
                min_above_one = m
            rows.append({"b": b, "c": c, "mahler_measure": C.s(m, 16),
                         "in_open_gap_1_phi": in_gap})
    return rows, min_above_one


def reciprocal_jump():
    """x^2 - b x + 1: M = 1 for |b|<=2, >= phi^2 for |b|>=3."""
    out = []
    for b in range(0, 13):
        m = C.mahler([1, -b, 1])
        out.append({"b": b, "mahler_measure": C.s(m, 16),
                    "equals_one": bool(abs(m - 1) < mp.mpf(10) ** (-25)),
                    "ge_phi_squared": bool(m >= PHI**2 - mp.mpf(10) ** (-20))})
    return out


def kron_sample():
    """Mahler of Kronecker-product spectra over the real quadratic seeds."""
    seeds = {k: C.CATALOG[k] for k in ("phi", "tau", "sqrt2", "sqrt3", "sqrt5")}
    vals, rows = [], []
    for na, ca in seeds.items():
        for nb, cb in seeds.items():
            rA, rB = C.mp_roots(ca), C.mp_roots(cb)
            prod = mp.mpf(1)
            for a in rA:
                for bb in rB:
                    v = abs(a * bb)
                    if v > 1:
                        prod *= v
            in_gap = bool(1 + mp.mpf("1e-9") < prod < PHI - mp.mpf("1e-9"))
            if prod > 1 + mp.mpf("1e-6"):
                vals.append(prod)
            rows.append({"A": na, "B": nb, "mahler_measure": C.s(prod, 16),
                         "in_open_gap_1_phi": in_gap})
    return rows, min(vals)


def main():
    scan_rows, min_above = degree2_scan()
    C.write_csv("degree2_quadratic_scan.csv",
                ["b", "c", "mahler_measure", "in_open_gap_1_phi"], scan_rows, SCRIPT)

    recip = reciprocal_jump()
    kron_rows, kron_min = kron_sample()

    # direct-sum multiplicativity and squaring examples (Lemma 6.2)
    m_phi = C.mahler(C.CATALOG["phi"])
    m_sqrt2 = C.mahler(C.CATALOG["sqrt2"])
    m_oplus = C.mahler([1, -1, -3, 2, 2])            # (x^2-x-1)(x^2-2), spectrum union
    m_phi_sq = C.mahler([1, -3, 1])                  # minpoly of phi^2

    summary = {
        "degree2_scan": {
            "count": len(scan_rows),
            "any_in_open_gap_1_phi": any(r["in_open_gap_1_phi"] for r in scan_rows),
            "min_measure_above_one": C.s(min_above, 20),
            "min_equals_phi": bool(abs(min_above - PHI) < mp.mpf(10) ** (-25)),
            "phi": C.s(PHI, 20), "mu_S": C.s(MU_S, 20),
        },
        "reciprocal_jump_x2_minus_bx_plus_1": recip,
        "oplus_multiplicativity": {
            "identity": "M(p (+) q) = M(p) M(q)",
            "M_phi": C.s(m_phi, 16), "M_sqrt2": C.s(m_sqrt2, 16),
            "M_phi_oplus_sqrt2": C.s(m_oplus, 16),
            "product_M_phi_times_M_sqrt2": C.s(m_phi * m_sqrt2, 16),
            "holds": bool(abs(m_oplus - m_phi * m_sqrt2) < mp.mpf(10) ** (-25)),
            "equals_3p236": bool(abs(m_oplus - PHI * 2) < mp.mpf(10) ** (-20)),
        },
        "squaring_squares_measure": {
            "identity": "M(theta -> theta^2) = M(theta)^2",
            "M_phi_squared_minpoly": "x^2 - 3x + 1",
            "M_of_that": C.s(m_phi_sq, 16), "M_phi_squared": C.s(m_phi**2, 16),
            "holds": bool(abs(m_phi_sq - m_phi**2) < mp.mpf(10) ** (-25)),
        },
        "kron_sample": {
            "count": len(kron_rows),
            "any_in_open_gap_1_phi": any(r["in_open_gap_1_phi"] for r in kron_rows),
            "min_nontrivial_measure": C.s(kron_min, 20),
            "min_equals_phi_squared": bool(abs(kron_min - PHI**2) < mp.mpf(10) ** (-18)),
        },
        "cyclotomic_quadratics_measure_one": {
            poly: C.s(C.mahler(coeffs), 6)
            for poly, coeffs in {"x^2-x+1": [1, -1, 1], "x^2+1": [1, 0, 1],
                                 "x-1": [1, -1]}.items()
        },
    }
    C.write_json("mahler_gap_summary.json", summary, SCRIPT)

    print(f"degree2_quadratic_scan.csv: {len(scan_rows)} quadratics, "
          f"min above 1 = {float(min_above)} (phi={float(PHI)}), "
          f"any in (1,phi): {summary['degree2_scan']['any_in_open_gap_1_phi']}")
    print(f"kron sample min nontrivial = {float(kron_min)} (phi^2={float(PHI**2)})")


if __name__ == "__main__":
    main()
