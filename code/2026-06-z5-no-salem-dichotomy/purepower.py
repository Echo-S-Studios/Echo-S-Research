"""
Producer: Theorem 5.1 (pure power x^5 - m) and Remark 6.2(1) (the psi^5
realification).

Source paper: papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex
Produces    : data/2026-06-z5-no-salem-dichotomy/purepower_family.csv
              data/2026-06-z5-no-salem-dichotomy/psi5_realification.json

Theorem 5.1: for m >= 2, x^5 - m has charge group Z/5Z, is non-reciprocal, and
M = m; the family realizes M in {1} cup [2, infty), attaining mu(5)=2 at x^5-2.
Remark 6.2(1): psi^5 (fifth powers of the roots) sends a charge-5 object to a
totally positive one with M(psi^5 O) = M(O)^5; for x^5-2 every fifth power is 2
(multiplicity 5), so psi^5 O = (x-2)^5, M = 32, and M(O) = 32^(1/5) = 2, while
the totally-positive (1,2) gap gives only M(O) >= 2^(1/5) = 1.1487.

Run: py code/2026-06-z5-no-salem-dichotomy/purepower.py
"""

import mpmath as mp
import sympy as sp

from z5_core import phi, charge_group, is_reciprocal, mahler, recognize_measure
from z5_io import write_csv, write_json

mp.mp.dps = 50


def family_rows(m_values=range(1, 13)):
    """x^5 - m for m in m_values: charge, M, reciprocity, gap membership."""
    rows = []
    for mm in m_values:
        coeffs = [1, 0, 0, 0, 0, -mm]
        M = mahler(coeffs)
        cg = charge_group(coeffs)
        rows.append({
            "m": mm,
            "polynomial": f"x^5 - {mm}",
            "charge_group": cg,
            "mahler": mp.nstr(M, 30),
            "mahler_closed_form": recognize_measure(M),
            "reciprocal": is_reciprocal(coeffs),
            "in_open_1_2": bool(1 + 1e-12 < float(M) < 2 - 1e-12),
        })
    return rows


def psi5_block():
    """Remark 6.2(1): the psi^5 realification, on x^5-2 and on the minimizer."""
    # x^5 - 2 : roots^5 all equal 2, psi^5 O = (x-2)^5
    roots = mp.polyroots([1, 0, 0, 0, 0, -2], maxsteps=200, extraprec=200)
    fifth = [r**5 for r in roots]
    all_two = all(mp.almosteq(f, mp.mpf(2), abs_eps=mp.mpf('1e-30')) for f in fifth)
    totally_positive = all(abs(mp.im(f)) < mp.mpf('1e-25') and mp.re(f) > 0
                           for f in fifth)
    M_O = mahler([1, 0, 0, 0, 0, -2])
    M_psi = mp.mpf(2) ** 5

    # the degree-4 minimizer: psi^5 image totally positive with M = phi^20
    mroots = mp.polyroots([1, -1, 6, 4, 1], maxsteps=200, extraprec=200)
    M_min = mp.mpf(1)
    for r in mroots:
        if abs(r) > 1:
            M_min *= abs(r)
    Mpsi_min = mp.mpf(1)
    for r in mroots:
        f = r ** 5
        if abs(f) > 1:
            Mpsi_min *= abs(f)
    phi20 = mp.mpf(str(sp.N(phi**20, 45)))

    two15 = mp.root(2, 5)
    return {
        "identity": "M(psi^5 O) = M(O)^5",
        "x5_minus_2": {
            "psi5_image": "(x - 2)^5",
            "fifth_powers_all_equal_2": bool(all_two),
            "totally_positive": bool(totally_positive),
            "M_psi5_O": mp.nstr(M_psi, 20),                  # 32
            "M_O": mp.nstr(M_O, 30),                         # 2
            "M_O_equals_32_pow_1_5": bool(
                mp.almosteq(M_O, mp.root(32, 5), abs_eps=mp.mpf('1e-30'))),
        },
        "minimizer_check": {
            "polynomial": "x^4 - x^3 + 6x^2 + 4x + 1",
            "M_O": mp.nstr(M_min, 30),
            "M_psi5_O": mp.nstr(Mpsi_min, 30),
            "M_psi5_equals_M_pow_5": bool(
                mp.almosteq(Mpsi_min, M_min**5, abs_eps=mp.mpf('1e-20'))),
            "equals_phi20": bool(mp.almosteq(Mpsi_min, phi20, abs_eps=mp.mpf('1e-15'))),
        },
        "realification_bound": {
            "value": "2^(1/5)",
            "value_dps": mp.nstr(two15, 30),
            "weaker_than_2": bool(two15 < 2),
            "reason": "the totally-positive (1,2) gap on psi^5 O yields only "
                      "M(O) >= 2^(1/5); the fifth power is unavoidable, so asking "
                      "psi^5 O to clear 32 is the original claim restated (circular).",
        },
    }


def main():
    fields = ["m", "polynomial", "charge_group", "mahler", "mahler_closed_form",
              "reciprocal", "in_open_1_2"]
    cpath = write_csv("purepower_family.csv", fields, family_rows(), __file__)

    payload = {
        "theorem": "Thm. 5.1 (pure power) + Rem. 6.2(1) (psi^5 realification)",
        "family": "x^5 - m realizes M in {1} cup [2, infty); mu(5)=2 at x^5-2",
        "mu_5": 2,
        "floor_object": "x^5 - 2",
        "psi5_realification": psi5_block(),
    }
    jpath = write_json("psi5_realification.json", payload, __file__)

    print(f"wrote {cpath}")
    print(f"wrote {jpath}")
    print("  x^5-m: charge 5, M=m, non-reciprocal; mu(5)=2 at x^5-2; "
          "M(psi^5 O)=M(O)^5")


if __name__ == "__main__":
    main()
