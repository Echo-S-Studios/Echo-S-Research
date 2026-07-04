"""
Producer: the consolidation ledger (Appendix A, entries A-M).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/consolidation_ledger.csv

Regenerates the paper's master numerical table in exact-ish 50-digit arithmetic.
For every named object it emits, in one row: charge group Z/n (or "none" =
charge-inadmissible), whether all n charges are attained, the Mahler measure
(value + recognised closed form), reciprocity, and admissibility.  Families
(x^n-2 for n=3..7; q_k for k=2..5) are expanded to one row per concrete member.
Tensor entries (D, M) are computed from the pairwise-product spectrum; the Adams
entry (C) reports the psi^k residual charge groups.

Run: py code/2026-06-charge-measure-coupling/consolidation_ledger.py
"""

import mpmath as mp
import sympy as sp

import cmc_core as core
from cmc_io import write_csv

mp.mp.dps = 50
_x = sp.symbols("x")
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))

# recognised closed forms for the Mahler column (value -> label)
_KNOWN = [
    ("1", mp.mpf(1)),
    ("2", mp.mpf(2)),
    ("phi", PHI),
    ("phi^2", PHI**2),
    ("phi^4", PHI**4),
    ("2+sqrt3", mp.mpf(2) + mp.sqrt(3)),
    ("128 = 2^7", mp.mpf(128)),
    ("tau (Lehmer)", mp.mpf("1.17628081825991750654")),
    ("beta_4 (Salem)", mp.mpf("1.72208380573904224502")),
]


def _closed_form(value):
    for label, v in _KNOWN:
        if abs(value - v) < mp.mpf(10) ** (-12):
            return label
    return ""


def _row(entry, obj, coeffs=None, roots=None, lead=1, note="", ref=""):
    """Assemble one ledger row from either a coefficient list or a root list."""
    if roots is None:
        roots = core.roots(coeffs)
        lead = coeffs[0]
    n = core.charge_group_from_roots(roots)
    m = core.mahler_from_roots(roots, lead=lead)
    if n is None:
        cg, allc, adm = "none", "", "no"
    else:
        cg = f"Z/{n}"
        charges = core.charges_from_roots(roots, n)
        allc = "yes" if set(charges) == set(range(n)) else "no"
        adm = "yes"
    recip = "yes" if (coeffs is not None and core.is_reciprocal(coeffs)) else (
        "" if coeffs is None else "no")
    return {
        "entry": entry,
        "object": obj,
        "degree": len(roots),
        "charge_group": cg,
        "all_charges_attained": allc,
        "mahler_value": mp.nstr(m, 20),
        "mahler_closed_form": _closed_form(m),
        "reciprocal": recip,
        "admissible": adm,
        "note": note,
        "paper_ref": ref,
    }


def build_rows():
    rows = []

    # A: x^n - 2, n = 3..7  ->  charge Z/n, all charges, M = 2
    for n in range(3, 8):
        coeffs = [1] + [0] * (n - 1) + [-2]
        rows.append(_row("A", f"x^{n}-2", coeffs,
                         note="realizes Z/n, all charges", ref="Thm 3.2"))

    # B: x^4 - 2  ->  cyclic Z/4 (single generator), not Z/2 x Z/2
    rows.append(_row("B", "x^4-2", [1, 0, 0, 0, -2],
                     note="cyclic Z/4 (charge 1 occurs), not Z/2 x Z/2",
                     ref="ledger B"))

    # C: x^6 - 2  ->  Z/6; psi^3 residual Z/2, psi^2 residual Z/3
    base6 = [1, 0, 0, 0, 0, 0, -2]
    res3 = core.charge_group_from_roots(core.adams_roots(base6, 3))
    res2 = core.charge_group_from_roots(core.adams_roots(base6, 2))
    rows.append(_row("C", "x^6-2", base6,
                     note=f"Z/6; psi^3 residual Z/{res3}, psi^2 residual Z/{res2}",
                     ref="Thm 3.4"))

    # D: (x^3-2) tensor (x^4-2)  ->  charpoly x^12-128, charge Z/12, M=128
    tr = core.tensor_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    cp = core.monic_charpoly_int(tr)
    cp_str = "x^12-128" if cp == [1] + [0] * 11 + [-128] else str(cp)
    rows.append(_row("D", "(x^3-2) (x) (x^4-2)", roots=tr,
                     note=f"charpoly {cp_str}, charge Z/lcm(3,4)=Z/12",
                     ref="Thm 3.5 / ledger D"))

    # E: q_k = x^{2k}+x^k-1, k=2..5  ->  charge Z/2k, all charges, M=phi
    for k in range(2, 6):
        coeffs = [0] * (2 * k + 1)
        coeffs[0], coeffs[k], coeffs[2 * k] = 1, 1, -1
        rows.append(_row("E", f"x^{2*k}+x^{k}-1 (q_{k})", coeffs,
                         note="even floor: charge Z/2k, all charges, M=phi",
                         ref="Thm 4.1"))

    # F: beta_4  ->  M in (phi,2), charge NONE (irrational angle; Salem)
    rows.append(_row("F", "x^4-x^3-x^2-x+1 (beta_4)", [1, -1, -1, -1, 1],
                     note="minimal deg-4 Salem; charge-inadmissible",
                     ref="ledger F / Lem 8.1"))

    # G: (x-1)(x-2)  ->  totally positive, M=2
    rows.append(_row("G", "(x-1)(x-2) = x^2-3x+2", [1, -3, 2],
                     note="totally positive; least M>1 is 2; next value phi^2",
                     ref="ledger G / Thm 4.6"))

    # I: pentagon quartic minimiser  ->  charge Z/5, M=phi^4
    rows.append(_row("I", "x^4-x^3+6x^2+4x+1", [1, -1, 6, 4, 1],
                     note="pentagon minimiser; irreducible, non-reciprocal",
                     ref="Thm 6.5 / ledger I"))

    # K: reciprocal Z/5 witnesses  ->  Phi_5*(x^2-3x+1) M=phi^2; also 2+sqrt3
    phi5_q3 = [int(a) for a in sp.Poly(sp.expand(
        (1 + _x + _x**2 + _x**3 + _x**4) * (_x**2 - 3 * _x + 1)), _x).all_coeffs()]
    rows.append(_row("K", "Phi_5*(x^2-3x+1)", phi5_q3,
                     note="reciprocal Z/5 witness, M=phi^2", ref="ledger K"))
    phi5_q4 = [int(a) for a in sp.Poly(sp.expand(
        (1 + _x + _x**2 + _x**3 + _x**4) * (_x**2 - 4 * _x + 1)), _x).all_coeffs()]
    rows.append(_row("K", "Phi_5*(x^2-4x+1)", phi5_q4,
                     note="reciprocal Z/5 witness, M=2+sqrt3 (erratum window)",
                     ref="ledger K / Thm 6.4 window"))
    rows.append(_row("K", "x^5-2", [1, 0, 0, 0, 0, -2],
                     note="non-reciprocal Z/5; M=2=mu(5)", ref="Thm 6.7 / ledger K"))

    # L: Lehmer commutator polynomial  ->  trace 0, M=tau in (1,phi), charge none
    Lx = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]
    rows.append(_row("L", "L(x)(x-1) = x^11-x^9-x^8+x^3+x^2-1", Lx,
                     note="trace 0 (commutator); M=tau in (1,phi); charge bottom",
                     ref="Prop 8.2 / ledger L"))

    # M: phi tensor phi  ->  charpoly (x+1)^2(x^2-3x+1), M=phi^2 (tropical)
    trM = core.tensor_roots([1, -1, -1], [1, -1, -1])
    cpM = core.monic_charpoly_int(trM)
    rows.append(_row("M", "(x^2-x-1) (x) (x^2-x-1)", roots=trM,
                     note=f"charpoly {cpM}; M=phi^2 tropical (NOT phi^4)",
                     ref="Prop 3.3 / ledger M"))

    return rows


def main():
    rows = build_rows()
    fields = ["entry", "object", "degree", "charge_group", "all_charges_attained",
              "mahler_value", "mahler_closed_form", "reciprocal", "admissible",
              "note", "paper_ref"]
    path = write_csv("consolidation_ledger.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} ledger rows (Appendix A entries A-M)")


if __name__ == "__main__":
    main()
