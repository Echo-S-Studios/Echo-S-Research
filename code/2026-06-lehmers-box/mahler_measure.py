"""
Producer: Mahler-measure preliminaries of *Lehmer's Box* (Section 2).

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Re-derives and EMITS the four foundational Mahler-measure facts:
  * Def. 2.1  -- the product form and the Jensen-integral form agree
                 (on polynomials with no root on the unit circle);
  * Lem. 2.3(ii) -- multiplicativity  Mah(pq) = Mah(p) Mah(q);
  * Lem. 2.3(iii)-- squaring roots     Mah(p^[2]) = Mah(p)^2;
  * Lem. 2.4  -- Kronecker: cyclotomic polynomials have Mahler measure 1,
                 non-cyclotomic integer polynomials exceed 1.

Output:
    data/2026-06-lehmers-box/mahler_measure_checks.csv
Backs: Def. 2.1, Lem. 2.3, Lem. 2.4 (Kronecker).
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

import box_core as C
from box_io import write_csv

mp.mp.dps = 50
_x = sp.symbols('x')


def _label(coeffs):
    """Human-readable monic polynomial string from highest-degree-first coeffs."""
    return sp.sstr(sp.Poly(coeffs, _x).as_expr())


def _squared_poly(coeffs):
    """Build p^[2] (roots squared) from (-1)^d p(x) p(-x) = p^[2](x^2)."""
    p = sp.Poly(coeffs, _x)
    d = p.degree()
    prod = sp.expand((-1) ** d * p.as_expr() * p.as_expr().subs(_x, -_x))
    prod_poly = sp.Poly(prod, _x)
    sq = {}
    for monom, coeff in prod_poly.terms():
        deg = monom[0]
        sq[deg // 2] = int(coeff)
    maxd = max(sq)
    return [sq.get(maxd - i, 0) for i in range(maxd + 1)]


def compute():
    rows = []

    # -- Def. 2.1: product form == Jensen-integral form (no on-circle roots) ---
    forms = [
        [1, -1, -1],          # x^2 - x - 1   (phi)
        [1, 0, -1, -1],       # x^3 - x - 1   (mu_S)
        [1, 0, -2],           # x^2 - 2
        [1, 0, -5],           # x^2 - 5
        [1, -7, 1],           # gap
    ]
    for coeffs in forms:
        prod = C.mahler_measure(coeffs)
        integ = C.mahler_measure_jensen(coeffs)
        rows.append(dict(
            check="def2.1_forms_agree", polynomial=_label(coeffs),
            lhs=mp.nstr(prod, 20), rhs=mp.nstr(integ, 20),
            abs_diff=mp.nstr(abs(prod - integ), 4),
            agree=bool(abs(prod - integ) < mp.mpf(10) ** -12)))

    # -- Lem. 2.3(ii): multiplicativity ---------------------------------------
    pairs = [([1, -1, -1], [1, 0, -2]),
             ([1, -7, 1], [1, -1, -1, -1, 1]),
             ([1, 0, -1, -1], [1, 1, -1])]
    for p, q in pairs:
        pq = sp.Poly(p, _x) * sp.Poly(q, _x)
        pq_coeffs = [int(c) for c in pq.all_coeffs()]
        lhs = C.mahler_measure(pq_coeffs)
        rhs = C.mahler_measure(p) * C.mahler_measure(q)
        rows.append(dict(
            check="lem2.3ii_multiplicative",
            polynomial=f"({_label(p)})*({_label(q)})",
            lhs=mp.nstr(lhs, 20), rhs=mp.nstr(rhs, 20),
            abs_diff=mp.nstr(abs(lhs - rhs), 4),
            agree=bool(abs(lhs - rhs) < mp.mpf(10) ** -12)))

    # -- Lem. 2.3(iii): squaring roots ----------------------------------------
    for coeffs in ([1, -1, -1], [1, 0, -2], [1, -1, -1, -1, 1], [1, 0, -1, -1]):
        squared = _squared_poly(coeffs)
        lhs = C.mahler_measure(squared)
        rhs = C.mahler_measure(coeffs) ** 2
        rows.append(dict(
            check="lem2.3iii_squaring", polynomial=_label(coeffs) + " -> p^[2]",
            lhs=mp.nstr(lhs, 20), rhs=mp.nstr(rhs, 20),
            abs_diff=mp.nstr(abs(lhs - rhs), 4),
            agree=bool(abs(lhs - rhs) < mp.mpf(10) ** -12)))

    # -- Lem. 2.4: Kronecker (cyclotomic => Mah = 1) --------------------------
    for n in range(1, 16):
        coeffs = [int(c) for c in sp.Poly(sp.cyclotomic_poly(n, _x), _x).all_coeffs()]
        m = C.mahler_measure(coeffs)
        rows.append(dict(
            check="lem2.4_kronecker_cyclotomic", polynomial=f"Phi_{n}",
            lhs=mp.nstr(m, 20), rhs="1",
            abs_diff=mp.nstr(abs(m - 1), 4),
            agree=bool(abs(m - 1) < mp.mpf(10) ** -12)))
    # non-cyclotomic sanity: x^2 - x - 1 has measure > 1
    m = C.mahler_measure([1, -1, -1])
    rows.append(dict(
        check="lem2.4_kronecker_noncyclotomic", polynomial="x^2 - x - 1",
        lhs=mp.nstr(m, 20), rhs=">1",
        abs_diff=mp.nstr(m - 1, 6), agree=bool(m > 1 + mp.mpf(10) ** -6)))

    return rows


def main():
    rows = compute()
    cols = ["check", "polynomial", "lhs", "rhs", "abs_diff", "agree"]
    path = write_csv("mahler_measure_checks.csv", cols, rows, __file__)
    n_ok = sum(1 for r in rows if r["agree"])
    print("wrote", path)
    print(f"  {n_ok}/{len(rows)} Mahler-measure identity checks agree")


if __name__ == "__main__":
    main()
