"""
Producer: the K-formation astride the fold (Section 5).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/kformation.json

Recomputes Proposition 5.1 and eqs. (6)-(8) for the seed quartic
Kf(x) = x^4 + 5x^2 - 5:
  * eq. (6): the inner quadratic y^2+5y-5 (y=x^2) has roots
    y_+ = (-5+3 sqrt5)/2 > 0 and y_- = (-5-3 sqrt5)/2 < 0 straddling zero;
  * eq. (7): the real roots are +/-K with K = sqrt(y_+) = 5^{1/4}/phi < 1
    (the terrain), confirmed two ways -- (5^{1/4}/phi)^2 = y_+ and
    minimal_polynomial(5^{1/4}/phi) = Kf;
  * beta = sqrt(|y_-|) > 1 gives the imaginary roots +/- i beta (the pi/2 turn);
  * eq. (8): Mah(Kf) = beta^2 = (5+3 sqrt5)/2 = phi^2 sqrt5, since |K|<1<beta,
    cross-checked by an independent numeric Mahler product over all four roots.

Run: py code/2026-07-helix-orthogonal-partner/kformation.py
"""

import sympy as sp

import helix_core as hc
from helix_io import write_json


def kformation():
    x = hc.x
    phi, sqrt5 = hc.phi, hc.sqrt5
    Kf = x**4 + 5 * x**2 - 5

    # eq. (6): inner quadratic
    y = sp.symbols("y")
    y_roots = sp.solve(y**2 + 5 * y - 5, y)
    yp = (-5 + 3 * sqrt5) / 2
    ym = (-5 - 3 * sqrt5) / 2

    # eq. (7): the real root K
    K = 5 ** sp.Rational(1, 4) / phi
    K_minpoly = sp.minimal_polynomial(K, x)

    # beta and Mahler
    beta = sp.sqrt(sp.Abs(ym))
    beta2 = (5 + 3 * sqrt5) / 2
    mah_exact = hc.mahler_of_multiset([K, -K, sp.I * beta, -sp.I * beta])
    mah_numeric = hc.mahler_numeric(Kf)
    mah_abs_err = hc.mp.nstr(abs(mah_numeric - hc.mp.mpf(str(sp.N(beta2, 50)))), 3)

    # explicit root split
    roots = sp.solve(Kf, x)
    reals = [r for r in roots if sp.im(r) == 0]
    imags = [r for r in roots if sp.re(r) == 0 and sp.im(r) != 0]

    return {
        "seed_quartic_Kf": "x^4 + 5x^2 - 5",
        "eq6_inner_quadratic": {
            "quadratic": "y^2 + 5y - 5  (y = x^2)",
            "roots_closed_form": ["(-5+3 sqrt5)/2", "(-5-3 sqrt5)/2"],
            "roots_match_solver_residual": str(
                sp.simplify(sp.prod([(y - r) for r in y_roots]) - (y**2 + 5 * y - 5))
            ),
            "y_plus_decimal": hc.dec(yp),
            "y_minus_decimal": hc.dec(ym),
            "straddle_zero": bool(sp.simplify(yp) > 0 and sp.simplify(ym) < 0),
        },
        "eq7_terrain_K": {
            "K_closed_form": "5^{1/4}/phi",
            "K_squared_equals_y_plus_residual": str(sp.simplify(K**2 - yp)),
            "sqrt5_over_phi2_rationalised_residual": str(
                sp.simplify(sqrt5 / phi**2 - (3 * sqrt5 - 5) / 2)
            ),
            "K_minimal_polynomial": str(K_minpoly.as_expr() if hasattr(K_minpoly, "as_expr") else K_minpoly),
            "K_minpoly_is_Kf_residual": str(sp.simplify(sp.Poly(K_minpoly, x).as_expr() - Kf)),
            "K_decimal": hc.dec(K),
            "K_less_than_1": bool(sp.simplify(K - 1) < 0),
        },
        "rotation_beta": {
            "beta_closed_form": "sqrt(|y_-|) = sqrt((5+3 sqrt5)/2)",
            "beta_squared_residual": str(sp.simplify(beta**2 - beta2)),
            "beta_decimal": hc.dec(beta),
            "beta_greater_than_1": bool(sp.simplify(beta - 1) > 0),
            "interpretation": "imaginary roots +/- i beta = a pi/2 turn about the origin",
        },
        "eq8_mahler": {
            "closed_form": "beta^2 = (5+3 sqrt5)/2 = phi^2 sqrt5",
            "beta2_equals_phi2_sqrt5_residual": str(sp.simplify(beta2 - phi**2 * sqrt5)),
            "exact_product_over_roots_residual": str(sp.simplify(mah_exact - beta2)),
            "mahler_decimal": hc.dec(beta2),
            "numeric_cross_check_abs_error": mah_abs_err,
        },
        "root_structure": {
            "num_real_roots": len(reals),
            "num_imaginary_roots": len(imags),
            "real_roots_sum_to_zero_residual": str(sp.simplify(sum(reals))),
            "imaginary_roots_sum_to_zero_residual": str(sp.simplify(sum(imags))),
            "K_lt_1_lt_beta": bool(sp.simplify(K - 1) < 0 and sp.simplify(beta - 1) > 0),
        },
    }


def main():
    payload = {
        "section": "5 -- The K-formation astride the fold",
        "results": {"prop_5_1_kformation": kformation()},
    }
    p_json = write_json("kformation.json", payload, __file__)
    print(f"wrote {p_json}")
    print("  Kf = x^4+5x^2-5 : real +/-K=5^{1/4}/phi, imaginary +/- i beta, Mah = phi^2 sqrt5")


if __name__ == "__main__":
    main()
