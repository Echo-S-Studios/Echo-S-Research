"""
Producer -- Section 12 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The K-formation seed f(x) = x^4 + 5x^2 - 5 straddles the fold:

  * prop:kform -- as g(y) = y^2 + 5y - 5 in y = x^2, the two y-roots straddle zero:
                  y_+ = (-5+3 sqrt5)/2 ~ 0.8541 > 0, y_- = (-5-3 sqrt5)/2 ~ -5.8541 < 0.
                  Real roots +-K with K = sqrt(y_+) = 5^{1/4}/phi ~ 0.9242 (inside the
                  unit circle); imaginary roots +- i beta with beta = sqrt|y_-| ~ 2.4195
                  (outside); Mah(f) = beta^2 = (5+3 sqrt5)/2 ~ 5.8541.
  * ssec:kformface -- the one Lorentzian field's complex place is off the circle:
                  |i beta| = 2.4195 and |5^{1/4} i| = 1.4953, both >> 1.

Emits:
    data/2026-06-lambda-2c/kformation.json
"""
import sympy as sp
import mpmath as mp
import lambda2c_common as cm

SCRIPT = "kformation.py"
x, y = sp.symbols('x y')


def main():
    ys = sp.solve(y**2 + 5 * y - 5, y)
    yp = (-5 + 3 * sp.sqrt(5)) / 2
    ym = (-5 - 3 * sp.sqrt(5)) / 2

    K = sp.sqrt(yp)
    K_alt = 5**sp.Rational(1, 4) / cm.PHI
    beta2 = (5 + 3 * sp.sqrt(5)) / 2
    beta = sp.sqrt(beta2)

    # independent Mahler measure of the quartic
    M = cm.mahler_mpmath([1, 0, 5, 0, -5])

    # complex place magnitudes
    fifth = mp.power(5, mp.mpf(1) / 4)
    beta_num = mp.sqrt(mp.mpf(str(sp.N(beta2, 40))))

    payload = {
        "quartic": "x^4 + 5x^2 - 5",
        "reduced_quadratic": "y^2 + 5y - 5  (y = x^2)",
        "y_roots": {
            "y_plus": str(sp.simplify(yp)), "y_plus_float": cm.approx(yp),   # 0.8541...
            "y_minus": str(sp.simplify(ym)), "y_minus_float": cm.approx(ym),  # -5.8541...
            "straddle_zero": bool(sp.N(yp) > 0 and sp.N(ym) < 0),
            "match_solved": all(
                any(sp.simplify(r - target) == 0 for r in ys) for target in (yp, ym)),
        },
        "real_roots": {
            "K": "5^(1/4)/phi", "K_float": cm.approx(K_alt),                 # 0.9242...
            "K_equals_sqrt_yplus": bool(sp.simplify(K**2 - K_alt**2) == 0),
            "inside_unit_circle": bool(cm.approx(K_alt) < 1),
        },
        "imaginary_roots": {
            "beta": str(sp.simplify(beta2)) + " under sqrt", "beta_float": cm.approx(beta),  # 2.4195
            "beta_equals_sqrt_absyminus": bool(sp.simplify(beta2 - (-ym)) == 0),
        },
        "mahler": {
            "mahler_exact": "(5+3 sqrt5)/2 = beta^2", "mahler_float": float(M),  # 5.8541
            "independent_mahler_matches_beta2": bool(
                abs(M - mp.mpf(str(sp.N(beta2, 40)))) < mp.mpf(10) ** (-20)),
        },
        "complex_place_off_circle": {
            "abs_i_beta": float(abs(mp.mpc(0, 1) * beta_num)),               # 2.4195
            "abs_fifthroot_i": float(abs(mp.mpc(0, 1) * fifth)),            # 1.4953
            "both_off_unit_circle": bool(beta_num > 1 and fifth > 1),
            "interpretation": "Lorentzian without Salem: complex embedding far from the circle",
        },
    }
    cm.write_json("kformation.json", payload, SCRIPT)
    print("wrote kformation.json")
    print(f"  K = {payload['real_roots']['K_float']:.6f} (inside), "
          f"beta = {payload['imaginary_roots']['beta_float']:.6f} (outside), "
          f"Mah = {payload['mahler']['mahler_float']:.6f}")


if __name__ == "__main__":
    main()
