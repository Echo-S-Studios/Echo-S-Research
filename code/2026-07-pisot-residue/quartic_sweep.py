"""Producer: the Pisot-quartic sweep and emission-gap probe (Section 7.2).

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the Section 7.2 sweep over the box [-3,3]^4 (2401 monic quartics).
It certifies 103 Pisot quartics = 102 complex-pair + 1 totally real.  The unique
totally-real hit is x^4-3x^3-2x^2+2x+1 (theta ~ 3.390 > phi), and the smallest
(minimal-theta) complex-pair instance is x^4-x^3-1 (theta ~ 1.3803).  The
totally-real subfamily satisfies the emission-gap prediction P4 theta >= phi.

Emits:
  data/2026-07-pisot-residue/quartic_sweep_summary.json
  data/2026-07-pisot-residue/quartic_pisots.csv       (the 103 certified Pisots)
"""
from sympy import symbols, Poly

from pisot_lib import is_pisot, dominant_root, phi, poly_str, write_json, write_csv

x = symbols('x')


def run_sweep(dps=40):
    """Certify Pisot quartics across [-3,3]^4 and record each hit."""
    pisots = []                        # (coeffs, theta, kind, n_real_in, n_pairs)
    totally_real = complex_pair = 0
    tr_hit = None
    min_cp = None
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                for d in range(-3, 4):
                    coeffs = [1, a, b, c, d]
                    if d == 0:
                        continue
                    if not Poly(coeffs, x).is_irreducible:
                        continue
                    ok, n_real_in, n_pairs = is_pisot(coeffs, dps)
                    if not ok:
                        continue
                    th = float(dominant_root(coeffs, dps))
                    if n_pairs == 0:
                        totally_real += 1
                        tr_hit = coeffs
                        kind = "totally_real"
                    else:
                        complex_pair += 1
                        kind = "complex_pair"
                        if min_cp is None or th < min_cp[0]:
                            min_cp = (th, coeffs)
                    pisots.append((coeffs, th, kind, n_real_in, n_pairs))
    return pisots, totally_real, complex_pair, tr_hit, min_cp


def main():
    pisots, tr, cp, tr_hit, min_cp = run_sweep()
    with __import__("mpmath").workdps(50):
        phi_val = float(phi(50))
        tr_theta = float(dominant_root(tr_hit, 50)) if tr_hit else None

    summary = {
        "description": "Section 7.2 Pisot-quartic sweep over [-3,3]^4: 103 certified "
                       "(102 complex-pair + 1 totally real) and the emission-gap probe P4.",
        "box": {"coeff_range": [-3, 3], "num_free_coeffs": 4, "size": 7 ** 4},
        "certified_pisot": len(pisots),
        "totally_real": tr,
        "complex_pair": cp,
        "split_sums": tr + cp == len(pisots),
        "all_scans_phi1_4_predicted": "{Phi_1^4} on all 103 by Theorem 3.1 (modulus pinning)",
        "unique_totally_real_hit": {
            "poly": poly_str(tr_hit),
            "coeffs_hi_to_lo": tr_hit,
            "theta_display": tr_theta,
            "phi_display": phi_val,
            "theta_greater_than_phi": tr_theta > phi_val,
            "emission_gap_not_falsified": tr_theta >= phi_val,
        },
        "smallest_complex_pair_instance": {
            "poly": poly_str(min_cp[1]),
            "coeffs_hi_to_lo": min_cp[1],
            "theta_display": float(min_cp[0]),
        },
    }
    jpath = write_json("quartic_sweep_summary.json", summary, "quartic_sweep.py")

    rows = [("[" + " ".join(str(int(v)) for v in c) + "]",
             poly_str(c).replace(",", ";"), f"{th:.6f}", kind, nri, npairs)
            for (c, th, kind, nri, npairs) in sorted(pisots, key=lambda t: t[1])]
    cpath = write_csv(
        "quartic_pisots.csv",
        ["coeffs_hi_to_lo", "poly", "theta_display", "kind", "n_real_inside",
         "n_nonreal_pairs"],
        rows, "quartic_sweep.py")

    print(f"wrote {jpath}")
    print(f"wrote {cpath}")
    print(f"  certified Pisot = {len(pisots)} = {tr} totally-real + {cp} complex-pair")
    print(f"  totally-real hit: {poly_str(tr_hit)} theta~{tr_theta:.4f} (> phi~{phi_val:.4f})")
    print(f"  smallest complex-pair: {poly_str(min_cp[1])} theta~{float(min_cp[0]):.4f}")


if __name__ == "__main__":
    main()
