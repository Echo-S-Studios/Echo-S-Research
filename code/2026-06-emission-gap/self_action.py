"""Producer: the self-action ad_M as a derivation and its difference spectrum.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces: data/2026-06-emission-gap/self_action.json
Paper results: sqrt5 = sqrt(disc(x^2-x-1)) (Lemma 10.1); the self-action
ad_R = [R, .] on 2x2 matrix space has spectrum {-sqrt5, 0, 0, sqrt5}, the
channel gap, off the unit circle (Lemma 10.1); for any M, ad_M has the
eigenvalue difference set {mu_i - mu_j} (Lemma 10.2); the K-formation's real
eigenvalue difference 2K has minimal polynomial x^4+20x^2-80, signature (2,1),
complex place at modulus 2*beta ~ 4.84 off the circle -- real, in K, not Salem.
"""
import mpmath as mp
import numpy as np
import sympy as sp

import emgap_core as C

SCRIPT = "self_action.py"
x = C.x


def adR_spectrum():
    R = C.companion(C.CATALOG["phi"])            # companion of x^2 - x - 1
    ev = sorted(float(v.real) for v in np.linalg.eigvals(C.ad_operator(R)))
    s5 = float(mp.sqrt(5))
    return {
        "host": "R = companion(x^2 - x - 1)",
        "operator": "ad_R = R (x) I - I (x) R^T on 2x2 matrix space (dim 4)",
        "spectrum": [round(v, 12) + 0.0 for v in ev],  # +0.0 normalizes -0.0 -> 0.0 (platform-stable)
        "expected": sorted([-s5, 0.0, 0.0, s5]),
        "matches_pm_sqrt5_and_0": bool(
            max(abs(a - b) for a, b in zip(ev, sorted([-s5, 0.0, 0.0, s5]))) < 1e-9),
        "off_unit_circle": bool(all(abs(v) < 1e-9 or abs(abs(v) - 1) > 0.5 for v in ev)),
    }


def difference_spectrum_general():
    """ad_M spectrum equals the eigenvalue difference set for a generic matrix."""
    M = np.array([[3.0, 1.0, 0.0], [0.0, -2.0, 4.0], [1.0, 0.0, 5.0]])
    mu = np.linalg.eigvals(M)
    # + 0.0 normalizes -0.0 -> 0.0 before sorting, so the output is platform-stable
    # (numpy/LAPACK yields -0.0 for the zero differences on Linux but +0.0 on Windows).
    diffs = sorted(round(float((a - b).real), 6) + 0.0 for a in mu for b in mu)
    ad_ev = sorted(round(float(v.real), 6) + 0.0 for v in np.linalg.eigvals(C.ad_operator(M)))
    return {
        "M": M.tolist(),
        "difference_set": diffs,
        "ad_M_spectrum": ad_ev,
        "matches": bool(max(abs(a - b) for a, b in zip(ad_ev, diffs)) < 1e-6),
    }


def K_difference_2K():
    """2K (real eigenvalue difference of the K-formation): min poly x^4+20x^2-80."""
    u = sp.symbols("u")
    Kpoly = sp.Poly(C.CATALOG["K"], x)
    sub = sp.expand(Kpoly.as_expr().subs(x, u / 2) * 16)     # x -> u/2, clear denom
    coeffs = [1, 0, 20, 0, -80]
    sig = C.signature_from_minpoly(coeffs)
    rts = C.mp_roots(coeffs)
    cplx = [r for r in rts if abs(r.real) < mp.mpf(10) ** (-20)]
    beta = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)
    return {
        "derivation": "substitute x -> u/2 into x^4+5x^2-5 and clear denominators",
        "min_poly_of_2K": str(sp.Poly(coeffs, x).as_expr()),
        "matches_substitution": bool(sp.expand(sub - (u**4 + 20 * u**2 - 80)) == 0),
        "signature_r1_r2": list(sig),
        "complex_place_modulus": C.s(abs(cplx[0]), 12) if cplx else None,
        "two_beta": C.s(2 * beta, 12),
        "complex_place_off_circle": bool(cplx and abs(abs(cplx[0]) - 2 * beta) < 1e-9),
        "is_salem": C.has_salem_factor(sp.Poly(coeffs, x)),
    }


def main():
    payload = {
        "golden_discriminant": {
            "poly": "x^2 - x - 1", "disc_b2_minus_4ac": (-1) ** 2 - 4 * 1 * (-1),
            "sqrt_disc": C.s(mp.sqrt(5), 20), "note": "channel threshold sqrt(D) = sqrt5",
        },
        "adR_spectrum": adR_spectrum(),
        "difference_spectrum_general": difference_spectrum_general(),
        "K_difference_2K": K_difference_2K(),
    }
    C.write_json("self_action.json", payload, SCRIPT)
    print(f"ad_R spectrum = {payload['adR_spectrum']['spectrum']} "
          f"(matches +/-sqrt5,0: {payload['adR_spectrum']['matches_pm_sqrt5_and_0']})")
    print(f"2K min poly = {payload['K_difference_2K']['min_poly_of_2K']}, "
          f"sig {tuple(payload['K_difference_2K']['signature_r1_r2'])}, "
          f"salem={payload['K_difference_2K']['is_salem']}")


if __name__ == "__main__":
    main()
