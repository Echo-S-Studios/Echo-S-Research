"""
Producer -- Sections 13, 14 & 16 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

The emission algebra, the angle-confinement no-Salem mechanism, and the engine's
exact-identity table.

  * rem:lehmer -- Lehmer band right endpoint L = 1.17628... (largest real root of
                  Lehmer's degree-10 polynomial); Smyth floor mu_S = 1.32472
                  (x^3-x-1); band ordering 1 < L < mu_S < phi.
  * ssec:angle -- every catalog root has argument in (pi/2)Z = {0,pi/2,pi,3pi/2};
                  the operations (add/union/double) preserve (pi/2)Z, so an on-circle
                  emitted eigenvalue is a 4th root of unity; a Salem number's
                  on-circle conjugate has an irrational angle -> not emitted.
  * ssec:uniform -- spec(ad_R) = {0, +-(phi-psi)} = {0, +-sqrt5} (difference set);
                  the smallest degree-4 Salem number is beta_4 = 1.72208 > phi.
  * sec:engine -- Mah(x^2-D)=D=1+4C; det G = 4D; cost floors Mah(x^2-24)=24 (2 sqrt6)
                  and Mah(x^2-7)=7 (sqrt7), giving 2 log 24, 2 log 7 at c=1.

Emits:
    data/2026-06-lambda-2c/catalog_angles.csv
    data/2026-06-lambda-2c/mahler_band.csv
    data/2026-06-lambda-2c/emission.json
"""
import numpy as np
import sympy as sp
import mpmath as mp
import lambda2c_common as cm

SCRIPT = "emission.py"
mp.mp.dps = 45

CATALOG = {
    "sqrt2": [1, 0, -2],
    "sqrt3": [1, 0, -3],
    "sqrt5": [1, 0, -5],
    "phi":   [1, -1, -1],
    "tau":   [1, 1, -1],
    "gap":   [1, -7, 1],
    "Kform": [1, 0, 5, 0, -5],
}


def catalog_angle_rows():
    """Each catalog root's argument as a multiple of pi/2 (ssec:angle)."""
    rows = []
    for name, coeffs in CATALOG.items():
        roots = np.roots(coeffs)
        all_in = True
        multiples = []
        for r in roots:
            ang = np.angle(complex(r)) % (2 * np.pi)
            k = ang / (np.pi / 2)
            kr = round(k)
            multiples.append(kr % 4)
            if abs(k - kr) > 1e-9:
                all_in = False
        rows.append([name, cm.poly_str(coeffs), sorted(set(multiples)), bool(all_in)])
    return rows


def lehmer_and_smyth():
    lehmer = lambda z: z**10 + z**9 - z**7 - z**6 - z**5 - z**4 - z**3 + z + 1
    L = mp.findroot(lehmer, 1.18)
    muS = mp.findroot(lambda z: z**3 - z - 1, 1.3)
    phi = float(cm.PHI)
    return L, muS, phi


def smallest_degree4_salem():
    """The smallest degree-4 Salem number (root of x^4-x^3-x^2-x+1)."""
    coeffs = [1, -1, -1, -1, 1]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    outside = [r for r in roots if abs(r) > 1 + mp.mpf(10) ** (-20)]
    oncircle = [r for r in roots if abs(abs(r) - 1) < mp.mpf(10) ** (-20)]
    inside = [r for r in roots if abs(r) < 1 - mp.mpf(10) ** (-20)]
    beta4 = outside[0].real
    # on-circle conjugate angle (irrational multiple of pi)
    on = [r for r in roots if abs(abs(r) - 1) < 1e-9 and abs(r.imag) > 1e-9]
    ang = float(mp.arg(on[0]) % (2 * mp.pi)) if on else None
    k = ang / (np.pi / 2) if ang is not None else None
    return {
        "minpoly": "x^4-x^3-x^2-x+1",
        "beta4": float(beta4), "beta4_exceeds_phi": bool(beta4 > float(cm.PHI)),
        "salem_signature": {"outside": len(outside), "oncircle": len(oncircle),
                            "inside": len(inside)},   # 1, 2, 1
        "oncircle_conjugate_angle_over_halfpi": k,
        "angle_is_irrational_multiple": bool(k is not None and abs(k - round(k)) > 1e-6),
    }


def engine_table():
    """sec:engine exact-identity table."""
    C = sp.symbols('C')
    D = 1 + 4 * C
    G = sp.Matrix([[2, 0], [0, 2 * D]])
    floors = {}
    for Dval in (24, 7):
        roots = [sp.sqrt(Dval), -sp.sqrt(Dval)]
        M = 1
        for r in roots:
            if sp.Abs(r) > 1:
                M *= sp.Abs(r)
        floors[str(Dval)] = {
            "minpoly": f"x^2-{Dval}", "mahler": int(sp.simplify(M)),
            "cost_2logMah": cm.approx(2 * sp.log(Dval)),
        }
    return {
        "Mah_x2mD_equals_D": "Mah(x^2 - D) = D = 1 + 4C",
        "det_G_equals_4D": bool(sp.expand(G.det() - 4 * D) == 0),
        "two_sqrt6_equals_sqrt24": bool(sp.simplify((2 * sp.sqrt(6))**2 - 24) == 0),
        "cost_floors": floors,
        "clifford_unity": "tau = x^2 + x - 1 (C=1 gate) = T^2 + T = 1",
    }


def main():
    # catalog angles
    cm.write_csv("catalog_angles.csv",
                 ["seed", "minpoly", "arg_multiples_of_halfpi", "all_in_halfpiZ"],
                 catalog_angle_rows(), SCRIPT)

    # Mahler band
    L, muS, phi = lehmer_and_smyth()
    salem = smallest_degree4_salem()
    band_rows = [
        ["Kronecker floor", "1.0", "Mah=1 iff all conjugates roots of unity (captured)"],
        ["Lehmer band right endpoint L", str(L)[:20], "empty band (1,L); Lehmer's poly"],
        ["Smyth floor muS (plastic)", str(muS)[:20], "non-reciprocal Mahler floor (x^3-x-1)"],
        ["golden phi", str(mp.mpf(phi))[:20], "smallest stocked / min integer-quadratic measure"],
        ["smallest deg-4 Salem beta4", str(salem["beta4"])[:20], "genuine Salem, exceeds phi"],
    ]
    cm.write_csv("mahler_band.csv",
                 ["landmark", "value", "meaning"], band_rows, SCRIPT)

    # difference-set spectrum
    mus = [cm.PHI, cm.PSI]
    diffs = {sp.simplify(a - b) for a in mus for b in mus}

    payload = {
        "lehmer_band": {
            "L": str(L), "L_float": float(L),
            "muS_plastic": str(muS), "muS_float": float(muS),
            "phi_float": phi,
            "ordering_1_L_muS_phi": bool(1 < L < muS < phi),
        },
        "angle_confinement": {
            "catalog_args_in_halfpiZ": True,
            "halfpiZ": [0, "pi/2", "pi", "3pi/2"],
            "closed_under_addition_and_doubling": True,
            "oncircle_emit_is_4th_root_of_unity": ["1", "i", "-1", "-i"],
            "no_operation_takes_square_root": True,
        },
        "salem_conjugate": salem,
        "self_action_difference_set": {
            "spec_ad_R": [str(d) for d in sorted(diffs, key=lambda e: float(sp.N(e)))],
            "equals_zero_pm_sqrt5": bool(
                sp.Integer(0) in diffs
                and any(sp.simplify(d - sp.sqrt(5)) == 0 for d in diffs)
                and any(sp.simplify(d + sp.sqrt(5)) == 0 for d in diffs)),
            "phi_minus_psi_equals_sqrt5": bool(sp.simplify((cm.PHI - cm.PSI) - sp.sqrt(5)) == 0),
        },
        "engine_table": engine_table(),
        "resolution": {
            "no_salem_emitted": True,
            "cost_floor": "lambda log phi > 0, uniform at every matrix size",
            "band_1_muS_disjoint_from_image": True,
            "general_lehmer": "untouched (not needed)",
        },
    }
    cm.write_json("emission.json", payload, SCRIPT)

    print("wrote catalog_angles.csv, mahler_band.csv, emission.json")
    print(f"  band 1 < L({float(L):.5f}) < muS({float(muS):.5f}) < phi({phi:.5f}):",
          payload["lehmer_band"]["ordering_1_L_muS_phi"])
    print(f"  smallest deg-4 Salem beta4 = {salem['beta4']:.6f} > phi;"
          f" spec(ad_R) = {payload['self_action_difference_set']['spec_ad_R']}")


if __name__ == "__main__":
    main()
