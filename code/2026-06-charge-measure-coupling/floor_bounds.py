"""
Producer: the floor / bound results (Sections 2, 4, 6).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/floor_bounds.json

  * Lem 2.6  real reciprocal-unit bound: an integer-trace pair {r,1/r} has
             trace >= 3 (trace 2 forces r=1), hence r >= phi^2.
  * Thm 6.4  realification bound 2^{1/5}=1.1487 is strictly weaker than the
             Smyth floor mu_S=1.3247 the theorem installs.
  * Prop 6.1 / Thm 6.5  pure-pentagon (degree-four Z/5) sector: the Galois-coupled
             construction expands to the minimiser x^4-x^3+6x^2+4x+1, whose
             +-72 pair has modulus phi^2 and +-144 pair modulus phi^-2, giving
             M = phi^4 (the degree-4 sector floor).
  * Prop 2.2 emission-gap sanity: every admissible construction of the paper has
             M in {1} U [phi, inf).

Run: py code/2026-06-charge-measure-coupling/floor_bounds.py
"""

import mpmath as mp
import sympy as sp

import cmc_core as core
from cmc_io import write_json

mp.mp.dps = 50
_x = sp.symbols("x")
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
TWO_PI = 2 * mp.pi


def reciprocal_unit_bound(tmax: int = 11):
    """Lem 2.6: larger root r of x^2 - t x + 1 for integer trace t."""
    rows = []
    for t in range(2, tmax + 1):
        r = max(sp.solve(_x**2 - t * _x + 1, _x))
        rv = mp.mpf(str(sp.N(r, 45)))
        rows.append({
            "trace": t,
            "larger_root": mp.nstr(rv, 18),
            "ge_phi_squared": bool(rv >= PHI**2 - mp.mpf(10) ** (-30)),
            "is_one": bool(abs(rv - 1) < mp.mpf(10) ** (-30)),
        })
    return {
        "claim": "integer-trace reciprocal unit: trace>=3 => r>=phi^2; trace 2 => r=1",
        "phi_squared": mp.nstr(PHI**2, 18),
        "rows": rows,
        "paper_ref": "Lem 2.6",
    }


def realification_bound():
    two_fifth = mp.mpf(2) ** (mp.mpf(1) / 5)
    muS = mp.findroot(lambda z: z**3 - z - 1, mp.mpf("1.3"))
    return {
        "realification_bound_2_pow_1_5": mp.nstr(two_fifth, 18),
        "smyth_floor_mu_S": mp.nstr(muS, 18),
        "smyth_strictly_stronger": bool(two_fifth < muS),
        "note": "psi^5 O totally positive with M(psi^5 O)=M(O)^5 gives only "
                "M(O)>=2^{1/5}; Smyth upgrades this to mu_S",
        "paper_ref": "Thm 6.4",
    }


def pentagon_sector():
    # Thm 6.5: Galois-coupled construction with s=phi^2 at +-72, t=phi^-2 at +-144
    phi = (1 + sp.sqrt(5)) / 2
    s, t = phi**2, phi**-2
    constructed = sp.expand((_x**2 - s * (phi - 1) * _x + s**2)
                            * (_x**2 + t * phi * _x + t**2))
    target = _x**4 - _x**3 + 6 * _x**2 + 4 * _x + 1
    expands_ok = sp.simplify(constructed - target) == 0

    # root geometry of the minimiser
    rts = core.roots([1, -1, 6, 4, 1])
    outer = [r for r in rts if abs(r) > 1]
    inner = [r for r in rts if abs(r) < 1]
    geom = {
        "outer_pair_modulus": mp.nstr(abs(outer[0]), 18),
        "outer_pair_is_phi2": all(abs(abs(r) - PHI**2) < mp.mpf(10) ** (-25) for r in outer),
        "outer_pair_angle_over_2pi": mp.nstr(abs(mp.arg(outer[0]) / TWO_PI), 12),
        "inner_pair_modulus": mp.nstr(abs(inner[0]), 18),
        "inner_pair_is_phi_minus2": all(abs(abs(r) - PHI**-2) < mp.mpf(10) ** (-25) for r in inner),
        "inner_pair_angle_over_2pi": mp.nstr(abs(mp.arg(inner[0]) / TWO_PI), 12),
    }
    return {
        "minimiser": "x^4-x^3+6x^2+4x+1",
        "construction_expands_to_minimiser": bool(expands_ok),
        "mahler": mp.nstr(core.mahler([1, -1, 6, 4, 1]), 18),
        "mahler_closed_form": "phi^4",
        "phi_fourth": mp.nstr(PHI**4, 18),
        "root_geometry": geom,
        "degree4_sector_floor": "M in {1} U [phi^4, inf)",
        "paper_ref": "Prop 6.1 / Thm 6.5",
    }


def emission_gap_sanity():
    objs = {
        "Phi_5 (x^4+x^3+x^2+x+1)": [1, 1, 1, 1, 1],
        "x^3-2": [1, 0, 0, -2],
        "q_2 = x^4+x^2-1": [1, 0, 1, 0, -1],
        "pentagon quartic": [1, -1, 6, 4, 1],
        "x^5-2": [1, 0, 0, 0, 0, -2],
        "x^2-3x+1": [1, -3, 1],
    }
    rows = []
    all_ok = True
    for name, c in objs.items():
        m = core.mahler(c)
        in_gap = (abs(m - 1) < mp.mpf(10) ** (-20)) or (m >= PHI - mp.mpf(10) ** (-20))
        all_ok = all_ok and in_gap
        rows.append({"object": name, "mahler": mp.nstr(m, 18),
                     "in_{1}U[phi,inf)": bool(in_gap)})
    return {"claim": "admissible => M in {1} U [phi, inf)",
            "all_pass": bool(all_ok), "objects": rows, "paper_ref": "Prop 2.2"}


def main():
    payload = {
        "reciprocal_unit_bound": reciprocal_unit_bound(),
        "realification_bound": realification_bound(),
        "pentagon_degree4_sector": pentagon_sector(),
        "emission_gap_sanity": emission_gap_sanity(),
    }
    path = write_json("floor_bounds.json", payload, __file__)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
