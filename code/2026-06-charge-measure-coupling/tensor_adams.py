"""
Producer: tensor and Adams operation laws (Section 2.1 backbone table,
Prop 3.3, Thm 3.4).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/tensor_adams.json

Computes, from pairwise-product / power spectra:
  * Prop 3.3 (ledger M): phi (x) phi -> eigenvalues {phi^2,-1,-1,phi'^2},
    charpoly (x+1)^2(x^2-3x+1), tropical Mahler = phi^2 (NOT the retracted
    factored value phi^4).
  * Thm 3.5 (ledger D): (x^3-2)(x)(x^4-2) -> charpoly x^12-128, charge Z/12,
    off-circle so the factored law happens to agree at 2^7.
  * Adams measure law M(psi^k A) = M(A)^k (Sec 2.1 table).
  * Thm 3.4 (ledger C): x^6-2 primary decomposition, psi^3 -> Z/2, psi^2 -> Z/3.
  * Sumset charge rule under tensor (Sec 2.1 table): composite charges equal
    { (k*a + m*b) mod lcm } for the coprime example.

Run: py code/2026-06-charge-measure-coupling/tensor_adams.py
"""

import mpmath as mp
import sympy as sp

import cmc_core as core
from cmc_io import write_json

mp.mp.dps = 50
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))


def phi_tensor_phi():
    rts = core.tensor_roots([1, -1, -1], [1, -1, -1])
    charpoly = core.monic_charpoly_int(rts)
    m_tropical = core.mahler_from_roots(rts)
    return {
        "object": "(x^2-x-1) (x) (x^2-x-1)",
        "eigenvalues_desc": "{phi^2, -1, -1, phi'^2}",
        "charpoly": charpoly,
        "charpoly_factored": "(x+1)^2 (x^2-3x+1)",
        "mahler_tropical": mp.nstr(m_tropical, 20),
        "equals_phi_squared": bool(abs(m_tropical - PHI**2) < mp.mpf(10) ** (-15)),
        "retracted_factored_value_phi4": mp.nstr(PHI**4, 20),
        "tropical_differs_from_factored":
            bool(abs(m_tropical - PHI**4) > mp.mpf("0.1")),
        "charge_group": core.charge_group_from_roots(rts),
        "paper_ref": "Prop 3.3 / ledger M",
    }


def offcircle_tensor():
    rts = core.tensor_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    charpoly = core.monic_charpoly_int(rts)
    m = core.mahler_from_roots(rts)
    all_off = all(abs(r) > 1 for r in rts)
    return {
        "object": "(x^3-2) (x) (x^4-2)",
        "charpoly": "x^12-128" if charpoly == [1] + [0] * 11 + [-128] else str(charpoly),
        "mahler": mp.nstr(m, 20),
        "mahler_equals_128": bool(abs(m - 128) < mp.mpf(10) ** (-15)),
        "all_products_off_unit_circle": bool(all_off),
        "factored_law_agrees_here": "yes (off-circle: (2^{7/12})^12 = 2^7)",
        "charge_group": core.charge_group_from_roots(rts),
        "paper_ref": "Thm 3.5 / ledger D",
    }


def adams_measure_law():
    out = []
    for n in (3, 5):
        base = [1] + [0] * (n - 1) + [-2]   # M = 2
        for k in (2, 3):
            m = core.mahler_from_roots(core.adams_roots(base, k))
            out.append({
                "base": f"x^{n}-2",
                "k": k,
                "mahler_psi_k": mp.nstr(m, 20),
                "predicted_M_pow_k": 2**k,
                "matches": bool(abs(m - 2**k) < mp.mpf(10) ** (-15)),
            })
    return {"law": "M(psi^k A) = M(A)^k", "checks": out, "paper_ref": "Sec 2.1 table"}


def adams_primary_decomposition():
    base = [1, 0, 0, 0, 0, 0, -2]   # x^6 - 2, charge Z/6
    g0 = core.charge_group_from_roots(core.roots(base))
    g_psi3 = core.charge_group_from_roots(core.adams_roots(base, 3))
    g_psi2 = core.charge_group_from_roots(core.adams_roots(base, 2))
    return {
        "object": "x^6-2",
        "base_charge_group": g0,
        "psi^3_residual_charge_group": g_psi3,
        "psi^2_residual_charge_group": g_psi2,
        "realizes": "CRT projectors: psi^3 -> Z/2 (2-primary), psi^2 -> Z/3 (3-primary)",
        "paper_ref": "Thm 3.4 / ledger C",
    }


def sumset_charge_rule():
    rts = core.tensor_roots([1, 0, 0, -2], [1, 0, 0, 0, -2])
    got = set(core.charges_from_roots(rts, 12))
    expected = {(4 * a + 3 * b) % 12 for a in range(3) for b in range(4)}
    return {
        "object": "(x^3-2) (x) (x^4-2)",
        "composite_charges_Z12": sorted(got),
        "sumset_prediction": sorted(expected),
        "matches": got == expected,
        "fills_full_Z12 (coprime)": got == set(range(12)),
        "rule": "charge(a (x) b) = charge(a) + charge(b) (sumset)",
        "paper_ref": "Sec 2.1 table / Thm 3.5",
    }


def main():
    payload = {
        "phi_tensor_phi_tropical": phi_tensor_phi(),
        "offcircle_tensor": offcircle_tensor(),
        "adams_measure_law": adams_measure_law(),
        "adams_primary_decomposition": adams_primary_decomposition(),
        "sumset_charge_rule": sumset_charge_rule(),
    }
    path = write_json("tensor_adams.json", payload, __file__)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
