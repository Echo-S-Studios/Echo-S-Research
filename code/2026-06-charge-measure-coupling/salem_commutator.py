"""
Producer: Salem exclusion (Lem 8.1) and the commutator escape (Prop 8.2).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/salem_commutator.json

Computational content only (the cited matrix theorems -- Shoda, Albert-Muckenhoupt,
Laffey-Reams -- are external inputs):
  * Lem 8.1: beta_4 and Lehmer's number are Salem-type (one root outside, its
    inverse inside, the rest on |z|=1 at irrational angles) and hence
    charge-inadmissible.
  * Prop 8.2 / ledger L: Lehmer's polynomial has trace -1, so L(x)(x-1) has
    trace 0; its integer companion matrix is a commutator whose spectrum carries
    Lehmer's number tau = 1.17628 in (1,phi) -- below the emission floor -- with
    charge group bottom.  The one non-abelian door through which the excluded
    Salem spectrum re-enters.

Run: py code/2026-06-charge-measure-coupling/salem_commutator.py
"""

import mpmath as mp
import numpy as np
import sympy as sp

import cmc_core as core
from cmc_io import write_json

mp.mp.dps = 50
_x = sp.symbols("x")
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
EPS = mp.mpf(10) ** (-20)

LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
COMMUTATOR = [1, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, -1]  # L(x)(x-1)


def _salem_structure(coeffs):
    rts = core.roots(coeffs)
    outside = [r for r in rts if abs(r) > 1 + EPS]
    inside = [r for r in rts if abs(r) < 1 - EPS]
    oncirc = [r for r in rts if abs(abs(r) - 1) <= EPS]
    return outside, oncirc, inside


def salem_witness(name, coeffs, ref):
    outside, oncirc, inside = _salem_structure(coeffs)
    m = core.mahler(coeffs)
    cg = core.charge_group(coeffs, nmax=400)
    return {
        "object": name,
        "reciprocal": core.is_reciprocal(coeffs),
        "roots_outside_unit_circle": len(outside),
        "roots_on_unit_circle": len(oncirc),
        "roots_inside_unit_circle": len(inside),
        "is_salem": len(outside) == 1 and len(inside) == 1 and len(oncirc) >= 2,
        "mahler": mp.nstr(m, 18),
        "mahler_equals_salem_number": bool(abs(m - abs(outside[0])) < EPS),
        "charge_group": cg,
        "charge_inadmissible": cg is None,
        "paper_ref": ref,
    }


def commutator_escape():
    # Lehmer trace = -(coeff of x^9)
    lehmer_trace = -LEHMER[1]
    # L(x)(x-1) via sympy, confirm the coeff list and trace 0
    Lpoly = sp.Poly(sum(c * _x**(10 - i) for i, c in enumerate(LEHMER)), _x)
    prod = sp.Poly(sp.expand(Lpoly.as_expr() * (_x - 1)), _x)
    prod_coeffs = [int(a) for a in prod.all_coeffs()]
    comm_trace = -prod_coeffs[1]

    # integer companion matrix of L(x)(x-1); check trace 0 and eigenvalues
    p = COMMUTATOR
    n = len(p) - 1
    C = np.zeros((n, n))
    C[0, :] = [-p[k] for k in range(1, n + 1)]
    for i in range(1, n):
        C[i, i - 1] = 1
    eig = np.linalg.eigvals(C)
    tau = float(core.mahler(LEHMER))
    m = core.mahler(COMMUTATOR)
    cg = core.charge_group(COMMUTATOR, nmax=400)

    return {
        "lehmer_polynomial": "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1",
        "lehmer_trace": lehmer_trace,
        "commutator_polynomial": "x^11-x^9-x^8+x^3+x^2-1",
        "commutator_poly_coeffs": prod_coeffs,
        "commutator_trace": comm_trace,
        "companion_matrix_trace": float(round(np.trace(C), 12)),
        "companion_is_integer": bool(np.allclose(C, np.round(C))),
        "carries_lehmer_number": bool(min(abs(eig - tau)) < 1e-6),
        "carries_eigenvalue_1": bool(min(abs(eig - 1.0)) < 1e-6),
        "mahler": mp.nstr(m, 18),
        "mahler_in_gap_1_phi": bool(1 < m < PHI),
        "charge_group": cg,
        "charge_bottom": cg is None,
        "interpretation": "trace-zero integer matrix (Shoda/Laffey-Reams commutator) "
                          "carrying a Salem number below the floor and off every "
                          "angle lattice: the abelian constraints do not extend",
        "paper_ref": "Prop 8.2 / ledger L",
    }


def main():
    payload = {
        "salem_witnesses": [
            salem_witness("beta_4 = x^4-x^3-x^2-x+1", [1, -1, -1, -1, 1],
                          "ledger F / Lem 8.1"),
            salem_witness("Lehmer = x^10+x^9-x^7-...-x^3+x+1", LEHMER,
                          "Prop 8.2 / Lem 8.1"),
        ],
        "commutator_escape": commutator_escape(),
    }
    path = write_json("salem_commutator.json", payload, __file__)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
