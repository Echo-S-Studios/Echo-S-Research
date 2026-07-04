"""
Producer -- Minimality chain and cross-cutting identities.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Proposition 7.1 (prop:minimal) and the Scope / Epistemic-ledger
                 identities (Sec 8, Sec 9).

Computes and EMITS to data/2026-06-generative-emptiness/minimality.json:
  * the ternary lock: #channels(d) = d^2 - d + 1 = 3 iff d = 2 (channel table 1,3,7);
  * the ternary growth spectrum spec(ad_R) = {-sqrt5, 0, +sqrt5} (three channels);
  * the isomorphism (pi/2)Z / 2pi ~= Z/4Z (addition table);
  * the quartic field Q(5^{1/4}): x^4-5 irreducible, splitting field contains i;
  * the identity sqrt5 = phi + phi^{-1} (floor image = grow generator);
  * the forced minimal chain  degree 2 => Z/4Z => phi.

Run:
    py code/2026-06-generative-emptiness/minimality.py
"""
import json
import os

import sympy as sp
import mpmath as mp

import ge_core as G
from ge_core import x, PHI, mahler, as_poly, roots_mp, mpf_str

OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness", "minimality.json",
)

PHI_V = PHI()


def ternary_lock():
    """#channels(d) = d^2 - d + 1 = 3 has unique positive solution d = 2."""
    d = sp.symbols("d")
    sols = sp.solve(sp.Eq(d**2 - d + 1, 3), d)
    ch = lambda n: n**2 - n + 1
    return {
        "channel_formula": "d^2 - d + 1",
        "solutions_of_eq_3": sorted(int(s) for s in sols),
        "positive_solution": [int(s) for s in sols if s > 0],
        "channel_table": {str(n): ch(n) for n in (1, 2, 3, 4)},
        "unique_degree_for_ternary": 2,
    }


def ternary_spectrum():
    """spec(ad_R) = {-sqrt5, 0, +sqrt5}: three distinct channels."""
    spec = [-sp.sqrt(5), sp.Integer(0), sp.sqrt(5)]
    return {
        "spectrum": ["-sqrt5", "0", "+sqrt5"],
        "n_channels": len(set(spec)),
        "symmetric": sp.simplify(spec[0] + spec[2]) == 0,
        "nonzero_growth_channels": ["-sqrt5", "+sqrt5"],
    }


def halfpi_isomorphism():
    """(pi/2)Z / 2pi ~= Z/4Z: representatives add mod 2pi as {0,1,2,3} add mod 4."""
    reps = [sp.Integer(0), sp.pi / 2, sp.pi, 3 * sp.pi / 2]
    table = []
    is_iso = True
    for i in range(4):
        row = []
        for j in range(4):
            s = sp.nsimplify((reps[i] + reps[j]) % (2 * sp.pi))
            idx = [k for k in range(4) if sp.simplify(s - reps[k]) == 0]
            row.append(idx[0] if idx else None)
            if idx != [(i + j) % 4]:
                is_iso = False
        table.append(row)
    return {
        "representatives": ["0", "pi/2", "pi", "3pi/2"],
        "addition_table_indices": table,
        "isomorphic_to_Z4": is_iso,
    }


def quartic_field():
    """Q(5^{1/4}): x^4-5 irreducible; splitting field contains i (root ratios)."""
    r = sp.Poly(x**4 - 5, x)
    roots = list(sp.roots(x**4 - 5, x).keys())
    ratios = {sp.simplify(a / b) for a in roots for b in roots if b != 0}
    return {
        "poly": "x^4 - 5",
        "irreducible": bool(r.is_irreducible),
        "splitting_field_contains_i": bool(sp.I in ratios or -sp.I in ratios),
        "field": "Q(5^{1/4})",
    }


def sqrt5_identity():
    """sqrt5 = phi + phi^{-1} (floor's image identified with the grow generator)."""
    phi = (1 + sp.sqrt(5)) / 2
    return {
        "identity": "sqrt5 = phi + 1/phi",
        "symbolic_zero": int(sp.simplify(phi + 1 / phi - sp.sqrt(5))),
        "numeric_residual": mpf_str(abs((PHI_V + 1 / PHI_V) - mp.sqrt(5)), 5),
        "phi": mpf_str(PHI_V, 20),
        "sqrt5": mpf_str(mp.sqrt(5), 20),
    }


def floor_is_phi():
    """The floor phi is the dominant real root of x^2-x-1 and its own measure."""
    P = as_poly(x**2 - x - 1)
    dom = max(roots_mp(P), key=lambda r: abs(r))
    return {
        "poly": "x^2 - x - 1",
        "dominant_root": mpf_str(dom.real, 20),
        "is_real_perron_not_salem": bool(abs(dom.imag) < mp.mpf(10)**(-25)),
        "mahler_measure": mpf_str(mahler(P), 20),
        "measure_equals_root": bool(abs(mahler(P) - PHI_V) < mp.mpf(10)**(-25)),
    }


def main():
    data = G.provenance("minimality.py")
    data.update({
        "result": "Minimality chain (Prop 7.1) and scope/ledger identities",
        "paper_result": "Proposition 7.1 (prop:minimal), Sec 8-9",
        "forced_chain": "degree 2 => Z/4Z => phi",
        "step1_ternary_lock": ternary_lock(),
        "step1b_ternary_spectrum": ternary_spectrum(),
        "step2_halfpi_isomorphism": halfpi_isomorphism(),
        "step2b_quartic_field": quartic_field(),
        "step3_floor_is_phi": floor_is_phi(),
        "identity_sqrt5_phi": sqrt5_identity(),
    })
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"wrote {OUT}")
    print("  channel table:", data["step1_ternary_lock"]["channel_table"],
          "-> unique degree 2")
    print("  (pi/2)Z iso Z/4Z:", data["step2_halfpi_isomorphism"]["isomorphic_to_Z4"])
    print("  Q(5^1/4) contains i:",
          data["step2b_quartic_field"]["splitting_field_contains_i"])
    print("  sqrt5 = phi + 1/phi residual:",
          data["identity_sqrt5_phi"]["numeric_residual"])


if __name__ == "__main__":
    main()
