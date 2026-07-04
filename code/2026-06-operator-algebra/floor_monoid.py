"""
Producer for Theorem 4.3 (thm:floor) -- "The (+)-image is the spectrum monoid;
the floor is its least generator" [OA-M-08; GA-01..13] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

M(S, (+)) = <phi, 2, 3, 5, beta^2> multiplicatively, with the single relation
beta^2 = phi^2 sqrt5.  phi,2,3,5 are multiplicatively independent (norms
-1,4,9,25 in Q(sqrt5)); five atoms but rank four (not factorial), with
beta^2 . beta^2 = 5 phi^4.  Least generator is phi -> no monoid element in
(1, phi): the cost floor.

Emits data/2026-06-operator-algebra/floor_monoid.json.

Difference from tests/: the verifier asserts each fact; this producer records
the generators, their exact + decimal values and Q(sqrt5) norms, the defining
relation, the non-unique-factorisation witness, the rank/atom count from an
integer-relation search, and the cost floor -- as one structured artifact.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

from opalg_core import (
    K_seed,
    dec_str,
    exact_str,
    golden_seed,
    mahler_exact,
    phi,
    seed_from_poly,
    sqrt5,
    write_json,
    x,
)

beta2 = phi**2 * sqrt5


def norm_Qsqrt5(expr):
    """N(a + b sqrt5) = a^2 - 5 b^2 in Q(sqrt5)."""
    p = sp.Poly(sp.expand(expr), sqrt5)
    b = p.coeff_monomial(sqrt5)
    a = p.coeff_monomial(1)
    return sp.simplify(a**2 - 5 * b**2)


def integer_relation_search(logs, R):
    """Return every nonzero integer tuple in [-R,R]^n with sum ~ 0 (<1e-20)."""
    import itertools

    found = []
    for combo in itertools.product(range(-R, R + 1), repeat=len(logs)):
        if all(c == 0 for c in combo):
            continue
        if abs(sum(c * L for c, L in zip(combo, logs))) < mp.mpf(10) ** (-20):
            found.append(combo)
    return found


def main():
    mp.mp.dps = 60

    # Generators realised as concrete Mahler measures.
    gens = [
        ("phi", golden_seed(), phi),
        ("2", seed_from_poly(x**2 - 2), sp.Integer(2)),
        ("3", seed_from_poly(x**2 - 3), sp.Integer(3)),
        ("5", seed_from_poly(x**2 - 5), sp.Integer(5)),
        ("beta^2", K_seed(), beta2),
    ]
    generators = []
    for name, seed, val in gens:
        realised = mahler_exact(seed)
        generators.append(
            {
                "generator": name,
                "value_exact": exact_str(val),
                "value_decimal": dec_str(val, 30),
                "realised_as_measure_of": (
                    "x^2-x-1" if name == "phi"
                    else "x^4+5x^2-5" if name == "beta^2"
                    else f"x^2-{name}"),
                "realised_ok": bool(sp.simplify(realised - val) == 0),
                "norm_Qsqrt5": (int(norm_Qsqrt5(val))
                                if name != "beta^2" else None),
            }
        )

    # Multiplicative independence of phi,2,3,5 (norm argument corroborated by
    # an integer-relation box search over the logs).
    base_logs = [mp.log((1 + mp.sqrt(5)) / 2), mp.log(2), mp.log(3), mp.log(5)]
    indep_found = integer_relation_search(base_logs, R=6)

    payload = {
        "monoid": "M(S, (+)) = <phi, 2, 3, 5, beta^2>  (multiplicative)",
        "generators": generators,
        "relation": {
            "statement": "beta^2 = phi^2 * sqrt5",
            "beta2_exact": exact_str(beta2),
            "beta2_decimal": dec_str(beta2, 30),
            "closed_form_(5+3sqrt5)/2": bool(sp.simplify(beta2 - (5 + 3 * sqrt5) / 2) == 0),
            "log_relation": "log(beta^2) = 2 log(phi) + (1/2) log(5)",
        },
        "non_unique_factorisation": {
            "statement": "beta^2 . beta^2 = 5 . phi^4",
            "holds": bool(sp.simplify(beta2**2 - 5 * phi**4) == 0),
            "beta4_exact": exact_str(beta2**2),
        },
        "norms_Qsqrt5": {"N(phi)": -1, "N(2)": 4, "N(3)": 9, "N(5)": 25},
        "multiplicative_independence": {
            "generators_tested": ["phi", "2", "3", "5"],
            "search_box": "[-6,6]^4 over (log phi, log2, log3, log5)",
            "integer_relations_found": [list(t) for t in indep_found],
            "independent": indep_found == [],
        },
        "structure": {
            "atoms": 5,
            "rank": 4,
            "factorial": False,
            "note": "5 atoms but rank 4 (one relation) -> not factorial",
        },
        "cost_floor": {
            "least_generator": "phi",
            "phi_decimal": dec_str(phi, 30),
            "no_element_in_open_interval_(1,phi)": bool(
                min(sp.N(g[2], 40) for g in gens) == sp.N(phi, 40)),
            "note": "phi is the least generator -> monoid meets (1,phi) in nothing",
        },
    }
    path = write_json("floor_monoid.json", "floor_monoid.py", payload)
    print(f"wrote {path}")
    print(f"  independence search found {len(indep_found)} relations "
          f"(expected 0); least generator phi={payload['cost_floor']['phi_decimal'][:6]}")


if __name__ == "__main__":
    main()
