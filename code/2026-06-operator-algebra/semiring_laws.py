"""
Producer for Theorem 2.2 (thm:semiring) -- "The emission semiring"
[OA-SR-01..08] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

Emits data/2026-06-operator-algebra/semiring_laws.csv : one row per semiring
axiom, each verified as a multiset identity over GENERIC eigenvalue symbols
(hence for every substitution), with a concrete degree-bookkeeping instance
drawn from the seed catalogue.

Difference from tests/: the verifier asserts each identity; this producer
EMITS a table of the axioms, their holds-flag, and the concrete degrees
deg(A(+)B)=degA+degB, deg(A(x)B)=degA*degB.
"""

from __future__ import annotations

import sympy as sp

from opalg_core import (
    ONE,
    ZERO,
    golden_seed,
    ms_equal,
    oplus,
    otimes,
    seed_from_poly,
    write_csv,
    x,
)

# Generic eigenvalue multisets (independent symbols): a law that holds for
# independent symbols holds for every numeric substitution.
a1, a2, a3 = sp.symbols("a1 a2 a3")
b1, b2 = sp.symbols("b1 b2")
c1, c2 = sp.symbols("c1 c2")
A, B, C = (a1, a2, a3), (b1, b2), (c1, c2)


def semiring_laws():
    """Return (law_id, name, statement, holds) tuples over generic symbols."""
    return [
        ("OA-SR-01", "oplus commutative", "A (+) B = B (+) A",
         ms_equal(oplus(A, B), oplus(B, A))),
        ("OA-SR-02", "oplus associative", "(A (+) B) (+) C = A (+) (B (+) C)",
         ms_equal(oplus(oplus(A, B), C), oplus(A, oplus(B, C)))),
        ("OA-SR-03", "otimes commutative", "A (x) B = B (x) A",
         ms_equal(otimes(A, B), otimes(B, A))),
        ("OA-SR-04", "otimes associative", "(A (x) B) (x) C = A (x) (B (x) C)",
         ms_equal(otimes(otimes(A, B), C), otimes(A, otimes(B, C)))),
        ("OA-SR-05", "left distributivity", "A (x) (B (+) C) = (A(x)B) (+) (A(x)C)",
         ms_equal(otimes(A, oplus(B, C)), oplus(otimes(A, B), otimes(A, C)))),
        ("OA-SR-06", "right distributivity", "(B (+) C) (x) A = (B(x)A) (+) (C(x)A)",
         ms_equal(otimes(oplus(B, C), A), oplus(otimes(B, A), otimes(C, A)))),
        ("OA-SR-07a", "0 additive identity", "A (+) 0 = A",
         ms_equal(oplus(A, ZERO), A) and ms_equal(oplus(ZERO, A), A)),
        ("OA-SR-07b", "0 multiplicative annihilator", "A (x) 0 = 0",
         ms_equal(otimes(A, ZERO), ZERO) and ms_equal(otimes(ZERO, A), ZERO)),
        ("OA-SR-08", "1={1} multiplicative identity", "A (x) 1 = A",
         ms_equal(otimes(A, ONE), A) and ms_equal(otimes(ONE, A), A)),
    ]


def main():
    # Concrete instance for degree bookkeeping: deg(A(+)B)=degA+degB,
    # deg(A(x)B)=degA*degB.  Use golden (deg 2) and x^2-2 (deg 2), K (deg 4).
    G = golden_seed()
    T = seed_from_poly(x**2 - 2)
    instance = (
        f"deg(golden(+)x^2-2)={len(oplus(G, T))}=2+2; "
        f"deg(golden(x)x^2-2)={len(otimes(G, T))}=2*2"
    )

    rows = []
    for law_id, name, statement, holds in semiring_laws():
        rows.append(
            {
                "law_id": law_id,
                "law": name,
                "statement": statement,
                "holds_over_generic_symbols": bool(holds),
                "concrete_instance": instance,
            }
        )
    # Explicit degree-bookkeeping rows across the catalogue.
    for (n1, s1) in [("golden", G), ("x^2-2", T)]:
        for (n2, s2) in [("golden", G), ("K", seed_from_poly(x**4 + 5 * x**2 - 5))]:
            rows.append(
                {
                    "law_id": "OA-SR-DEG",
                    "law": "degree bookkeeping",
                    "statement": f"deg({n1}(+){n2}) = deg+deg ; deg({n1}(x){n2}) = deg*deg",
                    "holds_over_generic_symbols": (
                        len(oplus(s1, s2)) == len(s1) + len(s2)
                        and len(otimes(s1, s2)) == len(s1) * len(s2)
                    ),
                    "concrete_instance": (
                        f"|{n1}(+){n2}|={len(oplus(s1, s2))}; "
                        f"|{n1}(x){n2}|={len(otimes(s1, s2))}"
                    ),
                }
            )

    path = write_csv(
        "semiring_laws.csv",
        "semiring_laws.py",
        ["law_id", "law", "statement", "holds_over_generic_symbols", "concrete_instance"],
        rows,
    )
    ok = all(r["holds_over_generic_symbols"] for r in rows)
    print(f"wrote {path} ({len(rows)} rows; all_hold={ok})")


if __name__ == "__main__":
    main()
