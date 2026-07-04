"""
Producer for Section 8 -- The unified picture, Theorem 8.1 (thm:unified) of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

The spectral operators form a commutative lambda-ring with two characters
(Mahler measure, angle charge).  The computationally load-bearing meeting
point is
    sqrt5 = phi + phi^{-1}
-- the (x)-logarithm of the floor and the grow generator of the self-action --
with companion identity phi - phi^{-1} = 1.  On the Lorentzian seed K both
characters are simultaneously defined: Character I gives M(K) = beta^2 =
phi^2 sqrt5, Character II gives the full Z/4Z.

Emits data/2026-06-operator-algebra/unified.json.

Difference from tests/: verifier asserts each identity; this producer emits the
two meeting-point identities (exact + decimal), the lambda-ring endomorphism
axiom evaluated on a concrete mixed object, and both characters read off K.
"""

from __future__ import annotations

import sympy as sp

from opalg_core import (
    K_seed,
    charge,
    charge_multiset_str,
    dec_str,
    exact_str,
    golden_seed,
    mahler_exact,
    ms_equal,
    oplus,
    otimes,
    phi,
    psi,
    sqrt5,
    write_json,
)


def main():
    G = golden_seed()
    K = K_seed()

    # the meeting-point identities
    id_sqrt5 = sp.simplify(sqrt5 - (phi + 1 / phi)) == 0
    id_diff = sp.simplify((phi - 1 / phi) - 1) == 0
    id_inv = sp.simplify((1 / phi) - (phi - 1)) == 0

    # phi as floor and as additive grow-generator
    floor_ok = sp.simplify(mahler_exact(G) - phi) == 0
    grow_ok = sp.simplify(mahler_exact(oplus(G, G)) - phi**2) == 0

    # lambda-ring endomorphism axiom on a concrete mixed object (golden + K)
    endo_add = ms_equal(psi(2, oplus(G, K)), oplus(psi(2, G), psi(2, K)))
    endo_mul = ms_equal(psi(2, otimes(G, K)), otimes(psi(2, G), psi(2, K)))

    # both characters on K
    MK = mahler_exact(K)
    charI_ok = sp.simplify(MK - phi**2 * sqrt5) == 0
    charII = sorted(charge(K).keys())

    payload = {
        "meeting_point_identities": {
            "sqrt5 = phi + phi^-1": {
                "holds": bool(id_sqrt5),
                "sqrt5_decimal": dec_str(sqrt5, 30),
                "phi_plus_inv_decimal": dec_str(phi + 1 / phi, 30),
            },
            "phi - phi^-1 = 1": {"holds": bool(id_diff)},
            "phi^-1 = phi - 1": {"holds": bool(id_inv),
                                 "phi_inv_decimal": dec_str(1 / phi, 30)},
        },
        "phi_is_floor_and_grow_generator": {
            "M(golden) = phi": bool(floor_ok),
            "M(golden (+) golden) = phi^2": bool(grow_ok),
            "note": "least Mahler generator (floor) = additive self-action base",
        },
        "lambda_ring_endomorphism_axiom_on_golden_plus_K": {
            "psi^2(A (+) B) = psi^2 A (+) psi^2 B": bool(endo_add),
            "psi^2(A (x) B) = psi^2 A (x) psi^2 B": bool(endo_mul),
        },
        "two_characters_on_K": {
            "Character_I_M(K)": exact_str(MK),
            "M(K) = phi^2 sqrt5 = beta^2": bool(charI_ok),
            "M(K)_decimal": dec_str(MK, 30),
            "Character_II_charge(K)": charII,
            "charge(K)_multiset": charge_multiset_str(charge(K)),
            "charge(K) = full Z/4Z": charII == [0, 1, 2, 3],
        },
    }
    path = write_json("unified.json", "unified.py", payload)
    print(f"wrote {path}")
    print(f"  sqrt5=phi+phi^-1: {bool(id_sqrt5)}; phi-phi^-1=1: {bool(id_diff)}; "
          f"M(K)=beta^2: {bool(charI_ok)}; charge(K)={charII}")


if __name__ == "__main__":
    main()
