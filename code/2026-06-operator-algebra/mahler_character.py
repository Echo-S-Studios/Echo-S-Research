"""
Producer for Character I -- the Mahler measure, Theorem 4.2 (thm:measure)
[OA-M-01..07] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

The Mahler measure is ADDITIVE on (+) and squares under psi^2, but TROPICAL
on (x): log M(A (x) B) = sum_{l in A} sum_{m in B} (log|l| + log|m|)^+ , a
(max,+) convolution -- never a product.

Emits:
  * data/2026-06-operator-algebra/mahler_pairs.csv -- for each ordered seed
    pair: the additive law M(A(+)B)=M(A)M(B), and the tensor side showing
    M(A(x)B) (exact + decimal) vs the multiplicative guess M(A)M(B), whether
    the tensor is multiplicative, and the tropical convolution match.
    Includes the paper's non-multiplicative witness golden(x)(x^2-2):
    M = 2 phi^2 = 3+sqrt5 = 5.2360..., while M(golden)M(x^2-2) = 2 phi = 3.2360...
  * data/2026-06-operator-algebra/mahler_squaring.csv -- M(psi^2 A) = M(A)^2.

Difference from tests/: the verifier asserts the equalities; this producer
tabulates the concrete measures on both sides so the additive/tropical split
is visible as data.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

from opalg_core import (
    dec_str,
    exact_str,
    log_mahler_mp,
    mahler_exact,
    oplus,
    otimes,
    psi,
    seed_catalogue,
    tropical_convolution,
    write_csv,
)

DPS = 60
TOL = mp.mpf(10) ** (-45)


def main():
    cat = seed_catalogue()
    # precompute exact measures
    M = {name: mahler_exact(A) for name, _, A in cat}

    pair_rows = []
    for na, _, A in cat:
        for nb, _, B in cat:
            M_oplus = mahler_exact(oplus(A, B))
            M_prod = sp.simplify(M[na] * M[nb])
            additive_ok = bool(sp.simplify(M_oplus - M_prod) == 0)

            M_otimes = mahler_exact(otimes(A, B))
            tensor_mult = bool(sp.simplify(M_otimes - M_prod) == 0)

            logM = log_mahler_mp(otimes(A, B), dps=DPS)
            trop = tropical_convolution(A, B, dps=DPS)
            tropical_ok = bool(abs(logM - trop) < TOL)

            pair_rows.append(
                {
                    "A": na,
                    "B": nb,
                    "M_A": exact_str(M[na]),
                    "M_B": exact_str(M[nb]),
                    "M_oplus": exact_str(M_oplus),
                    "M_A_times_M_B": exact_str(M_prod),
                    "additive_ok": additive_ok,
                    "M_otimes_exact": exact_str(M_otimes),
                    "M_otimes_decimal": dec_str(M_otimes, 24),
                    "tensor_is_multiplicative": tensor_mult,
                    "logM_otimes": mp.nstr(logM, 24),
                    "tropical_convolution": mp.nstr(trop, 24),
                    "tropical_ok": tropical_ok,
                }
            )

    p_pairs = write_csv(
        "mahler_pairs.csv",
        "mahler_character.py",
        ["A", "B", "M_A", "M_B", "M_oplus", "M_A_times_M_B", "additive_ok",
         "M_otimes_exact", "M_otimes_decimal", "tensor_is_multiplicative",
         "logM_otimes", "tropical_convolution", "tropical_ok"],
        pair_rows,
    )

    sq_rows = []
    for name, _, A in cat:
        M_psi2 = mahler_exact(psi(2, A))
        M_sq = sp.simplify(M[name] ** 2)
        sq_rows.append(
            {
                "seed": name,
                "M_A": exact_str(M[name]),
                "M_psi2_A": exact_str(M_psi2),
                "M_A_squared": exact_str(M_sq),
                "squares_ok": bool(sp.simplify(M_psi2 - M_sq) == 0),
            }
        )
    p_sq = write_csv(
        "mahler_squaring.csv",
        "mahler_character.py",
        ["seed", "M_A", "M_psi2_A", "M_A_squared", "squares_ok"],
        sq_rows,
    )

    # headline witness sanity print
    witness = next(r for r in pair_rows if r["A"] == "golden" and r["B"] == "x^2-2")
    print(f"wrote {p_pairs} ({len(pair_rows)} pairs)")
    print(f"wrote {p_sq} ({len(sq_rows)} seeds)")
    print(f"  witness golden(x)(x^2-2): M_otimes={witness['M_otimes_exact']} "
          f"({witness['M_otimes_decimal'][:6]}...), "
          f"mult-would-give={witness['M_A_times_M_B']}, "
          f"multiplicative={witness['tensor_is_multiplicative']}")


if __name__ == "__main__":
    main()
