"""Producer: the keystone, the Fibonacci power law, and the classical identities.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 1 -- the associative core):
  * Thm 1.3   power law  R^n = F_n R + F_{n-1} I = [[F_{n-1},F_n],[F_n,F_{n+1}]]
  * Cor 1.5   Lucas trace Tr(R^n)=L_n, det(R^n)=(-1)^n, Cassini, Binet
  * Prop 1.2  keystone R^2=R+I, Vieta, root5 identity; low powers R^2..R^5

Outputs:
  data/.../power_law.csv          -- one row per n in [-16,16]
  data/.../fib_lucas_identities.json

Unlike tests/test_core.py (which asserts each fact), this script builds the
closed forms from the recurrences, cross-checks them against the actual matrix
power R**n once, and EMITS the table plus a per-row verified flag.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, I2, phi, psi, sqrt5, fib, luc, power_law, cassini,
                      mat_eq, is_zero)
from eap_io import write_csv, write_json

N_LO, N_HI = -16, 16


def power_law_rows():
    rows = []
    for n in range(N_LO, N_HI + 1):
        M = power_law(n)                       # closed form F_n R + F_{n-1} I
        verified = mat_eq(R ** n, M)           # vs actual matrix exponentiation
        rows.append({
            "n": n,
            "F_n": int(fib(n)),
            "F_{n-1}": int(fib(n - 1)),
            "F_{n+1}": int(fib(n + 1)),
            "R^n_00": int(M[0, 0]), "R^n_01": int(M[0, 1]),
            "R^n_10": int(M[1, 0]), "R^n_11": int(M[1, 1]),
            "trace_Ln": int(sp.trace(M)),          # Cor 1.5a: = L_n
            "L_n": int(luc(n)),
            "det_(-1)^n": int(M.det()),            # Cor 1.5b
            "cassini": int(cassini(n)),            # Cor 1.5c: = (-1)^n
            "verified_matrix_power": bool(verified),
        })
    return rows


def identities_payload():
    # keystone and Vieta / root5 (Prop 1.2)
    keystone = mat_eq(R * R, R + I2)
    vieta_sum = is_zero(phi + psi - 1)
    vieta_prod = is_zero(phi * psi + 1)
    root5 = is_zero((phi - psi) - sqrt5) and is_zero((phi + 1 / phi) - sqrt5)

    # low powers R^2..R^5 in (R, I) coordinates (Sec. 1.3)
    low_powers = {f"R^{k}": [int(fib(k)), int(fib(k - 1))] for k in range(2, 6)}

    # Binet phi^n = F_n phi + F_{n-1}, and F_n = (phi^n - psi^n)/sqrt5
    binet_ok = all(
        is_zero(phi ** n - (fib(n) * phi + fib(n - 1))) and
        is_zero(fib(n) - (phi ** n - psi ** n) / sqrt5)
        for n in range(-8, 9)
    )
    return {
        "keystone_R2_eq_R_plus_I": bool(keystone),
        "trace_R": int(sp.trace(R)),
        "det_R": int(R.det()),
        "vieta_sum_phi_plus_psi_eq_1": bool(vieta_sum),
        "vieta_prod_phi_psi_eq_minus1": bool(vieta_prod),
        "root5_identity_phi_minus_psi_eq_sqrt5_eq_phi_plus_inv": bool(root5),
        "low_powers_in_R_I_coords": low_powers,
        "low_powers_note": "R^k = a R + b I with (a,b) = (F_k, F_{k-1}); "
                           "(1,1),(2,1),(3,2),(5,3) are consecutive Fibonacci pairs",
        "binet_verified_n_-8_8": bool(binet_ok),
        "n_range": [N_LO, N_HI],
    }


def main():
    rows = power_law_rows()
    fn = list(rows[0].keys())
    p1 = write_csv("power_law.csv", fn, rows, __file__)
    p2 = write_json("fib_lucas_identities.json", identities_payload(), __file__)
    allok = all(r["verified_matrix_power"] for r in rows)
    print(f"wrote {p1}  ({len(rows)} rows, matrix-power check all pass: {allok})")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
