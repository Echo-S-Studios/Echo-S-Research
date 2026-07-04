"""
Producer: the parity-graded Mahler floor mu(n) (Theorem 4.1, the coupling).

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    : data/2026-06-charge-measure-coupling/parity_floor.csv

For each charge-group order n it records the coupling result
    mu(n) = phi  (n even)   /   2  (n odd),
computed from the paper's explicit floor-attaining constructions:
  * even n = 2k : q_k = x^{2k}+x^k-1 has charge group Z/2k, attains all charges,
    and Mahler measure phi  (the golden seed sits on the lattice because pi is
    in it) -- realized attainment, [forced].
  * odd n       : x^n-2 has charge group Z/n and Mahler measure 2 -- the pi-ray
    obstruction (Lem 4.5) blocks the golden seed, so the realized floor is 2.
The row also carries the parity bit "lattice contains pi (argument of phi')",
which is the single fact driving the whole dichotomy.

Run: py code/2026-06-charge-measure-coupling/parity_floor.py
"""

import mpmath as mp
import sympy as sp

import cmc_core as core
from cmc_io import write_csv

mp.mp.dps = 50
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))


def q_k(k):
    coeffs = [0] * (2 * k + 1)
    coeffs[0], coeffs[k], coeffs[2 * k] = 1, 1, -1
    return coeffs


def build_rows(nmax: int = 12):
    rows = []
    for n in range(2, nmax + 1):
        even = (n % 2 == 0)
        # lattice (2pi/n)Z contains the ray pi  <=>  n even
        contains_pi = even
        if even:
            k = n // 2
            coeffs = q_k(k)
            witness = f"q_{k} = x^{n}+x^{k}-1"
            m = core.mahler(coeffs)
            cg = core.charge_group(coeffs)
            allc = core.attains_all_charges(coeffs, n)
            mu_closed = "phi"
        else:
            coeffs = [1] + [0] * (n - 1) + [-2]
            witness = f"x^{n}-2"
            m = core.mahler(coeffs)
            cg = core.charge_group(coeffs)
            allc = core.attains_all_charges(coeffs, n)
            mu_closed = "2"
        rows.append({
            "n": n,
            "parity": "even" if even else "odd",
            "lattice_contains_pi": "yes" if contains_pi else "no",
            "mu_n_closed_form": mu_closed,
            "mu_n_value": mp.nstr(m, 20),
            "witness": witness,
            "witness_charge_group": f"Z/{cg}",
            "all_charges_attained": "yes" if allc else "no",
            "golden_seed_hostable": "yes" if contains_pi else "no",
        })
    return rows


def main():
    rows = build_rows()
    fields = ["n", "parity", "lattice_contains_pi", "mu_n_closed_form", "mu_n_value",
              "witness", "witness_charge_group", "all_charges_attained",
              "golden_seed_hostable"]
    path = write_csv("parity_floor.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  mu(n) parity floor for n=2..{rows[-1]['n']} "
          f"(even -> phi, odd -> 2)")


if __name__ == "__main__":
    main()
