r"""
Producer: the completeness bound of the contact scan (Lemma 4.5, ledger N).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/completeness_bound.json

Refactors Lemma 4.5 (lem:complete): if Phi_M | P then phi(M) <= deg P, and for
every M >= 1, phi(M) >= sqrt(M/2), i.e. 2 phi(M)^2 >= M; hence scanning all
M <= 2 (deg P)^2 decides completely which cyclotomics divide P.  Ledger N: the
inequality 2 phi(M)^2 >= M is verified by exact integer sieve for all
1 <= M <= 2*10^5, with M=2 the unique tight case.

Emits: the sieve verdict (no violations, unique tight case), agreement with
sympy.totient on a range, the prime-power proof steps, and the derived scan
bounds 2 (deg Rat)^2 quoted in the paper for every instance.

Run: py code/2026-07-relational-charge/completeness_bound.py
"""

import sympy as sp

import relcharge_core as C
from relcharge_io import write_json

# (label, deg Rat, quoted bound) from the paper
QUOTED_BOUNDS = [
    ("Rat_beta4", 16, 512),        # ledger G
    ("Rat_Lehmer", 100, 20000),    # ledger H
    ("Rat_S6", 36, 2592),          # ledger P
    ("Rat_S8", 64, 8192),          # ledger P
    ("census_deg12 Rat", 144, 41472),   # ledger T
    ("mixed Rat_{L,beta4}", 40, 3200),   # ledger I
    ("nested Rat_{Rat_p} x^4-x+1", 256, 131072),  # ledger X
]


def main():
    N = 200000
    phi = C.totients_upto(N)
    violations = [M for M in range(1, N + 1) if 2 * phi[M] * phi[M] < M]
    tight = [M for M in range(1, N + 1) if 2 * phi[M] * phi[M] == M]

    # cross-check the sieve against sympy.totient on a range
    sieve_matches_sympy = all(phi[M] == int(sp.totient(M)) for M in range(1, 3001))

    # prime-power proof steps
    two_part = all(sp.totient(2**a) == 2 ** (a - 1) for a in range(1, 25))
    odd_prime_power = all(
        (p - 1) ** 2 >= p and all(int(sp.totient(p**b)) ** 2 >= p**b for b in range(1, 5))
        for p in sp.primerange(3, 200)
    )

    # every quoted scan bound equals 2 (deg Rat)^2
    bounds = []
    for label, d, quoted in QUOTED_BOUNDS:
        bounds.append({
            "object": label,
            "rat_degree": d,
            "scan_bound_2d2": 2 * d * d,
            "quoted_in_paper": quoted,
            "matches": (2 * d * d == quoted),
        })

    # phi(M) <= deg forces M <= 2 deg^2 (checked for several deg bounds)
    forcing = []
    for dcap in (1, 2, 5, 10, 16, 40, 100, 144, 256):
        beyond = [M for M in range(2 * dcap * dcap + 1, min(N, 4 * dcap * dcap) + 1)
                  if phi[M] <= dcap]
        forcing.append({
            "deg_cap": dcap,
            "max_M_with_phi_le_cap": max(M for M in range(1, 2 * dcap * dcap + 1)
                                         if phi[M] <= dcap),
            "any_M_beyond_2dcap2_with_phi_le_cap": len(beyond) > 0,
        })

    payload = {
        "lemma_4_5_totient_bound": {
            "inequality": "2 phi(M)^2 >= M for all M >= 1",
            "range_checked": N,
            "violations": violations,
            "tight_cases": tight,
            "unique_tight_case_is_2": tight == [2],
            "status": "[forced] (proof) + [computed] (sieve, ledger N)",
        },
        "sieve_matches_sympy_totient_to_3000": sieve_matches_sympy,
        "prime_power_steps": {
            "two_part_phi_2a_equals_2a_minus_1": bool(two_part),
            "odd_prime_power_phi_ge_sqrt": bool(odd_prime_power),
        },
        "phi_le_deg_forces_M_le_2deg2": forcing,
        "quoted_scan_bounds": bounds,
        "all_quoted_bounds_match": all(b["matches"] for b in bounds),
    }
    path = write_json("completeness_bound.json", payload, __file__)
    print(f"wrote {path}")
    print(f"  totient bound 2phi(M)^2>=M: {len(violations)} violations, "
          f"tight cases {tight} (Lemma 4.5 / ledger N)")


if __name__ == "__main__":
    main()
