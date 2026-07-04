"""Producer: the completeness bound and totient-filter candidate counts.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: Lemma 2.3 (completeness bound) -- 2 phi(M)^2 >= M with equality iff
M=2 (swept exactly to M <= 2e5), the derived finite scan bound M <= 2 n^4, and
the exact totient-filter candidate counts sizing the scans that appear in the
paper:
    quintic Rat scan  (deg 25):  53 candidates
    deg-12 Salem scan (deg 144): 290 candidates, largest m = 630
    C_2 scan          (deg 400): 790 candidates, largest m = 1680
plus the Lehmer scan sizing (deg 100, bound 20000).

Emits:
  data/2026-07-pisot-residue/totient_bound.json
"""
from sympy import totient

from pisot_lib import write_json


def _phi(m):
    return int(totient(m))


def completeness_sweep(limit=200001):
    """Sweep 2 phi(M)^2 >= M over 1..limit-1; record violations and equalities."""
    violations, equalities = [], []
    for M in range(1, limit):
        val = 2 * _phi(M) ** 2
        if val < M:
            violations.append(M)
        elif val == M:
            equalities.append(M)
    return violations, equalities


def scan_bounds():
    """The 2 n^4 bounds for the degrees that appear in the paper."""
    quoted = {5: 1250, 10: 20000, 12: 41472}      # quintic, Lehmer, deg-12 Salem
    out = []
    for n, bound in quoted.items():
        out.append({
            "n": n,
            "deg_Rat_p": n * n,
            "scan_bound_2n4": 2 * n ** 4,
            "matches_paper": 2 * n ** 4 == bound,
        })
    return out


def candidate_counts():
    """Exact #{m >= 1 : phi(m) <= K} and the largest such m for K in {25,144,400}."""
    def count_and_max(K):
        ms = [m for m in range(1, 4 * K * K + 10) if _phi(m) <= K]
        return len(ms), max(ms)

    out = {}
    for label, K, deg in [("quintic_Rat", 25, 25),
                          ("deg12_Salem", 144, 144),
                          ("C2", 400, 400)]:
        cnt, mx = count_and_max(K)
        out[label] = {"K_equals_deg": K, "num_candidates": cnt, "largest_m": mx}
    return out


def main():
    viol, eq = completeness_sweep()
    payload = {
        "description": "Lemma 2.3 completeness bound 2 phi(M)^2 >= M (equality iff M=2), "
                       "derived scan bound M <= 2 n^4, and totient-filter candidate counts.",
        "completeness_bound": {
            "swept_to": 200000,
            "violations": viol,
            "equality_cases": eq,
            "unique_tight_case_is_M2": eq == [2],
        },
        "derived_scan_bounds": scan_bounds(),
        "totient_filter_candidate_counts": candidate_counts(),
    }
    path = write_json("totient_bound.json", payload, "totient_bound.py")
    print(f"wrote {path}")
    print(f"  completeness bound: violations={viol}, equality cases={eq}")
    for label, rec in payload["totient_filter_candidate_counts"].items():
        print(f"  {label:<12} phi<= {rec['K_equals_deg']:>3}: "
              f"{rec['num_candidates']} candidates, largest m={rec['largest_m']}")


if __name__ == "__main__":
    main()
