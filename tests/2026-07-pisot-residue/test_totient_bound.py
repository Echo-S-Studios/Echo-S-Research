"""Lemma 2.3 (completeness bound) and the totient-filter candidate counts.

    2 phi(M)^2 >= M   with equality iff M=2,
so  Phi_M | P  forces  M <= 2 (deg P)^2,  giving the finite 2 n^4 scan bound.
The paper also quotes exact candidate counts #{m: phi(m) <= K} used to size the
several scans; we recompute them from scratch.
"""
import os
import sys

from sympy import totient

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _phi(m):
    return int(totient(m))


def test_completeness_bound_no_violation_and_unique_equality():
    """Lemma 2.3 / its proof: swept exactly, 2 phi(M)^2 >= M for all M, and the
    only tight case is M=2.  The paper states this sweep to M <= 2e5; we sweep
    the same range."""
    violations = []
    equalities = []
    for M in range(1, 200001):
        ph = _phi(M)
        val = 2 * ph * ph
        if val < M:
            violations.append(M)
        elif val == M:
            equalities.append(M)
    assert violations == []
    assert equalities == [2]


def test_derived_scan_bound_2n4():
    """Lemma 2.3: if Phi_M | Rat_p then phi(M) <= deg Rat_p = n^2, hence
    M <= 2 phi(M)^2 <= 2 n^4.  Check the derivation is self-consistent for the
    degrees that appear in the paper (n = 4,5,10,12) and that the quoted bounds
    match 2 n^4."""
    quoted = {5: 1250, 10: 20000, 12: 41472}    # quintic, Lehmer, deg-12 Salem
    for n, bound in quoted.items():
        assert 2 * n ** 4 == bound
        # any M with phi(M) <= n^2 indeed satisfies M <= 2 n^4
        for M in range(1, 2 * n ** 4 + 1):
            if _phi(M) <= n * n:
                assert M <= 2 * n ** 4


def test_totient_candidate_counts():
    """Exact candidate counts #{m >= 1 : phi(m) <= K} and the largest such m,
    as quoted for the three scans:
      quintic Rat scan  (deg 25):  53 candidates
      deg-12 Salem scan (deg 144): 290 candidates, largest m = 630
      C_2 scan          (deg 400): 790 candidates, largest m = 1680
    """
    def count_and_max(K):
        ms = [m for m in range(1, 4 * K * K + 10) if _phi(m) <= K]
        return len(ms), max(ms)

    c25, _ = count_and_max(25)
    assert c25 == 53                        # "53 totient-filtered candidates"

    c144, m144 = count_and_max(144)
    assert (c144, m144) == (290, 630)

    c400, m400 = count_and_max(400)
    assert (c400, m400) == (790, 1680)


def test_totient_lower_bound_proof_pieces():
    """Lemma 2.3 proof: odd prime powers use  p-1 >= sqrt(p) <=> p^2-3p+1 >= 0
    for p >= 3; and phi(M) >= sqrt(M/2) overall.  Spot-check both."""
    from sympy import isprime, sqrt, Rational
    for p in [3, 5, 7, 11, 13, 97]:
        assert isprime(p)
        assert (p - 1) ** 2 >= p                      # p-1 >= sqrt(p)
        assert p * p - 3 * p + 1 >= 0
    for M in range(1, 5000):
        assert 2 * _phi(M) ** 2 >= M                  # phi(M) >= sqrt(M/2)
