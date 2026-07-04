r"""Lemma 4.5 (lem:complete) and ledger N: the completeness of the contact scan.

Paper (Lemma 4.5): "If Phi_M | P then phi(M) <= deg P; and for every M >= 1,
phi(M) >= sqrt(M/2).  Hence scanning all M <= 2 (deg P)^2 decides ... which
cyclotomics divide P.  For P = Rat_p, of degree n^2, the bound is 2 n^4."
"The bound is tight at M=2 (phi(2)=1=sqrt(2/2)), so the constant 1/2 cannot be
removed."

Paper (ledger N): "2 phi(M)^2 >= M verified by exact integer sieve for all
1 <= M <= 2*10^5; unique tight case M=2."

We verify the inequality with an independent totient sieve (and cross-check it
against sympy.totient), confirm M=2 is the unique equality case, verify the
prime-power steps used in the proof, and confirm the derived scan bounds
2 n^4 quoted for the ledger instances.
"""
import mpmath as mp
import pytest
import sympy as sp

import _relcharge as R

x = R.x


def test_totient_bound_holds_and_is_tight_only_at_2():
    """2 phi(M)^2 >= M for all 1 <= M <= 2e5, with equality iff M = 2
    (ledger N).  Uses an independent sieve."""
    N = 200000
    phi = R.totient_sieve(N)
    violations = [M for M in range(1, N + 1) if 2 * phi[M] * phi[M] < M]
    tight = [M for M in range(1, N + 1) if 2 * phi[M] * phi[M] == M]
    assert violations == []
    assert tight == [2]


def test_sieve_agrees_with_sympy_totient():
    """Cross-engine: the hand-rolled sieve matches sympy.totient."""
    phi = R.totient_sieve(3000)
    for M in range(1, 3001):
        assert phi[M] == int(sp.totient(M))


def test_totient_lower_bound_symbolic_form():
    """phi(M) >= sqrt(M/2) is equivalent to 2 phi(M)^2 >= M; verify the
    equivalent integer inequality directly on a dense range."""
    phi = R.totient_sieve(20000)
    for M in range(1, 20001):
        assert mp.mpf(phi[M]) >= mp.sqrt(mp.mpf(M) / 2) - mp.mpf(10) ** -30


def test_prime_power_steps_of_the_proof():
    """The proof bounds phi on prime powers:
      * 2-part, a>=1: phi(2^a)=2^{a-1} >= sqrt(2^a/2)=2^{(a-1)/2};
      * odd prime power: phi(p^b)=p^{b-1}(p-1) >= p^{b/2}=sqrt(p^b),
        using p-1 >= sqrt(p) for p >= 3."""
    for a in range(1, 25):
        assert sp.totient(2**a) == 2 ** (a - 1)
        assert mp.mpf(2) ** (a - 1) >= mp.sqrt(mp.mpf(2) ** a / 2) - mp.mpf(10) ** -30
    for p in sp.primerange(3, 200):
        assert sp.Integer(p - 1) ** 2 >= p  # p-1 >= sqrt(p)
        for b in range(1, 5):
            val = sp.totient(p**b)
            assert val**2 >= p**b  # phi(p^b) >= sqrt(p^b)


def test_phi_M_at_most_degP_forces_M_bound():
    """Combining the two inequalities: phi(M) <= deg P forces M <= 2(deg P)^2.
    Check for a range of degree bounds d that every M with phi(M) <= d indeed
    satisfies M <= 2 d^2."""
    phi = R.totient_sieve(20000)
    for d in [1, 2, 5, 10, 16, 40]:
        Ms = [M for M in range(1, 2 * d * d + 1) if phi[M] <= d]
        # the claim: no M with phi(M) <= d exceeds 2 d^2
        beyond = [M for M in range(2 * d * d + 1, 20001) if phi[M] <= d]
        assert beyond == []
        assert max(Ms) <= 2 * d * d


@pytest.mark.parametrize(
    "degRat,expected_bound",
    [
        (16, 512),      # Rat_beta4 (ledger G)
        (100, 20000),   # Rat_Lehmer (ledger H)
        (36, 2592),     # Rat_S6 (ledger P)
        (64, 8192),     # Rat_S8 (ledger P)
        (144, 41472),   # census degree-12 (ledger T)
        (40, 3200),     # mixed Rat_{L,beta4} (ledger I)
        (256, 131072),  # nested Rat_{Rat_p} for x^4-x+1 (ledger X)
    ],
)
def test_quoted_scan_bounds_equal_2_deg_squared(degRat, expected_bound):
    """Every scan bound quoted in the paper equals 2 (deg Rat)^2."""
    assert 2 * degRat * degRat == expected_bound
