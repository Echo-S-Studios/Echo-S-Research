"""
Character I -- the Mahler measure.

Theorem 4.2 (thm:measure) [OA-M-01..07]:
    M(A (+) B) = M(A) M(B) ,   M(psi^2 A) = M(A)^2 ,
    log M(A (x) B) = sum_{l in A} sum_{m in B} (log|l| + log|m|)^+   (tropical).
    Non-multiplicativity witness: golden (x) (x^2-2) has M = 2 phi^2 = 5.2360...,
    while M(golden) M(x^2-2) = 2 phi = 3.2360...

Theorem 4.3 (thm:floor) [OA-M-08; GA-01..13]:
    (+)-image = <phi, 2, 3, 5, beta^2> with the single relation beta^2 = phi^2 sqrt5;
    generators phi,2,3,5 mult. independent (norms -1,4,9,25); beta^2 * beta^2 = 5 phi^4;
    five atoms, rank four (not factorial); least generator phi -> no element in (1, phi).
"""

import mpmath as mp
import sympy as sp

from _opalg_ops import (
    golden_seed,
    K_seed,
    log_mahler_mp,
    mahler_exact,
    mahler_mp,
    modulus_mp,
    oplus,
    otimes,
    phi,
    psi,
    seed_from_poly,
    x,
)

sqrt5 = sp.sqrt(5)


# --------------------------------------------------------------------------
# additive on (+) ,  squares under psi^2
# --------------------------------------------------------------------------
def test_mahler_additive_on_oplus():
    """M(A (+) B) = M(A) M(B)  (exact) for several concrete seeds."""
    seeds = [
        golden_seed(),
        seed_from_poly(x**2 - 2),
        seed_from_poly(x**2 - 3),
        K_seed(),
    ]
    for A in seeds:
        for B in seeds:
            lhs = mahler_exact(oplus(A, B))
            rhs = sp.simplify(mahler_exact(A) * mahler_exact(B))
            assert sp.simplify(lhs - rhs) == 0


def test_mahler_squares_under_psi2():
    """M(psi^2 A) = M(A)^2  (exact)."""
    for A in [golden_seed(), seed_from_poly(x**2 - 2), K_seed()]:
        lhs = mahler_exact(psi(2, A))
        rhs = sp.simplify(mahler_exact(A) ** 2)
        assert sp.simplify(lhs - rhs) == 0


# --------------------------------------------------------------------------
# tropical on (x)
# --------------------------------------------------------------------------
def test_mahler_tropical_convolution():
    """log M(A (x) B) equals the (max,+) convolution of the SEPARATE log-spectra.

    LHS is built from the actual product object A(x)B; RHS is assembled only
    from |lambda| (A) and |mu| (B) via sum (log|l|+log|m|)^+ -- the two agree.
    """
    mp.mp.dps = 60
    seeds = [
        golden_seed(),
        seed_from_poly(x**2 - 2),
        seed_from_poly(x**2 - 3),
        K_seed(),
    ]
    for A in seeds:
        for B in seeds:
            lhs = log_mahler_mp(otimes(A, B), dps=60)
            rhs = mp.mpf(0)
            for a in A:
                la = mp.log(modulus_mp(a, 60))
                for b in B:
                    lb = mp.log(modulus_mp(b, 60))
                    rhs += max(mp.mpf(0), la + lb)
            assert abs(lhs - rhs) < mp.mpf(10) ** (-45)


def test_tensor_non_multiplicative_witness():
    """golden (x) (x^2-2):  M = 2 phi^2 (exact) != M(golden) M(x^2-2) = 2 phi."""
    A = golden_seed()
    B = seed_from_poly(x**2 - 2)
    M_prod = mahler_exact(otimes(A, B))
    # independently: 2 phi^2 = 3 + sqrt5
    assert sp.simplify(M_prod - 2 * phi**2) == 0
    assert sp.simplify(M_prod - (3 + sqrt5)) == 0
    # multiplicative would give 2 phi ; the two genuinely differ
    mult = sp.simplify(mahler_exact(A) * mahler_exact(B))
    assert sp.simplify(mult - 2 * phi) == 0
    assert sp.simplify(M_prod - mult) != 0


def test_witness_decimal_values():
    """Paper's quoted decimals: 2 phi^2 = 5.2360..., 2 phi = 3.2360... .

    Built independently from phi, then checked (i) against the closed forms
    3 + sqrt5 and 1 + sqrt5 to 40 digits, and (ii) truncated to the 4 decimals
    the paper actually prints.
    """
    mp.mp.dps = 50
    val_phi = (1 + mp.sqrt(5)) / 2
    two_phi2 = 2 * val_phi**2
    two_phi = 2 * val_phi
    assert mp.almosteq(two_phi2, 3 + mp.sqrt(5), rel_eps=mp.mpf(10) ** -40)
    assert mp.almosteq(two_phi, 1 + mp.sqrt(5), rel_eps=mp.mpf(10) ** -40)
    # truncation to 4 decimals reproduces the printed "5.2360" and "3.2360"
    assert int(mp.floor(two_phi2 * 10000)) == 52360
    assert int(mp.floor(two_phi * 10000)) == 32360


# --------------------------------------------------------------------------
# the floor monoid  <phi, 2, 3, 5, beta^2>
# --------------------------------------------------------------------------
def test_generators_are_realised_measures():
    """Each named generator is an exact Mahler measure of a concrete object."""
    assert sp.simplify(mahler_exact(golden_seed()) - phi) == 0
    assert sp.simplify(mahler_exact(seed_from_poly(x**2 - 2)) - 2) == 0
    assert sp.simplify(mahler_exact(seed_from_poly(x**2 - 3)) - 3) == 0
    assert sp.simplify(mahler_exact(seed_from_poly(x**2 - 5)) - 5) == 0


def test_beta2_equals_phi2_sqrt5():
    """beta^2 = M(K) = phi^2 sqrt5  for K = x^4 + 5x^2 - 5  (the single relation)."""
    beta2 = mahler_exact(K_seed())
    assert sp.simplify(beta2 - phi**2 * sqrt5) == 0
    # closed form (5 + 3 sqrt5)/2
    assert sp.simplify(beta2 - (5 + 3 * sqrt5) / 2) == 0


def test_beta4_non_unique_factorisation():
    """beta^2 * beta^2 = 5 * phi^4  -> the non-factorial witness."""
    beta2 = phi**2 * sqrt5
    assert sp.simplify(beta2**2 - 5 * phi**4) == 0


def test_generator_norms_in_Qsqrt5():
    """Norms in Q(sqrt5): N(phi) = -1, N(2) = 4, N(3) = 9, N(5) = 25."""

    def norm(expr):
        # a + b sqrt5  ->  a^2 - 5 b^2
        p = sp.Poly(sp.expand(expr), sqrt5)
        b = p.coeff_monomial(sqrt5)
        a = p.coeff_monomial(1)
        return sp.simplify(a**2 - 5 * b**2)

    assert norm(phi) == -1
    assert norm(sp.Integer(2)) == 4
    assert norm(sp.Integer(3)) == 9
    assert norm(sp.Integer(5)) == 25


def test_multiplicative_independence_phi_2_3_5():
    """No nonzero integer relation phi^a 2^b 3^c 5^d = 1 in a search box.

    Corroborates the norm argument: taking N(.) of any relation gives
    (-1)^a 2^{2b} 3^{2c} 5^{2d} = 1, forcing b=c=d=0 by unique factorisation,
    then phi^a = 1 forces a=0 (phi is not a root of unity).
    """
    mp.mp.dps = 60
    lphi = mp.log((1 + mp.sqrt(5)) / 2)
    l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
    R = 6
    found = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            for c in range(-R, R + 1):
                for d in range(-R, R + 1):
                    if a == b == c == d == 0:
                        continue
                    val = a * lphi + b * l2 + c * l3 + d * l5
                    if abs(val) < mp.mpf(10) ** (-25):
                        found.append((a, b, c, d))
    assert found == []


def test_rank_four_five_atoms():
    """5 atoms, rank 4: exactly one relation, log(beta^2) = 2 log phi + (1/2) log 5."""
    # the relation, exact:
    assert sp.simplify((phi**2 * sqrt5) - phi**2 * 5 ** sp.Rational(1, 2)) == 0
    # the 4 base atoms are independent (numeric, high precision):
    mp.mp.dps = 50
    logs = [mp.log((1 + mp.sqrt(5)) / 2), mp.log(2), mp.log(3), mp.log(5)]
    # no rational dependence with small denominators among the 4
    R = 5
    dep = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            for c in range(-R, R + 1):
                for d in range(-R, R + 1):
                    if a == b == c == d == 0:
                        continue
                    if abs(a * logs[0] + b * logs[1] + c * logs[2] + d * logs[3]) < mp.mpf(10) ** -20:
                        dep.append((a, b, c, d))
    assert dep == []


def test_cost_floor_no_element_in_open_1_phi():
    """Least generator is phi: every generator is >= phi, so the monoid meets
    (1, phi) in nothing -- the cost floor."""
    gens = [phi, sp.Integer(2), sp.Integer(3), sp.Integer(5), phi**2 * sqrt5]
    vals = [sp.N(g, 40) for g in gens]
    assert min(vals) == sp.N(phi, 40)  # phi is the least generator
    for v in vals:
        assert v >= sp.N(phi, 40)  # nothing strictly between 1 and phi is a generator
