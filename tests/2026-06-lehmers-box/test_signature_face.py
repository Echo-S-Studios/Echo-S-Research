"""The signature face.  Section 5: Prop. 5.1 (Salem forces a Lorentzian
signature (2,m-1)), Lemma 5.2 ([K:Q]=16), Theorem 5.3 (the 27-subfield signature
census of K, re-derived by an independent Galois computation), and Cor. 5.5
(beta_4 > phi decided exactly by a sign in Q(sqrt5); complex place off circle).
"""
import random

import mpmath as mp
import numpy as np
import sympy as sp
import pytest

from _helpers import (is_salem, root_classification, to_ab_sqrt5,
                      sign_ab_sqrt5, PHI)

mp.mp.dps = 50
x = sp.symbols('x')


# ---------------------------------------------------------------------------
# Prop. 5.1  Salem -> signature (2, m-1); trace form (m+1, m-1)
# ---------------------------------------------------------------------------
def _integer_power_sums(int_coeffs, upto):
    """p_n = sum of roots^n (an integer for a monic integer polynomial)."""
    roots = mp.polyroots([mp.mpf(c) for c in int_coeffs],
                         maxsteps=500, extraprec=400)
    out = []
    for n in range(upto + 1):
        s = sum(r ** n for r in roots)
        val = mp.re(s)
        rounded = int(mp.nint(val))
        assert abs(val - rounded) < mp.mpf(10) ** -8   # is a rational integer
        out.append(rounded)
    return out


def _trace_form_signature(int_coeffs):
    """Signature (pos, neg) of the trace form Tr(theta^{i+j}) on 1..theta^{d-1}."""
    d = len(int_coeffs) - 1
    p = _integer_power_sums(int_coeffs, 2 * (d - 1))
    M = np.array([[p[i + j] for j in range(d)] for i in range(d)], dtype=float)
    ev = np.linalg.eigvalsh(M)
    pos = int(np.sum(ev > 1e-6))
    neg = int(np.sum(ev < -1e-6))
    return pos, neg


@pytest.mark.parametrize("coeffs,m", [
    ([1, -1, -1, -1, 1], 2),                          # beta_4, degree 4
    ([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], 5),      # Lehmer, degree 10
])
def test_salem_signature_and_trace_form(coeffs, m):
    """Prop. 5.1 (P2-NL-02/03): a Salem number of degree 2m has unit-rank
    signature (r1,r2)=(2,m-1) and trace-form signature (m+1,m-1), indefinite."""
    assert is_salem(coeffs)
    out, ins, onc = root_classification(coeffs)
    r1 = out + ins                # the two real embeddings beta, 1/beta
    r2 = onc // 2                 # on-circle conjugate pairs
    assert (r1, r2) == (2, m - 1)
    assert r2 >= 1                # indefinite for m >= 2
    assert _trace_form_signature(coeffs) == (m + 1, m - 1)


# ---------------------------------------------------------------------------
# Lemma 5.2  [K : Q] = 16
# ---------------------------------------------------------------------------
def test_K_degree_16():
    """Lemma 5.2 (P2-UNIF-02): K = Q(sqrt2, sqrt3, 5^(1/4)) has degree
    2*2*4 = 16; verified via the minimal polynomial of a primitive element."""
    assert 2 * 2 * 4 == 16
    theta = sp.sqrt(2) + sp.sqrt(3) + sp.root(5, 4)
    mpoly = sp.minimal_polynomial(theta, x)
    assert sp.degree(mpoly, x) == 16


def test_catalog_real_differences_lie_in_K():
    """Lemma 5.2: real differences of catalog eigenvalues lie in K.  Spot check:
    phi - psi = sqrt5 and (sqrt2) - (sqrt3) both live in K (verified by exhibiting
    them as elements of Q(sqrt2,sqrt3,5^(1/4)))."""
    # phi - psi where phi,psi are roots of x^2-x-1:  phi-psi = sqrt5 = (5^(1/4))^2
    assert sp.simplify((sp.root(5, 4)) ** 2 - sp.sqrt(5)) == 0
    # sqrt2 - sqrt3 is visibly in K; its square 5 - 2 sqrt6 is in K too
    val = (sp.sqrt(2) - sp.sqrt(3)) ** 2
    assert sp.simplify(val - (5 - 2 * sp.sqrt(6))) == 0


# ---------------------------------------------------------------------------
# Theorem 5.3  the 27-subfield signature census (independent Galois computation)
# ---------------------------------------------------------------------------
def _build_group():
    """G = Gal(K(i)/Q) ~ C2 x C2 x D4 as tuples g=(e1,e2,k,d):
    sqrt2->(-1)^e1 sqrt2; sqrt3->(-1)^e2 sqrt3; 5^(1/4)->i^k 5^(1/4); i->(-1)^d i.
    Composition (ga after gb): e's,d add mod2; k'' = kb*(1+2 da) + ka  (mod4)."""
    elems = [(e1, e2, k, d) for e1 in (0, 1) for e2 in (0, 1)
             for k in (0, 1, 2, 3) for d in (0, 1)]

    def mul(ga, gb):
        e1a, e2a, ka, da = ga
        e1b, e2b, kb, db = gb
        return ((e1a + e1b) % 2, (e2a + e2b) % 2,
                (kb * (1 + 2 * da) + ka) % 4, (da + db) % 2)

    return elems, mul


def _all_subgroups(elems, mul, E):
    def close(gens):
        S = set([E]) | set(gens)
        changed = True
        while changed:
            changed = False
            for a in list(S):
                for b in list(S):
                    p = mul(a, b)
                    if p not in S:
                        S.add(p)
                        changed = True
        return frozenset(S)

    subs = set([frozenset([E])])
    frontier = [frozenset([E])]
    while frontier:
        H = frontier.pop()
        for g in elems:
            if g in H:
                continue
            H2 = close(list(H) + [g])
            if H2 not in subs:
                subs.add(H2)
                frontier.append(H2)
    return subs


def _galois_census():
    elems, mul = _build_group()
    E = (0, 0, 0, 0)
    c = (0, 0, 0, 1)                       # complex conjugation
    # group axioms
    assert len(elems) == 32
    assert all(mul(E, g) == g and mul(g, E) == g for g in elems)
    random.seed(0)
    for _ in range(3000):
        a, b, cc = (random.choice(elems), random.choice(elems),
                    random.choice(elems))
        assert mul(mul(a, b), cc) == mul(a, mul(b, cc))

    def inv(g):
        for h in elems:
            if mul(g, h) == E:
                return h
        raise AssertionError

    assert mul(c, c) == E
    J = frozenset([E, c])                  # Gal(K(i)/K)

    subs = _all_subgroups(elems, mul, E)
    subK = [H for H in subs if J <= H]     # subfields of K

    def signature(H):
        Hs = set(H)
        deg = 32 // len(H)
        seen, reps = set(), []
        for g in elems:
            coset = frozenset(mul(g, h) for h in H)
            if coset not in seen:
                seen.add(coset)
                reps.append(g)
        r1 = sum(1 for g in reps if mul(mul(inv(g), c), g) in Hs)
        return deg, r1, (deg - r1) // 2

    census = {}
    for H in subK:
        sig = signature(H)
        census[sig] = census.get(sig, 0) + 1
    return census, len(subK), len(subs)


PAPER_CENSUS = {(1, 1, 0): 1, (2, 2, 0): 7, (4, 2, 1): 4, (4, 4, 0): 7,
                (8, 4, 2): 6, (8, 8, 0): 1, (16, 8, 4): 1}


def test_signature_census_matches():
    """Theorem 5.3 (P2-UNIF-03): the 27 subfields of K have exactly the census
    (deg,r1,r2): (1,1,0):1, (2,2,0):7, (4,2,1):4, (4,4,0):7, (8,4,2):6,
    (8,8,0):1, (16,8,4):1.  Re-derived by Galois correspondence on K(i)."""
    census, n_subfields, _ = _galois_census()
    assert census == PAPER_CENSUS
    assert n_subfields == 27
    assert sum(census.values()) == 27
    # the Salem shape (2,1) occurs for exactly four quartic subfields
    assert census[(4, 2, 1)] == 4
    # K itself has signature (8,4)
    assert census[(16, 8, 4)] == 1


def test_paper_27_is_the_subfield_count():
    """Theorem 5.3 backing note (corrected 2026-07-04): now reads '27 subgroups
    fixing a subfield (the 27 subfields)'.  The 27 is the number of subgroups of
    G = C2 x C2 x D4 fixing a subfield of K -- i.e. the subfield count, 27 --
    not the total number of subgroups (which is 158, see next test)."""
    _, n_subfields, n_total = _galois_census()
    assert n_subfields == 27
    assert n_total == 158


def test_total_subgroup_count_is_158():
    """Records the true total subgroup count of G = C2 x C2 x D4 (=158), so the
    '27' in the backing note is understood as the subfield count."""
    _, n_subfields, n_total = _galois_census()
    assert n_total == 158
    assert n_subfields == 27


# ---------------------------------------------------------------------------
# Cor. 5.5  beta_4 > phi by exact sign in Q(sqrt5); complex place off circle
# ---------------------------------------------------------------------------
def test_beta4_above_phi_by_exact_sign_rule():
    """Cor. 5.5 (P2-UNIF-04): beta_4 > phi  <=>  m_{beta_4}(phi) < 0, and
    m_{beta_4}(phi) = phi^4 - phi^3 - phi^2 - phi + 1 is a NEGATIVE element of
    Q(sqrt5), decided with no floating point."""
    val = PHI ** 4 - PHI ** 3 - PHI ** 2 - PHI + 1
    a, b = to_ab_sqrt5(val)
    assert (a, b) == (sp.Rational(1, 2), sp.Rational(-1, 2))   # = (1 - sqrt5)/2
    assert sign_ab_sqrt5(a, b) == -1                            # negative => beta_4 > phi


def test_four_salem_quartics_are_only_salem_subfields():
    """Cor. 5.5: the only Salem-bearing subfields of K are the four quartics of
    signature (2,1); so any Salem realisable in K has degree 4, hence
    Mah >= beta_4 > phi."""
    census, _, _ = _galois_census()
    salem_shaped = {sig: n for sig, n in census.items()
                    if sig[1] == 2 and sig[2] == (sig[0] // 2) - 1 and sig[0] >= 4}
    # signature (2, m-1) with r1 = 2: among subfields the only one is (4,2,1)
    assert salem_shaped == {(4, 2, 1): 4}


def test_complex_place_off_the_unit_circle():
    """Sec 5 closing: the Lorentzian quartic Q(5^(1/4)) has its complex place
    OFF the circle -- |i beta| = 2.4195.. and |5^(1/4) i| = 1.4953.., both != 1,
    so it is not a Salem on-circle conjugate."""
    beta = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)      # imag root of x^4+5x^2-5
    assert abs(beta - mp.mpf("2.4195251530")) < mp.mpf(10) ** -8
    assert abs(mp.mpf(5) ** mp.mpf("0.25") - mp.mpf("1.4953487812")) < mp.mpf(10) ** -8
    assert abs(beta - 1) > mp.mpf("0.1")          # off the unit circle
    assert abs(mp.mpf(5) ** mp.mpf("0.25") - 1) > mp.mpf("0.1")
