"""Entry-level operators: commutator, circulant, Cartan
(Sec. 8.4 Prop. 8.5, Prop. 8.6; Sec. 8.3 Prop. 12.3 circulant; Prop. 12.4 Cartan;
Sec. 11 Prop. 11.1; App. A).

Reproduces the paper's scans structurally: no emitted entry-level matrix carries
a Salem factor, and the size-2 commutator Mahler image omits (1, mu_S).
"""
import random

import mpmath as mp
import numpy as np
import sympy as sp

import emgap_util as U

mp.mp.dps = 40

x = U.x
PHI = (1 + mp.sqrt(5)) / 2


def test_size2_commutator_mahler_image_omits_gap():
    """Prop. 8.5 / App. A: a traceless 2x2 integer matrix has charpoly x^2+det,
    with Mahler image {1} U [2, infinity) -- no value in (1, 2) ⊃ (1, mu_S),
    hence no Salem."""
    seen = set()
    for a in range(-6, 7):
        for b in range(-6, 7):
            for c in range(-6, 7):
                det = -a * a - b * c            # det of [[a,b],[c,-a]]
                m = U.mahler([1, 0, det])
                seen.add(int(round(float(m))))
                assert not (1 + mp.mpf("1e-9") < m < mp.mpf(2) - mp.mpf("1e-9"))
    # image realised: 1 and integers >= 2, nothing in between
    assert 1 in seen
    assert all(v == 1 or v >= 2 for v in seen)


def test_circulant_emits_no_salem():
    """Prop. 12.3 / App. A: integer circulants (n = 4,5,6) have eigenvalues in
    the cyclotomic (abelian, totally real or CM) field Q(zeta_n); a Salem field
    is neither. A scan of random integer circulants finds zero Salem factors."""
    rng = random.Random(20260604)
    salem_hits = 0
    for _ in range(400):
        n = rng.choice([4, 5, 6])
        row = [rng.randint(-3, 3) for _ in range(n)]
        C = sp.Matrix([[row[(j - i) % n] for j in range(n)] for i in range(n)])
        cp = C.charpoly(x)
        salem_hits += len(U.salem_factors(sp.Poly(cp.as_expr(), x)))
    assert salem_hits == 0


def test_commutator_over_catalog_emits_no_salem():
    """Prop. 8.6 / App. A: commutators [A,B] with A,B ranging over catalog
    companions and their kron / direct-sum products (sizes 2 and 4) carry no
    Salem-eligible factor."""
    catalog = [[1, -1, -1], [1, 1, -1], [1, 0, -2], [1, 0, -3],
               [1, 0, -5], [1, -7, 1], [1, 0, 5, 0, -5]]
    comps = [sp.Matrix(U.companion(c).astype(int).tolist()) for c in catalog]

    def dsum(A, B):
        n, m = A.shape[0], B.shape[0]
        M = sp.zeros(n + m, n + m)
        M[:n, :n] = A
        M[n:, n:] = B
        return M

    # build a generating slice: companions + a few size-4 products
    mats = list(comps)
    quad = [c for c in comps if c.shape[0] == 4]
    twos = [c for c in comps if c.shape[0] == 2]
    for A in twos[:3]:
        for B in twos[:3]:
            mats.append(sp.Matrix(np.kron(np.array(A.tolist()),
                                          np.array(B.tolist())).astype(int).tolist()))
            mats.append(dsum(A, B))

    salem_hits = 0
    recip_ge4 = 0
    for A in mats:
        for B in mats:
            if A.shape != B.shape:
                continue
            comm = A * B - B * A
            cp = sp.Poly(comm.charpoly(x).as_expr(), x)
            salem_hits += len(U.salem_factors(cp))
            for fac, _ in sp.factor_list(cp.as_expr(), x)[1]:
                Q = sp.Poly(fac, x)
                if Q.degree() >= 4 and Q.degree() % 2 == 0 and \
                        U.is_palindromic(Q.all_coeffs()):
                    recip_ge4 += 1
    assert salem_hits == 0


def test_cartan_An_eigenvalues_totally_real_in_0_4():
    """Prop. 12.4 / App. A: Cartan matrices A_n (n = 3,5,8) have eigenvalues
    2 - 2cos(k pi/(n+1)) in [0, 4], all real -- totally real, no on-circle
    point, so no Salem."""
    for n in (3, 5, 8):
        M = 2 * np.eye(n)
        for i in range(n - 1):
            M[i, i + 1] = -1
            M[i + 1, i] = -1
        ev = np.linalg.eigvalsh(M)
        assert np.all(ev >= -1e-9) and np.all(ev <= 4 + 1e-9)
        expected = sorted(2 - 2 * np.cos(k * np.pi / (n + 1)) for k in range(1, n + 1))
        for got, exp in zip(sorted(ev), expected):
            assert abs(got - exp) < 1e-9
        # no Salem factor in the (integer) Cartan characteristic polynomial
        cp = sp.Matrix(M.astype(int).tolist()).charpoly(x)
        assert not U.has_salem_factor(sp.Poly(cp.as_expr(), x))
