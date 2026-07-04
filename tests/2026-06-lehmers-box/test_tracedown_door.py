"""The trace-down face and the one door.  Section 6 (Def. 6.1, Lemma 6.2 Salem =
flip-straddle, Prop. 6.3) and Section 7 (circulant Prop. 7.1, Shoda Lemma 7.2,
self-action difference spectrum Lemma 7.4, and the exact closure guard Prop. 7.5).
"""
import random

import mpmath as mp
import numpy as np
import sympy as sp

from _helpers import (is_salem, is_palindromic, to_ab_sqrt5, sign_ab_sqrt5)

mp.mp.dps = 50
x, t = sp.symbols('x t')
PHI_SYM = (1 + sp.sqrt(5)) / 2


def _tracedown(int_coeffs):
    """Def. 6.1: for reciprocal m_theta of degree 2m return T with
    m_theta(x) = x^m T(x + 1/x)."""
    p = sp.Poly(int_coeffs, x)
    d = p.degree()
    m = d // 2
    coeffs_hi = p.all_coeffs()
    a = {d - i: coeffs_hi[i] for i in range(len(coeffs_hi))}

    def s(j):                     # x^j + x^-j as a polynomial in t
        if j == 0:
            return sp.Integer(2)
        sm2, sm1 = sp.Integer(2), t
        for _ in range(2, j + 1):
            sm2, sm1 = sm1, sp.expand(t * sm1 - sm2)
        return sm1 if j >= 1 else sm2

    T = sp.Integer(a.get(m, 0))
    for j in range(1, m + 1):
        T += a.get(m + j, 0) * s(j)
    return sp.expand(T)


def test_tracedown_identity():
    """Def. 6.1: the trace-down T satisfies m_theta(x) = x^m * T(x + 1/x)."""
    for coeffs in ([1, -1, -1, -1, 1],                     # beta_4
                   [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]):  # Lehmer
        T = _tracedown(coeffs)
        m = (len(coeffs) - 1) // 2
        recon = sp.expand(x ** m * T.subs(t, x + 1 / x))
        recon = sp.expand(recon * 1)     # clear x^{-.} into polynomial form
        diff = sp.simplify(recon - sp.Poly(coeffs, x).as_expr())
        assert diff == 0, coeffs


def _tracedown_root_pattern(int_coeffs):
    T = _tracedown(int_coeffs)
    Tc = [complex(c) for c in sp.Poly(T, t).all_coeffs()]
    roots = np.roots(Tc)
    totally_real = all(abs(r.imag) < 1e-9 for r in roots)
    above2 = sum(1 for r in roots if r.real > 2 + 1e-9)
    inside = sum(1 for r in roots if -2 + 1e-9 < r.real < 2 - 1e-9)
    return totally_real, above2, inside, roots


def test_salem_is_flip_straddle():
    """Lemma 6.2 (P2-DELTA-02): a reciprocal integer theta is Salem iff its
    trace-down T is totally real with exactly one root in (2,inf) and the rest in
    (-2,2).  Verified on Lehmer (deg 10 -> T deg 5) and beta_4 (deg 4 -> T deg 2)."""
    for coeffs in ([1, -1, -1, -1, 1], [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]):
        assert is_salem(coeffs)
        m = (len(coeffs) - 1) // 2
        tot_real, above2, inside, _ = _tracedown_root_pattern(coeffs)
        assert tot_real
        assert above2 == 1
        assert inside == m - 1


def test_flip_boundary_discriminant():
    """Lemma 6.2 proof: for x^2 - t x + 1 the discriminant is D = t^2 - 4;
    D < 0 on the circle (|t|<2), D > 0 off it (|t|>2), D = 0 at t = +/-2."""
    for tv in (mp.mpf("1.3"), mp.mpf("0"), mp.mpf("-1.7")):
        assert tv ** 2 - 4 < 0
    for tv in (mp.mpf("2.3"), mp.mpf("-3")):
        assert tv ** 2 - 4 > 0
    assert mp.mpf(2) ** 2 - 4 == 0


def test_fourth_roots_map_to_cyclotomic_lattice():
    """Prop. 6.3 (P2-DELTA-01): under rho(z)=z+1/z the fourth roots of unity map
    to the cyclotomic lattice {2, 0, -2}; a genuine Salem straddle instead needs
    an interior value t in (-2,2) NOT in {2,0,-2} (irrational angle)."""
    rho = {z: sp.simplify(z + 1 / z) for z in (sp.Integer(1), sp.I,
                                               sp.Integer(-1), -sp.I)}
    assert set(rho.values()) == {2, sp.Integer(0), -2}
    # Lehmer's trace-down has interior roots strictly off the cyclotomic lattice
    _, _, _, roots = _tracedown_root_pattern([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
    interior = [r.real for r in roots if -2 < r.real < 2]
    assert all(min(abs(v - c) for c in (-2, 0, 2)) > 1e-3 for v in interior)


def test_circulant_eigenvalue_formula_and_no_salem():
    """Prop. 7.1 (P2-CIRC-01): eigenvalues of circ(c_0..c_{n-1}) are
    lambda_j = sum_k c_k omega^{jk} (omega = e^{2 pi i/n}); and no integer
    circulant eigenvalue is a Salem number (its char-poly factors carry no Salem
    minimal polynomial)."""
    random.seed(3)
    for _ in range(6):
        n = random.randint(3, 6)
        c = [random.randint(-3, 3) for _ in range(n)]
        circ = np.array([[c[(j - k) % n] for k in range(n)] for j in range(n)],
                        dtype=complex)
        ev_np = list(np.linalg.eigvals(circ))
        w = mp.e ** (2j * mp.pi / n)
        ev_f = [complex(sum(c[k] * w ** (j * k) for k in range(n)))
                for j in range(n)]
        # match the two eigenvalue multisets to within tolerance (avoids
        # fragile sort ties on complex-conjugate pairs)
        remaining = ev_np[:]
        for b in ev_f:
            j = min(range(len(remaining)), key=lambda idx: abs(remaining[idx] - b))
            assert abs(remaining[j] - b) < 1e-9
            remaining.pop(j)
        assert remaining == []
        # no Salem factor in the (exact integer) characteristic polynomial
        cp = sp.Matrix([[c[(j - k) % n] for k in range(n)]
                        for j in range(n)]).charpoly(x)
        for fac, _mult in sp.factor_list(cp.as_expr())[1]:
            fc = [int(cc) for cc in sp.Poly(fac, x).all_coeffs()]
            if fc and fc[0] < 0:
                fc = [-v for v in fc]
            assert not is_salem(fc)


def test_shoda_commutator_trace_zero():
    """Lemma 7.2 (Shoda), necessary direction: every commutator [X,Y]=XY-YX is
    traceless.  A traceless matrix can therefore carry any char poly, e.g. a
    'foreign' carrier of Lehmer (the one door)."""
    random.seed(5)
    for _ in range(30):
        n = random.randint(2, 5)
        X = np.random.randint(-3, 4, (n, n))
        Y = np.random.randint(-3, 4, (n, n))
        comm = X @ Y - Y @ X
        assert abs(np.trace(comm)) < 1e-9
    # a traceless integer carrier of Lehmer: companion(L) (trace -1) plus [1]
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    comp = sp.Matrix(sp.Poly(L, x).all_coeffs()[::-1][:-1])  # placeholder
    companion = sp.zeros(10, 10)
    for i in range(9):
        companion[i + 1, i] = 1
    for i, coeff in enumerate(L[1:][::-1]):
        companion[i, 9] = -coeff
    carrier = sp.zeros(11, 11)
    carrier[:10, :10] = companion
    carrier[10, 10] = 1                     # add [1] to zero the trace
    assert carrier.trace() == 0             # eligible commutator by Shoda
    cp = carrier.charpoly(x).as_expr()
    quo, rem = sp.div(cp, sp.Poly(L, x).as_expr(), x)
    assert sp.simplify(rem) == 0            # Lehmer divides its char poly


def test_selfaction_difference_spectrum():
    """Lemma 7.4 (P2-UNIF-01): spectrum of ad_R = [R,.] is the difference set
    {lambda_i - lambda_j}; for the golden seed R=C(x^2-x-1) it is {0, +sqrt5,
    -sqrt5}."""
    R = np.array([[0, 1], [1, 1]], dtype=float)     # companion of x^2 - x - 1
    I2 = np.eye(2)
    adR = np.kron(R, I2) - np.kron(I2, R.T)
    ev = sorted(np.linalg.eigvals(adR).real)
    s5 = float(mp.sqrt(5))
    assert abs(ev[0] + s5) < 1e-9
    assert abs(ev[-1] - s5) < 1e-9
    assert sum(1 for v in ev if abs(v) < 1e-9) == 2   # two zeros
    # general: eigenvalues are exactly pairwise differences of R's eigenvalues
    lam = np.linalg.eigvals(R)
    diffs = sorted(round(a - b, 8) for a in lam for b in lam)
    assert sorted(round(v, 8) for v in ev) == diffs


# -- Prop. 7.5  the exact closure guard --------------------------------------
def _validate(int_coeffs):
    """Prop. 7.5 guard: FORCED (no Salem factor), FORCED_ABOVE_FLOOR (Salem
    factor with beta >= phi), INVALID_CLOSURE (Salem factor with beta < phi).
    The beta vs phi decision is the exact sign of m_beta(phi) in Q(sqrt5)."""
    verdict = "FORCED"
    poly = sp.Poly(int_coeffs, x)
    for fac, _m in sp.factor_list(poly.as_expr())[1]:
        fc = [int(c) for c in sp.Poly(fac, x).all_coeffs()]
        if fc and fc[0] < 0:
            fc = [-v for v in fc]
        if is_salem(fc):
            m_beta_phi = sp.Poly(fc, x).as_expr().subs(x, PHI_SYM)
            a, b = to_ab_sqrt5(m_beta_phi)
            s = sign_ab_sqrt5(a, b)
            if s > 0:                       # m_beta(phi) > 0  <=>  beta < phi
                return "INVALID_CLOSURE"
            else:                           # beta >= phi
                verdict = "FORCED_ABOVE_FLOOR"
    return verdict


def test_guard_ladder_three_verdicts():
    """Prop. 7.5 (P2-GUARD-01..05): framework objects read FORCED; a planted
    minimal degree-four Salem beta_4 reads FORCED_ABOVE_FLOOR; a carrier of
    Lehmer's number reads INVALID_CLOSURE."""
    # framework / catalog objects: no Salem factor
    assert _validate([1, -1, -1]) == "FORCED"                 # phi seed
    assert _validate([1, -7, 1]) == "FORCED"                  # gap
    assert _validate([1, 0, 5, 0, -5]) == "FORCED"            # K (non-reciprocal)
    # beta_4: Salem, above the floor
    assert _validate([1, -1, -1, -1, 1]) == "FORCED_ABOVE_FLOOR"
    # Lehmer: Salem, below the floor
    assert _validate([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]) == "INVALID_CLOSURE"


def test_guard_sign_is_exact_in_Q_sqrt5():
    """Prop. 7.5 mechanism: the sign of a + b sqrt5 is decided with no floating
    point -- shared signs give that sign; opposite signs compare a^2 with 5 b^2."""
    # same-sign short-circuit
    assert sign_ab_sqrt5(sp.Rational(3), sp.Rational(1)) == 1
    assert sign_ab_sqrt5(sp.Rational(-2), sp.Rational(-5)) == -1
    # opposite signs decided by a^2 vs 5 b^2
    assert sign_ab_sqrt5(sp.Rational(3), sp.Rational(-1)) == 1     # 9 > 5
    assert sign_ab_sqrt5(sp.Rational(2), sp.Rational(-1)) == -1    # 4 < 5
    assert sign_ab_sqrt5(sp.Rational(1), sp.Rational(-2)) == -1    # 1 < 20
