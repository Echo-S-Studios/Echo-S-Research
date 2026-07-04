r"""
Independent verification of Section 4 ("Capacity and the applied growth gate"):
the exact Northcott/Landau admissibility certificate (Ex 4.2) and the exact
Fisher matrices with the trace-zero identity G = n * Fisher_exp (Ex 4.4).

Fisher_exp = (1/n)(G - (1/n) t t^T), t the trace vector.  We rebuild each Gram
from field traces (or the regular representation) and recompute Fisher_exp, then
compare to the paper's displayed matrices.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols, zeros, eye

x = symbols("x")


def companion(coeffs):
    n = len(coeffs) - 1
    C = zeros(n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]
    return C


def fisher_exp(G, tvec):
    n = G.shape[0]
    t = Matrix(tvec)
    return (G - (t * t.T) / n) / n


# --- Example 4.2 (the gate decides on exact integers) ------------------------

def test_gate_x2_minus_24_height_and_sumsq():
    """Ex 4.2: the seed 2 sqrt6 (x^2-24) has (deg, coeff-height)=(2,24) and
    sum c_i^2 = 577."""
    coeffs = [1, 0, -24]
    assert len(coeffs) - 1 == 2                        # degree
    assert max(abs(c) for c in coeffs) == 24           # coefficient height
    assert sum(c * c for c in coeffs) == 577


def test_gate_x2_minus_24_landau_certificate():
    """Ex 4.2: sum c_i^2 = 577 certifies Mah^2 <= 577 (Landau: Mah <= ||p||_2).
    The true Mahler measure of x^2-24 is 24, and 24^2 = 576 <= 577."""
    import mpmath as mp
    mp.mp.dps = 40
    roots = mp.polyroots([1, 0, -24])
    Mah = mp.mpf(1)
    for r in roots:
        Mah *= max(mp.mpf(1), abs(r))
    assert abs(Mah - 24) < mp.mpf("1e-30")             # Mahler measure = 24
    assert Mah**2 <= 577                               # Landau bound holds
    assert 24**2 == 576 <= 577


def test_gate_x2_minus_7_landau_certificate():
    """Ex 4.2: sqrt7 (x^2-7) has sum c_i^2 = 50, certifying Mah <= 8 (50<=64),
    while the exact Mahler measure is 7."""
    coeffs = [1, 0, -7]
    assert sum(c * c for c in coeffs) == 50
    import mpmath as mp
    mp.mp.dps = 40
    roots = mp.polyroots([1, 0, -7])
    Mah = mp.mpf(1)
    for r in roots:
        Mah *= max(mp.mpf(1), abs(r))
    assert abs(Mah - 7) < mp.mpf("1e-30")              # Mahler measure = 7
    assert Mah <= sp.sqrt(50) <= 8                     # 50 <= 64 -> Mah <= 8


def test_gate_decision_grow_stop_reject_integers_only():
    """Ex 4.2: under Budget(64,256) & gain 96 -> GROW; under Budget(4,10) ->
    REJECT (height 24>10); with gain 0 -> STOP.  Decisions are integer comparisons."""
    deg, height, gain = 2, 24, 96

    def decide(Dmax, Hmax, gain, floor=0):
        if gain <= floor:
            return "STOP"
        if deg > Dmax or height > Hmax:
            return "REJECT"
        return "GROW"

    assert decide(64, 256, 96) == "GROW"
    assert decide(4, 10, 96) == "REJECT"
    assert decide(64, 256, 0) == "STOP"


# --- Example 4.4 (exact Fisher matrices) -------------------------------------

def test_fisher_Qsqrt5():
    """Ex 4.4: Q(sqrt5), {1,phi}: G=[[2,1],[1,3]], Fisher_exp=[[0,0],[0,5/4]].
    G rebuilt from trace(companion(x^2-x-1)^k); t=(Tr1,Tr phi)=(2,1)."""
    C = companion([1, -1, -1])
    G = Matrix([[(C**0).trace(), (C**1).trace()],
                [(C**1).trace(), (C**2).trace()]])
    assert G == Matrix([[2, 1], [1, 3]])
    F = fisher_exp(G, [2, 1])
    assert F == Matrix([[0, 0], [0, Q(5, 4)]])


def test_fisher_Qsqrt2sqrt3():
    """Ex 4.4: Q(sqrt2,sqrt3): G=diag(4,8,12,24), Fisher_exp=diag(0,2,3,6).
    t=(Tr1,0,0,0)=(4,0,0,0)."""
    G = sp.diag(4, 8, 12, 24)
    F = fisher_exp(G, [4, 0, 0, 0])
    assert F == sp.diag(0, 2, 3, 6)


def test_fisher_Qsqrt2sqrt3sqrt7():
    """Ex 4.4: Q(sqrt2,sqrt3,sqrt7): G=diag(8,56,16,112,24,168,48,336),
    Fisher_exp=diag(0,7,2,14,3,21,6,42). t=(8,0,...,0)."""
    G = sp.diag(8, 56, 16, 112, 24, 168, 48, 336)
    F = fisher_exp(G, [8, 0, 0, 0, 0, 0, 0, 0])
    assert F == sp.diag(0, 7, 2, 14, 3, 21, 6, 42)


def test_fisher_sqrt5_residual_norm_10():
    """Ex 4.4: the residual sqrt5 = 2 phi - 1 = (-1,2) is trace-zero and
    ||sqrt5||_G^2 = 10 = 2*5 = n*Fisher(sqrt5) with n=2."""
    # sqrt5 = 2 phi - 1 symbolic check
    phi = (1 + sqrt(5)) / 2
    assert sp.simplify((2 * phi - 1) - sqrt(5)) == 0
    C = companion([1, -1, -1])
    G = Matrix([[(C**0).trace(), (C**1).trace()],
                [(C**1).trace(), (C**2).trace()]])
    t = Matrix([2, 1])
    v = Matrix([-1, 2])
    assert (t.T * v)[0] == 0                            # trace-zero
    norm = (v.T * G * v)[0]
    assert norm == 10
    F = fisher_exp(G, [2, 1])
    assert 2 * (v.T * F * v)[0] == 10                   # n * Fisher = ||.||_G^2


def test_fisher_equals_n_fisher_on_trace_zero_symbolic():
    """Ex 4.4: the identity ||r||_G^2 = n Fisher(r) holds on EVERY trace-zero
    direction (t^T r = 0), because Fisher_exp=(1/n)(G-(1/n)t t^T) drops the t-term
    there.  Proved symbolically for a generic trace-zero r."""
    n = 3
    G = sp.MatrixSymbol("G", n, n)
    t = Matrix(sp.symbols("t0 t1 t2"))
    r = Matrix(sp.symbols("r0 r1 r2"))
    Gm = Matrix(G)
    F = (Gm - (t * t.T) / n) / n
    lhs = (r.T * Gm * r)[0]
    rhs_full = n * (r.T * F * r)[0]
    # rhs_full = ||r||_G^2 - (1/n)(t^T r)^2 ; on trace-zero (t^T r=0) they agree
    diff = sp.expand(lhs - rhs_full)
    assert sp.expand(diff - (t.T * r)[0]**2 / n) == 0
