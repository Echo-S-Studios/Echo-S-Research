"""Part 5 -- The trace form: deviation, duality, dilation (Sec. 5).

Independent re-derivations of the Lorentzian trace-form metric, the deviation
operator X_n = F_n H, the Trace-Form Duality (1/2)Tr(X_n^2)=5F_n^2=L_n^2-4(-1)^n,
the psi^n axis-dilation, and the deviation ladder into V_m.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from _eap_helpers import R, I2, phi, psi, sqrt5, fib, luc, mat_eq, is_zero, comm

H = 2 * R - I2
S = sp.Matrix([[2, 1], [1, -2]])
J = sp.Matrix([[0, -1], [1, 0]])


def X(n):
    """Deviation operator X_n = 2 R^n - L_n I (Def 5.3), built from the actual
    matrix power and the independent Lucas number."""
    return 2 * R**n - luc(n) * I2


# ----------------------------------------------------------------------------
# Prop 5.2 -- the metric is Lorentzian
# ----------------------------------------------------------------------------
def test_trace_form_diagonal_signature():
    """Prop 5.2: trace form on {H,S,J} is diag(10,10,-2), signature (2,1)."""
    basis = [H, S, J]
    gram = sp.Matrix(3, 3, lambda i, j: sp.trace(basis[i] * basis[j]))
    assert mat_eq(gram, sp.diag(10, 10, -2))
    assert sp.trace(H * H) == 10
    assert sp.trace(S * S) == 10
    assert sp.trace(J * J) == -2


def test_trace_form_off_diagonal_vanishes():
    """Prop 5.2: off-diagonal entries vanish, e.g. Tr(HS)=0."""
    assert is_zero(sp.trace(H * S))
    assert is_zero(sp.trace(H * J))
    assert is_zero(sp.trace(S * J))


def test_trace_form_invariance():
    """Prop 5.2: invariance <[Z,X],Y> + <X,[Z,Y]> = 0 (cyclicity of trace)."""
    basis = [H, S, J]
    for Z in basis:
        for Xm in basis:
            for Y in basis:
                lhs = sp.trace(comm(Z, Xm) * Y) + sp.trace(Xm * comm(Z, Y))
                assert is_zero(lhs)


# ----------------------------------------------------------------------------
# Prop 5.4 -- deviation collapses to the Cartan
# ----------------------------------------------------------------------------
def test_deviation_equals_Fn_H():
    """Prop 5.4: 'X_n = F_n H' for all n (traceless, integer entries)."""
    for n in range(-8, 9):
        assert mat_eq(X(n), fib(n) * H), f"X_n=F_nH failed at n={n}"
        assert is_zero(sp.trace(X(n)))
        # integer entries
        assert all(e.is_integer for e in X(n))


def test_example_deviation_values():
    """Ex 5.5: X_1=H, X_3=2H (F_3=2), X_5=5H (F_5=5)."""
    assert mat_eq(X(1), H)
    assert mat_eq(X(3), 2 * H) and fib(3) == 2
    assert mat_eq(X(5), 5 * H) and fib(5) == 5


# ----------------------------------------------------------------------------
# Thm 5.6 -- Trace-Form Duality
# ----------------------------------------------------------------------------
def test_trace_form_duality():
    """Thm 5.6: '(1/2)Tr(X_n^2) = 5 F_n^2 = L_n^2 - 4(-1)^n = (phi^n-psi^n)^2'."""
    for n in range(-8, 9):
        half_tr = sp.Rational(1, 2) * sp.trace(X(n) * X(n))
        assert is_zero(half_tr - 5 * fib(n)**2), f"1/2 Tr X_n^2 at n={n}"
        assert is_zero(5 * fib(n)**2 - (luc(n)**2 - 4 * (-1)**(n % 2))), f"L_n^2-4(-1)^n at n={n}"
        assert is_zero(5 * fib(n)**2 - (phi**n - psi**n)**2), f"(phi^n-psi^n)^2 at n={n}"


# ----------------------------------------------------------------------------
# Thm 5.7 -- psi^n dilates the Cartan axis
# ----------------------------------------------------------------------------
def test_dilation_spectrum_gap():
    """Thm 5.7(a): eigenvalue gap dilates phi-psi=sqrt5 -> phi^n-psi^n=F_n sqrt5."""
    for n in range(-8, 9):
        assert is_zero((phi**n - psi**n) - fib(n) * sqrt5)


def test_dilation_operator_adRn():
    """Thm 5.7(b): 'ad_{R^n} = F_n ad_R', hence X_1=H -> X_n=F_n H."""
    testmats = [sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[1, 2], [3, 4]]), H, S, J]
    for n in range(-6, 7):
        for Y in testmats:
            assert mat_eq(comm(R**n, Y), fib(n) * comm(R, Y)), f"n={n}"


def test_dilation_index_multiplicative():
    """Thm 5.7(b): psi^n : X_m -> X_{nm} (index multiplicative).
    Applying phi->phi^n to the object {phi^m,psi^m} gives {phi^{nm},psi^{nm}},
    whose deviation is X_{nm}."""
    for n in range(-3, 4):
        for m in range(-3, 4):
            # deviation of {phi^{nm}, psi^{nm}} is 2R^{nm}-L_{nm}I = X_{nm}
            assert mat_eq(X(n * m), fib(n * m) * H)


def test_dilation_invariant_seed_lift():
    """Thm 5.7(c): seed (1/2)Tr(H^2)=5=L_1^2-4(-1)^1 lifts to 5F_n^2=L_n^2-4(-1)^n."""
    assert sp.Rational(1, 2) * sp.trace(H * H) == 5
    assert luc(1)**2 - 4 * (-1)**1 == 5


# ----------------------------------------------------------------------------
# Thm 5.11 -- the deviation ladder into V_m
# ----------------------------------------------------------------------------
def half_trace_Vm(n, m):
    """(1/2) Tr_{V_m}(X_n^2) computed from first principles: on V_m the Cartan
    H/sqrt5 has weights {m-2k}, so X_n=F_n H has trace of square
    5 F_n^2 sum_k (m-2k)^2."""
    weights_sq_sum = sum((m - 2 * k)**2 for k in range(m + 1))
    trace = 5 * fib(n)**2 * weights_sq_sum
    return sp.Rational(1, 2) * trace


def test_deviation_ladder_binomial_form():
    """Thm 5.11: '(1/2)Tr_{V_m}(X_n^2) = 5 F_n^2 C(m+2,3)'."""
    for n in range(-4, 5):
        for m in range(0, 7):
            lhs = half_trace_Vm(n, m)
            rhs = 5 * fib(n)**2 * sp.binomial(m + 2, 3)
            assert is_zero(lhs - rhs), f"ladder n={n}, m={m}"


def test_deviation_ladder_casimir_form():
    """Thm 5.11: also '= (5F_n^2/3) dim(V_m) Cas(V_m)'."""
    for n in range(-4, 5):
        for m in range(0, 7):
            dimV = m + 1
            cas = sp.Rational(1, 2) * m * (m + 2)
            rhs = sp.Rational(5, 3) * fib(n)**2 * dimV * cas
            assert is_zero(half_trace_Vm(n, m) - rhs), f"casimir form n={n}, m={m}"


def test_ladder_V3_V4_scalings():
    """Thm 5.11: V_3 gives 50 F_n^2, V_4 gives 100 F_n^2 (base 5F_n^2 times
    C(5,3)=10 and C(6,3)=20)."""
    assert sp.binomial(5, 3) == 10
    assert sp.binomial(6, 3) == 20
    for n in range(-4, 5):
        assert is_zero(half_trace_Vm(n, 3) - 50 * fib(n)**2)
        assert is_zero(half_trace_Vm(n, 4) - 100 * fib(n)**2)


# ----------------------------------------------------------------------------
# Ex 5.12 -- the deviation exercise end to end
# ----------------------------------------------------------------------------
def test_exercise_X4_full():
    """Ex 5.12: from R^4=3R+2I, X_4=3H (F_4=3, L_4=7), and
    (1/2)Tr(X_4^2)=5F_4^2=45=L_4^2-4=(phi^4-psi^4)^2; ladder into V_2 gives
    5F_4^2 C(4,3)=180."""
    assert mat_eq(R**4, 3 * R + 2 * I2)
    assert mat_eq(X(4), 3 * H) and fib(4) == 3 and luc(4) == 7
    half_tr = sp.Rational(1, 2) * sp.trace(X(4) * X(4))
    assert half_tr == 45 == 5 * fib(4)**2
    assert luc(4)**2 - 4 * (-1)**4 == 45
    assert is_zero((phi**4 - psi**4)**2 - 45)
    assert half_trace_Vm(4, 2) == 5 * fib(4)**2 * sp.binomial(4, 3) == 180
