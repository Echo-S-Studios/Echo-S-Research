"""Part 1 -- The associative core (Sec. 1 of the primer).

Independent re-derivations of the keystone relation, the Fibonacci power law,
the classical Fibonacci/Lucas identities, and the maximal-order claim.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from _eap_helpers import (R, I2, phi, psi, sqrt5, fib, luc, mat_eq, is_zero, comm)


# ----------------------------------------------------------------------------
# Prop 1.2 (the keystone) and eq. (2) (root5)
# ----------------------------------------------------------------------------
def test_char_poly_is_x2_minus_x_minus_1():
    """Def 1.1 / Prop 1.2: det(xI - R) = x^2 - x - 1."""
    x = sp.symbols('x')
    cp = sp.expand((x * I2 - R).det())
    assert is_zero(cp - (x**2 - x - 1))


def test_keystone_relation():
    """Prop 1.2, eq. (1): 'R^2 = R + I'.  Computed by direct multiplication."""
    assert mat_eq(R * R, R + I2)


def test_trace_and_det_of_R():
    """Prop 1.2: 'Tr R = 1, det R = -1'."""
    assert sp.trace(R) == 1
    assert R.det() == -1


def test_eigenvalues_are_phi_psi():
    """Prop 1.2: eigenvalues of R are phi, psi (roots of x^2-x-1)."""
    eigs = set(R.eigenvals().keys())
    assert eigs == {sp.nsimplify(phi), sp.nsimplify(psi)}


def test_vieta_relations():
    """Prop 1.2: 'phi + psi = 1, phi*psi = -1'."""
    assert is_zero(phi + psi - 1)
    assert is_zero(phi * psi + 1)


def test_root5_identity():
    """Eq. (2): 'phi - psi = sqrt5 = phi + phi^{-1}'."""
    assert is_zero((phi - psi) - sqrt5)
    assert is_zero((phi + 1 / phi) - sqrt5)


# ----------------------------------------------------------------------------
# Sec 1.3 hand computation of low powers
# ----------------------------------------------------------------------------
def test_low_powers_reduction():
    """Sec 1.3: R^2=R+I, R^3=2R+I, R^4=3R+2I, R^5=5R+3I.
    Each RHS is built from R,I with the stated integer coefficients and
    compared to the actual matrix power."""
    assert mat_eq(R**2, R + I2)
    assert mat_eq(R**3, 2 * R + I2)
    assert mat_eq(R**4, 3 * R + 2 * I2)
    assert mat_eq(R**5, 5 * R + 3 * I2)


# ----------------------------------------------------------------------------
# Thm 1.3 (power law) and eq. (3)
# ----------------------------------------------------------------------------
def test_power_law_scalar_form():
    """Thm 1.3, eq. (3): 'R^n = F_n R + F_{n-1} I' for all integer n.
    R**n computed by actual matrix exponentiation; RHS from the bare
    Fibonacci recurrence."""
    for n in range(-16, 17):
        lhs = R**n
        rhs = fib(n) * R + fib(n - 1) * I2
        assert mat_eq(lhs, rhs), f"power law failed at n={n}"


def test_power_law_matrix_entry_form():
    """Thm 1.3, eq. (3): 'R^n = [[F_{n-1}, F_n],[F_n, F_{n+1}]]'."""
    for n in range(-12, 13):
        entryform = sp.Matrix([[fib(n - 1), fib(n)], [fib(n), fib(n + 1)]])
        assert mat_eq(R**n, entryform), f"entry form failed at n={n}"


def test_example_R5_and_Rinv():
    """Ex 1.4: R^5 = [[3,5],[5,8]]; R^{-1}=[[-1,1],[1,0]] with R R^{-1}=I."""
    assert mat_eq(R**5, sp.Matrix([[3, 5], [5, 8]]))
    Rinv = fib(-1) * R + fib(-2) * I2       # F_{-1}=1, F_{-2}=-1
    assert mat_eq(Rinv, sp.Matrix([[-1, 1], [1, 0]]))
    assert mat_eq(R * Rinv, I2)


# ----------------------------------------------------------------------------
# Cor 1.5 -- the three classical identities
# ----------------------------------------------------------------------------
def test_lucas_trace():
    """Cor 1.5(a): 'Tr(R^n) = L_n = phi^n + psi^n'.
    Trace of the actual matrix power vs. independent Lucas recurrence,
    and vs. the eigenvalue power sum."""
    for n in range(-14, 15):
        tr = sp.expand(sp.trace(R**n))
        assert is_zero(tr - luc(n)), f"Lucas trace failed at n={n}"
        assert is_zero(luc(n) - (phi**n + psi**n)), f"L_n=phi^n+psi^n at n={n}"


def test_determinant_power():
    """Cor 1.5(b): 'det(R^n) = (-1)^n'."""
    for n in range(-14, 15):
        assert sp.expand((R**n).det()) == (-1)**(n % 2)


def test_cassini():
    """Cor 1.5(c): 'F_{n-1} F_{n+1} - F_n^2 = (-1)^n'."""
    for n in range(-14, 15):
        val = fib(n - 1) * fib(n + 1) - fib(n)**2
        assert val == (-1)**(n % 2), f"Cassini failed at n={n}"


def test_binet_eigenvalue_shadow():
    """Cor 1.5(d): 'phi^n = F_n phi + F_{n-1}' and 'F_n=(phi^n-psi^n)/sqrt5'."""
    for n in range(-14, 15):
        assert is_zero(phi**n - (fib(n) * phi + fib(n - 1)))
        assert is_zero(psi**n - (fib(n) * psi + fib(n - 1)))
        assert is_zero(fib(n) - (phi**n - psi**n) / sqrt5)


def test_example_cassini_n5():
    """Ex 1.6: at n=5, Tr=11=L_5; det=-1; F_4 F_6 - F_5^2 = -1; phi^5=5phi+3."""
    M = R**5
    assert sp.trace(M) == 11 == luc(5)
    assert M.det() == -1
    assert fib(4) * fib(6) - fib(5)**2 == -1
    assert is_zero(phi**5 - (5 * phi + 3))
    # decimal illustration phi^5 ~ 11.09
    assert abs(float(phi**5) - 11.09) < 0.01


def test_exercise_R_minus_2():
    """Ex 1.9: R^{-2}=[[2,-1],[-1,1]] two ways; det=1=(-1)^{-2}; Tr=3=L_{-2}."""
    way1 = fib(-2) * R + fib(-3) * I2       # F_{-2}=-1, F_{-3}=2
    way2 = (R**(-1))**2
    target = sp.Matrix([[2, -1], [-1, 1]])
    assert mat_eq(way1, target)
    assert mat_eq(way2, target)
    assert target.det() == 1 == sp.Integer(-1)**(-2)
    assert sp.trace(target) == 3 == luc(-2)


# ----------------------------------------------------------------------------
# Thm 1.7 (the golden order is the maximal order)
# ----------------------------------------------------------------------------
def test_I_and_R_independent_rank2():
    """Thm 1.7: {I,R} is Q-linearly independent, so Z[R] is free rank 2."""
    a, b = sp.symbols('a b')
    sol = sp.solve([sp.Eq(a * 1 + b * 0, 0),    # (1,1) entry: a
                    sp.Eq(a * 0 + b * 1, 0)],   # (1,2) entry: b
                   [a, b])
    assert sol == {a: 0, b: 0}


def test_order_discriminant_is_fundamental():
    """Thm 1.7: disc(x^2-x-1)=5, field disc(Q(sqrt5))=5 (since 5=1 mod 4),
    hence [O_K : Z[phi]]^2 = 5/5 = 1, index 1, so Z[R]=O_{Q(sqrt5)}."""
    x = sp.symbols('x')
    disc_poly = sp.discriminant(x**2 - x - 1, x)
    assert disc_poly == 5
    # field discriminant of Q(sqrt d) with d=5 squarefree, 5 = 1 (mod 4) -> d
    d = 5
    field_disc = d if d % 4 == 1 else 4 * d
    assert field_disc == 5
    index_sq = sp.Rational(disc_poly, field_disc)
    assert index_sq == 1
    assert sp.sqrt(index_sq) == 1


def test_clumsy_seed_2phi_has_index_2():
    """Intuition after Thm 1.7: the seed 2phi has defining polynomial of
    discriminant 20 = 2^2 * 5, giving a suborder of index 2."""
    x = sp.symbols('x')
    # minimal polynomial of 2*phi = 1 + sqrt5:  (x-1)^2 = 5  ->  x^2-2x-4
    twophi = sp.expand(2 * phi)
    minpoly = sp.minimal_polynomial(twophi, x)
    assert sp.expand(minpoly) == x**2 - 2 * x - 4
    disc = sp.discriminant(minpoly, x)
    assert disc == 20
    field_disc = 5
    assert sp.Rational(disc, field_disc) == 4      # index^2
    assert sp.sqrt(sp.Rational(disc, field_disc)) == 2
