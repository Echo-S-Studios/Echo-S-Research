"""eap_core -- the Emission Algebra A, as a clean reusable library.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
  "The Emission Algebra A -- one matrix, one relation, three layers."

This module packages the paper's computational procedures as importable
functions (no assertions, no I/O).  The producer scripts d1..d7 and
helix_orthogonal_partner import from here and EMIT machine-readable data.
It is an independent construction from the paper's premises:

  * golden field  phi, psi = (1 +/- sqrt5)/2, the seed R = [[0,1],[1,1]];
  * Fibonacci / Lucas from their bare recurrences (extended to n < 0);
  * the Lie-layer elements H, S, J, N_+, N_-, the deviation X_n = 2R^n - L_n I;
  * the semiring characters -- Mahler measure M and charge chi in Z/4Z.

Everything that gates a decision is settled in exact arithmetic over Q,
Q(sqrt5) or Q(5^{1/4}); mpmath is used only to decide |lambda| vs 1 (with a
guard band) and for the cited Lehmer/Kronecker numerics.
"""
from __future__ import annotations

import sympy as sp
import mpmath as mp

# --- the golden field (Prop 1.2, eq. 2) ------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2            # golden ratio, root of x^2 - x - 1
psi = (1 - sqrt5) / 2            # its conjugate

# --- the seed (Definition 1.1) ---------------------------------------------
R = sp.Matrix([[0, 1], [1, 1]])
I2 = sp.eye(2)

# high-precision decimal working precision for the modulus/charge decisions
_DPS = 60
_UNIT_TOL = mp.mpf("1e-9")       # guard band: every object sits >= phi or == 1


# --- Fibonacci / Lucas from the bare recurrences ---------------------------
def fib(n: int):
    """F_n with F_0=0, F_1=1, F_{n+1}=F_n+F_{n-1}; F_{-n}=(-1)^{n+1} F_n."""
    n = int(n)
    if n == 0:
        return sp.Integer(0)
    if n > 0:
        a, b = sp.Integer(0), sp.Integer(1)
        for _ in range(n):
            a, b = b, a + b
        return a
    return sp.Integer((-1) ** (n + 1)) * fib(-n)


def luc(n: int):
    """L_n with L_0=2, L_1=1, same recurrence; L_{-n}=(-1)^n L_n."""
    n = int(n)
    if n == 0:
        return sp.Integer(2)
    if n > 0:
        a, b = sp.Integer(2), sp.Integer(1)
        for _ in range(n):
            a, b = b, a + b
        return a
    return sp.Integer((-1) ** n) * luc(-n)


# --- exact-arithmetic conveniences -----------------------------------------
def comm(A, B):
    """Matrix commutator [A,B] = AB - BA."""
    return A * B - B * A


def is_zero(expr) -> bool:
    """True iff a scalar sympy expression is identically zero (exact)."""
    return sp.simplify(sp.expand(expr)) == 0


def mat_eq(A, B) -> bool:
    """True iff matrices A, B are symbolically equal (entrywise)."""
    A, B = sp.Matrix(A), sp.Matrix(B)
    if A.shape != B.shape:
        return False
    return all(sp.simplify(e) == 0 for e in (A - B))


def mat_to_rows(M):
    """Nested list of canonical string forms of a sympy matrix's entries."""
    M = sp.Matrix(M)
    return [[sp.sstr(sp.nsimplify(M[i, j])) for j in range(M.cols)]
            for i in range(M.rows)]


# --- the associative core (Sec. 1) -----------------------------------------
def power_law(n: int):
    """R^n = F_n R + F_{n-1} I = [[F_{n-1},F_n],[F_n,F_{n+1}]]  (Thm 1.3)."""
    Fn, Fn1, Fnp1 = fib(n), fib(n - 1), fib(n + 1)
    return sp.Matrix([[Fn1, Fn], [Fn, Fnp1]])


def cassini(n: int):
    """F_{n-1} F_{n+1} - F_n^2  (= (-1)^n, Cor 1.5c)."""
    return fib(n - 1) * fib(n + 1) - fib(n) ** 2


# --- the Lie layer (Sec. 2) ------------------------------------------------
vphi = sp.Matrix([1, phi])
vpsi = sp.Matrix([1, psi])

H = 2 * R - I2                       # Cartan, = [[-1,2],[2,1]]
Np = vphi * vpsi.T                   # N_+ = v_phi v_psi^T (Thm 2.5)
Nm = vpsi * vphi.T                   # N_- = v_psi v_phi^T
S = sp.Matrix([[2, 1], [1, -2]])     # S = N_+ + N_-  (a second boost)
J = sp.Matrix([[0, -1], [1, 0]])     # the rotation generator

# eigenvector conjugator that diagonalises H (Thm 2.10)
V = sp.Matrix([[1, 1], [phi, psi]])


def adR_matrix_on_M2():
    """4x4 matrix of ad_R = [R, .] on M_2 in the basis E11,E12,E21,E22."""
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]),
             sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    cols = []
    for E in basis:
        img = comm(R, E)
        cols.append([img[0, 0], img[0, 1], img[1, 0], img[1, 1]])
    return sp.Matrix(cols).T


# --- representation theory (Sec. 4) ----------------------------------------
def casimir(m: int):
    """Quadratic Casimir on V_m: (1/2) m (m+2)  (Thm 4.1)."""
    return sp.Rational(1, 2) * m * (m + 2)


def dim_V(m: int) -> int:
    return m + 1


def weights_V(m: int):
    """h-weights of V_m: {m, m-2, ..., -m}  (Thm 4.1)."""
    return [m - 2 * k for k in range(m + 1)]


def sympow_tower(n: int):
    """spec(R|V_n) = {phi^{n-k} psi^k : k=0..n}  (Prop 4.3)."""
    return [phi ** (n - k) * psi ** k for k in range(n + 1)]


def cg_decomp(a: int, b: int):
    """Clebsch-Gordan highest weights j from |a-b| to a+b in steps of 2."""
    return list(range(abs(a - b), a + b + 1, 2))


# --- the trace form and deviation (Sec. 5) ---------------------------------
def trace_form(X, Y):
    """<X,Y> = Tr(XY)  (Def 5.1)."""
    return sp.trace(sp.Matrix(X) * sp.Matrix(Y))


def X(n: int):
    """Deviation operator X_n = 2 R^n - L_n I  (Def 5.3); collapses to F_n H."""
    return 2 * R ** n - luc(n) * I2


def half_trace_Vm(n: int, m: int):
    """(1/2) Tr_{V_m}(X_n^2) from first principles: X_n = F_n H, and on V_m the
    normalised Cartan H/sqrt5 carries weights {m-2k}, so Tr = 5 F_n^2 sum(m-2k)^2."""
    weights_sq = sum((m - 2 * k) ** 2 for k in range(m + 1))
    return sp.Rational(1, 2) * 5 * fib(n) ** 2 * weights_sq


# --- the semiring characters (Definition 3.1) ------------------------------
def _mpc(z):
    """High-precision (dps=_DPS) complex value of a sympy expression, routed
    through decimal strings so no double-precision truncation occurs."""
    mp.mp.dps = _DPS
    v = sp.N(z, _DPS)
    return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))


def modulus(z):
    """|z| as a high-precision mpmath real."""
    return abs(_mpc(z))


def outside_unit(vals):
    """Sub-list of vals with |v| > 1 (guard band _UNIT_TOL)."""
    return [v for v in vals if modulus(v) > 1 + _UNIT_TOL]


def mahler_numeric(vals):
    """M(A) = prod max(1,|lambda|) as an mpmath real  (Definition 3.1)."""
    mp.mp.dps = _DPS
    prod = mp.mpf(1)
    for v in vals:
        m = modulus(v)
        if m > 1 + _UNIT_TOL:
            prod *= m
    return prod


def mahler_exact(vals):
    """Exact symbolic Mahler measure: the guard band selects the |v|>1 values,
    then their moduli are multiplied SYMBOLICALLY (algebraic values)."""
    out = outside_unit(vals)
    if not out:
        return sp.Integer(1)
    return sp.nsimplify(sp.simplify(sp.prod([sp.Abs(v) for v in out])))


def charge_one(z) -> int:
    """chi of a single number: round(2 arg(z)/pi) mod 4  (Definition 3.1)."""
    mp.mp.dps = _DPS
    a = mp.arg(_mpc(z))
    return int(mp.nint(2 * a / mp.pi)) % 4


def charge_set(vals):
    """chi of a multiset: the SET of per-element charges."""
    return sorted({charge_one(v) for v in vals})


def roots_numeric(coeffs, dps: int = 50):
    """Numeric roots (mpmath.mpc) of an integer/real polynomial given by its
    coefficient list, leading coefficient first.  Fast for higher degree,
    where symbolic all_roots() of a CRootOf becomes prohibitively slow."""
    mp.mp.dps = dps
    return mp.polyroots(coeffs, maxsteps=500, extraprec=200)


def mahler_of_roots(roots):
    """M(A) = prod max(1, |r|) over already-numeric roots (mpmath)."""
    mp.mp.dps = _DPS
    prod = mp.mpf(1)
    for r in roots:
        if abs(r) > 1 + _UNIT_TOL:
            prod *= abs(r)
    return prod


def charge_of_roots(roots):
    """chi = {round(2 arg(r)/pi) mod 4} over already-numeric roots (mpmath)."""
    return sorted({int(mp.nint(2 * mp.arg(r) / mp.pi)) % 4 for r in roots})


def lehmer_mahler(dps: int = 50):
    """Mahler measure of Lehmer's degree-10 polynomial
    x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1, from its roots at high precision."""
    mp.mp.dps = dps
    coeffs = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    roots = mp.polyroots(coeffs, maxsteps=500, extraprec=300)
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    return M
