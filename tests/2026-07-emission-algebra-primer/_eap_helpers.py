"""Independent helpers for verifying *The Emission Algebra A* (2026-07 primer).

Nothing here is copied from the paper's final displayed numbers.  These are
standard, from-scratch constructions built only from the paper's premises:

  * the golden companion matrix R = [[0,1],[1,1]]  (Definition 1.1),
  * the golden-field elements phi, psi as (1 +/- sqrt5)/2,
  * Fibonacci / Lucas numbers from their bare recurrences (not from any
    closed form in the paper), valid for negative indices too,
  * the semiring characters (Mahler measure, charge) from Definition 3.1.

Every test re-derives a paper value from these premises and only THEN compares
to the paper's stated result.
"""
import sympy as sp
import mpmath as mp

# --- the golden field ------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2          # the golden ratio
psi = (1 - sqrt5) / 2          # its conjugate

# --- the seed (Definition 1.1) ---------------------------------------------
R = sp.Matrix([[0, 1], [1, 1]])
I2 = sp.eye(2)

_DPS = 60


# --- Fibonacci / Lucas from the bare recurrence ----------------------------
def fib(n):
    """F_n from F_0=0, F_1=1 and F_{n+1}=F_n+F_{n-1}; extended to n<0 by
    the standard reflection F_{-n} = (-1)^{n+1} F_n.  Independent of any
    closed form asserted in the paper."""
    n = int(n)
    if n == 0:
        return sp.Integer(0)
    if n > 0:
        a, b = sp.Integer(0), sp.Integer(1)
        for _ in range(n):
            a, b = b, a + b
        return a
    return sp.Integer((-1) ** (n + 1)) * fib(-n)


def luc(n):
    """L_n from L_0=2, L_1=1 and the same recurrence; L_{-n}=(-1)^n L_n."""
    n = int(n)
    if n == 0:
        return sp.Integer(2)
    if n > 0:
        a, b = sp.Integer(2), sp.Integer(1)
        for _ in range(n):
            a, b = b, a + b
        return a
    return sp.Integer((-1) ** n) * luc(-n)


# --- symbolic matrix / scalar equality -------------------------------------
def mat_eq(A, B):
    """True iff matrices A,B are symbolically equal (entrywise simplify)."""
    A = sp.Matrix(A)
    B = sp.Matrix(B)
    if A.shape != B.shape:
        return False
    return all(sp.simplify(e) == 0 for e in (A - B))


def is_zero(expr):
    """True iff a scalar sympy expression is identically zero."""
    return sp.simplify(sp.expand(expr)) == 0


def comm(A, B):
    """Matrix commutator [A,B] = AB - BA."""
    return A * B - B * A


# --- semiring characters (Definition 3.1) ----------------------------------
# Branch tolerance for the |lambda| vs 1 decision.  Every object in the primer
# has magnitude either exactly 1 (on the unit circle) or >= phi ~ 1.618, so a
# 1e-9 guard band safely separates the two without masking any real value.
_UNIT_TOL = mp.mpf('1e-9')


def _mpc(z):
    """High-precision (dps=_DPS) complex value of a sympy expression, routed
    through decimal strings so no double-precision truncation occurs."""
    mp.mp.dps = _DPS
    v = sp.N(z, _DPS)
    re = mp.mpf(str(sp.re(v)))
    im = mp.mpf(str(sp.im(v)))
    return mp.mpc(re, im)


def modulus(z):
    """|z| as a high-precision mpmath real."""
    return abs(_mpc(z))


def outside_unit(vals):
    """Sub-list of vals with |v| > 1 (branch decided numerically with a guard
    band; all values in the primer sit well clear of the unit circle)."""
    return [v for v in vals if modulus(v) > 1 + _UNIT_TOL]


def mahler(vals):
    """Mahler measure M(A)=prod max(1,|lambda|), as an mpmath real
    (Definition 3.1).  Unit-circle values contribute exactly 1."""
    mp.mp.dps = _DPS
    prod = mp.mpf(1)
    for v in vals:
        m = modulus(v)
        if m > 1 + _UNIT_TOL:
            prod *= m
    return prod


def mahler_exact(vals):
    """Exact symbolic Mahler measure: the numeric guard band decides which
    values lie outside the unit circle, then their moduli are multiplied
    SYMBOLICALLY (for algebraic values in Q(sqrt5), Q(5^{1/4}), ...)."""
    out = outside_unit(vals)
    if not out:
        return sp.Integer(1)
    return sp.nsimplify(sp.simplify(sp.prod([sp.Abs(v) for v in out])))


def charge_one(z):
    """chi of a single number: round(2*arg(z)/pi) mod 4  (Definition 3.1)."""
    mp.mp.dps = _DPS
    a = mp.arg(_mpc(z))              # principal arg in (-pi, pi]
    k = int(mp.nint(2 * a / mp.pi))
    return k % 4


def charge_set(vals):
    """chi of a multiset: the SET of per-element charges (Definition 3.1)."""
    return set(charge_one(v) for v in vals)
