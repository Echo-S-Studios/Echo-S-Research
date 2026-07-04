"""
Independent re-derivation toolkit for the whitepaper

    "The Operator Algebra of the Emission Semiring"
    (papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex)

An *object* (Def. 2.1) is a finite multiset of nonzero algebraic numbers,
realised here as a plain tuple of sympy expressions.  Every operator below is
coded directly from the paper's Definitions -- NOT lifted from any companion
code -- so tests importing this module constitute an independent check.

Conventions
-----------
oplus  :  A (+) B   = multiset union                 (Def. 2.1)
otimes :  A (x) B   = { a*b : a in A, b in B }        (Def. 2.1, Kronecker)
psi    :  psi^n(A)  = { a**n : a in A }               (Def. 3.1, Adams op)
Mahler :  M(A)      = prod_a max(1, |a|)              (Def. 4.1)
charge :  q(a)      = round(2*arg(a)/pi) mod 4        (Def. 5.1)
"""

from collections import Counter

import mpmath as mp
import sympy as sp

x = sp.Symbol("x")
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2  # golden ratio, root of x^2 - x - 1


# --------------------------------------------------------------------------
# seeds (eigenvalue multisets)
# --------------------------------------------------------------------------
def golden_seed():
    """Roots of x^2 - x - 1 : {phi, 1 - phi}.  Mahler measure phi."""
    return (phi, 1 - phi)


def seed_from_poly(poly_expr):
    """Eigenvalue multiset = roots (with multiplicity) of an integer poly."""
    rd = sp.roots(sp.Poly(poly_expr, x))
    out = []
    for r, m in rd.items():
        out.extend([r] * m)
    return tuple(out)


# The Lorentzian seed K = x^4 + 5 x^2 - 5 used throughout the paper.
def K_seed():
    return seed_from_poly(x**4 + 5 * x**2 - 5)


# --------------------------------------------------------------------------
# semiring / power operators
# --------------------------------------------------------------------------
ZERO = ()  # empty multiset  (additive identity, Def. 2.1)
ONE = (sp.Integer(1),)  # {1}  (multiplicative identity, Def. 2.1)


def oplus(A, B):
    """Direct sum = multiset union."""
    return tuple(A) + tuple(B)


def otimes(A, B):
    """Tensor product = elementwise products with multiplicity."""
    return tuple(sp.expand(sp.sympify(a) * sp.sympify(b)) for a in A for b in B)


def psi(n, A):
    """Adams / power operation lambda -> lambda**n."""
    return tuple(sp.expand(sp.sympify(a) ** n) for a in A)


# --------------------------------------------------------------------------
# multiset equality (canonical, exact)
# --------------------------------------------------------------------------
def canon(A):
    return Counter(sp.expand(sp.sympify(a)) for a in A)


def ms_equal(A, B):
    return canon(A) == canon(B)


# --------------------------------------------------------------------------
# symmetric functions of an eigenvalue multiset
# --------------------------------------------------------------------------
def power_sum(A, k):
    """p_k = sum lambda_i**k."""
    return sp.expand(sum(sp.sympify(a) ** k for a in A))


def elem_sym(A, k):
    """e_k = k-th elementary symmetric polynomial of the multiset."""
    return sp.expand(sp.symmetric_poly(k, list(A)))


def sym2(A):
    """Sym^2 A = { lambda_i * lambda_j : i <= j }."""
    L = list(A)
    n = len(L)
    return tuple(sp.expand(L[i] * L[j]) for i in range(n) for j in range(i, n))


def wedge2(A):
    """Lambda^2 A = { lambda_i * lambda_j : i < j }."""
    L = list(A)
    n = len(L)
    return tuple(sp.expand(L[i] * L[j]) for i in range(n) for j in range(i + 1, n))


# --------------------------------------------------------------------------
# Character I : the Mahler measure
# --------------------------------------------------------------------------
def mahler_exact(A):
    """M(A) = prod_a max(1, |a|), simplified to closed form."""
    p = sp.Integer(1)
    for a in A:
        p *= sp.Max(1, sp.Abs(sp.sympify(a)))
    return sp.nsimplify(sp.radsimp(sp.simplify(p)))


def _to_mp_re_im(a, dps):
    a = sp.sympify(a)
    re = mp.mpf(str(sp.N(sp.re(a), dps + 15)))
    im = mp.mpf(str(sp.N(sp.im(a), dps + 15)))
    return re, im


def modulus_mp(a, dps=60):
    mp.mp.dps = dps
    re, im = _to_mp_re_im(a, dps)
    return mp.sqrt(re * re + im * im)


def mahler_mp(A, dps=60):
    mp.mp.dps = dps
    p = mp.mpf(1)
    for a in A:
        m = modulus_mp(a, dps)
        if m > 1:
            p *= m
    return p


def log_mahler_mp(A, dps=60):
    mp.mp.dps = dps
    s = mp.mpf(0)
    for a in A:
        m = modulus_mp(a, dps)
        if m > 1:
            s += mp.log(m)
    return s


# --------------------------------------------------------------------------
# Character II : the angle charge  (Def. 5.1)
# --------------------------------------------------------------------------
def charge_one(a, dps=60):
    """q = round(2 * arg(a) / pi)  mod 4."""
    mp.mp.dps = dps
    re, im = _to_mp_re_im(a, dps)
    theta = mp.atan2(im, re)
    return int(mp.nint(2 * theta / mp.pi)) % 4


def charge(A, dps=60):
    """Multiset of charges of an object."""
    return Counter(charge_one(a, dps) for a in A)


def sumset_mod4(cA, cB):
    """Sumset (mod 4) of two charge multisets: {qa + qb} with multiplicity."""
    out = Counter()
    for qa, na in cA.items():
        for qb, nb in cB.items():
            out[(qa + qb) % 4] += na * nb
    return out
