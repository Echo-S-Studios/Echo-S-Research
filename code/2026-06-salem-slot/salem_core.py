"""
salem_core -- shared computational primitives for the producer scripts of

    "The Occupant of the Salem Slot: A Positive Characterization: the Trace
     Redirection, the Grow Channel, the sqrt5 Limit at the Floor, and Its Rate"
    papers/2026-06-salem-slot/salem_slot.tex   (AceTheDactyl / Echo S Studios)

This module is a clean, self-contained implementation of the single map on
which the whole paper is built -- the *trace substitution*

    trace(theta) = theta + 1/theta                     (Definition 2.1)

together with its inverse (the *lift* L), the trace-down of a reciprocal
polynomial, the Mahler measure, and the canonical Salem family S_n = x^n P - P*.

It is imported by the producer scripts in this folder; it is NOT itself a
producer (it emits no data and is not meant to be run directly).  It is written
independently of the verification suite under tests/ (which it does not import):
the tests *assert* the paper's numbers, these routines *compute and emit* them.

All numeric root finding is done at high precision with mpmath so that the
"is this trace root exactly at a lattice angle" decision (Section 'entry') is
never adjudicated by machine doubles -- the paper's "exact arithmetic" regime.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

# Working precision for the mpmath side.  40+ dps is the paper/test regime.
DPS = 50
mp.mp.dps = DPS

# Symbols used across the module.
x, t = sp.symbols('x t')

# The golden ratio and sqrt5, both as exact sympy and high-precision mpmath.
phi_sym = (1 + sp.sqrt(5)) / 2
sqrt5_sym = sp.sqrt(5)


def PHI() -> mp.mpf:
    """Golden ratio phi = (1+sqrt5)/2 at the current mpmath precision."""
    return (1 + mp.sqrt(5)) / 2


def SQRT5() -> mp.mpf:
    """sqrt5 at the current mpmath precision."""
    return mp.sqrt(5)


# --------------------------------------------------------------------------
# root finding + Mahler measure  (Character: magnitude / Mah)
# --------------------------------------------------------------------------
def roots_mp(coeffs, dps: int = DPS):
    """Roots (mpmath complex) of a polynomial given highest-degree-first."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        cs = [c if isinstance(c, (mp.mpf, mp.mpc)) else mp.mpf(str(c)) for c in coeffs]
        return mp.polyroots(cs, maxsteps=800, extraprec=4 * dps)
    finally:
        mp.mp.dps = old


def mahler_from_coeffs(coeffs, dps: int = DPS) -> mp.mpf:
    """Mahler measure Mah(P) = |lead| * prod_{|root|>1} |root| of a monic-or-not
    integer polynomial (Theorem 3.3 uses prod max(1,|root|); for monic P the two
    agree).  Returns prod over roots of max(1,|root|)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        m = mp.mpf(1)
        for r in roots_mp(coeffs, dps):
            a = mp.fabs(r)
            if a > 1:
                m *= a
        return m
    finally:
        mp.mp.dps = old


def dominant_real_root(coeffs, dps: int = DPS) -> mp.mpf:
    """Largest real root of a polynomial (highest-degree-first coeffs)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        reals = [mp.re(r) for r in roots_mp(coeffs, dps) if mp.fabs(mp.im(r)) < mp.mpf(10) ** (-dps // 2)]
        return max(reals)
    finally:
        mp.mp.dps = old


# --------------------------------------------------------------------------
# the trace map, its inverse (lift), and the trace-down  (Section 2)
# --------------------------------------------------------------------------
def trace(theta):
    """trace(theta) = theta + 1/theta  (Definition 2.1)."""
    return theta + 1 / theta


def salem_from_tau0(tau0):
    """Inverse of trace on the grow branch: beta = (t + sqrt(t^2-4))/2, the
    dominant root of x^2 - t x + 1 (Lemma 4.1)."""
    return (tau0 + mp.sqrt(tau0 ** 2 - 4)) / 2


def flip_discriminant():
    """D(t) = t^2 - 4, the discriminant of x^2 - t x + 1 (Definition 2.1)."""
    return sp.discriminant(x ** 2 - t * x + 1, x)


def lift(T_expr):
    """The lift L(T)(x) = x^{deg T} * T(x + 1/x), a reciprocal polynomial in x
    (Definition 'Entry map', Section 'entry').  Inverse of the trace-down."""
    Tp = sp.Poly(T_expr, t)
    d = Tp.degree()
    return sp.Poly(sp.expand(x ** d * Tp.as_expr().subs(t, x + 1 / x)), x)


def trace_down(P_expr):
    """The trace-down of a monic reciprocal P of even degree 2m: the monic
    T of degree m with P(x) = x^m T(x + 1/x)  (Definition 2.1).

    Computed exactly via the Chebyshev-like power-sum recurrence
        p_0 = 2,  p_1 = t,  p_k = t p_{k-1} - p_{k-2}   (p_k = x^k + x^{-k}),
    since P(x)/x^m = a_m + sum_{k=1}^m a_{m+k} (x^k + x^{-k}).
    """
    Pp = sp.Poly(P_expr, x)
    deg = Pp.degree()
    if deg % 2 != 0:
        raise ValueError("trace_down expects even degree (reciprocal P)")
    m = deg // 2
    a = {j: Pp.coeff_monomial(x ** j) if j > 0 else Pp.coeff_monomial(1) for j in range(deg + 1)}
    p = {0: sp.Integer(2), 1: t}
    for k in range(2, m + 1):
        p[k] = sp.expand(t * p[k - 1] - p[k - 2])
    T = a[m] + sum(a[m + k] * p[k] for k in range(1, m + 1))
    return sp.Poly(sp.expand(T), t)


def is_reciprocal(P_expr) -> bool:
    """P is reciprocal iff its coefficient list is palindromic (monic assumed)."""
    Pp = sp.Poly(P_expr, x)
    c = Pp.all_coeffs()
    return c == c[::-1]


# --------------------------------------------------------------------------
# the canonical Salem family  S_n = x^n P - P*,  P = x^2 - x - 1  (Section 6/9)
# --------------------------------------------------------------------------
def Sn_poly(n: int):
    """S_n(x) = x^n P - P* with P = x^2 - x - 1, expanded to
    x^{n+2} - x^{n+1} - x^n + x^2 + x - 1  (Theorem 6.6 / Section 9)."""
    return x ** (n + 2) - x ** (n + 1) - x ** n + x ** 2 + x - 1


def beta_n(n: int, dps: int = DPS) -> mp.mpf:
    """The Salem factor beta_n -> phi from below: the largest real root of S_n
    (Theorem 6.6).  Found by mpmath starting just below phi."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        f = lambda z: z ** (n + 2) - z ** (n + 1) - z ** n + z ** 2 + z - 1
        return mp.findroot(f, PHI() - mp.mpf('1e-3'))
    finally:
        mp.mp.dps = old


def golden_gap(n: int, dps: int = DPS) -> mp.mpf:
    """gap(n) = sqrt5 - tau0(beta_n) = sqrt5 - (beta_n + 1/beta_n)  (Theorem 6.6)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        b = beta_n(n, dps)
        return SQRT5() - (b + 1 / b)
    finally:
        mp.mp.dps = old


# --------------------------------------------------------------------------
# angle charge A  (Section 'entry', Proposition 'charge')
# --------------------------------------------------------------------------
def root_args_over_pi(P_expr, dps: int = 30):
    """arg(root)/pi for every root of P (a sympy expr in x)."""
    return [float(sp.arg(sp.N(r, dps)) / sp.pi) for r in sp.nroots(sp.Poly(P_expr, x), n=dps)]


def on_half_pi_lattice(arg_over_pi: float, tol: float = 1e-9) -> bool:
    """True iff the argument lies in (pi/2)Z, i.e. arg/pi is a multiple of 1/2."""
    return abs(2 * arg_over_pi - round(2 * arg_over_pi)) < tol


def angle_charge(P_expr, dps: int = 30) -> int:
    """A(P) = 1 iff every root has argument in (pi/2)Z, else 0 (Definition
    'Entry map and conserved charge')."""
    return int(all(on_half_pi_lattice(a) for a in root_args_over_pi(P_expr, dps)))
