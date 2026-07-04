"""Independent verification engine for

    "The Z/5Z Case of the No-Salem Dichotomy"  (2026-06-30, AceTheDactyl).

NOTHING here is copied from the paper's displayed answers.  Everything is a
from-scratch construction built only from the paper's *definitions*:

  * Mahler measure of a monic integer polynomial (Def. 2.3):
        M(p) = prod_j max(1, |alpha_j|)          (leading coeff = 1)
    computed independently, either by high-precision root finding (mpmath)
    or by numpy screening.
  * charge group (Def. 2.1): for a squarefree monic p in Z[x] with roots
        alpha_j = r_j exp(2 pi i theta_j),
    the charge group is Z/nZ with n = lcm_j d_j, where theta_j = a_j/d_j is the
    reduced fraction of the (rational) argument; if some theta_j is irrational
    the polynomial is charge-inadmissible.
  * reciprocal test: p reciprocal iff its coefficient vector is (anti)palindromic
        x^deg p(1/x) = +- p(x).

The window enumeration reproduces the paper's Sec. 7 verification table using
these definitions only.
"""
import os
import sys
from math import gcd
from functools import reduce
from fractions import Fraction
import itertools

import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

# ---------------------------------------------------------------------------
# Golden field constants -- built from the bare definition phi = (1+sqrt5)/2,
# NOT from any decimal quoted in the paper.
# ---------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2            # golden ratio (Prop. 3.1)


def sigma(expr):
    """Galois conjugation of Q(sqrt5): sqrt5 -> -sqrt5."""
    return expr.subs(sqrt5, -sqrt5)


# ---------------------------------------------------------------------------
# Mahler measure -- high precision, independent of the paper.
# ---------------------------------------------------------------------------
def mahler_hp(coeffs):
    """M(p) for a monic integer poly given highest-degree-first `coeffs`,
    via mpmath high-precision roots.  Returns an mpf."""
    coeffs = [mp.mpf(int(c)) for c in coeffs]
    roots = mp.polyroots(coeffs, maxsteps=800, extraprec=200)
    prod = mp.mpf(1)
    for r in roots:
        a = abs(r)
        if a > 1:
            prod *= a
    return prod


def mahler_np(roots):
    """M from a numpy root array (screening precision)."""
    return float(np.prod(np.maximum(1.0, np.abs(roots))))


# ---------------------------------------------------------------------------
# charge group (Definition 2.1), general version.
# Returns:
#   int n         -> charge group Z/nZ
#   None          -> charge-inadmissible (an irrational argument)
#   'nonsqfree'   -> repeated root (definition requires squarefree)
# ---------------------------------------------------------------------------
def charge_group(coeffs, tol=1e-7, maxden=60):
    roots = np.roots(coeffs)
    if len(roots) > 1:
        d = np.abs(roots[:, None] - roots[None, :])
        np.fill_diagonal(d, np.inf)
        if d.min() < 1e-6:
            return 'nonsqfree'
    dens = []
    for a in roots:
        theta = (np.angle(a) / (2 * np.pi)) % 1.0
        fr = Fraction(theta).limit_denominator(maxden)
        if abs(float(fr) - theta) < tol or abs(float(fr) - theta - 1.0) < tol:
            dens.append(fr.denominator if fr.numerator != 0 else 1)
        else:
            return None
    return reduce(lambda a, b: a * b // gcd(a, b), dens, 1)


def is_charge5(coeffs, tol=1e-6):
    """Fast specialization: True iff charge group is exactly Z/5Z.

    charge = 5  <=>  every root angle is a multiple of 1/5 (so every d_j in
    {1,5})  AND  at least one root has a non-integer angle (so lcm = 5, not 1).
    Squarefree required.  Returns (True, M) or (False, None).
    """
    roots = np.roots(coeffs)
    if len(roots) > 1:
        d = np.abs(roots[:, None] - roots[None, :])
        np.fill_diagonal(d, np.inf)
        if d.min() < 1e-6:
            return (False, None)
    theta = (np.angle(roots) / (2 * np.pi)) % 1.0
    t5 = (5 * theta) % 1.0
    if not np.all(np.minimum(t5, 1 - t5) < tol):      # some angle not a fifth
        return (False, None)
    if np.all(np.minimum(theta, 1 - theta) < tol):    # all positive real -> Z/1Z
        return (False, None)
    return (True, mahler_np(roots))


def is_reciprocal(coeffs):
    """p reciprocal iff coeff vector is palindromic or anti-palindromic."""
    c = list(int(v) for v in coeffs)
    return c == c[::-1] or c == [-v for v in c[::-1]]


# ---------------------------------------------------------------------------
# The paper's verification window (Sec. 7 / Prop. 6.1):
#   quartics |c| <= 10, quintics |c| <= 4, sextics |c| <= 3.
# We enumerate every monic integer poly in that box, keep the charge-exactly-5
# squarefree ones, and record (coeffs, M, reciprocal?).  Cached once.
# ---------------------------------------------------------------------------
_WINDOW = [(4, 10), (5, 4), (6, 3)]
_cache = None


def window_objects():
    global _cache
    if _cache is not None:
        return _cache
    found = []
    for deg, b in _WINDOW:
        for tail in itertools.product(range(-b, b + 1), repeat=deg):
            if tail[-1] == 0:          # x | p  ->  root at 0, angle undefined
                continue
            coeffs = (1,) + tail
            ok, M = is_charge5(coeffs)
            if ok:
                found.append((coeffs, M, is_reciprocal(coeffs)))
    _cache = found
    return found
