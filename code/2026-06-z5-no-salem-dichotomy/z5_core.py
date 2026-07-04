"""
z5_core -- shared computational primitives for the producer scripts of

    "The Z/5Z Case of the No-Salem Dichotomy"  (AceTheDactyl, 2026-06-30)
    papers/2026-06-z5-no-salem-dichotomy/Z5-no-salem-dichotomy-whitepaper.tex

Clean, self-contained implementation of the paper's definitions.  Imported by
the producer scripts in this folder; NOT a producer itself (emits no data, not
meant to be run directly).  Written independently of the verification suite in
tests/2026-06-z5-no-salem-dichotomy/ (which it does not import).

Definitions implemented (paper Sec. 2 / Sec. 7 engine):
  * Mahler measure (Def. 2.3):  M(p) = |lead| * prod_{|alpha|>1} |alpha|,
    computed at high precision with mpmath (the paper's dps=45 exact regime).
  * charge group (Def. 2.1): for squarefree monic p with roots
    alpha_j = r_j exp(2 pi i theta_j), charge = Z/nZ, n = lcm_j d_j where
    theta_j = a_j/d_j is the reduced argument fraction; irrational angle =>
    charge-inadmissible.  A fast numpy screen ("floats only screen", Sec. 7)
    isolates the charge-exactly-5 objects; the retained handful then get an
    exact/high-precision Mahler value and an exact Z-factorization.
  * reciprocity: p reciprocal iff its coefficient vector is (anti)palindromic,
    x^deg p(1/x) = +- p(x).

The golden-field constants are built from the bare definition phi=(1+sqrt5)/2,
never transcribed from a decimal in the paper.
"""

from __future__ import annotations

from math import gcd
from functools import reduce
from fractions import Fraction

import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

# ---------------------------------------------------------------------------
# Golden field Q(sqrt5) -- exact symbolic constants (Prop. 3.1).
# ---------------------------------------------------------------------------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2                     # golden ratio, phi^2 - phi - 1 = 0
_x = sp.symbols('x')


def sigma(expr):
    """Galois conjugation of Q(sqrt5): sqrt5 -> -sqrt5 (Thm. 4.1, Prop. 3.1)."""
    return expr.subs(sqrt5, -sqrt5)


# ---------------------------------------------------------------------------
# Character I : Mahler measure (Def. 2.3), high precision.
# ---------------------------------------------------------------------------
def mahler(coeffs, dps: int = 45):
    """M(p) for a monic integer poly (highest-degree-first `coeffs`).

    M(p) = |lead| * prod over roots with |root| > 1 of |root|, via mpmath
    high-precision roots.  Returns an mpf.
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        cs = [mp.mpf(int(c)) for c in coeffs]
        rts = mp.polyroots(cs, maxsteps=1200, extraprec=8 * dps)
        m = mp.fabs(cs[0])
        for r in rts:
            a = mp.fabs(r)
            if a > 1:
                m *= a
    finally:
        mp.mp.dps = old
    return m


def mahler_np(roots) -> float:
    """M from a numpy root array (screening precision only)."""
    return float(np.prod(np.maximum(1.0, np.abs(roots))))


# ---------------------------------------------------------------------------
# Character II : charge group (Def. 2.1).
# ---------------------------------------------------------------------------
def charge_group(coeffs, tol: float = 1e-7, maxden: int = 60):
    """Charge group Z/nZ of a monic integer poly (Def. 2.1).

    Returns:
        int n        -> charge group Z/nZ (n = lcm of the argument denominators)
        None         -> charge-inadmissible (some root has irrational argument)
        'nonsqfree'  -> repeated root (the definition requires squarefree)
    """
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


def is_charge5(coeffs, tol: float = 1e-6):
    """Fast specialization: is the charge group exactly Z/5Z?

    charge = 5  <=>  every root angle is a multiple of 1/5 (each d_j in {1,5})
    AND at least one root has a non-integer angle (so the lcm is 5, not 1).
    Squarefree required.  Returns (True, M_screen) or (False, None), where
    M_screen is the numpy-precision Mahler value (screening only).
    """
    roots = np.roots(coeffs)
    if len(roots) > 1:
        d = np.abs(roots[:, None] - roots[None, :])
        np.fill_diagonal(d, np.inf)
        if d.min() < 1e-6:
            return (False, None)
    theta = (np.angle(roots) / (2 * np.pi)) % 1.0
    t5 = (5 * theta) % 1.0
    if not np.all(np.minimum(t5, 1 - t5) < tol):       # some angle not a fifth
        return (False, None)
    if np.all(np.minimum(theta, 1 - theta) < tol):     # all positive real -> Z/1Z
        return (False, None)
    return (True, mahler_np(roots))


# ---------------------------------------------------------------------------
# reciprocity (Sec. 2).
# ---------------------------------------------------------------------------
def is_reciprocal(coeffs) -> bool:
    """p reciprocal iff its coefficient vector is palindromic or anti-palindromic."""
    c = [int(v) for v in coeffs]
    return c == c[::-1] or c == [-v for v in c[::-1]]


# ---------------------------------------------------------------------------
# exact Z-factorization + human-readable factor string.
# ---------------------------------------------------------------------------
def _poly_to_str(coeffs) -> str:
    """A degree-first integer coeff list -> compact polynomial string."""
    poly = sp.Poly([int(c) for c in coeffs], _x)
    return str(poly.as_expr())


def factor_over_Z(coeffs):
    """Factor a monic integer poly over Z.

    Returns (factor_string, [(factor_coeffs, multiplicity), ...]) with the
    cyclotomic factor Phi_5 named explicitly when present.
    """
    P = sp.Poly([int(c) for c in coeffs], _x)
    const, facs = P.factor_list()
    parts = []
    factors = []
    phi5 = sp.Poly([1, 1, 1, 1, 1], _x)
    for f, mult in facs:
        fc = [int(c) for c in f.all_coeffs()]
        factors.append((fc, int(mult)))
        if f == phi5:
            name = "Phi_5"
        else:
            name = "(" + str(f.as_expr()) + ")"
        parts.append(name if mult == 1 else f"{name}^{mult}")
    prefix = "" if const == 1 else f"{const}*"
    return prefix + "*".join(parts) if parts else str(const), factors


# ---------------------------------------------------------------------------
# closed-form recognition for the Mahler measures that occur in this paper.
# ---------------------------------------------------------------------------
def _known_measures():
    """Named algebraic values that arise as charge-5 Mahler measures."""
    tbl = {
        "1": mp.mpf(1),
        "phi^2": mp.mpf(str(sp.N(phi**2, 45))),          # (3+sqrt5)/2 ~ 2.618
        "2+sqrt3": mp.mpf(2) + mp.sqrt(3),               # ~ 3.732
        "phi^4": mp.mpf(str(sp.N(phi**4, 45))),          # (7+3sqrt5)/2 ~ 6.854
    }
    return tbl


def recognize_measure(M, tol=mp.mpf('1e-9')) -> str:
    """Map a high-precision Mahler value to a closed form when recognizable
    ({1, phi^2, 2+sqrt3, phi^4} or a small positive integer); else 'numeric'."""
    M = mp.mpf(M)
    for name, val in _known_measures().items():
        if mp.almosteq(M, val, abs_eps=tol):
            return name
    n = int(mp.nint(M))
    if n >= 1 and mp.almosteq(M, mp.mpf(n), abs_eps=tol):
        return str(n)
    return "numeric"


# ---------------------------------------------------------------------------
# The paper's verification window (Sec. 7 / Prop. 6.1):
#     quartics |c| <= 10, quintics |c| <= 4, sextics |c| <= 3.
# Enumerate every monic integer poly in the box, keep the charge-exactly-5
# squarefree ones, and record the object with an EXACT (high-precision) Mahler
# measure, its closed form, reciprocity, and its Z-factorization.
# ---------------------------------------------------------------------------
WINDOW = [(4, 10), (5, 4), (6, 3)]


def _iter_window():
    import itertools
    for deg, b in WINDOW:
        for tail in itertools.product(range(-b, b + 1), repeat=deg):
            if tail[-1] == 0:                  # x | p -> root at 0, angle undefined
                continue
            yield deg, b, (1,) + tail


def window_objects(dps: int = 45):
    """Enumerate the Sec. 7 window and return the charge-exactly-5 objects.

    Each item is a dict:
        degree, coeffs (tuple), poly (str), reciprocal (bool),
        mahler (mpf, high precision), mahler_str (40-digit str),
        measure_closed_form (str), factorization (str), window_bound (int).
    """
    found = []
    for deg, b, coeffs in _iter_window():
        ok, _ = is_charge5(coeffs)
        if not ok:
            continue
        M = mahler(coeffs, dps=dps)
        fac_str, _ = factor_over_Z(coeffs)
        found.append({
            "degree": deg,
            "coeffs": coeffs,
            "poly": _poly_to_str(coeffs),
            "reciprocal": is_reciprocal(coeffs),
            "mahler": M,
            "mahler_str": mp.nstr(M, 40),
            "measure_closed_form": recognize_measure(M),
            "factorization": fac_str,
            "window_bound": b,
        })
    return found
