"""
Independent numerical/algebraic primitives for validating the paper
"The Charge-Measure Coupling on a Spectral Semiring" (v4).

These re-implement, from scratch, the paper's two characters:
  * Character I : the Mahler measure  M(p) = |lead| * prod max(1,|root|)
  * Character II: the charge group    Z/n  (n = least common denominator of
                  the conjugate arguments arg(alpha)/(2*pi)), or None if some
                  conjugate angle is irrational (charge-inadmissible).

Uniquely-named module (`_cmc_helpers`) so concurrent per-paper test folders do
not collide in sys.modules.
"""

import mpmath as mp
import numpy as np

mp.mp.dps = 50


# --------------------------------------------------------------------------
# root finding
# --------------------------------------------------------------------------
def roots_mp(coeffs, dps=50):
    """Roots of a polynomial given by coeffs (highest degree first) at
    `dps` decimal digits.  Coeffs may be ints, floats, sympy numbers, mpf."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        cs = [mp.mpf(str(c)) if not isinstance(c, (mp.mpf, mp.mpc)) else c
              for c in coeffs]
        r = mp.polyroots(cs, maxsteps=500, extraprec=4 * dps)
    finally:
        mp.mp.dps = old
    return r


# --------------------------------------------------------------------------
# Character I : Mahler measure
# --------------------------------------------------------------------------
def mahler(coeffs, dps=50):
    """M(p) = |leading| * prod_{|root|>1} |root|.  Independent of the paper."""
    rs = roots_mp(coeffs, dps)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        m = mp.fabs(mp.mpf(str(coeffs[0])))
        for r in rs:
            a = mp.fabs(r)
            if a > 1:
                m *= a
    finally:
        mp.mp.dps = old
    return m


# --------------------------------------------------------------------------
# Character II : charge group
# --------------------------------------------------------------------------
def _t_values(coeffs, dps=50):
    """For each root, t = arg(root)/(2*pi) in (-1/2, 1/2]."""
    rs = roots_mp(coeffs, dps)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        two_pi = 2 * mp.pi
        ts = [mp.arg(r) / two_pi for r in rs]
    finally:
        mp.mp.dps = old
    return ts


def charge_group(coeffs, Nmax=240, tol=mp.mpf(10) ** (-12), dps=50):
    """
    Least n>=1 with alpha^n in R_{>0} for every root alpha  (== least common
    denominator of the t = arg/(2*pi)).  Returns int n, or None if no such n
    exists up to Nmax (charge-inadmissible: some conjugate angle irrational).
    """
    ts = _t_values(coeffs, dps)
    for n in range(1, Nmax + 1):
        ok = True
        for t in ts:
            nt = n * t
            if mp.fabs(nt - mp.nint(nt)) > tol:
                ok = False
                break
        if ok:
            return n
    return None


def charges(coeffs, n=None, dps=50):
    """Multiset of charges chi_n(alpha)=round(n*t) mod n over the roots.
    If n is None it is taken to be the charge group order."""
    if n is None:
        n = charge_group(coeffs, dps=dps)
        if n is None:
            raise ValueError("charge-inadmissible: no finite charge group")
    ts = _t_values(coeffs, dps)
    out = []
    for t in ts:
        out.append(int(mp.nint(n * t)) % n)
    return sorted(out)


def is_reciprocal(coeffs):
    """p reciprocal iff coefficient list is palindromic (x^deg p(1/x)=p(x)),
    up to an overall sign.  coeffs highest-first, monic assumed."""
    c = list(coeffs)
    rev = c[::-1]
    return c == rev or c == [-a for a in rev]


# --------------------------------------------------------------------------
# tensor of two spectra : polynomial whose roots are the pairwise products
# --------------------------------------------------------------------------
def tensor_charpoly_roots(cA, cB, dps=50):
    """Roots (as mpc) of A (x) B : the products alpha_i * beta_j."""
    rA = roots_mp(cA, dps)
    rB = roots_mp(cB, dps)
    return [a * b for a in rA for b in rB]


def mahler_from_roots(rs, lead=1, dps=50):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        m = mp.fabs(mp.mpf(str(lead)))
        for r in rs:
            a = mp.fabs(r)
            if a > 1:
                m *= a
    finally:
        mp.mp.dps = old
    return m


def charge_group_from_roots(rs, Nmax=240, tol=mp.mpf(10) ** (-12), dps=50):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        two_pi = 2 * mp.pi
        ts = [mp.arg(r) / two_pi for r in rs]
        for n in range(1, Nmax + 1):
            if all(mp.fabs(n * t - mp.nint(n * t)) <= tol for t in ts):
                return n
        return None
    finally:
        mp.mp.dps = old


# --------------------------------------------------------------------------
# fast numpy screen (for finite-window "computed" scans only)
# --------------------------------------------------------------------------
def mahler_np(coeffs):
    """Double-precision Mahler measure via numpy.roots (screening only)."""
    rs = np.roots(coeffs)
    m = abs(coeffs[0])
    for r in rs:
        a = abs(r)
        if a > 1.0:
            m *= a
    return m


def all_real_positive_np(coeffs, tol=1e-9):
    """True if every root is real (|Im| tiny) and positive."""
    rs = np.roots(coeffs)
    for r in rs:
        if abs(r.imag) > tol * max(1.0, abs(r)):
            return False
        if r.real <= tol:
            return False
    return True
