"""
cmc_core -- shared computational primitives for the producer scripts of

    "The Charge-Measure Coupling on a Spectral Semiring: Cyclicity, Composition,
     and a Parity-Graded Mahler Floor" (v4, revised)
    papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex

This module is a clean, self-contained implementation of the paper's two
characters and the semiring operations acting on them.  It is imported by the
producer scripts in this folder; it is NOT a producer itself (it emits no data
and is not meant to be run directly).  It is written independently of the
verification suite under tests/ (which it does not import).

    * Character I  (magnitude) : the Mahler measure   M(p) = |lead| * prod_{|a|>1} |a|
    * Character II (phase)      : the charge group      Z/n  (least n with a^n in R_{>0}
                                  for every root a), or None when a conjugate angle
                                  is irrational (charge-inadmissible).

All root finding is done in exact-ish high precision with mpmath so that the
"is this argument exactly rational" decision is never adjudicated by machine
doubles (Section 7, "Exact arithmetic").
"""

from __future__ import annotations

from math import gcd
from functools import reduce

import mpmath as mp

# Working precision.  30-35 digits is the paper's stated regime; we use 50 so
# the rational-angle margin (~1e-46) sits far below the membership tolerance.
DPS = 50
mp.mp.dps = DPS

# membership tolerance for "n*t is an integer" (Section 7: ~1e-18 at 30 digits;
# at 50 digits the true-zero residual is ~1e-46, so 1e-12 is safely between).
ANGLE_TOL = mp.mpf(10) ** (-12)


# --------------------------------------------------------------------------
# root finding
# --------------------------------------------------------------------------
def roots(coeffs, dps: int = DPS):
    """Roots of a monic-or-not integer/rational polynomial.

    coeffs are highest-degree-first.  Returned as a list of mpmath complex
    numbers at `dps` decimal digits.
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        cs = [c if isinstance(c, (mp.mpf, mp.mpc)) else mp.mpf(str(c)) for c in coeffs]
        r = mp.polyroots(cs, maxsteps=800, extraprec=4 * dps)
    finally:
        mp.mp.dps = old
    return r


# --------------------------------------------------------------------------
# Character I : Mahler measure
# --------------------------------------------------------------------------
def mahler(coeffs, dps: int = DPS):
    """M(p) = |leading coeff| * prod over roots with |root| > 1 of |root|."""
    return mahler_from_roots(roots(coeffs, dps), lead=coeffs[0], dps=dps)


def mahler_from_roots(rts, lead=1, dps: int = DPS):
    """Mahler measure from a precomputed root multiset (used for tensor/Adams
    spectra whose polynomial we never form explicitly)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        m = mp.fabs(mp.mpf(str(lead)))
        for r in rts:
            a = mp.fabs(r)
            if a > 1:
                m *= a
    finally:
        mp.mp.dps = old
    return m


# --------------------------------------------------------------------------
# Character II : charge group and charges
# --------------------------------------------------------------------------
def _t_values(rts, dps: int = DPS):
    """Normalised arguments t = arg(root)/(2*pi) in (-1/2, 1/2] for each root."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        two_pi = 2 * mp.pi
        ts = [mp.arg(r) / two_pi for r in rts]
    finally:
        mp.mp.dps = old
    return ts


def charge_group(coeffs, nmax: int = 240, tol=ANGLE_TOL, dps: int = DPS):
    """Least n >= 1 with a^n in R_{>0} for every root a (= least common
    denominator of the t = arg/(2*pi)).  Returns int, or None if no such n
    <= nmax exists (charge-inadmissible: an irrational conjugate angle)."""
    return charge_group_from_roots(roots(coeffs, dps), nmax=nmax, tol=tol, dps=dps)


def charge_group_from_roots(rts, nmax: int = 240, tol=ANGLE_TOL, dps: int = DPS):
    ts = _t_values(rts, dps)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        for n in range(1, nmax + 1):
            if all(mp.fabs(n * t - mp.nint(n * t)) <= tol for t in ts):
                return n
        return None
    finally:
        mp.mp.dps = old


def charges(coeffs, n=None, dps: int = DPS):
    """Sorted multiset of charges chi_n(a) = round(n*t) mod n over the roots.
    If n is None the charge-group order is used."""
    rts = roots(coeffs, dps)
    if n is None:
        n = charge_group_from_roots(rts, dps=dps)
        if n is None:
            raise ValueError("charge-inadmissible: no finite charge group")
    ts = _t_values(rts, dps)
    return sorted(int(mp.nint(n * t)) % n for t in ts)


def charges_from_roots(rts, n, dps: int = DPS):
    ts = _t_values(rts, dps)
    return sorted(int(mp.nint(n * t)) % n for t in ts)


def attains_all_charges(coeffs, n, dps: int = DPS) -> bool:
    """True iff every residue 0..n-1 occurs among the charges."""
    return set(charges(coeffs, n, dps)) == set(range(n))


# --------------------------------------------------------------------------
# reciprocity
# --------------------------------------------------------------------------
def is_reciprocal(coeffs) -> bool:
    """p is reciprocal iff x^deg p(1/x) = +-p(x), i.e. its coefficient list is
    palindromic up to an overall sign (monic assumed)."""
    c = list(coeffs)
    rev = c[::-1]
    return c == rev or c == [-a for a in rev]


# --------------------------------------------------------------------------
# semiring operations on spectra
# --------------------------------------------------------------------------
def tensor_roots(cA, cB, dps: int = DPS):
    """Roots of A (x) B : the deg(A)*deg(B) pairwise products a_i * b_j."""
    rA = roots(cA, dps)
    rB = roots(cB, dps)
    return [a * b for a in rA for b in rB]


def adams_roots(coeffs, k: int, dps: int = DPS):
    """Roots of psi^k(A) : each eigenvalue raised to the k-th power."""
    return [r ** k for r in roots(coeffs, dps)]


def monic_charpoly_int(rts, dps: int = DPS):
    """Reconstruct integer coefficients (highest-first) of the monic polynomial
    with the given root multiset, rounding to the nearest integer.  Used to
    recover e.g. x^12 - 128 from a tensor spectrum."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        coeffs = [mp.mpc(1)]
        for r in rts:
            new = [mp.mpc(0)] * (len(coeffs) + 1)
            for i, c in enumerate(coeffs):
                new[i] += c
                new[i + 1] -= c * r
            coeffs = new
        out = []
        for c in coeffs:
            assert mp.fabs(mp.im(c)) < mp.mpf(10) ** (-15), "non-real coefficient"
            out.append(int(mp.nint(mp.re(c))))
    finally:
        mp.mp.dps = old
    return out


# --------------------------------------------------------------------------
# exact finite-group helpers (Z/L arithmetic; no floating point)
# --------------------------------------------------------------------------
def lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def subgroup_closure(gens, L: int):
    """Additive closure of {gens} inside Z/L (always contains 0)."""
    S = {0}
    frontier = {g % L for g in gens} | {0}
    while frontier:
        added = set()
        for s in S:
            for g in frontier:
                v = (s + g) % L
                if v not in S:
                    added.add(v)
        if not added:
            break
        S |= added
        frontier = added
    return S


def cyclic_generator_mod(gens, L: int) -> int:
    """The gcd generating <gens> as a subgroup of Z/L (the subgroup is cyclic,
    equal to <g> with g = gcd(gens..., L))."""
    return reduce(gcd, [g % L for g in gens] + [L])
