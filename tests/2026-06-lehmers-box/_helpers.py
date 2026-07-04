"""Shared independent re-derivation helpers for the Lehmer's Box tests.

Nothing here restates a paper value; these are generic tools (Mahler measure by
the product formula, Salem classification, exact Q(sqrt5) sign) used by the
tests to *build up* results from first principles before comparing to the paper.
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 50


def mahler_product(int_coeffs, dps=60):
    """Mahler measure of a monic integer polynomial via the PRODUCT formula
    M(p) = |a_d| * prod max(1,|root|).  int_coeffs highest-degree first."""
    with mp.workdps(dps):
        lead = abs(mp.mpf(int(int_coeffs[0])))
        roots = mp.polyroots([mp.mpf(int(c)) for c in int_coeffs],
                             maxsteps=400, extraprec=400)
        m = lead
        for r in roots:
            a = abs(r)
            if a > 1:
                m *= a
        return m


def root_classification(int_coeffs, tol=mp.mpf(10) ** -18, dps=80):
    """Return (n_outside, n_inside, n_oncircle) for roots of an integer poly."""
    with mp.workdps(dps):
        roots = mp.polyroots([mp.mpf(int(c)) for c in int_coeffs],
                             maxsteps=600, extraprec=600)
        out = ins = onc = 0
        for r in roots:
            a = abs(r)
            if a > 1 + tol:
                out += 1
            elif a < 1 - tol:
                ins += 1
            else:
                onc += 1
        return out, ins, onc


def is_palindromic(int_coeffs):
    """Reciprocal test for a monic integer polynomial: coefficient palindrome."""
    c = list(int_coeffs)
    return c == c[::-1]


def is_salem(int_coeffs):
    """True iff int_coeffs is the (irreducible, reciprocal) minimal polynomial
    of a Salem number: even degree >= 4, reciprocal, exactly one real root > 1,
    one real root in (0,1), the rest a positive number of on-circle pairs."""
    d = len(int_coeffs) - 1
    if d < 4 or d % 2 != 0:
        return False
    if not is_palindromic(int_coeffs):
        return False
    out, ins, onc = root_classification(int_coeffs)
    # one real reciprocal pair off circle, the remaining d-2 on the circle
    return out == 1 and ins == 1 and onc == d - 2 and onc >= 2


# ---- exact arithmetic in Q(sqrt5) ------------------------------------------
_S5 = sp.sqrt(5)


def to_ab_sqrt5(expr):
    """Write an element of Q(sqrt5) as (a, b) with expr = a + b*sqrt5, using the
    field automorphism sqrt5 -> -sqrt5 (exact, no floating point)."""
    expr = sp.expand(expr)
    conj = expr.subs(_S5, -_S5)
    a = sp.nsimplify(sp.simplify((expr + conj) / 2))
    b = sp.nsimplify(sp.simplify((expr - conj) / (2 * _S5)))
    return sp.Rational(a), sp.Rational(b)


def sign_ab_sqrt5(a, b):
    """Exact sign of a + b*sqrt5 (a,b rational) by the paper's Prop. 7.5 rule:
    if a,b share a sign the sign is theirs; otherwise compare a^2 with 5 b^2."""
    a = sp.Rational(a)
    b = sp.Rational(b)
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    # opposite signs: decide by |a| vs |b|sqrt5  <=>  a^2 vs 5 b^2
    mag = 1 if a * a > 5 * b * b else (-1 if a * a < 5 * b * b else 0)
    return mag if a > 0 else -mag


PHI = (1 + _S5) / 2  # exact symbolic golden ratio in Q(sqrt5)
