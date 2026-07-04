"""
helix_core -- shared exact-arithmetic kernel for the producer scripts of
"The Dissolved Helix and Its Orthogonal Partner".

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex

This module is imported by the section producers; it is not run directly. It
provides the golden data (R, phi, psi, sqrt5; Def. 2.1), the two characters of
Def. 2.4 -- Character I (Mahler measure M) and Character II (the Z/4Z phase
charge chi) -- and the ad-operator builder used both for spec(ad_R) (Prop. 2.3)
and for the D<0 rotation channel (Thm. 4.2). All computation is exact over
Q, Q(sqrt5), or Q(5^{1/4}); mpmath is used only for decimal display and for
independent high-precision Mahler cross-checks (no float crosses a decision).
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
DIGITS = 40

x, lam = sp.symbols("x lambda")
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2                 # golden ratio, root of x^2 - x - 1
psi = (1 - sqrt5) / 2                 # conjugate = 1 - phi = -1/phi
R = sp.Matrix([[0, 1], [1, 1]])       # Fibonacci companion (Def. 2.1)
I2 = sp.eye(2)


def dec(expr, digits: int = DIGITS) -> str:
    """Exact sympy expression -> plain decimal string (display only, mpmath)."""
    return mp.nstr(mp.mpf(str(sp.N(expr, digits + 12))), digits)


def chi_of(z) -> int:
    """Character II (Def. 2.4): chi(z) = round(2*arg(z)/pi) mod 4 in Z/4Z.

    Implemented as floor(2*arg/pi + 1/2) (round-half-up) so the four cardinal
    directions map 1 -> 0, i -> 1, -1 -> 2, -i -> 3.
    """
    a = sp.arg(z)
    return int(sp.floor(2 * a / sp.pi + sp.Rational(1, 2)) % 4)


def mahler_of_multiset(roots) -> sp.Expr:
    """Character I (Def. 2.4): M(A) = prod max(1, |lambda|) over a multiset (exact)."""
    m = sp.Integer(1)
    for r in roots:
        m *= sp.Max(1, sp.Abs(r))
    return sp.simplify(m)


def ad_matrix(M: sp.Matrix) -> sp.Matrix:
    """ad_M = [M, .] as a 4x4 matrix on M_2 ~ C^4 (row-major vec of the 2x2 basis)."""
    basis = [
        sp.Matrix([[1, 0], [0, 0]]),
        sp.Matrix([[0, 1], [0, 0]]),
        sp.Matrix([[0, 0], [1, 0]]),
        sp.Matrix([[0, 0], [0, 1]]),
    ]
    vec = lambda A: [A[0, 0], A[0, 1], A[1, 0], A[1, 1]]
    cols = [sp.Matrix(4, 1, vec(M * E - E * M)) for E in basis]
    return sp.Matrix.hstack(*cols)


def eig_multiset(matrix: sp.Matrix):
    """Return {simplified eigenvalue: multiplicity} for an exact matrix."""
    return {sp.simplify(sp.radsimp(k)): v for k, v in matrix.eigenvals().items()}


def mahler_numeric(poly, var=x) -> mp.mpf:
    """Independent high-precision Mahler measure |lead|*prod max(1,|root|)."""
    P = sp.Poly(poly, var)
    coeffs = [mp.mpf(str(c)) for c in P.all_coeffs()]
    rts = mp.polyroots(coeffs, maxsteps=400, extraprec=300)
    m = abs(coeffs[0])
    for r in rts:
        if abs(r) > 1:
            m *= abs(r)
    return m


def is_reciprocal(poly, var=x) -> bool:
    """True iff the coefficient vector is palindromic (up to an overall sign)."""
    c = sp.Poly(poly, var).all_coeffs()
    rc = c[::-1]
    return c == rc or c == [-v for v in rc]


def is_even_poly(poly, var=x) -> bool:
    """True iff poly(-x) == poly(x) (a polynomial in x^2)."""
    return sp.expand(poly.subs(var, -var) - poly) == 0
