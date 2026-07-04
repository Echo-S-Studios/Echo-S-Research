"""Independent re-derivation helpers for the Pisot Cross-Shell Residue whitepaper.

Uniquely-named module (imported by the test files via a local sys.path insert) so
it never collides with helpers in the sibling test folders that run concurrently.

Everything here builds the paper's objects FROM SCRATCH (companion matrices,
resultants, high-precision roots) so the tests genuinely re-derive rather than
restate.  Paper: pisot_residue_whitepaper.tex (AceTheDactyl, Echo S Studios).
"""
import sympy as sp
from sympy import symbols, Poly, resultant, cyclotomic_poly, totient
import mpmath as mp

x, y = symbols('x y')


# --------------------------------------------------------------------------
# Ratio object  Rat_p = prim  Res_y( p(y), p(x*y) )      (Prop. 2.1)
# --------------------------------------------------------------------------
def rat_object(coeffs):
    """coeffs = monic integer poly, high->low degree, e.g. [1,a,b,c,d,e].
    Returns the primitive part of Res_y(p(y), p(x*y)) as a sympy Poly in x."""
    p_y = Poly(coeffs, y)
    n = p_y.degree()
    expr = sum(coeffs[k] * (x * y) ** (n - k) for k in range(len(coeffs)))
    p_xy = Poly(sp.expand(expr), y)
    R = resultant(p_y.as_expr(), p_xy.as_expr(), y)
    return Poly(sp.expand(R), x).primitive()[1]


# --------------------------------------------------------------------------
# Cyclotomic "scan": multiset {m: multiplicity} of Phi_m factors of a poly.
# Exact: factor over Z, then match each irreducible factor against Phi_m for
# the finitely many m with phi(m)=deg(factor).   (Lemma 2.2 machinery.)
# --------------------------------------------------------------------------
def cyclotomic_scan(poly_in_x):
    P = Poly(poly_in_x, x)
    out = {}
    for fac, mult in P.factor_list()[1]:
        d = fac.degree()
        for m in range(1, 2 * d * d + 11):           # phi(m)=d  =>  m <= 2 d^2
            if int(totient(m)) == d and Poly(cyclotomic_poly(m, x), x) == fac:
                out[m] = out.get(m, 0) + mult
                break
    return out


def multiplicity_of_factor(poly_in_x, lin):
    """Exact multiplicity of the linear factor `lin` (e.g. x-1) in poly."""
    q = Poly(poly_in_x, x)
    L = Poly(lin, x)
    m = 0
    while q.rem(L).is_zero:
        q = q.quo(L)
        m += 1
    return m


# --------------------------------------------------------------------------
# Companion matrix (last-column form): eigenvalues = roots of the poly.
# --------------------------------------------------------------------------
def companion(coeffs):
    n = len(coeffs) - 1
    C = sp.zeros(n, n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]          # -a0 .. -a_{n-1}
    return C


# --------------------------------------------------------------------------
# High-precision roots (mpmath) with generous working precision.
# --------------------------------------------------------------------------
def hp_roots(coeffs, dps=60):
    with mp.workdps(dps):
        return mp.polyroots([mp.mpf(int(c)) for c in coeffs],
                            maxsteps=600, extraprec=4 * dps)


def is_pisot(coeffs, dps=50):
    """Exactly one root outside the closed unit disk, that root real > 1,
    and NO root on the unit circle.  Returns (is_pisot, n_realstrict_inside,
    n_nonreal_pairs) with the counts taken over the roots other than theta."""
    with mp.workdps(dps):
        tol = mp.mpf(10) ** (-dps // 4)
        rts = hp_roots(coeffs, dps)
        if any(abs(abs(r) - 1) < tol for r in rts):
            return (False, None, None)
        outside = [r for r in rts if abs(r) > 1]
        if len(outside) != 1:
            return (False, None, None)
        th = outside[0]
        if abs(th.imag) > tol or th.real <= 1:
            return (False, None, None)
        inside = [r for r in rts if abs(r) < 1]
        nonreal = [r for r in inside if abs(r.imag) > tol]
        real_in = [r for r in inside if abs(r.imag) <= tol]
        return (True, len(real_in), len(nonreal) // 2)


def dominant_root(coeffs, dps=50):
    with mp.workdps(dps):
        rts = hp_roots(coeffs, dps)
        return max((r.real for r in rts
                    if abs(r.imag) < mp.mpf(10) ** (-dps // 4)), default=None)


# --------------------------------------------------------------------------
# Root-of-unity test at high precision.  z is a root of unity of order <= Mmax
# iff |z|=1 and z**m == 1 for some 1<=m<=Mmax.  A unimodular z at an irrational
# angle returns False (its powers never come near 1).  Used for the C_2 census.
# --------------------------------------------------------------------------
def is_root_of_unity(z, Mmax=1680, dps=80):
    with mp.workdps(dps):
        eps = mp.mpf(10) ** (-dps // 3)
        if abs(abs(z) - 1) > eps:
            return False
        p = mp.mpf(1) * z
        for m in range(1, Mmax + 1):
            p = p * z
            if abs(p - 1) < eps:
                return True
        return False


# Golden ratio
def phi(dps=50):
    with mp.workdps(dps):
        return (1 + mp.sqrt(5)) / 2


# reciprocal / anti-reciprocal test on a coefficient vector
def is_pm_reciprocal(coeffs):
    rev = coeffs[::-1]
    return coeffs == rev or coeffs == [-c for c in rev]
