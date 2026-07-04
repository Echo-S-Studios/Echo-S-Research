r"""Independent re-derivation utilities for the relational-charge paper.

Everything here is written from the *definitions* in the paper, not copied from
its Appendix-B probe.  The two engines used are:

  * sympy  -- exact integer/rational polynomial arithmetic (resultants,
    factorisation, Sturm root counts, cyclotomic polynomials).
  * mpmath -- high-precision numerics, used ONLY as an independent cross-check
    of the exact results (never as the primary decision).

The paper's decision procedure trial-divides a ratio polynomial by every
cyclotomic Phi_M up to a totient bound.  We instead *factor* the ratio
polynomial over Q and identify which irreducible factors are cyclotomic --- a
genuinely different algorithm that is complete by construction (factorisation
finds every irreducible factor), so it also independently corroborates the
completeness claim of Lemma 4.5 (lem:complete).
"""
import math
from fractions import Fraction

import mpmath as mp
import sympy as sp

mp.mp.dps = 50

x, y, z = sp.symbols("x y z")


# --------------------------------------------------------------------------
# ratio objects (Definition 6.3 def:rat, and the mixed object of Sec. 8)
# --------------------------------------------------------------------------
def ratio_poly(p_expr):
    r"""Rat_p = primitive part of Res_y(p(y), p(xy)); root multiset = all
    ordered ratios alpha_j/alpha_i (Definition 6.3)."""
    P = sp.expand(p_expr)
    R = sp.resultant(P.subs(x, y), sp.expand(P.subs(x, x * y)), y)
    return sp.Poly(sp.expand(R), x).primitive()[1]


def mixed_ratio_poly(p_expr, q_expr):
    r"""Rat_{p,q} = primitive part of Res_y(q(y), p(xy)); root multiset =
    {alpha/beta : p(alpha)=0, q(beta)=0} (Definition 8.1)."""
    P, Q = sp.expand(p_expr), sp.expand(q_expr)
    R = sp.resultant(Q.subs(x, y), sp.expand(P.subs(x, x * y)), y)
    return sp.Poly(sp.expand(R), x).primitive()[1]


def companion_matrix(p_expr):
    """Companion matrix of a monic integer polynomial (rational entries)."""
    p = sp.Poly(p_expr, x)
    n = p.degree()
    c = p.all_coeffs()  # leading first
    lead = c[0]
    C = sp.zeros(n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        # last column holds -a_i / a_n where a_i is coeff of x^i
        C[i, n - 1] = -sp.Rational(c[n - i], lead)
    return C


def ratio_poly_via_kronecker(p_expr):
    r"""Independent construction (ledger G): for p with |p(0)|=1 the matrix
    C_p (x) C_p^{-1} has characteristic polynomial with the same root multiset
    as Rat_p (the ordered ratios).  Uses no resultant."""
    C = companion_matrix(p_expr)
    K = sp.Matrix(sp.kronecker_product(C, C.inv()))
    F = sp.Poly(K.charpoly(x).as_expr(), x)
    return F.primitive()[1]


# --------------------------------------------------------------------------
# cyclotomic contact signatures
# --------------------------------------------------------------------------
def _identify_cyclotomic(fac_poly):
    """Return M if the (monic, irreducible) factor equals Phi_M, else None.

    Independent of the paper's scan: we already have an irreducible factor and
    only need its *name*.  A factor of degree e can be Phi_M only for M with
    totient(M)=e, and any such M satisfies M <= 2 e**2 (the totient bound we
    separately verify in test_completeness_bound.py); we compare the actual
    polynomials, so identification is exact and unambiguous.
    """
    e = fac_poly.degree()
    fm = fac_poly.monic().as_expr()
    for M in range(1, 2 * e * e + 3):
        if sp.totient(M) == e and sp.expand(sp.cyclotomic_poly(M, x) - fm) == 0:
            return M
    return None


def contact_signature(P):
    """Exact cyclotomic-contact signature {M: multiplicity} of a polynomial,
    obtained by factoring over Q and naming the cyclotomic factors."""
    P = sp.Poly(P, x)
    sig = {}
    for fac, mult in P.factor_list()[1]:
        M = _identify_cyclotomic(fac)
        if M is not None:
            sig[M] = sig.get(M, 0) + mult
    return dict(sorted(sig.items()))


def phi1_multiplicity(P):
    """Multiplicity of Phi_1 = (x-1) in P (exact)."""
    P = sp.Poly(P, x)
    m = 0
    xm1 = sp.Poly(x - 1, x)
    while True:
        q, r = sp.div(P, xm1)
        if r.is_zero:
            m += 1
            P = q
        else:
            return m


# --------------------------------------------------------------------------
# numeric cross-check of a contact signature (independent engine: mpmath)
# --------------------------------------------------------------------------
def roots_mp(p_expr):
    """High-precision complex roots of an integer polynomial."""
    p = sp.Poly(p_expr, x)
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    return mp.polyroots(coeffs, maxsteps=400, extraprec=300)


def contact_signature_numeric(p_expr, Mmax=None, tol=mp.mpf(10) ** -22):
    """Reproduce the contact signature of Rat_p purely numerically.

    Form every ordered ratio of the roots of p; a ratio r is a root of unity of
    order M iff |r|=1 and r**M=1 for the least such M.  Returns {M: count}.
    This never uses sympy factorisation, so it is engine-independent.
    """
    rts = roots_mp(p_expr)
    n = len(rts)
    if Mmax is None:
        Mmax = 2 * (n * n) * (n * n) + 2  # 2 (deg Rat_p)^2 ; generous
    root_orders = {}  # M -> number of ratios that are primitive M-th roots of 1
    for a in rts:
        for b in rts:
            r = a / b
            if abs(abs(r) - 1) > tol:
                continue  # not unimodular -> never a root of unity
            found = None
            rp = mp.mpc(1, 0)
            for M in range(1, Mmax + 1):
                rp = rp * r
                if abs(rp - 1) < tol:
                    found = M  # least M with r**M = 1 = primitive order
                    break
            if found is not None:
                root_orders[found] = root_orders.get(found, 0) + 1
    # Phi_M contributes phi(M) roots per copy, so its multiplicity is
    # (#primitive M-th roots) / phi(M).
    sig = {}
    for M, cnt in root_orders.items():
        e = int(sp.totient(M))
        assert cnt % e == 0, f"root count {cnt} for order {M} not divisible by phi={e}"
        sig[M] = cnt // e
    return dict(sorted(sig.items()))


# --------------------------------------------------------------------------
# angle coordinate and absolute / relational groups (Sec. 2--3)
# --------------------------------------------------------------------------
def angles_of(p_expr, maxden=5000):
    """t(alpha)=Arg(alpha)/2pi in [0,1) for every root, as exact Fractions
    recovered from high-precision numerics (rational-angle objects only)."""
    out = []
    for r in roots_mp(p_expr):
        a = mp.arg(r) / (2 * mp.pi)
        a = a - mp.floor(a)  # into [0,1)
        fr = Fraction(mp.nstr(a, 40)).limit_denominator(maxden)
        out.append(fr)
    return out


def absolute_lcd(p_expr, maxden=5000):
    """Least common denominator of the angles = absolute charge order n."""
    L = 1
    for fr in angles_of(p_expr, maxden):
        L = L * fr.denominator // math.gcd(L, fr.denominator)
    return L


def relational_order(p_expr, maxden=5000):
    """Order m of the relational group Delta = <t(alpha)-t(beta)> <= Q/Z:
    lcm of denominators of all pairwise angle differences."""
    ang = angles_of(p_expr, maxden)
    L = 1
    for i in range(len(ang)):
        for j in range(len(ang)):
            d = (ang[i] - ang[j]) % 1
            L = L * d.denominator // math.gcd(L, d.denominator)
    return L


# --------------------------------------------------------------------------
# Mahler measure (magnitude character)
# --------------------------------------------------------------------------
def mahler_measure(p_expr):
    """M(p) = |lead| * prod max(1,|root|), high precision."""
    p = sp.Poly(p_expr, x)
    lead = abs(mp.mpf(int(p.LC())))
    prod = lead
    for r in roots_mp(p_expr):
        if abs(r) > 1:
            prod *= abs(r)
    return prod


# --------------------------------------------------------------------------
# trace polynomial and Salem certification (Thm 7.10, ledger O)
# --------------------------------------------------------------------------
def trace_poly(p_expr):
    """For reciprocal p of even degree 2d, the T in Z[z] with
    p(x) = x^d T(x + 1/x), built from x^j + x^-j = P_j(z), P_0=2, P_1=z,
    P_{j+1}=z P_j - P_{j-1}."""
    p = sp.Poly(p_expr, x)
    n = p.degree()
    assert n % 2 == 0
    d = n // 2
    c = p.all_coeffs()[::-1]  # low -> high, c[k] = coeff of x^k
    Pj = [sp.Integer(2), z]
    for j in range(2, d + 1):
        Pj.append(sp.expand(z * Pj[j - 1] - Pj[j - 2]))
    T = sp.Integer(c[d])
    for j in range(1, d + 1):
        T += c[d - j] * Pj[j]
    return sp.Poly(sp.expand(T), z)


def is_reciprocal(p_expr):
    p = sp.Poly(p_expr, x)
    c = p.all_coeffs()
    n = len(c)
    return all(c[i] == c[n - 1 - i] for i in range(n))


def trace_sturm_pattern(p_expr):
    """(a, b, mid) = number of real roots of the trace polynomial strictly in
    (2,inf), strictly in (-inf,-2), and in the open (-2,2), with roots at the
    endpoints +-2 reported separately as at2, atm2."""
    T = trace_poly(p_expr)
    at2 = 1 if T.eval(2) == 0 else 0
    atm2 = 1 if T.eval(-2) == 0 else 0
    a = T.count_roots(2, sp.oo) - at2
    b = T.count_roots(-sp.oo, -2) - atm2
    mid = T.count_roots(-2, 2)
    return a, b, mid, at2, atm2


def is_salem_polynomial(p_expr):
    """Exact Salem certificate: reciprocal + irreducible + trace-Sturm pattern
    (1, 0, d-1) with no trace root at +-2 (Theorem 7.10)."""
    p = sp.Poly(p_expr, x)
    if p.degree() % 2 != 0:
        return False
    if not is_reciprocal(p_expr):
        return False
    if not p.is_irreducible:
        return False
    d = p.degree() // 2
    a, b, mid, at2, atm2 = trace_sturm_pattern(p_expr)
    return a == 1 and b == 0 and mid == d - 1 and at2 == 0 and atm2 == 0


# --------------------------------------------------------------------------
# totient sieve (independent of sympy.totient; used for the completeness bound)
# --------------------------------------------------------------------------
def totient_sieve(N):
    phi = list(range(N + 1))
    for i in range(2, N + 1):
        if phi[i] == i:  # i is prime
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi


# convenient named polynomials -------------------------------------------------
B4 = x**4 - x**3 - x**2 - x + 1                                   # beta_4
LEHMER = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1  # Lehmer
S6 = x**6 - x**4 - x**3 - x**2 + 1
S8 = x**8 - x**5 - x**4 - x**3 + 1
PLASTIC = x**3 - x - 1
