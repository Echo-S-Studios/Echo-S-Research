r"""
relcharge_core -- shared computational primitives for the producer scripts of

    "Relational Charge on the Spectral Semiring: Gauge Rigidity, Coherence
     Types, and a Reference-Free Parity Floor" (AceTheDactyl / Echo S Studios)
    papers/2026-07-relational-charge/relational_charge_paper.tex

This module is a clean, self-contained implementation of the paper's decision
procedures, written directly from the definitions in the manuscript.  It is
imported by the producer scripts in this folder; it is NOT a producer itself
(it emits no data and is not run directly).  It does NOT import from tests/.

Method (faithful to the paper).  The contact routine here is the paper's OWN
probe of Appendix B: it trial-divides the ratio polynomial by every cyclotomic
Phi_M up to the totient bound M <= 2 (deg P)^2 of Lemma 4.5, in exact integer
/ rational polynomial arithmetic (no floats in any decision).  This is the
procedure the paper describes and ran.  The independent verification suite under
tests/ deliberately uses a *different* algorithm (factor the ratio polynomial
over Q and name the cyclotomic factors); the two agree by construction.  mpmath
is used here only for the numeric-valued readouts (angles, Mahler measure),
never to decide a torsion / cyclotomic-divisibility question.

Everything in the torsion path -- resultants over Z[x], primitive parts, exact
trial division by cyclotomics, Sturm root counts -- is exact.
"""

from __future__ import annotations

import math
from fractions import Fraction

import mpmath as mp
import sympy as sp

mp.mp.dps = 50

x, y = sp.symbols("x y")


# --------------------------------------------------------------------------
# ratio objects  (Definition 6.3 def:rat; the mixed object of Section 8)
# --------------------------------------------------------------------------
def ratio_poly(p_expr):
    r"""Rat_p = primitive part of Res_y(p(y), p(xy)) in Z[x] (Definition 6.3).

    Because p is monic the resultant has degree exactly n^2 and its root
    multiset is all n^2 ordered ratios alpha_j / alpha_i.
    """
    P = sp.expand(p_expr)
    R = sp.resultant(P.subs(x, y), sp.expand(P.subs(x, x * y)), y)
    return sp.Poly(sp.expand(R), x).primitive()[1]


def mixed_ratio_poly(p_expr, q_expr):
    r"""Rat_{p,q} = primitive part of Res_y(q(y), p(xy)) (Definition 8.1).

    Root multiset = {alpha/beta : p(alpha)=0, q(beta)=0}, degree (deg p)(deg q).
    """
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
        C[i, n - 1] = -sp.Rational(c[n - i], lead)
    return C


def ratio_poly_via_kronecker(p_expr):
    r"""Independent construction (ledger G): for |p(0)|=1 the integer matrix
    C_p (x) C_p^{-1} has characteristic polynomial with the same root multiset
    as Rat_p (all ordered ratios).  Uses no resultant.
    """
    C = companion_matrix(p_expr)
    K = sp.Matrix(sp.kronecker_product(C, C.inv()))
    F = sp.Poly(K.charpoly(x).as_expr(), x)
    return F.primitive()[1]


# --------------------------------------------------------------------------
# totient sieve and the complete cyclotomic-contact scan (paper's Appendix B)
# --------------------------------------------------------------------------
def totients_upto(N):
    """Euler totient for 1..N by linear sieve (independent of sympy.totient)."""
    phi = list(range(N + 1))
    for i in range(2, N + 1):
        if phi[i] == i:  # i is prime
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi


def cyclotomic_contacts(P):
    r"""Complete cyclotomic-contact signature {M: multiplicity} of a polynomial,
    computed by the paper's probe (Appendix B): trial-divide by every Phi_M with
    phi(M) <= deg P, scanning M up to the complete bound 2 (deg P)^2 of
    Lemma 4.5.  Exact rational polynomial arithmetic throughout.
    """
    P = sp.Poly(P, x)
    d = P.degree()
    if d <= 0:
        return {}
    N = 2 * d * d  # complete: Phi_M | P forces M <= 2 d^2
    phi = totients_upto(N)
    hits = {}
    for M in range(1, N + 1):
        if phi[M] > d:
            continue
        cm = sp.Poly(sp.cyclotomic_poly(M, x), x)
        _, r = P.div(cm)
        if r.is_zero:
            mult, T = 0, P
            while True:
                q, r = T.div(cm)
                if r.is_zero:
                    mult, T = mult + 1, q
                else:
                    break
            hits[M] = mult
    return dict(sorted(hits.items()))


def phi1_multiplicity(P):
    """Multiplicity of Phi_1 = (x - 1) in P (exact (x-1)-adic valuation)."""
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


def signature_str(sig):
    r"""Format a contact signature dict as '{Phi_1^4, Phi_2^4}' (or 'empty')."""
    if not sig:
        return "empty"
    parts = [f"Phi_{M}^{mult}" for M, mult in sorted(sig.items())]
    return "{" + ", ".join(parts) + "}"


# --------------------------------------------------------------------------
# high-precision roots and the numeric-valued readouts (Sec. 2-3)
# --------------------------------------------------------------------------
def roots_mp(p_expr):
    """High-precision complex roots of an integer polynomial (mpmath)."""
    p = sp.Poly(p_expr, x)
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    return mp.polyroots(coeffs, maxsteps=400, extraprec=300)


def angles_of(p_expr, maxden=5000):
    """t(alpha) = Arg(alpha)/2pi in [0,1) for every root, as exact Fractions
    (meaningful only for rational-angle / charge-admissible objects)."""
    out = []
    for r in roots_mp(p_expr):
        a = mp.arg(r) / (2 * mp.pi)
        a = a - mp.floor(a)
        out.append(Fraction(mp.nstr(a, 40)).limit_denominator(maxden))
    return out


def absolute_lcd(p_expr, maxden=5000):
    """Absolute charge order n = least common denominator of the angles."""
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


def mahler_measure(p_expr):
    """M(p) = |lead| * prod_{|root|>1} |root|  (high precision)."""
    p = sp.Poly(p_expr, x)
    prod = abs(mp.mpf(int(p.LC())))
    for r in roots_mp(p_expr):
        if abs(r) > 1:
            prod *= abs(r)
    return prod


# --------------------------------------------------------------------------
# trace polynomial and exact Salem certification (Thm 7.10, ledger O)
# --------------------------------------------------------------------------
def is_reciprocal(p_expr):
    p = sp.Poly(p_expr, x)
    c = p.all_coeffs()
    n = len(c)
    return all(c[i] == c[n - 1 - i] for i in range(n))


def trace_poly(p_expr):
    r"""For reciprocal p of even degree 2d, the T in Z[z] with
    p(x) = x^d T(x + 1/x), built from x^j + x^-j = P_j(z): P_0=2, P_1=z,
    P_{j+1} = z P_j - P_{j-1}."""
    z = sp.symbols("z")
    p = sp.Poly(p_expr, x)
    n = p.degree()
    if n % 2 != 0:
        raise ValueError("trace_poly requires even degree")
    d = n // 2
    c = p.all_coeffs()[::-1]  # low -> high, c[k] = coeff of x^k
    Pj = [sp.Integer(2), z]
    for j in range(2, d + 1):
        Pj.append(sp.expand(z * Pj[j - 1] - Pj[j - 2]))
    T = sp.Integer(c[d])
    for j in range(1, d + 1):
        T += c[d - j] * Pj[j]
    return sp.Poly(sp.expand(T), z)


def trace_sturm_pattern(p_expr):
    """(a, b, mid, at2, atm2): count of real roots of the trace polynomial
    strictly in (2, inf), strictly in (-inf, -2), in the open (-2, 2), and
    exactly at +2 / -2.  Exact Sturm counts."""
    z = sp.symbols("z")
    T = trace_poly(p_expr)
    at2 = 1 if T.eval(2) == 0 else 0
    atm2 = 1 if T.eval(-2) == 0 else 0
    a = T.count_roots(2, sp.oo) - at2
    b = T.count_roots(-sp.oo, -2) - atm2
    mid = T.count_roots(-2, 2)
    return a, b, mid, at2, atm2


def is_salem_polynomial(p_expr):
    """Exact Salem certificate (Theorem 7.10): reciprocal + irreducible +
    trace-Sturm pattern (1, 0, d-1) with no trace root at +-2."""
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
# named polynomials referenced throughout the paper
# --------------------------------------------------------------------------
B4 = x**4 - x**3 - x**2 - x + 1                                    # beta_4 (Salem)
LEHMER = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1   # Lehmer (Salem)
S6 = x**6 - x**4 - x**3 - x**2 + 1                                 # Salem sextic
S8 = x**8 - x**5 - x**4 - x**3 + 1                                 # Salem octic
PLASTIC = x**3 - x - 1                                             # plastic / smallest Pisot
Q2 = x**4 + x**2 - 1                                              # q_2 (golden floor)
KSEED = x**4 + 5 * x**2 - 5                                       # catalog seed K
GROUPDROP = x**4 + 5 * x**2 + 5                                   # group-drop quartic
X4MX1 = x**4 - x + 1                                              # fully-rigid quartic
TWISTSHELL = x**4 + x**2 + 2                                      # twisted-shell non-inert witness


def qk(k):
    """q_k = x^{2k} + x^k - 1 (Lemma 5.2): even parity-floor family, M = phi."""
    return x ** (2 * k) + x**k - 1
