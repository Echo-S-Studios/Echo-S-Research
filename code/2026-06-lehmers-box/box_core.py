"""
box_core -- computational core for the Lehmer's Box producer scripts.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex
  "Lehmer's Box: A Golden Floor and an Angle Lattice that Confine Spectral
   Emission Away from Salem Numbers -- Without Resolving Lehmer's Problem".

This module holds the reusable engines the paper's procedures are built from:
Mahler measure by the PRODUCT formula (Def. 2.1), root classification and Salem
detection (Def. 2.6), exact sign arithmetic in Q(sqrt5) (Prop. 7.5), the seed
CATALOG (Def. 2.9), the trace-down substitution (Def. 6.1), the closure-guard
ladder (Prop. 7.5), the self-action difference spectrum (Lem. 7.4), and the
Galois-correspondence signature census of K (Thm. 5.3).

These are PRODUCERS' building blocks: every function RETURNS a computed value so
the calling scripts can EMIT it to data/.  Nothing here asserts against a paper
number and nothing here imports from tests/ -- the two derivations are kept
independent on purpose.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

mp.mp.dps = 60

# --------------------------------------------------------------------------- #
#  Named catalog and numbers (Def. 2.9, Sec. 1.2, Cor. 5.5)
# --------------------------------------------------------------------------- #
# Catalog minimal polynomials, highest-degree-first coefficients (Def. 2.9).
CATALOG = {
    "phi":   [1, -1, -1],        # x^2 - x - 1    golden seed
    "tau":   [1, 1, -1],         # x^2 + x - 1    tau = phi^{-1}
    "sqrt2": [1, 0, -2],         # x^2 - 2
    "sqrt3": [1, 0, -3],         # x^2 - 3
    "sqrt5": [1, 0, -5],         # x^2 - 5
    "gap":   [1, -7, 1],         # x^2 - 7x + 1   roots phi^4, phi^-4
    "K":     [1, 0, 5, 0, -5],   # x^4 + 5x^2 - 5  Lorentzian quartic
}

# Lehmer's number L(x) = x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1 (Sec.1.2)
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]

# Minimal degree-four Salem beta_4, x^4 - x^3 - x^2 - x + 1 (Cor. 5.5)
BETA4 = [1, -1, -1, -1, 1]


# --------------------------------------------------------------------------- #
#  Mahler measure and root geometry (Def. 2.1, Def. 2.6)
# --------------------------------------------------------------------------- #
def mahler_measure(int_coeffs, dps=60):
    """Mahler measure via the PRODUCT formula  M(p) = |a_d| * prod max(1,|root|)
    (Def. 2.1).  int_coeffs are highest-degree first.  Returns an mpf."""
    with mp.workdps(dps):
        lead = abs(mp.mpf(int(int_coeffs[0])))
        roots = mp.polyroots([mp.mpf(int(c)) for c in int_coeffs],
                             maxsteps=600, extraprec=400)
        m = lead
        for r in roots:
            if abs(r) > 1:
                m *= abs(r)
        return m


def mahler_measure_jensen(int_coeffs, dps=40):
    """Mahler measure via the Jensen-integral form
    exp(int_0^1 log|p(e^{2 pi i theta})| d theta) (Def. 2.1, second expression).
    Finite only when p has no root on the unit circle."""
    with mp.workdps(dps):
        cs = [mp.mpf(int(c)) for c in int_coeffs]
        d = len(cs) - 1

        def integrand(theta):
            z = mp.e ** (2j * mp.pi * theta)
            val = sum(cs[i] * z ** (d - i) for i in range(len(cs)))
            return mp.log(abs(val))

        integral = mp.quad(integrand, [0, mp.mpf(1) / 4, mp.mpf(1) / 2,
                                       mp.mpf(3) / 4, 1])
        return mp.e ** integral


def roots_of(int_coeffs, dps=80):
    """All complex roots (mpc) of an integer polynomial."""
    with mp.workdps(dps):
        return mp.polyroots([mp.mpf(int(c)) for c in int_coeffs],
                            maxsteps=800, extraprec=600)


def root_classification(int_coeffs, tol=mp.mpf(10) ** -18, dps=80):
    """Return (n_outside, n_inside, n_oncircle) for the roots (Def. 2.6 pattern)."""
    out = ins = onc = 0
    for r in roots_of(int_coeffs, dps=dps):
        a = abs(r)
        if a > 1 + tol:
            out += 1
        elif a < 1 - tol:
            ins += 1
        else:
            onc += 1
    return out, ins, onc


def dominant_real_root(int_coeffs, dps=80):
    """Largest real root (the algebraic number a polynomial names)."""
    with mp.workdps(dps):
        reals = [mp.re(r) for r in roots_of(int_coeffs, dps=dps)
                 if abs(mp.im(r)) < mp.mpf(10) ** -20]
        return max(reals)


def is_palindromic(int_coeffs):
    """Reciprocal test: the coefficient list is a palindrome."""
    c = list(int_coeffs)
    return c == c[::-1]


def is_salem(int_coeffs):
    """True iff int_coeffs is the minimal polynomial of a Salem number (Def. 2.6):
    even degree >= 4, reciprocal, one real root > 1, one in (0,1), and the rest a
    positive number of conjugate pairs ON the unit circle."""
    d = len(int_coeffs) - 1
    if d < 4 or d % 2 != 0:
        return False
    if not is_palindromic(int_coeffs):
        return False
    out, ins, onc = root_classification(int_coeffs)
    return out == 1 and ins == 1 and onc == d - 2 and onc >= 2


# --------------------------------------------------------------------------- #
#  Exact sign arithmetic in Q(sqrt5) (Prop. 7.5)
# --------------------------------------------------------------------------- #
_S5 = sp.sqrt(5)
PHI_SYM = (1 + _S5) / 2          # exact symbolic golden ratio


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
    if a,b share a sign the sign is theirs; otherwise compare a^2 with 5 b^2 as
    integers.  Returns -1, 0, or +1 with no irrational evaluation."""
    a = sp.Rational(a)
    b = sp.Rational(b)
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    mag = 1 if a * a > 5 * b * b else (-1 if a * a < 5 * b * b else 0)
    return mag if a > 0 else -mag


# --------------------------------------------------------------------------- #
#  Trace-down and the closure guard (Def. 6.1, Lem. 6.2, Prop. 7.5)
# --------------------------------------------------------------------------- #
_x, _t = sp.symbols('x t')


def tracedown(int_coeffs):
    """Def. 6.1: for a reciprocal m_theta of degree 2m return the degree-m
    polynomial T (as a sympy expr in t) with m_theta(x) = x^m * T(x + 1/x)."""
    p = sp.Poly(int_coeffs, _x)
    d = p.degree()
    m = d // 2
    coeffs_hi = p.all_coeffs()
    a = {d - i: coeffs_hi[i] for i in range(len(coeffs_hi))}

    def s(j):                     # x^j + x^-j as a Chebyshev-like poly in t
        if j == 0:
            return sp.Integer(2)
        sm2, sm1 = sp.Integer(2), _t
        for _ in range(2, j + 1):
            sm2, sm1 = sm1, sp.expand(_t * sm1 - sm2)
        return sm1

    T = sp.Integer(a.get(m, 0))
    for j in range(1, m + 1):
        T += a.get(m + j, 0) * s(j)
    return sp.expand(T)


def tracedown_root_pattern(int_coeffs):
    """Return (totally_real, n_roots_above_2, n_roots_in_open_(-2,2), roots) for
    the trace-down T (Lem. 6.2: Salem iff totally real, one root > 2, rest in
    (-2,2))."""
    import numpy as np
    T = tracedown(int_coeffs)
    Tc = [complex(c) for c in sp.Poly(T, _t).all_coeffs()]
    roots = list(np.roots(Tc))
    totally_real = all(abs(r.imag) < 1e-9 for r in roots)
    above2 = sum(1 for r in roots if r.real > 2 + 1e-9)
    inside = sum(1 for r in roots if -2 + 1e-9 < r.real < 2 - 1e-9)
    return totally_real, above2, inside, roots


def validate(int_coeffs):
    """Prop. 7.5 closure guard.  Returns a dict with the verdict and, when a
    Salem factor is present, the exact test element m_beta(phi) = a + b*sqrt5 and
    its sign.  Verdict ladder:
        FORCED               -- no Salem factor
        FORCED_ABOVE_FLOOR   -- a Salem factor with beta >= phi (m_beta(phi) <= 0)
        INVALID_CLOSURE      -- a Salem factor with beta <  phi (m_beta(phi) >  0)
    The beta-vs-phi decision is the exact sign of m_beta(phi) in Q(sqrt5)."""
    verdict = "FORCED"
    salem_factor = None
    ab = None
    sign = None
    poly = sp.Poly(int_coeffs, _x)
    for fac, _mult in sp.factor_list(poly.as_expr())[1]:
        fc = [int(c) for c in sp.Poly(fac, _x).all_coeffs()]
        if fc and fc[0] < 0:
            fc = [-v for v in fc]
        if is_salem(fc):
            salem_factor = fc
            m_beta_phi = sp.Poly(fc, _x).as_expr().subs(_x, PHI_SYM)
            a, b = to_ab_sqrt5(m_beta_phi)
            ab = (a, b)
            sign = sign_ab_sqrt5(a, b)
            if sign > 0:                       # m_beta(phi) > 0  <=>  beta < phi
                verdict = "INVALID_CLOSURE"
                break
            else:                              # beta >= phi
                verdict = "FORCED_ABOVE_FLOOR"
    return {
        "verdict": verdict,
        "salem_factor": salem_factor,
        "m_beta_phi_ab": ab,        # (a, b) with m_beta(phi) = a + b sqrt5, or None
        "sign": sign,               # sign of m_beta(phi), or None
    }


# --------------------------------------------------------------------------- #
#  Self-action difference spectrum (Lem. 7.4)
# --------------------------------------------------------------------------- #
def selfaction_spectrum(int_coeffs):
    """Lem. 7.4: spectrum of ad_R = [R, .] equals the pairwise difference set of
    the eigenvalues of the companion R.  Built as the Kronecker form
    R (x) I - I (x) R^T.  Returns the sorted real spectrum (floats)."""
    import numpy as np
    d = len(int_coeffs) - 1
    # companion matrix of the monic integer polynomial (highest-degree first)
    comp = np.zeros((d, d))
    for i in range(d - 1):
        comp[i + 1, i] = 1.0
    for i, coeff in enumerate(int_coeffs[1:][::-1]):
        comp[i, d - 1] = -float(coeff)
    Id = np.eye(d)
    adR = np.kron(comp, Id) - np.kron(Id, comp.T)
    ev = sorted(np.linalg.eigvals(adR).real)
    return ev


# --------------------------------------------------------------------------- #
#  Galois-correspondence signature census of K (Thm. 5.3)
# --------------------------------------------------------------------------- #
def build_group():
    """G = Gal(K(i)/Q) ~ C2 x C2 x D4, order 32, as tuples g=(e1,e2,k,d):
        sqrt2   -> (-1)^e1 sqrt2
        sqrt3   -> (-1)^e2 sqrt3
        5^(1/4) -> i^k  5^(1/4)
        i       -> (-1)^d i
    Composition (ga after gb): e's and d add mod 2;
        k'' = kb*(1 + 2 da) + ka   (mod 4)
    (the twist 1+2da is the semidirect action of conjugation on 5^(1/4))."""
    elems = [(e1, e2, k, d) for e1 in (0, 1) for e2 in (0, 1)
             for k in (0, 1, 2, 3) for d in (0, 1)]

    def mul(ga, gb):
        e1a, e2a, ka, da = ga
        e1b, e2b, kb, db = gb
        return ((e1a + e1b) % 2, (e2a + e2b) % 2,
                (kb * (1 + 2 * da) + ka) % 4, (da + db) % 2)

    return elems, mul


def all_subgroups(elems, mul, E):
    """Enumerate every subgroup of (elems, mul) with identity E."""
    def close(gens):
        S = set([E]) | set(gens)
        changed = True
        while changed:
            changed = False
            for a in list(S):
                for b in list(S):
                    p = mul(a, b)
                    if p not in S:
                        S.add(p)
                        changed = True
        return frozenset(S)

    subs = set([frozenset([E])])
    frontier = [frozenset([E])]
    while frontier:
        H = frontier.pop()
        for g in elems:
            if g in H:
                continue
            H2 = close(list(H) + [g])
            if H2 not in subs:
                subs.add(H2)
                frontier.append(H2)
    return subs


def galois_census():
    """Thm. 5.3: enumerate the subfields of K via the Galois correspondence on
    K(i), and compute the signature (deg, r1, r2) of each.

    Returns a dict:
        census            -- {(deg, r1, r2): count} over the subfields of K
        n_subfields       -- number of subfields of K (= subgroups containing c)
        n_total_subgroups -- total subgroups of G ~ C2 x C2 x D4
        order_breakdown   -- {subgroup_order: count} over ALL subgroups
    Subfields of K correspond to subgroups H >= J = <c> (c = complex
    conjugation); r1 counts the real embeddings by the fixed-coset criterion
    g^{-1} c g in H."""
    elems, mul = build_group()
    E = (0, 0, 0, 0)
    c = (0, 0, 0, 1)                        # complex conjugation

    def inv(g):
        for h in elems:
            if mul(g, h) == E:
                return h
        raise AssertionError("no inverse")

    J = frozenset([E, c])                   # Gal(K(i)/K)
    subs = all_subgroups(elems, mul, E)
    subK = [H for H in subs if J <= H]      # subfields of K

    def signature(H):
        Hs = set(H)
        deg = 32 // len(H)
        seen, reps = set(), []
        for g in elems:
            coset = frozenset(mul(g, h) for h in H)
            if coset not in seen:
                seen.add(coset)
                reps.append(g)
        r1 = sum(1 for g in reps if mul(mul(inv(g), c), g) in Hs)
        return deg, r1, (deg - r1) // 2

    census = {}
    for H in subK:
        sig = signature(H)
        census[sig] = census.get(sig, 0) + 1

    order_breakdown = {}
    for H in subs:
        order_breakdown[len(H)] = order_breakdown.get(len(H), 0) + 1

    return {
        "census": census,
        "n_subfields": len(subK),
        "n_total_subgroups": len(subs),
        "order_breakdown": order_breakdown,
    }
