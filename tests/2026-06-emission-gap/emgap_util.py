"""Independent re-derivation helpers for the Emission-Gap paper tests.

Everything here rebuilds quantities from first principles (roots, Mahler
measures, trace-downs, field signatures, self-actions) so the test files can
assert the paper's stated values against a genuinely independent computation.

Precision policy: transcendental / root checks use mpmath at >= 40 dps.
Exact structural decisions (Salem flip-straddle, root counting) use sympy Sturm.
"""
import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 50

x, t = sp.symbols("x t")

PHI = (1 + mp.sqrt(5)) / 2                    # golden ratio
PSI = (1 - mp.sqrt(5)) / 2                    # conjugate
MU_S = mp.findroot(lambda z: z**3 - z - 1, 1.324)   # plastic number (Smyth)

# ---------------------------------------------------------------------------
# roots / Mahler measure
# ---------------------------------------------------------------------------

def mp_roots(coeffs):
    """High-precision roots of a polynomial given highest-degree-first coeffs."""
    cc = [mp.mpf(str(c)) if not isinstance(c, mp.mpc) else c for c in coeffs]
    return mp.polyroots(cc, maxsteps=2000, extraprec=400)


def mahler(coeffs):
    """Mahler measure M(p) = |a_n| * prod max(1,|root|), coeffs highest-first."""
    cc = [mp.mpf(str(c)) for c in coeffs]
    lead = abs(cc[0])
    rts = mp.polyroots(cc, maxsteps=2000, extraprec=400)
    prod = mp.mpf(1)
    for r in rts:
        m = abs(r)
        if m > 1:
            prod *= m
    return lead * prod


def mahler_poly(poly):
    """Mahler measure of a sympy Poly / expr in x."""
    P = as_poly(poly)
    return mahler([int(c) if c == int(c) else float(c) for c in P.all_coeffs()])


# ---------------------------------------------------------------------------
# sympy poly plumbing
# ---------------------------------------------------------------------------

def as_poly(obj, var=x):
    if isinstance(obj, sp.Poly):
        return obj
    if isinstance(obj, (list, tuple)):
        return sp.Poly(list(obj), var)
    return sp.Poly(obj, var)


def is_palindromic(coeffs):
    """Self-reciprocal (palindromic) coefficient list (highest-first)."""
    c = [sp.nsimplify(v) for v in coeffs]
    return c == c[::-1]


# ---------------------------------------------------------------------------
# trace-down  R(x) = x^m T(x + 1/x)  and the Salem flip-straddle
# ---------------------------------------------------------------------------

def trace_down(Rpoly):
    """Return T(t) (sympy Poly) for a monic, self-reciprocal R of even degree 2m,
    using  x^k + x^{-k} = c_k(t)  with c_0=2, c_1=t, c_{k+1}=t c_k - c_{k-1}."""
    R = as_poly(Rpoly)
    deg = R.degree()
    assert deg % 2 == 0, "trace-down needs even degree"
    m = deg // 2
    coeffs = R.all_coeffs()            # highest-first, length 2m+1
    a = list(reversed(coeffs))          # a[k] = coeff of x^k
    c = [sp.Integer(2), t]              # c_0, c_1
    for k in range(1, m):
        c.append(sp.expand(t * c[k] - c[k - 1]))
    T = sp.sympify(a[m])
    for k in range(1, m + 1):
        T += a[m + k] * c[k]
    return sp.Poly(sp.expand(T), t)


def flip_straddle(Tpoly):
    """Salem flip-straddle (Lemma 8.2): T totally real, exactly one root in
    (2, inf), the remaining m-1 roots in (-2, 2)."""
    T = as_poly(Tpoly, t)
    m = T.degree()
    # no root exactly on the flip boundary for a genuine Salem trace-down
    if T.eval(2) == 0 or T.eval(-2) == 0:
        return False
    n_real = T.count_roots()
    if n_real != m:
        return False
    n_out = T.count_roots(sp.Integer(2), sp.oo)
    n_in = T.count_roots(sp.Integer(-2), sp.Integer(2))
    return n_out == 1 and n_in == m - 1


# ---------------------------------------------------------------------------
# Salem detection
# ---------------------------------------------------------------------------

def is_salem_numeric(coeffs, tol=1e-7):
    """Numeric Salem test on a polynomial's root geometry (highest-first)."""
    deg = len(coeffs) - 1
    if deg < 4 or deg % 2:
        return False
    rts = np.roots([complex(c) for c in coeffs])
    mods = np.abs(rts)
    outside = int(np.sum(mods > 1 + tol))
    inside = int(np.sum(mods < 1 - tol))
    oncirc = int(np.sum(np.abs(mods - 1) <= tol))
    if outside != 1 or inside != 1 or oncirc < 2:
        return False
    beta = rts[int(np.argmax(mods))]
    return abs(beta.imag) <= tol and beta.real > 1


def salem_factors(poly):
    """Return list of irreducible integer factors that are Salem polynomials,
    decided exactly via palindromy + trace-down flip-straddle (the paper's guard)."""
    P = as_poly(poly)
    out = []
    _, facs = sp.factor_list(P.as_expr(), x)
    for fac, _mult in facs:
        Q = sp.Poly(fac, x)
        d = Q.degree()
        if d < 4 or d % 2:
            continue
        coeffs = Q.all_coeffs()
        lead = coeffs[0]
        if lead == -1:
            coeffs = [-c for c in coeffs]
            Q = sp.Poly(coeffs, x)
            lead = coeffs[0]
        if lead != 1:                     # Salem numbers are algebraic integers
            continue
        if not is_palindromic(coeffs):
            continue
        T = trace_down(Q)
        if flip_straddle(T):
            out.append(Q)
    return out


def has_salem_factor(poly):
    return len(salem_factors(poly)) > 0


# ---------------------------------------------------------------------------
# field signatures & the trace form
# ---------------------------------------------------------------------------

def signature_from_minpoly(coeffs):
    """(r1, r2) of Q(theta) from an irreducible min poly (highest-first)."""
    P = sp.Poly(coeffs, x)
    deg = P.degree()
    r1 = P.count_roots()
    r2 = (deg - r1) // 2
    return (r1, r2)


def field_signature(gen_expr):
    """(r1, r2) and degree of Q(gen_expr) via its minimal polynomial (exact)."""
    mpoly = sp.minimal_polynomial(gen_expr, x)
    P = sp.Poly(mpoly, x)
    deg = P.degree()
    r1 = P.count_roots()
    r2 = (deg - r1) // 2
    return (r1, r2), deg


def power_sums(coeffs, K):
    """Newton's identities: power sums p_0..p_K of the roots of a monic
    integer polynomial (coeffs highest-first, leading == 1)."""
    c = [sp.Rational(v) for v in coeffs]
    assert c[0] == 1, "power_sums needs a monic polynomial"
    n = len(c) - 1
    cs = c[1:]                             # cs[i-1] = coeff c_i of x^{n-i}
    p = [sp.Integer(n)]                    # p_0 = n
    for k in range(1, K + 1):
        if k <= n:
            s = sum(cs[i - 1] * p[k - i] for i in range(1, k))
            pk = -s - k * cs[k - 1]
        else:
            s = sum(cs[i - 1] * p[k - i] for i in range(1, n + 1))
            pk = -s
        p.append(pk)
    return p


def trace_form_signature(coeffs):
    """Signature (n_pos, n_neg) of the trace form Tr_{K/Q}(z^2) on Q(theta),
    theta a root of the given monic irreducible poly. Built from exact power
    sums G_{ij} = Tr(theta^{i+j}) = p_{i+j}; signs of eigenvalues counted."""
    P = sp.Poly(coeffs, x)
    n = P.degree()
    p = power_sums(coeffs, 2 * n - 2)
    G = np.array([[float(p[i + j]) for j in range(n)] for i in range(n)],
                 dtype=float)
    G = (G + G.T) / 2
    ev = np.linalg.eigvalsh(G)
    scale = max(1.0, float(np.max(np.abs(ev))))
    pos = int(np.sum(ev > 1e-9 * scale))
    neg = int(np.sum(ev < -1e-9 * scale))
    return pos, neg


# ---------------------------------------------------------------------------
# companions & self-action (ad_M)
# ---------------------------------------------------------------------------

def companion(coeffs):
    """Companion matrix (numpy) of a monic polynomial, coeffs highest-first."""
    P = sp.Poly(coeffs, x)
    assert P.LC() == 1
    c = [float(v) for v in P.all_coeffs()]   # [1, c1, ..., cn]
    n = P.degree()
    M = np.zeros((n, n))
    M[1:, :-1] = np.eye(n - 1)
    M[:, -1] = [-c[n - k] for k in range(n)]   # last column: -c_n,...,-c_1
    return M


def ad_operator(M):
    """Matrix of ad_M : X -> M X - X M on the space of matrices (as a flattened
    linear operator), using the Kronecker identity ad_M = M (x) I - I (x) M^T."""
    n = M.shape[0]
    I = np.eye(n)
    return np.kron(M, I) - np.kron(I, M.T)
