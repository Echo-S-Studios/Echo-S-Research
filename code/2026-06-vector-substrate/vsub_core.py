"""
vsub_core -- exact number-field engine for the PRODUCER scripts of

    "The Vector Substrate: Number Fields as Exact Learning Geometry"
    papers/2026-06-vector-substrate/vector_substrate.tex
    (AceTheDactyl, Echo S Studios Research Developments, June 2026)

This module is the reusable computational core shared by the producer scripts
in this folder.  It re-implements, from first principles and exactly over Q,
the machinery the paper describes: the regular representation / companion
matrix, the trace-form Gram, invariant factors via determinantal divisors,
heights / Mahler measure, the trace-form projector and residual, the Minkowski
embedding, and the two Fisher metrics.

It is INDEPENDENT of tests/2026-06-vector-substrate/ (nothing is imported from
there).  The difference in intent: the test modules ASSERT the paper's values;
these producers COMPUTE them and EMIT machine-readable data (see the provenance
writers `write_json` / `write_csv`).
"""
import csv
import json
import os
from itertools import combinations

import sympy as sp

x = sp.symbols("x")

# Provenance stamped onto every emitted artifact.
SOURCE_PAPER = "papers/2026-06-vector-substrate/vector_substrate.tex"
DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data",
    "2026-06-vector-substrate",
)


# ----------------------------------------------------------------------
# Regular representation (Sec. 2.1-2.2, Prop. 2.2)
# ----------------------------------------------------------------------
def companion(nonlead_coeffs):
    """Companion matrix C(m) of m = x^n + c_{n-1}x^{n-1} + ... + c_0, in the
    paper's convention (multiplication by theta on the power basis
    {1, theta, ..., theta^{n-1}}): theta^k -> theta^{k+1} for k<n-1 (ones on the
    sub-diagonal), theta^{n-1} -> -sum_j c_j theta^j (last column = -c).

    `nonlead_coeffs` = [c_0, c_1, ..., c_{n-1}] (constant term first, length n).
    """
    c = [sp.sympify(v) for v in nonlead_coeffs]
    n = len(c)
    C = sp.zeros(n, n)
    for k in range(n - 1):
        C[k + 1, k] = 1
    for j in range(n):
        C[j, n - 1] = -c[j]
    return C


def companion_from_poly(poly_expr):
    """Companion matrix from a monic sympy polynomial in x."""
    p = sp.Poly(poly_expr, x)
    if p.LC() != 1:
        raise ValueError("minimal polynomial must be monic")
    high_to_low = p.all_coeffs()               # [1, c_{n-1}, ..., c_0]
    nonlead = list(reversed(high_to_low[1:]))  # [c_0, ..., c_{n-1}]
    return companion(nonlead)


def rho(coord, C):
    """Regular representation rho(alpha) = sum_i a_i C^i of the element with
    power-basis coordinate vector `coord` = (a_0, ..., a_{n-1}) (alpha = sum a_i
    theta^i), using the algebra-homomorphism property rho(theta^i) = C^i."""
    coord = [sp.sympify(v) for v in coord]
    n = C.shape[0]
    R = sp.zeros(n, n)
    P = sp.eye(n)
    for i in range(len(coord)):
        R += coord[i] * P
        if i < n - 1:
            P = P * C
    return R


def field_trace(coord, C):
    """Tr_{K/Q}(alpha) = trace of the regular representation."""
    return sp.expand(rho(coord, C).trace())


def field_norm(coord, C):
    """N_{K/Q}(alpha) = det of the regular representation."""
    return sp.expand(rho(coord, C).det())


# ----------------------------------------------------------------------
# Trace form / Gram (Sec. 2.4-2.6, Thm. 2.14)
# ----------------------------------------------------------------------
def gram(basis_coords, C):
    """Trace-form Gram G_ij = Tr(w_i w_j) = trace(rho(w_i) rho(w_j)) for a basis
    given by power-basis coordinate vectors `basis_coords`."""
    R = [rho(b, C) for b in basis_coords]
    m = len(basis_coords)
    G = sp.zeros(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = sp.expand((R[i] * R[j]).trace())
    return G


def power_gram(C):
    """Trace-form Gram in the power basis: G_ij = Tr(theta^{i+j}) = trace(C^{i+j})."""
    n = C.shape[0]
    powers = [sp.eye(n)]
    for _ in range(2 * n):
        powers.append(powers[-1] * C)
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = powers[i + j].trace()
    return G


def trace_vector(basis_coords, C):
    """Trace vector t_i = Tr(w_i)."""
    return sp.Matrix([field_trace(b, C) for b in basis_coords])


# ----------------------------------------------------------------------
# Invariant factors / Smith normal form (Sec. 2.3, Thm. 4.2 bridge)
# ----------------------------------------------------------------------
def invariant_factors(A):
    """Non-constant invariant factors (monic, ascending degree) of xI - A,
    computed from determinantal divisors: D_i = monic gcd of all i x i minors of
    (xI - A), invariant factors s_i = D_i / D_{i-1}.  Largest = minimal
    polynomial; product = characteristic polynomial.  Independent of any library
    Smith-normal-form routine."""
    n = A.shape[0]
    XA = x * sp.eye(n) - A
    D = [sp.Poly(1, x)]
    for i in range(1, n + 1):
        g = None
        for rows in combinations(range(n), i):
            for cols in combinations(range(n), i):
                minor = sp.Poly(XA[list(rows), list(cols)].det(), x)
                g = minor if g is None else g.gcd(minor)
        D.append(g.monic() if (g is not None and g.as_expr() != 0) else sp.Poly(0, x))
    facs = []
    for i in range(1, n + 1):
        q = sp.Poly(sp.div(D[i].as_expr(), D[i - 1].as_expr(), x)[0], x)
        if q.degree() >= 1:
            facs.append(q.monic().as_expr())
    return facs


def matrix_minpoly(A):
    """Minimal polynomial (monic) of a rational square matrix A, as the
    least-degree monic annihilator found by vectorized power dependence."""
    n = A.shape[0]
    cols = [sp.eye(n).reshape(n * n, 1)]
    P = sp.eye(n)
    for _ in range(1, n + 1):
        P = P * A
        cols.append(P.reshape(n * n, 1))
        ns = sp.Matrix.hstack(*cols).nullspace()
        if ns:
            v = ns[0]
            lead = v[-1]
            return sp.expand(sum((v[k] / lead) * x ** k for k in range(len(v))))
    return A.charpoly(x).as_expr()


# ----------------------------------------------------------------------
# Heights / Mahler measure (Sec. 5, Thm. 5.2)
# ----------------------------------------------------------------------
def coeff_height(poly_expr):
    """ch(m) = max_i |c_i| over the NON-leading coefficients of a monic poly."""
    p = sp.Poly(poly_expr, x)
    nonlead = p.all_coeffs()[1:]
    return max((abs(int(c)) for c in nonlead), default=0)


def two_norm_sq(poly_expr):
    """||p||_2^2 = sum of squares of ALL integer coefficients (incl. leading)."""
    p = sp.Poly(poly_expr, x)
    return sum(int(c) ** 2 for c in p.all_coeffs())


def mahler_measure_mp(poly_expr, dps=50):
    """Mahler measure |a_n| * prod max(1, |root|), high precision via mpmath."""
    import mpmath as mp

    mp.mp.dps = dps
    p = sp.Poly(poly_expr, x)
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    lead = coeffs[0]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    prod = abs(lead)
    for r in roots:
        if abs(r) > 1:
            prod *= abs(r)
    return prod


def admissible(poly_expr, Dmax, Hmax):
    """Northcott admissibility (Def. 5.4): deg <= Dmax and ch <= Hmax; a pure
    integer decision that never touches the (float) Mahler measure."""
    p = sp.Poly(poly_expr, x)
    return (p.degree() <= Dmax) and (coeff_height(poly_expr) <= Hmax)


def log_interval(value, places=5, dps=60):
    """Rigorous rational enclosure [lo, hi] of log(value), rounded OUTWARD to
    `places` decimals, so `hi` is a certified rational UPPER bound on the cost
    side.  Used for the certified-GROW decisions (Ex. 8.7): an exact rational
    gain is compared against `lambda * hi`."""
    import mpmath as mp

    mp.mp.dps = dps
    v = sp.sympify(value)
    L = mp.log(mp.mpf(str(sp.N(v, dps))))
    Lr = sp.Rational(mp.nstr(L, dps - 5))      # exact rational within ~1e-50 of L
    scale = sp.Integer(10) ** places
    lo = sp.Rational(sp.floor(Lr * scale), scale)
    hi = sp.Rational(sp.ceiling(Lr * scale), scale)
    return lo, hi


# ----------------------------------------------------------------------
# Projector / residual (Sec. 3)
# ----------------------------------------------------------------------
def projector(B, G):
    """Trace-form projector P = B (B^T G B)^{-1} B^T G onto col(B)."""
    B = sp.Matrix(B)
    G = sp.Matrix(G)
    inner = B.T * G * B
    return B * inner.inv() * B.T * G


def residual(xcoord, B, G):
    """Residual r = x - P x (coordinate vector)."""
    P = projector(B, G)
    xv = sp.Matrix(xcoord)
    return xv - P * xv


def gnorm2(v, G):
    """Squared G-norm v^T G v."""
    v = sp.Matrix(v)
    return sp.expand((v.T * sp.Matrix(G) * v)[0])


# ----------------------------------------------------------------------
# Fisher metrics (Sec. 7)
# ----------------------------------------------------------------------
def fisher_exp(G, t, n):
    """Conjugate-family (max-entropy) Fisher metric at the uniform point,
    Fisher_exp = (1/n) G - (1/n^2) t t^T  (Thm. 7.7)."""
    G = sp.Matrix(G)
    t = sp.Matrix(t)
    return sp.simplify(G / n - (t * t.T) / n ** 2)


# ----------------------------------------------------------------------
# Numeric Minkowski embeddings (cross-check G = M^T M, Thm. 2.14)
# ----------------------------------------------------------------------
def embedding_matrix_power(min_poly_expr, dps=50):
    """Numeric embedding matrix M in the power basis: M[k, i] = sigma_k(theta)^i,
    sigma_k over the complex roots of the minimal polynomial."""
    import mpmath as mp

    mp.mp.dps = dps
    p = sp.Poly(min_poly_expr, x)
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=300, extraprec=300)
    n = len(roots)
    M = mp.matrix(n, n)
    for k in range(n):
        for i in range(n):
            M[k, i] = roots[k] ** i
    return M


# ----------------------------------------------------------------------
# Provenance-stamped writers + serialization helpers
# ----------------------------------------------------------------------
def sval(v):
    """Serialize an exact sympy scalar: int -> int, rational -> 'a/b',
    otherwise a compact string (e.g. 'sqrt(5)')."""
    v = sp.sympify(v)
    if v.is_Integer:
        return int(v)
    if v.is_Rational:
        return str(v)
    return sp.sstr(v)


def mat_to_list(M):
    """A sympy matrix as a nested list of serialized scalars (row-major)."""
    M = sp.Matrix(M)
    return [[sval(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]


def write_json(filename, payload, script):
    """Write `payload` (a dict) to data/2026-06-vector-substrate/<filename> with
    provenance fields _source_paper and _generated_by prepended."""
    os.makedirs(DATA_DIR, exist_ok=True)
    out = {
        "_source_paper": SOURCE_PAPER,
        "_generated_by": f"code/2026-06-vector-substrate/{script}",
    }
    out.update(payload)
    path = os.path.join(DATA_DIR, filename)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return path


def write_csv(filename, header, rows, script):
    """Write a CSV to data/2026-06-vector-substrate/<filename> with a leading
    provenance comment line."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, filename)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(
            f"# source: {SOURCE_PAPER}; "
            f"generated by: code/2026-06-vector-substrate/{script}\n"
        )
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    return path
