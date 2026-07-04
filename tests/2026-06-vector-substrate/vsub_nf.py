"""
Shared, exact number-field helpers for INDEPENDENT verification of
"The Vector Substrate: Number Fields as Exact Learning Geometry".

Nothing here restates the paper's numbers: every routine re-derives the
regular representation, trace form, invariant factors, heights, etc. from
first principles (sympy, exact over Q). Tests import this module and compare
the freshly-derived objects against the paper's stated values.
"""
import sympy as sp
from itertools import combinations

x = sp.symbols('x')


# ----------------------------------------------------------------------
# Regular representation
# ----------------------------------------------------------------------
def companion(nonlead_coeffs):
    """Companion matrix C(m), paper convention, for
        m = x^n + c_{n-1} x^{n-1} + ... + c_1 x + c_0.
    `nonlead_coeffs` = [c_0, c_1, ..., c_{n-1}] (constant term first, length n).

    This is multiplication-by-theta on the power basis {1,theta,...,theta^{n-1}}:
    theta^k -> theta^{k+1} for k<n-1 (sub-diagonal ones), and
    theta^{n-1} -> -sum_j c_j theta^j (last column = -c).
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
    """Companion matrix from a monic sympy polynomial expression in x."""
    p = sp.Poly(poly_expr, x)
    assert p.LC() == 1, "minimal polynomial must be monic"
    high_to_low = p.all_coeffs()          # [1, c_{n-1}, ..., c_0]
    nonlead = list(reversed(high_to_low[1:]))  # [c_0, ..., c_{n-1}]
    return companion(nonlead)


def rho(coord, C):
    """Regular representation of the element with power-basis coordinate vector
    `coord` = (a_0,...,a_{n-1}) i.e. sum_i a_i theta^i, as rho = sum_i a_i C^i.
    Uses the algebra-homomorphism property rho(theta^i)=C^i."""
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
    """Field trace Tr_{K/Q}(alpha) = trace of the regular representation."""
    return sp.expand(rho(coord, C).trace())


def field_norm(coord, C):
    """Field norm N_{K/Q}(alpha) = det of the regular representation."""
    return sp.expand(rho(coord, C).det())


# ----------------------------------------------------------------------
# Trace form / Gram
# ----------------------------------------------------------------------
def gram(basis_coords, C):
    """Trace-form Gram G_ij = Tr(w_i w_j) = trace(rho(w_i) rho(w_j)) for an
    arbitrary basis given by power-basis coordinate vectors `basis_coords`."""
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
    """t_i = Tr(w_i)."""
    return sp.Matrix([field_trace(b, C) for b in basis_coords])


# ----------------------------------------------------------------------
# Invariant factors via determinantal divisors (independent of any library SNF)
# ----------------------------------------------------------------------
def invariant_factors(A):
    """Invariant factors of xI - A for a rational square matrix A, computed from
    determinantal divisors: D_i = monic gcd of all i x i minors of (xI - A);
    the invariant factors are s_i = D_i / D_{i-1}. Returns the list of
    non-constant invariant factors (monic), ascending by degree.
    (Largest = minimal polynomial; product = characteristic polynomial.)"""
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
    """Minimal polynomial (monic) of a rational square matrix A, found as the
    least-degree monic annihilator via vectorized power dependence."""
    n = A.shape[0]
    cols = [sp.eye(n).reshape(n * n, 1)]
    P = sp.eye(n)
    for d in range(1, n + 1):
        P = P * A
        cols.append(P.reshape(n * n, 1))
        ns = sp.Matrix.hstack(*cols).nullspace()
        if ns:
            v = ns[0]
            lead = v[-1]
            return sp.expand(sum((v[k] / lead) * x**k for k in range(len(v))))
    return A.charpoly(x).as_expr()


# ----------------------------------------------------------------------
# Heights / Mahler measure
# ----------------------------------------------------------------------
def coeff_height(poly_expr):
    """ch(m) = max_i |c_i| over NON-leading coefficients of a monic poly."""
    p = sp.Poly(poly_expr, x)
    nonlead = p.all_coeffs()[1:]
    return max((abs(int(c)) for c in nonlead), default=0)


def two_norm_sq(poly_expr):
    """||p||_2^2 = sum of squares of ALL integer coefficients (incl. leading)."""
    p = sp.Poly(poly_expr, x)
    return sum(int(c)**2 for c in p.all_coeffs())


def mahler_measure_mp(poly_expr, dps=50):
    """Mahler measure |a_n| * prod max(1,|root|), high precision via mpmath."""
    import mpmath as mp
    mp.mp.dps = dps
    p = sp.Poly(poly_expr, x)
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    lead = coeffs[0]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    prod = abs(lead)
    for r in roots:
        m = abs(r)
        if m > 1:
            prod *= m
    return prod


def admissible(poly_expr, Dmax, Hmax):
    """Northcott admissibility: deg <= Dmax and ch <= Hmax (pure integer decision)."""
    p = sp.Poly(poly_expr, x)
    return (p.degree() <= Dmax) and (coeff_height(poly_expr) <= Hmax)


# ----------------------------------------------------------------------
# Projector / residual
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
# Numeric embeddings (for independent cross-check of G = M^T M)
# ----------------------------------------------------------------------
def embedding_matrix_power(min_poly_expr, dps=50):
    """Numeric embedding matrix M in the POWER basis: M[k,i] = sigma_k(theta)^i,
    where sigma_k runs over the complex roots of the minimal polynomial."""
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
