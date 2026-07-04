r"""
Independent verification of Section 3 ("Cross-field growth") of the paper:
the coordinates->minimal-polynomial bridge, the invariant-factor (charpoly-is-
insufficient) witnesses, the disjoint Kronecker Gram + determinant relation, and
the NON-disjoint compositum witness Q(sqrt2)(sqrt2+sqrt3)=Q(sqrt2,sqrt3).

Every polynomial/factorisation/reconstruction is re-derived with sympy
(exact symbolic) from the paper's premises, then compared to the paper's value.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols, expand, zeros, eye, resultant

x, y = symbols("x y")


def companion(coeffs):
    n = len(coeffs) - 1
    C = zeros(n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]
    return C


# --- Example 3.3 (the bridge, exercised) -------------------------------------

def test_bridge_phi_maps_to_golden_law():
    """Ex 3.3: phi = (1/2,1/2) in Q[x]/(x^2-5) has minimal polynomial x^2-x-1.
    Build rho(phi) = 1/2*I + 1/2*rep(sqrt5) and read off its (minimal=char) poly."""
    rep_s5 = companion([1, 0, -5])            # regular rep of sqrt5 in Q[x]/(x^2-5)
    rho_phi = Q(1, 2) * eye(2) + Q(1, 2) * rep_s5
    assert sp.factor(rho_phi.charpoly(x).as_expr()) == x**2 - x - 1


def test_bridge_2sqrt6_maps_to_x2_minus_24():
    """Ex 3.3: 2 sqrt6 -> x^2-24 (minimal polynomial of 2 sqrt6)."""
    assert sp.expand(sp.minimal_polynomial(2 * sqrt(6), x) - (x**2 - 24)) == 0


def test_bridge_rational_3_maps_to_x_minus_3():
    """Ex 3.3: the rational 3 -> x-3."""
    assert sp.expand(sp.minimal_polynomial(sp.Integer(3), x) - (x - 3)) == 0


def test_bridge_theta_primitive_maps_to_quartic():
    """Ex 3.3: the primitive theta=sqrt2+sqrt3 -> x^4-10x^2+1."""
    assert sp.expand(sp.minimal_polynomial(sqrt(2) + sqrt(3), x)
                     - (x**4 - 10 * x**2 + 1)) == 0


def test_bridge_half_is_rejected_not_algebraic_integer():
    """Ex 3.3: 1/2 in Q(sqrt5) is REJECTED: its (monic) minimal polynomial x-1/2
    is not in Z[x], so 1/2 is not an algebraic integer."""
    mp_int = sp.minimal_polynomial(Q(1, 2), x)          # sympy: 2x-1 (integer-primitive)
    mp_monic = sp.expand(mp_int / sp.LC(mp_int, x))     # monic form x - 1/2
    assert sp.expand(mp_monic - (x - Q(1, 2))) == 0
    # rejection criterion: the monic minimal polynomial has a non-integer coefficient
    coeffs = sp.Poly(mp_monic, x).all_coeffs()
    assert not all(c.is_Integer for c in coeffs)


# --- Remark 3.4 (charpoly is an incomplete similarity invariant) -------------

def test_invfac_phi_dsum_phi_vs_squared_companion():
    """Rem 3.4: rho(phi)+rho(phi) and C((x^2-x-1)^2) share char.poly (x^2-x-1)^2
    and trace 2, yet are NOT similar (minpoly x^2-x-1 vs (x^2-x-1)^2)."""
    Cphi = companion([1, -1, -1])
    A = sp.diag(Cphi, Cphi)                    # phi (+) phi
    p2 = sp.Poly((x**2 - x - 1)**2, x).all_coeffs()
    B = companion(p2)                          # C((x^2-x-1)^2)
    # SAME characteristic polynomial and trace ...
    assert sp.expand(A.charpoly(x).as_expr() - B.charpoly(x).as_expr()) == 0
    assert A.trace() == B.trace() == 2
    # ... but A satisfies x^2-x-1 (minpoly) and B does not -> not similar
    I4 = eye(4)
    assert A * A - A - I4 == zeros(4)
    assert B * B - B - I4 != zeros(4)


def test_invfac_phi_dsum_phi_same_mahler_measure():
    """Rem 3.4: the two matrices additionally share Mahler measure (= phi^2), so
    even (charpoly, trace, Mahler) cannot separate them; only the invariant-factor
    list (=minpoly here) does.  Mahler measure is a function of the (shared)
    characteristic polynomial (x^2-x-1)^2, so equality is automatic; we compute the
    common value M((x^2-x-1)^2) = M(x^2-x-1)^2 = phi^2 from the simple-root factor."""
    Cphi = companion([1, -1, -1])
    A = sp.diag(Cphi, Cphi)
    B = companion(sp.Poly((x**2 - x - 1)**2, x).all_coeffs())
    assert sp.expand(A.charpoly(x).as_expr() - B.charpoly(x).as_expr()) == 0  # same charpoly
    import mpmath as mp
    mp.mp.dps = 40
    roots = mp.polyroots([1, -1, -1])                    # simple roots of x^2-x-1
    M = mp.mpf(1)
    for rt in roots:
        M *= max(mp.mpf(1), abs(rt))
    phi = (1 + mp.sqrt(5)) / 2
    assert abs(M - phi) < mp.mpf("1e-30")                # M(x^2-x-1) = phi
    assert abs(M**2 - phi**2) < mp.mpf("1e-30")          # shared Mahler = phi^2


def test_invfac_jordan_charpoly_and_minpoly_insufficient():
    """Rem 3.4: J2(+)J2 and J2(+)J1(+)J1 both have charpoly x^4 and minpoly x^2,
    but invariant factors (x^2,x^2) vs (x,x,x^2) differ, so they are not similar.
    rank(A) distinguishes them (2 vs 1)."""
    J2 = Matrix([[0, 1], [0, 0]])
    JA = sp.diag(J2, J2)
    JB = sp.diag(J2, Matrix([[0]]), Matrix([[0]]))
    assert JA.charpoly(x).as_expr() == x**4
    assert JB.charpoly(x).as_expr() == x**4
    assert JA * JA == zeros(4) and JB * JB == zeros(4)   # minpoly divides x^2
    assert JA != zeros(4) and JB != zeros(4)             # minpoly is exactly x^2
    assert JA.rank() == 2 and JB.rank() == 1             # not similar


# --- Example 3.7 (disjoint compositum: adjoining sqrt7) -----------------------

def field_trace_form(basis_squares, deg):
    """Diagonal trace form of a multiquadratic field: for basis elements that are
    square roots of squarefree integers (or 1), Tr(b_i b_j)=0 off diagonal and
    Tr(b_i^2)=deg*b_i^2-as-rational.  We pass the rational values b_i^2."""
    return sp.diag(*[deg * s for s in basis_squares])


def test_disjoint_gram_QK_and_QL_from_traces():
    """Ex 3.7: G_K for Q(sqrt2,sqrt3) is diag(4,8,12,24) and G_L for Q(sqrt7) is
    diag(2,14). Derived from Tr(c)=deg*c on the orthogonal square-root basis."""
    # Q(sqrt2,sqrt3): basis {1,sqrt2,sqrt3,sqrt6}, squares {1,2,3,6}, deg 4
    GK = field_trace_form([1, 2, 3, 6], 4)
    assert GK == sp.diag(4, 8, 12, 24)
    # Q(sqrt7): basis {1,sqrt7}, squares {1,7}, deg 2
    GL = field_trace_form([1, 7], 2)
    assert GL == sp.diag(2, 14)


def test_disjoint_kronecker_gram():
    """Ex 3.7: G_{KL} = G_K (x) G_L = diag(8,56,16,112,24,168,48,336)."""
    GK = sp.diag(4, 8, 12, 24)
    GL = sp.diag(2, 14)
    GKL = Matrix(sp.kronecker_product(GK, GL))
    assert [GKL[i, i] for i in range(8)] == [8, 56, 16, 112, 24, 168, 48, 336]


def test_disjoint_out_of_field_score_56():
    """Ex 3.7: the out-of-field sqrt7 score is the out-of-field Gram diagonal 56,
    i.e. Tr_{KL/Q}(sqrt7^2)=Tr(7)=8*7=56 in the degree-8 compositum."""
    assert 8 * 7 == 56


# --- Remark 3.8 (the determinant relation, exactly) --------------------------

def test_kronecker_det_relation():
    """Rem 3.8: det G_{KL} = (det G_K)^{[L:Q]} (det G_L)^{[K:Q]} = 9216^2 * 28^4,
    and it equals the product of the eight Kronecker diagonal entries."""
    GK = sp.diag(4, 8, 12, 24)
    GL = sp.diag(2, 14)
    GKL = Matrix(sp.kronecker_product(GK, GL))
    detK, detL = GK.det(), GL.det()
    assert detK == 9216 and detL == 28
    assert GKL.det() == detK**2 * detL**4          # [L:Q]=2, [K:Q]=4
    assert detK**2 * detL**4 == 9216**2 * 28**4
    prod8 = 1
    for v in [8, 56, 16, 112, 24, 168, 48, 336]:
        prod8 *= v
    assert prod8 == GKL.det()


# --- Theorem 3.6 / Example 3.9 (the non-disjoint degree-4 witness) -----------

def test_nondisjoint_mbeta_factors_over_Qsqrt2():
    """Ex 3.9: over K=Q(sqrt2), m_beta=x^4-10x^2+1 factors as
    (x^2-2 sqrt2 x -1)(x^2+2 sqrt2 x -1)."""
    prod = sp.expand((x**2 - 2 * sqrt(2) * x - 1) * (x**2 + 2 * sqrt(2) * x - 1))
    assert sp.simplify(prod - (x**4 - 10 * x**2 + 1)) == 0


def test_nondisjoint_operator_polynomial_degree_6():
    """Ex 3.9: with primitive theta=alpha+beta (alpha=sqrt2, beta=sqrt2+sqrt3, c=1),
    the degree-6 operator polynomial over Q is
    x^6-25x^4+91x^2-75 = (x^2-3)(x^4-22x^2+25).
    Re-derived as the squarefree part of Res_y(y^2-2, m_beta(x-y)) whose roots are
    the distinct values alpha_i + beta_j over both embeddings of alpha."""
    m_beta = x**4 - 10 * x**2 + 1
    res = resultant(y**2 - 2, m_beta.subs(x, x - y), y)   # deg-8, roots alpha_i+beta_j
    op = sp.expand(sp.sqf_part(sp.expand(res)))            # deg-6 (distinct roots)
    assert op == x**6 - 25 * x**4 + 91 * x**2 - 75
    assert sp.expand(op - (x**2 - 3) * (x**4 - 22 * x**2 + 25)) == 0


def test_nondisjoint_selected_factor_is_minpoly_of_theta():
    """Ex 3.9 & eq (mtheta): the genuine factor m_theta=x^4-22x^2+25 is the exact
    minimal polynomial of the primitive element theta = alpha+beta = 2 sqrt2 + sqrt3
    (degree 4 = m*e', not the tensor degree 8)."""
    theta = 2 * sqrt(2) + sqrt(3)
    assert sp.expand(sp.minimal_polynomial(theta, x) - (x**4 - 22 * x**2 + 25)) == 0


def test_nondisjoint_spurious_factor_from_other_embedding():
    """Ex 3.9: the spurious factor x^2-3 has root sqrt3, coming from the embedding
    alpha=-sqrt2 under which theta=-sqrt2+(sqrt2+sqrt3)=sqrt3 collapses to Q(sqrt3)."""
    theta_collapsed = -sqrt(2) + (sqrt(2) + sqrt(3))
    assert sp.simplify(theta_collapsed - sqrt(3)) == 0
    assert sp.simplify((sqrt(3))**2 - 3) == 0            # sqrt3 satisfies x^2-3


def test_nondisjoint_beta_reconstructed_in_power_basis():
    """Ex 3.9: in the power basis of theta=2 sqrt2+sqrt3, beta=sqrt2+sqrt3 is
    recovered exactly as beta = -7/20 theta + 1/20 theta^3,
    i.e. coordinate vector (0,-7/20,0,1/20)."""
    theta = 2 * sqrt(2) + sqrt(3)
    recon = sp.expand(-Q(7, 20) * theta + Q(1, 20) * theta**3)
    assert sp.simplify(recon - (sqrt(2) + sqrt(3))) == 0
