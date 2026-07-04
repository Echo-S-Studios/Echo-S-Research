"""Part 2 -- The Lie layer: an sl2 with golden root length (Sec. 2).

Independent re-derivations of the self-action spectrum, root spaces, the
sl2-triple, Lorentzian signature (2,1), the rational structure constants,
and the null-frame transition (incl. the conjugator V).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from _eap_helpers import R, I2, phi, psi, sqrt5, mat_eq, is_zero, comm

# eigenvectors of R (Sec 2.2): v_phi=(1,phi), v_psi=(1,psi)
vphi = sp.Matrix([1, phi])
vpsi = sp.Matrix([1, psi])

# canonical Lie-layer elements
H = 2 * R - I2                      # Sec 2.3
Np = vphi * vpsi.T                  # Thm 2.5:  N_+ = v_phi v_psi^T
Nm = vpsi * vphi.T                  # Thm 2.5:  N_- = v_psi v_phi^T
S = sp.Matrix([[2, 1], [1, -2]])    # Prop 2.7:  S = N_+ + N_-
J = sp.Matrix([[0, -1], [1, 0]])    # Prop 2.7


def _adR_matrix_on_M2():
    """Matrix of ad_R = [R, .] on M2 in the basis E11,E12,E21,E22."""
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]),
             sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    cols = []
    for E in basis:
        img = comm(R, E)
        cols.append([img[0, 0], img[0, 1], img[1, 0], img[1, 1]])
    return sp.Matrix(cols).T


# ----------------------------------------------------------------------------
# Prop 2.2 -- spectrum of the self-action
# ----------------------------------------------------------------------------
def test_adR_spectrum():
    """Prop 2.2: 'spec(ad_R) = {0,0,+sqrt5,-sqrt5}' on M2."""
    M = _adR_matrix_on_M2()
    eigs = M.eigenvals()             # dict value:multiplicity
    # normalise keys
    got = {}
    for val, mult in eigs.items():
        got[sp.simplify(val)] = mult
    assert got.get(sp.Integer(0)) == 2
    assert got.get(sp.simplify(sqrt5)) == 1
    assert got.get(sp.simplify(-sqrt5)) == 1


# ----------------------------------------------------------------------------
# Thm 2.5 -- root spaces
# ----------------------------------------------------------------------------
def test_root_vectors_explicit_entries():
    """Thm 2.5: 'N_+=[[1,psi],[phi,-1]]', 'N_-=[[1,phi],[psi,-1]]'."""
    assert mat_eq(Np, sp.Matrix([[1, psi], [phi, -1]]))
    assert mat_eq(Nm, sp.Matrix([[1, phi], [psi, -1]]))


def test_eigenvectors_orthogonal():
    """Sec 2.2: 'v_phi^T v_psi = 1 + phi*psi = 0'."""
    assert is_zero((vphi.T * vpsi)[0])


def test_root_eigen_equations():
    """Thm 2.5: '[R,N_+]=+sqrt5 N_+', '[R,N_-]=-sqrt5 N_-'."""
    assert mat_eq(comm(R, Np), sqrt5 * Np)
    assert mat_eq(comm(R, Nm), -sqrt5 * Nm)


def test_zero_space_is_span_I_R():
    """Thm 2.5: the 0-eigenspace of ad_R (matrices commuting with R) is
    exactly span{I,R}.  Solve RX=XR and check the solution space is 2-dim
    and spanned by I,R."""
    p, q, r, s = sp.symbols('p q r s')
    X = sp.Matrix([[p, q], [r, s]])
    eqs = list(comm(R, X))
    sol = sp.linsolve(eqs, [p, q, r, s])
    # parametrise: free params -> dimension 2
    (solset,) = sol
    free = sorted(sp.Matrix(solset).free_symbols, key=lambda t: t.name)
    assert len(free) == 2
    # every solution is a*I + b*R for some a,b
    a, b = sp.symbols('a b')
    Xsol = X.subs({p: solset[0], q: solset[1], r: solset[2], s: solset[3]})
    # match against a*I+b*R
    target = a * I2 + b * R
    match = sp.solve([sp.Eq(Xsol[i], target[i]) for i in range(4)], [a, b], dict=True)
    assert match, "commuting matrix is not in span{I,R}"


# ----------------------------------------------------------------------------
# Thm 2.6 -- the sl2 triple
# ----------------------------------------------------------------------------
def test_H_explicit():
    """Sec 2.3: 'H := 2R - I = [[-1,2],[2,1]]'."""
    assert mat_eq(H, sp.Matrix([[-1, 2], [2, 1]]))


def test_sl2_triple_brackets():
    """Thm 2.6: '[H,N_+]=2sqrt5 N_+', '[H,N_-]=-2sqrt5 N_-', '[N_+,N_-]=sqrt5 H'."""
    assert mat_eq(comm(H, Np), 2 * sqrt5 * Np)
    assert mat_eq(comm(H, Nm), -2 * sqrt5 * Nm)
    assert mat_eq(comm(Np, Nm), sqrt5 * H)


def test_sl2_normalised_standard_form():
    """Thm 2.6: h=H/sqrt5, e=N_+, f=N_-/5 give '[h,e]=2e,[h,f]=-2f,[e,f]=h'."""
    h = H / sqrt5
    e = Np
    f = Nm / 5
    assert mat_eq(comm(h, e), 2 * e)
    assert mat_eq(comm(h, f), -2 * f)
    assert mat_eq(comm(e, f), h)


def test_all_triple_elements_traceless():
    """Thm 2.6 setup: H,N_+,N_- are traceless (live in sl2)."""
    for M in (H, Np, Nm):
        assert is_zero(sp.trace(M))


# ----------------------------------------------------------------------------
# Prop 2.7 -- Lorentzian signature (2,1)
# ----------------------------------------------------------------------------
def test_squares_and_signature():
    """Prop 2.7: 'H^2 = S^2 = 5I', 'J^2 = -I', 'N_+^2 = N_-^2 = 0'."""
    assert mat_eq(H * H, 5 * I2)
    assert mat_eq(S * S, 5 * I2)
    assert mat_eq(J * J, -I2)
    assert mat_eq(Np * Np, sp.zeros(2))
    assert mat_eq(Nm * Nm, sp.zeros(2))


def test_S_equals_Nplus_plus_Nminus():
    """Prop 2.7 / Thm 2.8: 'S = N_+ + N_-'."""
    assert mat_eq(S, Np + Nm)


def test_signature_from_eigenvalues():
    """Prop 2.7: matrix squaring to +5I has eigenvalues +/-sqrt5 (hyperbolic);
    squaring to -I has eigenvalues +/-i (elliptic): signature (2,1)."""
    assert set(H.eigenvals().keys()) == {sqrt5, -sqrt5}
    assert set(S.eigenvals().keys()) == {sqrt5, -sqrt5}
    assert set(J.eigenvals().keys()) == {sp.I, -sp.I}


# ----------------------------------------------------------------------------
# Prop 2.9 -- rational structure constants & splitting obstruction
# ----------------------------------------------------------------------------
def test_rational_structure_constants():
    """Prop 2.9 / Ex 2.11: '[H,S]=10J', '[H,J]=2S', '[S,J]=-2H'."""
    assert mat_eq(comm(H, S), 10 * J)
    assert mat_eq(comm(H, J), 2 * S)
    assert mat_eq(comm(S, J), -2 * H)


def test_HS_and_SH_products():
    """Prop 2.9 proof: 'HS=[[0,-5],[5,0]]', 'SH=[[0,5],[-5,0]]'."""
    assert mat_eq(H * S, sp.Matrix([[0, -5], [5, 0]]))
    assert mat_eq(S * H, sp.Matrix([[0, 5], [-5, 0]]))


def test_adH_on_SJ_plane_char_poly():
    """Prop 2.9: ad_H on span{S,J} is [[0,2],[10,0]], char poly lambda^2-20,
    roots +/-2sqrt5 not in Q."""
    # columns are coords of [H,S]=10J and [H,J]=2S in basis (S,J)
    adH = sp.Matrix([[0, 2], [10, 0]])
    lam = sp.symbols('lambda')
    cp = sp.expand((lam * sp.eye(2) - adH).det())
    assert is_zero(cp - (lam**2 - 20))
    roots = sp.solve(cp, lam)
    assert set(roots) == {2 * sqrt5, -2 * sqrt5}
    assert sp.Rational(20).is_integer and not sp.sqrt(20).is_rational


# ----------------------------------------------------------------------------
# Thm 2.10 -- the null-frame transition and the conjugator V
# ----------------------------------------------------------------------------
def test_null_frame_change():
    """Thm 2.10: 'N_+ = (S+sqrt5 J)/2', 'N_- = (S-sqrt5 J)/2', and inversely
    'S=N_++N_-', 'sqrt5 J = N_+ - N_-'."""
    assert mat_eq(Np, (S + sqrt5 * J) / 2)
    assert mat_eq(Nm, (S - sqrt5 * J) / 2)
    assert mat_eq(S, Np + Nm)
    assert mat_eq(sqrt5 * J, Np - Nm)


def test_transition_matrices_inverse_and_det():
    """Thm 2.10: M_{R->E} and M_{E->R} are mutually inverse and
    det(M_{R->E}) = -2/sqrt5 (not rational)."""
    MRE = sp.Matrix([[1, 0, 0],
                     [0, 1, 1 / sqrt5],
                     [0, 1, -1 / sqrt5]])
    MER = sp.Matrix([[1, 0, 0],
                     [0, sp.Rational(1, 2), sp.Rational(1, 2)],
                     [0, sqrt5 / 2, -sqrt5 / 2]])
    assert mat_eq(MRE * MER, sp.eye(3))
    assert is_zero(MRE.det() - (-2 / sqrt5))
    assert not sp.nsimplify(MRE.det()).is_rational


def test_transition_matrices_realise_basis_change():
    """Thm 2.10: M_{R->E} sends R-coords (a,b,c) of aH+bS+cJ to the
    E-coords (p,q,r) of pH+qN_+ +rN_-.  Verified on the actual matrices."""
    a, b, c = sp.symbols('a b c')
    MRE = sp.Matrix([[1, 0, 0],
                     [0, 1, 1 / sqrt5],
                     [0, 1, -1 / sqrt5]])
    p, q, r = MRE * sp.Matrix([a, b, c])
    lhs = a * H + b * S + c * J
    rhs = p * H + q * Np + r * Nm
    assert mat_eq(sp.expand(lhs - rhs), sp.zeros(2))


def test_conjugator_V():
    """Thm 2.10: V=[[1,1],[phi,psi]] diagonalises H and triangularises N_+,N_-:
    'V^{-1}HV=diag(sqrt5,-sqrt5)', 'V^{-1}N_+V=[[0,(5-sqrt5)/2],[0,0]]',
    'V^{-1}N_-V=[[0,0],[(5+sqrt5)/2,0]]'."""
    V = sp.Matrix([[1, 1], [phi, psi]])
    Vi = V.inv()
    assert mat_eq(Vi * H * V, sp.diag(sqrt5, -sqrt5))
    assert mat_eq(Vi * Np * V, sp.Matrix([[0, (5 - sqrt5) / 2], [0, 0]]))
    assert mat_eq(Vi * Nm * V, sp.Matrix([[0, 0], [(5 + sqrt5) / 2, 0]]))


# ----------------------------------------------------------------------------
# Ex 2.12 -- the triple by hand
# ----------------------------------------------------------------------------
def test_example_triple_products():
    """Ex 2.12: 'N_+ N_- = [[1+psi^2, phi-psi],[phi-psi, phi^2+1]]' and
    'N_- N_+ = [[1+phi^2, psi-phi],[psi-phi, psi^2+1]]', giving [N_+,N_-]=sqrt5 H."""
    assert mat_eq(Np * Nm, sp.Matrix([[1 + psi**2, phi - psi],
                                      [phi - psi, phi**2 + 1]]))
    assert mat_eq(Nm * Np, sp.Matrix([[1 + phi**2, psi - phi],
                                      [psi - phi, psi**2 + 1]]))
    assert mat_eq(comm(Np, Nm), sqrt5 * H)
    # phi^2 - psi^2 = sqrt5 and 2(phi-psi)=2sqrt5
    assert is_zero((phi**2 - psi**2) - sqrt5)
    assert is_zero(2 * (phi - psi) - 2 * sqrt5)
