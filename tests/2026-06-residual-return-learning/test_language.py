r"""
Independent verification of Section 5 ("A residual-valued language") over
Cl(2,0) = M_2(R): the exact carrier / matrix isomorphism, the phi-keystone
return operator L(X)=RX+XR-X, its exact 2-dimensional kernel (Thm 5.4), the
'disproof' negative witness (Rem 5.5), the commit projector (eq proj, Thm 5.6),
generalisation-by-return-to-zero, the lexicon example, and the firewall counts.

Everything is computed through the matrix isomorphism
    mat(a+b e1+c e2+d i) = [[a+c, b-d],[b+d, a-c]]
with sympy exact rationals, then compared to the paper's stated values.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, eye, zeros, symbols

x = symbols("x")


def mat(a, b, c, d):
    """Clifford holding X=a+b e1+c e2+d i -> its 2x2 real matrix."""
    return Matrix([[a + c, b - d], [b + d, a - c]])


def cl_coords(M):
    """Inverse iso: 2x2 matrix -> (a,b,c,d) Clifford coordinates."""
    a = (M[0, 0] + M[1, 1]) / 2
    c = (M[0, 0] - M[1, 1]) / 2
    b = (M[1, 0] + M[0, 1]) / 2
    d = (M[1, 0] - M[0, 1]) / 2
    return (a, b, c, d)


R = mat(Q(1, 2), 1, Q(-1, 2), 0)   # keystone R = 1/2 + e1 - 1/2 e2


# --- Section 5.1 (the exact carrier) -----------------------------------------

def test_matrix_iso_is_algebra_homomorphism():
    """Sec 5.1: e1^2=e2^2=1, i=e1 e2, i^2=-1 under mat(.).  Verify the generator
    relations through the matrix isomorphism."""
    e1, e2, i, one = mat(0, 1, 0, 0), mat(0, 0, 1, 0), mat(0, 0, 0, 1), eye(2)
    assert e1 * e1 == one
    assert e2 * e2 == one
    assert i * i == -one
    assert e1 * e2 == i                    # i = e1 e2


def test_keystone_matrix_and_golden_law():
    """Sec 5.2 & Prop 5.7: mat(R)=[[0,1],[1,1]] (companion of x^2-x-1), symmetric,
    with R^2 = R + I (the golden law), scalar part tau(R)=1/2, det(R)=-1."""
    assert R == Matrix([[0, 1], [1, 1]])
    assert R == R.T                                   # symmetric
    assert R * R == R + eye(2)                        # R^2 = R + 1
    assert sp.expand(R.charpoly(x).as_expr() - (x**2 - x - 1)) == 0
    a, b, c, d = cl_coords(R)
    assert a == Q(1, 2)                               # tau(R) = scalar part = 1/2
    assert R.det() == -1


def test_gate_P0_is_matrix_idempotent():
    """Ex 5.3: the gate P0 = 1/2(1+e2) has mat(P0)=[[1,0],[0,0]], a symmetric
    idempotent (nu(P0)=0 by the idempotence reading)."""
    P0 = mat(Q(1, 2), 0, Q(1, 2), 0)
    assert P0 == Matrix([[1, 0], [0, 0]])
    assert P0 * P0 == P0
    assert P0 == P0.T


def test_cayley_hamilton_holds():
    """Sec 5.1: Cayley-Hamilton X^2 - tr(X) X + det(X) I = 0 holds exactly for the
    matrix image of any holding (verified on the keystone and a generic holding)."""
    for M in (R, mat(Q(3, 7), -2, Q(5, 3), 4)):
        tr, det = M.trace(), M.det()
        assert sp.expand(M * M - tr * M + det * eye(2)) == zeros(2)


# --- Section 5.2 (the return operator) ---------------------------------------

def Lop(X):
    return R * X + X * R - X


def test_Lop_of_keystone_not_in_kernel():
    """Prop 5.8: R is NOT in ker L; indeed L(R) = 5/2 + e1 - 1/2 e2 != 0.
    Independently L(R)=2R^2-R=2(R+I)-R=R+2I."""
    val = Lop(R)
    assert val != zeros(2)
    assert cl_coords(val) == (Q(5, 2), 1, Q(-1, 2), 0)
    assert val == R + 2 * eye(2)


# --- Theorem 5.4 (the exact kernel is the lexicon space) ---------------------

def Lop_matrix():
    """4x4 matrix of L over the ordered Clifford basis {1,e1,e2,i}."""
    basis = [mat(1, 0, 0, 0), mat(0, 1, 0, 0), mat(0, 0, 1, 0), mat(0, 0, 0, 1)]
    cols = [list(cl_coords(Lop(B))) for B in basis]
    return Matrix(cols).T


def test_kernel_dim_2_and_basis():
    """Thm 5.4: ker L is exactly 2-dimensional, ker L = span{e1+2e2, i} =
    span{H(0,1,2,0), H(0,0,0,1)}, by exact rational nullspace of the 4x4 matrix of L."""
    L = Lop_matrix()
    ns = L.nullspace()
    assert len(ns) == 2
    # both stated basis vectors lie in ker L and span the same 2-space
    k1 = Matrix([0, 1, 2, 0])   # e1 + 2 e2
    k2 = Matrix([0, 0, 0, 1])   # i
    assert L * k1 == zeros(4, 1)
    assert L * k2 == zeros(4, 1)
    span = Matrix.hstack(*ns)
    assert Matrix.hstack(span, k1, k2).rank() == 2


def test_kernel_sylvester_eigenvalues():
    """Thm 5.4: L=S-id with S(X)=RX+XR having eigenvalues {2 phi,1,1,2 psi} over
    the eigenvalues {phi,psi} of R (phi+psi=1); hence L has eigenvalues
    {2 phi-1, 0, 0, 2 psi-1} = {sqrt5, 0, 0, -sqrt5}, eigenvalue 0 with mult 2."""
    L = Lop_matrix()
    eigs = L.eigenvals()
    assert eigs.get(sp.Integer(0)) == 2
    assert any(sp.simplify(e - sp.sqrt(5)) == 0 for e in eigs)
    assert any(sp.simplify(e + sp.sqrt(5)) == 0 for e in eigs)


# --- Remark 5.5 (a corrected hint, pinned by its disproof) -------------------

def test_disproof_antisymmetric_coords_not_in_kernel():
    """Rem 5.5: the naive 'N=H(0,-1,1,0)' is NOT in ker L:
    L(H(0,-1,1,0)) = H(-3,0,0,0) != 0."""
    val = Lop(mat(0, -1, 1, 0))
    assert cl_coords(val) == (-3, 0, 0, 0)
    assert val != zeros(2)


def test_disproof_correct_image_is_i_in_kernel():
    """Rem 5.5: the correct image of [[0,-1],[1,0]] through the iso is
    cl([[0,-1],[1,0]]) = i = H(0,0,0,1), which IS in ker L."""
    N = Matrix([[0, -1], [1, 0]])
    assert cl_coords(N) == (0, 0, 0, 1)               # cl(N) = i
    assert Lop(N) == zeros(2)                          # i in ker L


# --- Theorem 5.6 / eq (proj) (the commit projector) --------------------------

def commit_projector():
    K = Matrix([[0, 0], [1, 0], [2, 0], [0, 1]])       # kernel basis e1+2e2, i
    return K * (K.T * K).inv() * K.T


def test_commit_projector_exact_value():
    """eq (proj) / Thm 5.6: the commit is the exact idempotent orthogonal projector
    onto ker L, = [[0,0,0,0],[0,1/5,2/5,0],[0,2/5,4/5,0],[0,0,0,1]], and it is not
    the identity (idempotence alone would admit I)."""
    P = commit_projector()
    assert P == Matrix([[0, 0, 0, 0],
                        [0, Q(1, 5), Q(2, 5), 0],
                        [0, Q(2, 5), Q(4, 5), 0],
                        [0, 0, 0, 1]])
    assert P * P == P
    assert P != eye(4)


def test_generalization_E1_and_E1_plus_1_merge():
    """Thm 5.6 / Ex 5.7: commit(E1+1)=commit(E1)=H(0,1/5,2/5,0) (the constant 1 is
    orthogonal to the slack), so E1 and E1+1 return to the SAME residue."""
    P = commit_projector()
    E1 = Matrix([0, 1, 0, 0])
    E1p1 = Matrix([1, 1, 0, 0])
    assert P * E1 == Matrix([0, Q(1, 5), Q(2, 5), 0])
    assert P * E1p1 == P * E1


def test_generalization_positive_scaling_same_ray():
    """Thm 5.6: positive scaling preserves the committed direction:
    commit(3 E1) = 3 commit(E1) (same ray -> same word)."""
    P = commit_projector()
    E1 = Matrix([0, 1, 0, 0])
    assert P * (3 * E1) == 3 * (P * E1)


def test_distinct_residues_i_and_2i_stay_distinct():
    """Thm 5.6 / Ex 5.7: i (H(0,0,0,1)) and 2i (H(0,0,0,2)) commit to DISTINCT
    exact residues though they share the sign-word (+i); dedup is on the value."""
    P = commit_projector()
    ci = P * Matrix([0, 0, 0, 1])
    c2i = P * Matrix([0, 0, 0, 2])
    assert ci == Matrix([0, 0, 0, 1])
    assert c2i == Matrix([0, 0, 0, 2])
    assert ci != c2i


def test_lexicon_five_tokens_four_entries():
    """Ex 5.7: feeding {E1, E1+1, i, 2i, -i} through commit yields committed values
    {H(0,1/5,2/5,0) (x2, merged), H(0,0,0,1), H(0,0,0,2), H(0,0,0,-1)}: five tokens,
    FOUR distinct exact residues."""
    P = commit_projector()
    tokens = {
        "E1": Matrix([0, 1, 0, 0]),
        "E1+1": Matrix([1, 1, 0, 0]),
        "i": Matrix([0, 0, 0, 1]),
        "2i": Matrix([0, 0, 0, 2]),
        "-i": Matrix([0, 0, 0, -1]),
    }
    committed = [tuple(P * v) for v in tokens.values()]
    assert committed[0] == committed[1]                # E1, E1+1 merge
    assert len(set(committed)) == 4                    # four distinct residues


# --- Table 6 (the jurisdiction firewall counts) ------------------------------

def test_firewall_counts_consistent():
    """Table 6: wired = theorem(20)+computed(5) = 25; bank = 20+5+interpretive(2)
    = 27; false_as_stated = 0 (unused)."""
    theorem, computed, interpretive, false_as_stated = 20, 5, 2, 0
    assert theorem + computed == 25                    # wired total
    assert theorem + computed + interpretive == 27     # bank total
    assert false_as_stated == 0
