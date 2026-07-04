r"""
Independent verification of Section 2 ("Learning as exact basis growth") of
    "Residual Return: Exact Learning Dynamics and Language over the Vector Substrate"
    (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Every value is re-derived from the paper's stated premises (the field
K = Q(sqrt2+sqrt3) = Q[x]/(x^4-10x^2+1), power basis {1,th,th^2,th^3}) using
exact sympy/stdlib arithmetic, then compared to the paper's stated result.

The trace-form Gram is built INDEPENDENTLY from the regular representation
(companion matrix) of theta: Tr_{K/Q}(theta^k) is exactly trace(C(m_theta)^k),
so we never transcribe the paper's matrix -- we reconstruct it.
"""
import hashlib
import json

import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols, zeros, eye

x = symbols("x")


# --- independent helpers -----------------------------------------------------

def companion(coeffs):
    """Companion matrix of a monic poly given coeffs=[1,c_{n-1},...,c_1,c_0].

    Its char. polynomial is the given poly and its eigenvalues are the roots,
    so trace(C^k) = sum(root_i^k) = Tr_{K/Q}(theta^k) for a primitive theta.
    """
    n = len(coeffs) - 1
    C = zeros(n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]
    return C


def trace_form_powerbasis(min_coeffs):
    """Exact trace-form Gram G_{ij}=Tr(theta^{i-1}theta^{j-1}) in the power basis,
    built from the regular representation (companion) of theta -- independent of
    the paper's displayed matrix."""
    C = companion(min_coeffs)
    n = C.shape[0]
    # power sums p_k = Tr(theta^k) = trace(C^k)
    powers = [eye(n)]
    for _ in range(2 * n):
        powers.append(powers[-1] * C)
    p = [powers[k].trace() for k in range(2 * n)]
    return Matrix(n, n, lambda i, j: p[i + j])


M_THETA = [1, 0, -10, 0, 1]  # theta = sqrt2+sqrt3, m_theta = x^4-10x^2+1


# --- eq:Gpower ---------------------------------------------------------------

def test_gpower_matrix_from_traces():
    """Eq (Gpower): the trace-form Gram of Q(sqrt2+sqrt3) in the power basis is
    [[4,0,20,0],[0,20,0,196],[20,0,196,0],[0,196,0,1940]], G_{ij}=Tr(th^{i-1}th^{j-1})."""
    G = trace_form_powerbasis(M_THETA)
    expected = Matrix([[4, 0, 20, 0],
                       [0, 20, 0, 196],
                       [20, 0, 196, 0],
                       [0, 196, 0, 1940]])
    assert G == expected


def test_gpower_sanity_companion_charpoly():
    """The companion used to build the Gram really is the regular rep of theta:
    its characteristic polynomial is m_theta = x^4-10x^2+1."""
    C = companion(M_THETA)
    assert sp.expand(C.charpoly(x).as_expr() - (x**4 - 10 * x**2 + 1)) == 0


def test_gpower_positive_definite():
    """Sec 2: 'positive definite (the field is totally real)'. Check all leading
    principal minors are > 0 (Sylvester's criterion) and all eigenvalues > 0."""
    G = trace_form_powerbasis(M_THETA)
    for k in range(1, 5):
        assert G[:k, :k].det() > 0
    assert all(ev > 0 for ev in G.eigenvals())


# --- exact projector P = B (B^T G B)^{-1} B^T G ------------------------------

def projector(B, G):
    return B * (B.T * G * B).inv() * B.T * G


def test_projector_idempotent_exact():
    """Prop 2.1 / Sec 2: P=B(B^T G B)^{-1}B^T G is an exact rational projector,
    P^2 = P exactly (B = {1,theta})."""
    G = trace_form_powerbasis(M_THETA)
    B = Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])  # columns 1, theta
    P = projector(B, G)
    assert P * P == P
    # projector fixes its column space
    assert P * B == B


# --- Example 2.4 (a worked episode): w = 2 sqrt6 = theta^2 - 5 ----------------

def test_episode_w_equals_theta2_minus_5():
    """Ex 2.4: 'the off-axis quantity w = 2 sqrt6 = theta^2 - 5', theta=sqrt2+sqrt3.
    Check symbolically that (sqrt2+sqrt3)^2 - 5 = 2 sqrt6."""
    th = sqrt(2) + sqrt(3)
    assert sp.simplify(th**2 - 5 - 2 * sqrt(6)) == 0


def test_episode_coords_and_orthogonality():
    """Ex 2.4: w has coordinate vector (-5,0,1,0) and is G-orthogonal to
    col(B)=span{1,theta}."""
    G = trace_form_powerbasis(M_THETA)
    w = Matrix([-5, 0, 1, 0])  # -5*1 + 1*theta^2 = theta^2-5
    for col in (Matrix([1, 0, 0, 0]), Matrix([0, 1, 0, 0])):
        assert (col.T * G * w)[0] == 0


def test_episode_residual_is_w_and_norm_96():
    """Ex 2.4: for x = 1 + w the residual r = w exactly and
    ||r||_G^2 = Tr((2 sqrt6)^2) = Tr(24) = 4*24 = 96."""
    G = trace_form_powerbasis(M_THETA)
    B = Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])
    P = projector(B, G)
    xvec = Matrix([-4, 0, 1, 0])  # x = 1 + w = theta^2 - 4
    r = xvec - P * xvec
    assert r == Matrix([-5, 0, 1, 0])  # residual equals w
    score = (r.T * G * r)[0]
    assert score == 96
    # independent trace route: Tr_{K/Q}(24) = trace of regular rep (24*I_4) = 96
    assert (24 * eye(4)).trace() == 96


def test_episode_seed_minpoly_x2_minus_24():
    """Ex 2.4 / Ex 2.9: the seed 2 sqrt6 has exact minimal polynomial x^2 - 24."""
    w = 2 * sqrt(6)
    mp = sp.minimal_polynomial(w, x)
    assert sp.expand(mp - (x**2 - 24)) == 0


# --- Example 2.6 (a hand-checkable witness link): SHA-256 digest --------------

def test_witness_digest_exact_31f1f1e05ac9a35a():
    """Ex 2.6 & Sec 8: the genesis witness body, canonically JSON-encoded (sorted
    keys, no whitespace), hashed as SHA-256('genesis'+json)[0:16], equals
    '31f1f1e05ac9a35a'. Recomputed from scratch from the FULL body given in the paper."""
    body = {
        "coords": ["-5", "0", "1", "0"],
        "event": "basis_growth",
        "index": 0,
        "min_poly": [1, 0, -24],
        "num_seeds": 3,
        "prev_hash": "genesis",
        "snap": "exact",
        "streak": 4,
    }
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"))
    # the paper displays the exact canonical string; confirm our encoder matches it
    assert canonical == (
        '{"coords":["-5","0","1","0"],"event":"basis_growth","index":0,'
        '"min_poly":[1,0,-24],"num_seeds":3,"prev_hash":"genesis",'
        '"snap":"exact","streak":4}'
    )
    digest = hashlib.sha256(("genesis" + canonical).encode()).hexdigest()[:16]
    assert digest == "31f1f1e05ac9a35a"
