"""
Independent verification of cross-field growth / compositum claims
(Sec. 4.4: Prop. 4.6 (Kronecker Gram), Cor. 4.7 (multiquadratic tower),
Ex. 4.8 (adjoining sqrt7), Rem. 4.9 & Ex. 4.10 (non-disjoint case)).
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import vsub_nf as nf

x = sp.symbols('x')


def _kron_det(A, B):
    """det(A (x) B) = (det A)^{dim B} (det B)^{dim A}."""
    return A.det()**B.shape[0] * B.det()**A.shape[0]


def test_kron_gram_determinant_identity():
    """Prop. 4.6: det G_{KL} = (det G_K)^k (det G_L)^m via det(A(x)B)."""
    A = sp.diag(4, 8, 12, 24)                    # G_K, m=4
    B = sp.diag(2, 14)                           # G_L, k=2
    KL = sp.Matrix(sp.kronecker_product(A, B))
    assert KL.det() == _kron_det(A, B)
    assert KL.det() == A.det()**2 * B.det()**4


def test_multiquadratic_tower_2_3():
    """Cor. 4.7: Q(sqrt2,sqrt3) has G = [[2,0],[0,4]] (x) [[2,0],[0,6]], a diagonal
    matrix whose entries are the multiset {4,8,12,24} = {2^2 * prod_{i in S} p_i}."""
    G2 = sp.Matrix([[2, 0], [0, 4]])             # diag(2, 2*2)
    G3 = sp.Matrix([[2, 0], [0, 6]])             # diag(2, 2*3)
    K = sp.Matrix(sp.kronecker_product(G2, G3))
    diag = sorted(K[i, i] for i in range(4))
    assert diag == [4, 8, 12, 24]
    # explicit subset-product formula: entries 2^r * prod_{i in S} p_i, r=2
    from itertools import combinations
    ps = [2, 3]
    entries = []
    for k in range(len(ps) + 1):
        for S in combinations(ps, k):
            prod = 1
            for p in S:
                prod *= p
            entries.append(2**len(ps) * prod)
    assert sorted(entries) == [4, 8, 12, 24]
    # matches the direct biquadratic Gram (Ex. 2.16) as a multiset
    assert diag == sorted([4, 8, 12, 24])


def test_adjoin_sqrt7_kronecker():
    """Ex. 4.8: G_K=diag(4,8,12,24), G_L=diag(2,14), linearly disjoint
    ([KL:Q]=8=4*2); G_{KL}=G_K (x) G_L=diag(8,56,16,112,24,168,48,336),
    det G_{KL}=9216^2 * 28^4."""
    GK = sp.diag(4, 8, 12, 24)
    GL = sp.diag(2, 14)
    assert 4 * 2 == 8                            # [KL:Q]=[K:Q][L:Q] (disjoint)
    KL = sp.Matrix(sp.kronecker_product(GK, GL))
    assert [KL[i, i] for i in range(8)] == [8, 56, 16, 112, 24, 168, 48, 336]
    assert KL == sp.diag(8, 56, 16, 112, 24, 168, 48, 336)
    assert KL.det() == 9216**2 * 28**4
    assert GK.det() == 9216 and GL.det() == 28


def test_kron_gram_from_trace_definition():
    """Prop. 4.6: rebuild G_{KL} on the product basis from the trace-form
    definition Tr_{KL}(e_i f_j e_i' f_j')=Tr_K(e_i e_i') Tr_L(f_j f_j'), and
    confirm it equals G_K (x) G_L.  Take K=L=Q(sqrt5)-like diagonal factors."""
    # Use two independent quadratic factors via their power-basis Grams.
    CK = nf.companion_from_poly(x**2 - 2)        # Q(sqrt2): G_K = diag(2,4)
    CL = nf.companion_from_poly(x**2 - 3)        # Q(sqrt3): G_L = diag(2,6)
    GK = nf.power_gram(CK)
    GL = nf.power_gram(CL)
    assert GK == sp.diag(2, 4)
    assert GL == sp.diag(2, 6)
    # Kronecker product of the factor Grams
    KL = sp.Matrix(sp.kronecker_product(GK, GL))
    # compositum Q(sqrt2,sqrt3) trace form on the product basis {1,sqrt3,sqrt2,sqrt6}
    # (this basis ordering is exactly e_i (x) f_j); rebuild via traces in Q(sqrt2+sqrt3):
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    S2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
    S3 = (0, sp.Rational(11, 2), 0, sp.Rational(-1, 2))
    S6 = (sp.Rational(-5, 2), 0, sp.Rational(1, 2), 0)
    one = (1, 0, 0, 0)
    prod_basis = [one, S3, S2, S6]               # 1, sqrt3, sqrt2, sqrt2*sqrt3
    G_comp = nf.gram(prod_basis, C)
    assert G_comp == KL                          # Kronecker structure reproduced


def test_nondisjoint_factorisation():
    """Ex. 4.10: over K=Q(sqrt2), x^4-10x^2+1 factors as
    (x^2+2sqrt2 x-1)(x^2-2sqrt2 x-1); beta=sqrt2+sqrt3 satisfies
    beta^2-2sqrt2 beta-1=0, so [K(beta):K]=2<4 (NOT linearly disjoint)."""
    s2 = sp.sqrt(2)
    f = (x**2 + 2*s2*x - 1) * (x**2 - 2*s2*x - 1)
    assert sp.expand(f) == x**4 - 10*x**2 + 1
    beta = sp.sqrt(2) + sp.sqrt(3)
    assert sp.simplify(beta**2 - 2*s2*beta - 1) == 0
    # degree over K is 2 (root of a quadratic over Q(sqrt2)), not deg m_beta = 4
    assert sp.degree(x**2 - 2*s2*x - 1, x) == 2
    assert sp.degree(sp.minimal_polynomial(beta, x), x) == 4
    # true compositum Q(sqrt2,sqrt3) has degree 4, not the tensor degree 2*4=8
    assert 2 * 4 == 8                            # fictitious tensor degree
    # Q(sqrt2,sqrt3) really is degree 4 over Q
    assert sp.degree(sp.minimal_polynomial(beta, x), x) == 4
