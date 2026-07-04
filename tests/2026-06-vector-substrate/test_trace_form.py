"""
Independent verification of the trace-form / Minkowski-embedding claims
(Sec. 2.4-2.6, 2.7: Thm. 2.14 (G=M^T M, det G=d_K), Ex. 2.15 golden,
Ex. 2.16 biquadratic, Ex. 2.18 Q(i), Ex. 2.19 Q(cbrt2), Table 2 catalog,
Ex. 3.9 trace duality / different).

Gram matrices are rebuilt from traces of products in the regular representation;
G=M^T M (resp. G_2=M^* M) is cross-checked against numeric Minkowski embeddings.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import sympy as sp
import mpmath as mp
import vsub_nf as nf

x = sp.symbols('x')

# theta = sqrt2 + sqrt3 ; coordinate vectors of sqrt2, sqrt3, sqrt6 in the
# power basis {1, theta, theta^2, theta^3} (derived, not restated):
#   sqrt2 = (theta^3 - 9 theta)/2, sqrt3 = (11 theta - theta^3)/2, sqrt6=(theta^2-5)/2
S2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
S3 = (0, sp.Rational(11, 2), 0, sp.Rational(-1, 2))
S6 = (sp.Rational(-5, 2), 0, sp.Rational(1, 2), 0)


def _quadratic_disc(m):
    """d_K for Q(sqrt m), m squarefree: m if m=1 mod 4 else 4m."""
    return m if m % 4 == 1 else 4 * m


def test_golden_gram_and_disc():
    """Ex. 2.15: basis {1,phi}, G=[[2,1],[1,3]], det G=5=d_{Q(sqrt5)};
    Tr(1)=2, Tr(phi)=1, Tr(phi^2)=3."""
    C = nf.companion_from_poly(x**2 - x - 1)
    G = nf.gram([(1, 0), (0, 1)], C)
    assert G == sp.Matrix([[2, 1], [1, 3]])
    assert G.det() == 5 == _quadratic_disc(5)
    assert nf.field_trace((1, 0), C) == 2
    assert nf.field_trace((0, 1), C) == 1
    # Tr(phi^2): phi^2 has coords (1,1) since phi^2 = phi + 1
    assert nf.field_trace((1, 1), C) == 3


def test_G_equals_MtM_golden_symbolic():
    """Thm. 2.14: G = M^T M exactly (golden field, symbolic)."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    M = sp.Matrix([[1, phi], [1, phic]])
    G = sp.simplify(M.T * M)
    assert G == sp.Matrix([[2, 1], [1, 3]])
    assert sp.simplify(G.det() - 5) == 0


def test_G_equals_MtM_numeric_powerbasis():
    """Thm. 2.14: G=M^T M in the power basis, cross-checked numerically against
    the exact trace-form Gram for several fields."""
    mp.mp.dps = 40
    cases = [x**2 - x - 1, x**2 - 2, x**2 - 7, x**3 - 2, x**2 + 1, x**4 - 10*x**2 + 1]
    for poly in cases:
        C = nf.companion_from_poly(poly)
        G_exact = nf.power_gram(C)               # exact rational trace form
        M = nf.embedding_matrix_power(poly, dps=40)
        MtM = M.T * M                            # (un-conjugated) => trace form
        n = C.shape[0]
        for i in range(n):
            for j in range(n):
                assert abs(mp.mpf(MtM[i, j].real) - mp.mpf(int(G_exact[i, j]))) < mp.mpf("1e-25")
                assert abs(mp.mpf(MtM[i, j].imag)) < mp.mpf("1e-25")


def test_biquadratic_gram_and_index():
    """Ex. 2.16: {1,sqrt2,sqrt3,sqrt6} gives G=diag(4,8,12,24), det=9216=2^10*3^2;
    d_K=2304=2^8*3^2 (product of the three subfield discs), index 2."""
    C = nf.companion_from_poly(x**4 - 10*x**2 + 1)
    G = nf.gram([(1, 0, 0, 0), S2, S3, S6], C)
    assert G == sp.diag(4, 8, 12, 24)
    assert G.det() == 9216 == 2**10 * 3**2
    dK = 8 * 12 * 24                             # product of Q(sqrt2),Q(sqrt3),Q(sqrt6) discs
    assert dK == 2304 == 2**8 * 3**2
    index = sp.sqrt(sp.Rational(int(G.det()), dK))
    assert index == 2


def test_gaussian_field_forms():
    """Ex. 2.18: Q(i): trace form G=diag(2,-2) (signature (1,1)); Hermitian
    G_2=M^* M=2I positive definite, det G_2=4=|d_{Q(i)}|; covol=sqrt4=2."""
    C = nf.companion_from_poly(x**2 + 1)
    G = nf.power_gram(C)
    assert G == sp.diag(2, -2)
    assert sorted(G.eigenvals().keys()) == [-2, 2]          # signature (1,1)
    # Hermitian G_2 = M^* M via exact embeddings i -> +-i
    M = sp.Matrix([[1, sp.I], [1, -sp.I]])
    G2 = sp.simplify(M.conjugate().T * M)
    assert G2 == sp.diag(2, 2)
    assert G2.det() == 4
    assert sp.sqrt(G2.det()) == 2                            # covolume sqrt|d_K|
    assert sp.Rational(1, 2) * sp.sqrt(G2.det()) == 1        # geometric rescaling 2^{-r2}


def test_cubic_complex_place_disc():
    """Ex. 2.19: Q(cbrt2), basis {1,cbrt2,cbrt4}: G=[[3,0,0],[0,0,6],[0,6,0]],
    det=-108=d_K, indefinite signature (2,1) with eigenvalues 3,+-6."""
    C = nf.companion_from_poly(x**3 - 2)
    G = nf.power_gram(C)
    assert G == sp.Matrix([[3, 0, 0], [0, 0, 6], [0, 6, 0]])
    assert G.det() == -108
    eig = sorted(G.eigenvals().keys())
    assert eig == [-6, 3, 6]
    pos = sum(1 for e in eig if e > 0)
    neg = sum(1 for e in eig if e < 0)
    assert (pos, neg) == (2, 1)                              # (r1+r2, r2)
    # Hermitian companion positive definite with det |d_K|=108
    M = nf.embedding_matrix_power(x**3 - 2, dps=40)
    G2 = M.conjugate().T * M
    assert abs(mp.mpf(G2[0, 0].real) - 3) < mp.mpf("1e-25")
    assert abs(mp.det(G2).real - 108) < mp.mpf("1e-20")


def test_catalog_discriminants_and_different_norms():
    """Table 2: det G=d_K and N(different)=N(m'(theta))=|d_K| for the catalog,
    covolume=sqrt|d_K|.  Different generator = m'(theta)."""
    # (poly, generator power basis, |d_K|, covol^2=|d_K|)
    catalog = {
        "Q(sqrt5)": (x**2 - x - 1, 5),
        "Q(sqrt2)": (x**2 - 2, 8),
        "Q(sqrt3)": (x**2 - 3, 12),
        "Q(sqrt7)": (x**2 - 7, 28),
        "Q(i)":     (x**2 + 1, 4),
        "Q(cbrt2)": (x**3 - 2, 108),
    }
    for name, (poly, absdK) in catalog.items():
        C = nf.companion_from_poly(poly)
        G = nf.power_gram(C)
        assert abs(int(G.det())) == absdK, name
        # different = (m'(theta)); N(different)=|N(m'(theta))| = |disc| = |d_K|
        mprime = sp.Poly(poly, x).diff(x)
        deriv_coords = [mprime.coeff_monomial(x**i) for i in range(C.shape[0])]
        Nd = nf.field_norm(deriv_coords, C)
        assert abs(int(Nd)) == absdK, name
        # covolume squared equals |d_K|
        covol = sp.sqrt(absdK)
        assert sp.simplify(covol**2 - absdK) == 0


def test_trace_duality_golden():
    """Ex. 3.9: G^{-1}=(1/5)[[3,-1],[-1,2]]; dual basis 1^v=(3-phi)/5,
    phi^v=(2phi-1)/5; Tr(1*1^v)=1, Tr(phi*1^v)=0; different (2phi-1)=(sqrt5),
    N=5."""
    C = nf.companion_from_poly(x**2 - x - 1)
    G = nf.gram([(1, 0), (0, 1)], C)
    Ginv = G.inv()
    assert Ginv == sp.Rational(1, 5) * sp.Matrix([[3, -1], [-1, 2]])
    # dual basis coordinates omega_i^v = sum_j (G^{-1})_{ij} omega_j (in {1,phi})
    one_dual = [Ginv[0, 0], Ginv[0, 1]]          # (3/5, -1/5) -> (3-phi)/5
    phi_dual = [Ginv[1, 0], Ginv[1, 1]]          # (-1/5, 2/5) -> (2phi-1)/5
    assert one_dual == [sp.Rational(3, 5), sp.Rational(-1, 5)]
    assert phi_dual == [sp.Rational(-1, 5), sp.Rational(2, 5)]
    # Tr(1 * 1^v) = 1, Tr(phi * 1^v) = 0  (dual-basis defining property)
    # 1 * 1^v has coords one_dual ; phi * 1^v = phi*(3-phi)/5, coords via rho
    assert nf.field_trace(one_dual, C) == 1
    prod = nf.rho((0, 1), C) * sp.Matrix(one_dual)   # coords of phi * 1^v
    assert nf.field_trace(list(prod), C) == 0
    # different (2phi-1) = (sqrt5): N(2phi-1) = -5, |.|=5=|d_K|
    assert nf.field_norm((-1, 2), C) == -5
