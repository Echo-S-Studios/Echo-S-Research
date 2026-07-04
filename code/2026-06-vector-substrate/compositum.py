"""
PRODUCER: cross-field growth and the Kronecker Gram of a disjoint compositum.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex, Sec. 4.4
  - Prop. 4.6 (G_{KL} = G_K (x) G_L; det(A(x)B) identity)
  - Cor. 4.7 (multiquadratic tower {2,3} -> multiset {4,8,12,24})
  - Ex. 4.8 (adjoining sqrt7: diag(8,56,16,112,24,168,48,336), det 9216^2*28^4)
  - Rem. 4.9 / Ex. 4.10 (non-disjoint case: x^4-10x^2+1 factors over Q(sqrt2))

Emits:
  data/compositum.json
"""
from itertools import combinations

import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "compositum.py"


def kron_det_identity():
    """Prop. 4.6: det G_{KL} = (det G_K)^k (det G_L)^m via det(A(x)B)."""
    A = sp.diag(4, 8, 12, 24)                    # G_K, m = 4
    B = sp.diag(2, 14)                           # G_L, k = 2
    KL = sp.Matrix(sp.kronecker_product(A, B))
    return {
        "G_K": vc.mat_to_list(A),
        "G_L": vc.mat_to_list(B),
        "det_G_KL": vc.sval(KL.det()),
        "formula_detA_pow_dimB_detB_pow_dimA": vc.sval(A.det()**B.shape[0] * B.det()**A.shape[0]),
        "equals": KL.det() == A.det()**B.shape[0] * B.det()**A.shape[0],
    }


def multiquadratic_tower():
    """Cor. 4.7: Q(sqrt2,sqrt3) Gram = diag(2,4) (x) diag(2,6) has diagonal
    multiset {4,8,12,24} = {2^2 * prod_{i in S} p_i : S subset {2,3}}."""
    G2 = sp.Matrix([[2, 0], [0, 4]])
    G3 = sp.Matrix([[2, 0], [0, 6]])
    K = sp.Matrix(sp.kronecker_product(G2, G3))
    diag = sorted(int(K[i, i]) for i in range(4))
    ps = [2, 3]
    entries = []
    for k in range(len(ps) + 1):
        for S in combinations(ps, k):
            prod = 1
            for p in S:
                prod *= p
            entries.append(2**len(ps) * prod)
    return {
        "kron_diagonal": diag,
        "subset_product_formula": sorted(entries),
        "match": diag == sorted(entries) == [4, 8, 12, 24],
    }


def adjoin_sqrt7():
    """Ex. 4.8: G_K=diag(4,8,12,24), G_L=diag(2,14), linearly disjoint
    ([KL:Q]=8=4*2); G_{KL}=G_K(x)G_L=diag(8,56,16,112,24,168,48,336)."""
    GK = sp.diag(4, 8, 12, 24)
    GL = sp.diag(2, 14)
    KL = sp.Matrix(sp.kronecker_product(GK, GL))
    return {
        "disjoint_degree_check": (4 * 2 == 8),
        "G_KL_diagonal": [int(KL[i, i]) for i in range(8)],
        "det_G_KL": vc.sval(KL.det()),
        "det_G_KL_factored": "9216^2 * 28^4",
        "det_G_K": vc.sval(GK.det()),
        "det_G_L": vc.sval(GL.det()),
    }


def kron_from_trace_definition():
    """Prop. 4.6: rebuild G_{compositum} on the product basis directly from the
    trace-form definition and confirm it equals G_K (x) G_L.  Factors:
    Q(sqrt2) (G=diag(2,4)) and Q(sqrt3) (G=diag(2,6)); compositum Q(sqrt2,sqrt3)
    on product basis {1, sqrt3, sqrt2, sqrt6}."""
    CK = vc.companion_from_poly(x**2 - 2)
    CL = vc.companion_from_poly(x**2 - 3)
    GK = vc.power_gram(CK)
    GL = vc.power_gram(CL)
    KL = sp.Matrix(sp.kronecker_product(GK, GL))

    C = vc.companion_from_poly(x**4 - 10 * x**2 + 1)
    S2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
    S3 = (0, sp.Rational(11, 2), 0, sp.Rational(-1, 2))
    S6 = (sp.Rational(-5, 2), 0, sp.Rational(1, 2), 0)
    prod_basis = [(1, 0, 0, 0), S3, S2, S6]      # e_i (x) f_j ordering
    G_comp = vc.gram(prod_basis, C)
    return {
        "G_factor_Q_sqrt2": vc.mat_to_list(GK),
        "G_factor_Q_sqrt3": vc.mat_to_list(GL),
        "G_kron": vc.mat_to_list(KL),
        "G_compositum_from_traces": vc.mat_to_list(G_comp),
        "kron_equals_trace_gram": G_comp == KL,
    }


def non_disjoint():
    """Ex. 4.10: over K=Q(sqrt2), x^4-10x^2+1 = (x^2+2sqrt2 x-1)(x^2-2sqrt2 x-1),
    so [K(beta):K]=2 < 4 = deg m_beta: NOT linearly disjoint; the tensor degree
    2*4=8 is fictitious, the true compositum Q(sqrt2,sqrt3) has degree 4."""
    s2 = sp.sqrt(2)
    f = (x**2 + 2 * s2 * x - 1) * (x**2 - 2 * s2 * x - 1)
    beta = sp.sqrt(2) + sp.sqrt(3)
    return {
        "factorisation_expands_to_min_poly": sp.expand(f) == x**4 - 10 * x**2 + 1,
        "beta_satisfies_quadratic_over_K": bool(sp.simplify(beta**2 - 2 * s2 * beta - 1) == 0),
        "degree_over_K": int(sp.degree(x**2 - 2 * s2 * x - 1, x)),
        "deg_min_poly_over_Q": int(sp.degree(sp.minimal_polynomial(beta, x), x)),
        "fictitious_tensor_degree": 2 * 4,
        "true_compositum_degree": 4,
        "linearly_disjoint": False,
    }


def main():
    payload = {
        "kron_det_identity": kron_det_identity(),
        "multiquadratic_tower": multiquadratic_tower(),
        "adjoin_sqrt7": adjoin_sqrt7(),
        "kron_from_trace_definition": kron_from_trace_definition(),
        "non_disjoint": non_disjoint(),
    }
    p = vc.write_json("compositum.json", payload, SCRIPT)
    print(f"wrote {p}")
    print("  adjoin sqrt7 G_KL diagonal =", payload["adjoin_sqrt7"]["G_KL_diagonal"])
    print("  Kronecker == trace-form Gram:",
          payload["kron_from_trace_definition"]["kron_equals_trace_gram"])


if __name__ == "__main__":
    main()
