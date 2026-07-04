r"""
Producer: Section 3 "Cross-field growth" of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits:
  * Ex 3.5 (the bridge, exercised): phi->x^2-x-1, 2sqrt6->x^2-24, 3->x-3,
      theta->x^4-10x^2+1, and 1/2 REJECTED (not an algebraic integer)
      -> data/bridge_minpolys.csv
  * Rem 3.7 (charpoly is an incomplete similarity invariant): the phi(+)phi vs
      C((x^2-x-1)^2) witness and the Jordan J2(+)J2 vs J2(+)J1(+)J1 witness,
      with full invariant-factor lists
      -> data/invariant_factor_witnesses.json
  * Ex 3.9 / Rem 3.10 (disjoint compositum, adjoining sqrt7): G_K, G_L, the
      Kronecker Gram, the out-of-field score 56, the det relation 9216^2*28^4
      -> data/disjoint_kronecker.json
  * Thm 3.6 / Ex 3.13 (the non-disjoint degree-4 witness
      Q(sqrt2)(sqrt2+sqrt3)=Q(sqrt2,sqrt3)): the factorisation over Q(sqrt2),
      the degree-6 operator polynomial x^6-25x^4+91x^2-75=(x^2-3)(x^4-22x^2+25),
      the selected m_theta = x^4-22x^2+25, and beta = -7/20 theta + 1/20 theta^3
      -> data/nondisjoint_compositum.json

Run:  py code/2026-06-residual-return-learning/crossfield.py

The bridge minimal polynomials are computed by the paper's LARGEST-INVARIANT-FACTOR
route (rrl_core.coords_to_minpoly), the non-disjoint operator polynomial by a
resultant + squarefree-part construction -- computed, not asserted.
"""
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols, resultant

import rrl_core as core

SCRIPT = "crossfield.py"
x = core.X
y = symbols("y")


def produce_bridge():
    """Ex 3.5: the coordinates -> minimal-polynomial bridge on five inputs."""
    cases = [
        ("phi",     [Q(1, 2), Q(1, 2)], [1, 0, -5],       "Q(sqrt5)"),
        ("2*sqrt6", [-5, 0, 1, 0],      [1, 0, -10, 0, 1], "Q(sqrt2+sqrt3)"),
        ("3",       [3, 0],             [1, 0, -5],        "Q(sqrt5)"),
        ("theta",   [0, 1, 0, 0],       [1, 0, -10, 0, 1], "Q(sqrt2+sqrt3)"),
        ("1/2",     [Q(1, 2), 0],       [1, 0, -5],        "Q(sqrt5)"),
    ]
    rows = []
    for name, coords, mc, field in cases:
        mp_expr, is_int = core.coords_to_minpoly(coords, mc)
        rows.append([name, field, core.poly_str(mp_expr).replace(" ", ""),
                     "accepted" if is_int else "rejected(non-integer)"])
    core.write_csv(SCRIPT, "bridge_minpolys.csv",
                   ["element", "field", "min_polynomial", "status"], rows)
    return rows


def _witness(A, label):
    facs = core.invariant_factors(A)
    nontrivial = [str(sp.factor(f)) for f in facs if sp.expand(f) != 1]
    charpoly = A.charpoly(x).as_expr()
    minpoly = facs[-1]
    return {
        "label": label,
        "charpoly": str(sp.factor(charpoly)),
        "minpoly_largest_invariant_factor": str(sp.factor(minpoly)),
        "invariant_factors_nontrivial": nontrivial,
        "trace": int(A.trace()),
        "rank": int(A.rank()),
    }


def produce_invariant_factor_witnesses():
    """Rem 3.7: charpoly (and even charpoly+minpoly) cannot separate similarity."""
    Cphi = core.companion([1, -1, -1])
    A = sp.diag(Cphi, Cphi)                                    # phi (+) phi
    B = core.companion(sp.Poly((x**2 - x - 1)**2, x).all_coeffs())  # C((x^2-x-1)^2)
    wA, wB = _witness(A, "rho(phi) (+) rho(phi)"), _witness(B, "C((x^2-x-1)^2)")
    # shared Mahler measure phi^2 (M((x^2-x-1)^2) = M(x^2-x-1)^2 = phi^2, from the
    # simple-root factor to avoid multiple-root instability)
    mah = core.mahler_measure([1, -1, -1]) ** 2

    J2 = Matrix([[0, 1], [0, 0]])
    JA = sp.diag(J2, J2)                                       # J2 (+) J2
    JB = sp.diag(J2, Matrix([[0]]), Matrix([[0]]))             # J2 (+) J1 (+) J1
    wJA, wJB = _witness(JA, "J2 (+) J2"), _witness(JB, "J2 (+) J1 (+) J1")

    payload = {
        "phi_witness": {
            "A": wA, "B": wB,
            "same_charpoly": bool(sp.expand(A.charpoly(x).as_expr()
                                            - B.charpoly(x).as_expr()) == 0),
            "same_trace": bool(A.trace() == B.trace() == 2),
            "shared_mahler_measure_phi_squared": str(mah),
            "similar": wA["invariant_factors_nontrivial"]
            == wB["invariant_factors_nontrivial"],
        },
        "jordan_witness": {
            "A": wJA, "B": wJB,
            "same_charpoly_x4": bool(JA.charpoly(x).as_expr()
                                     == JB.charpoly(x).as_expr() == x**4),
            "same_minpoly_x2": bool(wJA["minpoly_largest_invariant_factor"]
                                    == wJB["minpoly_largest_invariant_factor"] == "x**2"),
            "distinguished_by_rank": [wJA["rank"], wJB["rank"]],   # 2 vs 1
            "similar": wJA["invariant_factors_nontrivial"]
            == wJB["invariant_factors_nontrivial"],
        },
    }
    core.write_json(SCRIPT, "invariant_factor_witnesses.json", payload)
    return payload


def produce_disjoint():
    """Ex 3.9 / Rem 3.10: disjoint compositum K=Q(sqrt2,sqrt3), L=Q(sqrt7)."""
    GK = sp.diag(4, 8, 12, 24)          # Tr(c)=deg*c on {1,sqrt2,sqrt3,sqrt6}, deg 4
    GL = sp.diag(2, 14)                 # {1,sqrt7}, deg 2
    GKL = Matrix(sp.kronecker_product(GK, GL))
    diag = [int(GKL[i, i]) for i in range(8)]
    detK, detL = int(GK.det()), int(GL.det())
    prod8 = 1
    for v in diag:
        prod8 *= v
    payload = {
        "K": "Q(sqrt2,sqrt3)", "L": "Q(sqrt7)",
        "G_K_diag": [int(GK[i, i]) for i in range(4)],
        "G_L_diag": [int(GL[i, i]) for i in range(2)],
        "kronecker_gram_diag": diag,                     # 8,56,16,112,24,168,48,336
        "out_of_field_sqrt7_score": 8 * 7,               # 56
        "det_G_K": detK, "det_G_L": detL,                # 9216, 28
        "det_relation": "det(G_KL) = (det G_K)^[L:Q] (det G_L)^[K:Q] = 9216^2 * 28^4",
        "det_G_KL": int(GKL.det()),
        "det_via_formula": detK**2 * detL**4,
        "det_via_product_of_diag": prod8,
        "all_three_agree": bool(int(GKL.det()) == detK**2 * detL**4 == prod8),
    }
    core.write_json(SCRIPT, "disjoint_kronecker.json", payload)
    return payload


def produce_nondisjoint():
    """Thm 3.6 / Ex 3.13: the non-disjoint witness, resolving open problem O4
    in the bounded-degree case."""
    m_beta = x**4 - 10 * x**2 + 1                          # minpoly of beta=sqrt2+sqrt3
    # factor over Q(sqrt2):  (x^2 - 2 sqrt2 x - 1)(x^2 + 2 sqrt2 x - 1)
    factored = sp.expand((x**2 - 2 * sqrt(2) * x - 1) * (x**2 + 2 * sqrt(2) * x - 1))
    factors_over_K = bool(sp.simplify(factored - m_beta) == 0)
    # operator polynomial of theta = alpha + beta (alpha=sqrt2, c=1) over Q:
    #   squarefree part of Res_y(y^2-2, m_beta(x-y)), degree 6
    res = resultant(y**2 - 2, m_beta.subs(x, x - y), y)
    op = sp.expand(sp.sqf_part(sp.expand(res)))
    op_factors = sp.factor(op)
    m_theta = sp.minimal_polynomial(2 * sqrt(2) + sqrt(3), x)  # selected genuine factor
    # beta reconstruction in the power basis of theta = 2 sqrt2 + sqrt3
    theta = 2 * sqrt(2) + sqrt(3)
    beta_recon = -Q(7, 20) * theta + Q(1, 20) * theta**3
    recon_ok = bool(sp.simplify(beta_recon - (sqrt(2) + sqrt(3))) == 0)
    payload = {
        "problem": "non-disjoint compositum Q(sqrt2)(sqrt2+sqrt3) = Q(sqrt2,sqrt3)",
        "K": "Q(sqrt2) = Q[x]/(x^2-2)", "m": 2,
        "beta": "sqrt2+sqrt3", "m_beta": core.poly_str(m_beta), "deg_m_beta": 4,
        "m_beta_factors_over_K": "(x^2-2*sqrt2*x-1)*(x^2+2*sqrt2*x-1)",
        "factorisation_verified": factors_over_K,
        "effective_degree_e_prime": 2,
        "true_compositum_degree": 4,
        "tensor_degree_does_not_exist": 8,
        "operator_polynomial_deg6": core.poly_str(op),         # x^6-25x^4+91x^2-75
        "operator_polynomial_factored": str(op_factors),       # (x^2-3)(x^4-22x^2+25)
        "selected_factor_m_theta": core.poly_str(m_theta),     # x^4-22x^2+25
        "selected_m_theta_coeffs": core.poly_coeffs(m_theta),
        "spurious_factor": "x^2-3 (root sqrt3 from embedding alpha=-sqrt2)",
        "beta_reconstruction": "beta = -7/20*theta + 1/20*theta^3",
        "beta_coords_in_power_basis": ["0", "-7/20", "0", "1/20"],
        "reconstruction_verified": recon_ok,
    }
    core.write_json(SCRIPT, "nondisjoint_compositum.json", payload)
    return payload


def main():
    br = produce_bridge()
    ifw = produce_invariant_factor_witnesses()
    dj = produce_disjoint()
    nd = produce_nondisjoint()
    print(f"[{SCRIPT}] bridge: {[(r[0], r[2], r[3].split('(')[0]) for r in br]}")
    print(f"[{SCRIPT}] phi(+)phi similar to C(sq)? "
          f"{ifw['phi_witness']['similar']} (should be False)")
    print(f"[{SCRIPT}] Kronecker det agree three ways: {dj['all_three_agree']}")
    print(f"[{SCRIPT}] non-disjoint m_theta = {nd['selected_factor_m_theta']}; "
          f"beta recon verified = {nd['reconstruction_verified']}")
    print(f"[{SCRIPT}] wrote 4 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
