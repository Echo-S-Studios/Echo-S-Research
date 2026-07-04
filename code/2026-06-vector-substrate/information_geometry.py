"""
PRODUCER: information geometry of the trace form (both Fisher metrics).

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex, Sec. 7
  - Thm. 7.1 (Gaussian location family: Fisher = M^T M = G)
  - Thm. 7.7 (conjugate/max-entropy family: Fisher_exp = (1/n)G - (1/n^2)tt^T)
  - Sec. 7.3 evidence table (Fisher_exp for Q(sqrt5), Q(sqrt2,sqrt3), +sqrt7)
  - Lemma 7.8 (trace-zero = residual subspace), Thm. 7.9 (||r||^2 = n Fisher(r))
  - Prop. 7.2 / Ex. 7.5 (discriminant = Jeffreys volume sqrt|d_K|; ~4.6x ratio)
  - Rem. 7.4 (MDL: (1/2)log det G = log covol)

The two Fisher metrics are RE-DERIVED here from first principles (Hessian of the
Gaussian log-likelihood; Hessian of the exponential-family log-partition) and
then emitted alongside the closed forms.

Emits:
  data/fisher_matrices.json
  data/jeffreys_volumes.csv
"""
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "information_geometry.py"


def gaussian_fisher_from_hessian():
    """Thm. 7.1: Fisher of N(Ma, I) = M^T M = G, derived as the constant Hessian
    of the negative log-likelihood (golden field, symbolic)."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    M = sp.Matrix([[1, phi], [1, phic]])
    a1, a2, y1, y2 = sp.symbols("a1 a2 y1 y2", real=True)
    a = sp.Matrix([a1, a2])
    y = sp.Matrix([y1, y2])
    mu = M * a
    ll = -sp.Rational(1, 2) * ((y - mu).T * (y - mu))[0]
    fisher = sp.simplify(-sp.hessian(ll, (a1, a2)))
    G = sp.simplify(M.T * M)
    return {
        "M": vc.mat_to_list(M),
        "fisher_gauss": vc.mat_to_list(fisher),
        "G_equals_MtM": vc.mat_to_list(G),
        "fisher_equals_G": fisher == G,
    }


def exp_fisher_from_logpartition():
    """Thm. 7.7: for p(k;a) ~ exp((Ma)_k), Fisher at a=0 = Hessian of the
    log-partition A(a) = log sum_k exp((Ma)_k) = (1/n)G - (1/n^2)tt^T (golden)."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    M = sp.Matrix([[1, phi], [1, phic]])
    a1, a2 = sp.symbols("a1 a2", real=True)
    eta = M * sp.Matrix([a1, a2])
    A = sp.log(sp.exp(eta[0]) + sp.exp(eta[1]))
    fisher = sp.simplify(sp.hessian(A, (a1, a2)).subs({a1: 0, a2: 0}))
    G = sp.Matrix([[2, 1], [1, 3]])
    t = sp.Matrix([2, 1])
    formula = vc.fisher_exp(G, t, 2)
    return {
        "fisher_exp_from_hessian": vc.mat_to_list(fisher),
        "fisher_exp_from_formula": vc.mat_to_list(formula),
        "match": fisher == formula,
        "value": vc.mat_to_list(sp.Matrix([[0, 0], [0, sp.Rational(5, 4)]])),
    }


def general_identity():
    """Thm. 7.7: M^T(I/n - 11^T/n^2)M = (1/n)G - (1/n^2)tt^T for a general M
    (symbolic 3x2)."""
    ms = sp.symbols("m11 m12 m21 m22 m31 m32")
    M = sp.Matrix([[ms[0], ms[1]], [ms[2], ms[3]], [ms[4], ms[5]]])
    n = 3
    one = sp.ones(n, 1)
    Cov = sp.eye(n) / n - (one * one.T) / n**2
    lhs = sp.simplify(M.T * Cov * M)
    rhs = sp.simplify(M.T * M / n - (M.T * one) * (M.T * one).T / n**2)
    return {"identity_holds": sp.simplify(lhs - rhs) == sp.zeros(2, 2)}


def evidence_table():
    """Sec. 7.3: Fisher_exp = (1/n)G - (1/n^2)tt^T for the three worked fields."""
    # Q(sqrt5), {1,phi}
    C = vc.companion_from_poly(x**2 - x - 1)
    G2 = vc.gram([(1, 0), (0, 1)], C)
    t2 = vc.trace_vector([(1, 0), (0, 1)], C)
    F2 = vc.fisher_exp(G2, t2, 2)
    # Q(sqrt2,sqrt3)
    G4 = sp.diag(4, 8, 12, 24)
    F4 = vc.fisher_exp(G4, sp.Matrix([4, 0, 0, 0]), 4)
    # Q(sqrt2,sqrt3,sqrt7)
    G8 = sp.diag(8, 56, 16, 112, 24, 168, 48, 336)
    F8 = vc.fisher_exp(G8, sp.Matrix([8, 0, 0, 0, 0, 0, 0, 0]), 8)
    return {
        "Q_sqrt5": {"n": 2, "G": vc.mat_to_list(G2), "fisher_exp": vc.mat_to_list(F2)},
        "Q_sqrt2_sqrt3": {"n": 4, "G": vc.mat_to_list(G4), "fisher_exp": vc.mat_to_list(F4)},
        "Q_sqrt2_sqrt3_sqrt7": {"n": 8, "G": vc.mat_to_list(G8), "fisher_exp": vc.mat_to_list(F8)},
    }


def residual_norm_is_fisher():
    """Lemma 7.8 / Thm. 7.9: on the trace-zero subspace G = n Fisher_exp; for
    sqrt5=(-1,2) in Q(sqrt5), ||sqrt5||_G^2 = 10 = 2*5 = n*Fisher(sqrt5).  Also
    the both-ways 1 in col(B) check (Rem. 7.11)."""
    C = vc.companion_from_poly(x**2 - x - 1)
    G = vc.gram([(1, 0), (0, 1)], C)
    t = vc.trace_vector([(1, 0), (0, 1)], C)
    F = vc.fisher_exp(G, t, 2)
    v = sp.Matrix([-1, 2])                        # sqrt5, trace zero
    # residual of the constant 1 against <phi> has trace 5/3 != 0
    r_of_1 = vc.residual((1, 0), sp.Matrix([0, 1]), G)
    return {
        "trace_sqrt5": vc.sval((t.T * v)[0]),                 # 0
        "sqrt5_norm2_G": vc.sval(vc.gnorm2(v, G)),            # 10
        "fisher_sqrt5": vc.sval((v.T * F * v)[0]),            # 5
        "identity_norm2_eq_n_fisher": vc.gnorm2(v, G) == 2 * (v.T * F * v)[0],
        "residual_of_1_off_phi_trace": vc.sval(vc.field_trace(list(r_of_1), C)),  # 5/3
    }


def jeffreys_rows():
    """Prop. 7.2 / Ex. 7.5 / Rem. 7.4: Jeffreys volume = sqrt|d_K|; ratio; and
    (1/2)log|d_K| = log covol."""
    cases = [("Q(sqrt5)", 5, "sqrt(5)"), ("Q(sqrt2)", 8, "2*sqrt(2)"),
             ("Q(cbrt2)", 108, "6*sqrt(3)")]
    rows = []
    v5 = sp.sqrt(5)
    for name, absdK, covol in cases:
        vol = sp.sqrt(absdK)
        rows.append([
            name,
            absdK,
            covol,
            round(float(vol), 6),
            round(float(vol / v5), 6),               # ratio vs Q(sqrt5)
            round(float(sp.Rational(1, 2) * sp.log(absdK)), 6),   # = log covol
        ])
    return rows


def main():
    payload = {
        "gaussian_fisher_from_hessian": gaussian_fisher_from_hessian(),
        "exp_fisher_from_logpartition": exp_fisher_from_logpartition(),
        "general_identity": general_identity(),
        "evidence_table": evidence_table(),
        "residual_norm_is_fisher": residual_norm_is_fisher(),
    }
    p1 = vc.write_json("fisher_matrices.json", payload, SCRIPT)
    p2 = vc.write_csv(
        "jeffreys_volumes.csv",
        ["field", "abs_d_K", "covol_exact", "jeffreys_volume", "ratio_vs_Q_sqrt5",
         "half_log_dK_eq_log_covol"],
        jeffreys_rows(),
        SCRIPT,
    )
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print("  Gaussian Fisher == G:", payload["gaussian_fisher_from_hessian"]["fisher_equals_G"])
    print("  exp Fisher (golden) =", payload["exp_fisher_from_logpartition"]["value"],
          "; ||sqrt5||^2 =", payload["residual_norm_is_fisher"]["sqrt5_norm2_G"], "= 2*5")


if __name__ == "__main__":
    main()
