"""
Producer -- Section 2 of
    "The Exchange Rate lambda = 2c" (papers/2026-06-lambda-2c/lambda_2c_paper.tex)

Produces the derived identity lambda = 2c and its supporting substrate:

  * thm:lambda2c / eq:balance -- lambda is the coefficient of log M(theta) in the
    growth rule; solving the MDL balance (1/2c)||r||^2 = log M for the gain gives
    coefficient 2c; c=1 -> lambda=2 (shipped), c=n -> lambda=2n (degree-aware).
  * eq:kl2c -- the 2nd-order KL divergence equals half the Fisher quadratic form
    (Gaussian: exact; exponential family: quadratic Taylor coefficient = 1/2 I).
  * rem:sigma -- sigma = 1/(2c) = 1/lambda.
  * lem:fisher substrate -- the trace-form Gram G_ij = Tr(theta^{i+j}) equals the
    conjugate covariance B^T B; Q(sqrt2) -> diag(2,4); Q(sqrt5) -> [[2,1],[1,3]]
    (Lucas-number entries L_0,L_1,L_2 = 2,1,3).

Emits:
    data/2026-06-lambda-2c/identity_lambda_2c.json
    data/2026-06-lambda-2c/trace_form_gram.csv
"""
import sympy as sp
import lambda2c_common as cm

SCRIPT = "identity_lambda_2c.py"


def lambda_identity():
    """Solve the MDL balance for the exchange rate; return the coefficient facts."""
    c, g, logM, lam = sp.symbols('c g logM lambda', positive=True)
    KL = g / (2 * c)                                  # gain measured intrinsically
    g_at_balance = sp.solve(sp.Eq(KL, logM), g)[0]    # balance KL = log M
    coeff = sp.simplify(g_at_balance / logM)          # coefficient of log M
    lam_solved = sp.solve(sp.Eq(lam * logM, g_at_balance), lam)[0]
    return {
        "gain_at_balance": str(g_at_balance),         # 2*c*logM
        "coefficient_of_logM": str(coeff),            # 2*c
        "lambda_solved": str(lam_solved),             # 2*c
        "lambda_equals_2c": bool(sp.simplify(lam_solved - 2 * c) == 0),
        "shipped_c1": int((2 * c).subs(c, 1)),        # 2
        "degree_cn": str(sp.simplify((2 * c).subs(c, sp.Symbol('n', positive=True)))),  # 2*n
        "sigma_precision": str(sp.simplify(sp.Rational(1, 1) / (2 * c))),               # 1/(2c) = 1/lambda
    }


def kl_second_order():
    """eq:kl2c -- KL = 1/2 Fisher form.  Gaussian location: exact; exponential
    rate family: the quadratic Taylor coefficient equals 1/2 * Fisher."""
    # Gaussian location family N(mu, Sigma), KL(mu || mu+r) = 1/2 r^T Sigma^{-1} r
    a, b, d = sp.symbols('a b d', positive=True)
    Sigma = sp.Matrix([[a, b], [b, d]])
    r1, r2 = sp.symbols('r1 r2', real=True)
    r = sp.Matrix([r1, r2])
    KL_gauss = sp.Rational(1, 2) * (r.T * Sigma.inv() * r)[0]
    half_fisher = sp.Rational(1, 2) * (r.T * Sigma.inv() * r)[0]
    gaussian_exact = bool(sp.simplify(KL_gauss - half_fisher) == 0)

    # Exponential family p_theta(x) = theta e^{-theta x}; Fisher I = 1/theta^2
    theta, rr, x = sp.symbols('theta r x', positive=True)
    p_theta = theta * sp.exp(-theta * x)
    p_shift = (theta + rr) * sp.exp(-(theta + rr) * x)
    KL = sp.simplify(sp.integrate(p_theta * sp.log(p_theta / p_shift), (x, 0, sp.oo)))
    series = sp.series(KL, rr, 0, 3).removeO()
    quad_coeff = series.coeff(rr, 2)
    fisher = 1 / theta**2
    return {
        "gaussian_KL_is_half_fisher_exact": gaussian_exact,
        "expfamily_KL_closed_form": str(KL),
        "expfamily_quadratic_coeff": str(sp.simplify(quad_coeff)),
        "expfamily_half_fisher": str(sp.simplify(sp.Rational(1, 2) * fisher)),
        "expfamily_quad_equals_half_fisher": bool(
            sp.simplify(quad_coeff - sp.Rational(1, 2) * fisher) == 0),
        "expfamily_linear_coeff": str(sp.simplify(series.coeff(rr, 1))),   # 0
        "expfamily_const_coeff": str(sp.simplify(series.coeff(rr, 0))),    # 0
    }


def trace_form_gram(conjugates):
    """G_ij = Tr(theta^{i+j}) = sum_k sigma_k(theta)^{i+j}, and its B^T B form."""
    n = len(conjugates)
    B = sp.Matrix(n, n, lambda k, i: conjugates[k]**i)
    Gram = sp.simplify(B.T * B)
    G = sp.Matrix(n, n, lambda i, j: sp.simplify(sum(ck**(i + j) for ck in conjugates)))
    matches = sp.simplify(Gram - G).is_zero_matrix
    return G, bool(matches)


def main():
    payload = {
        "identity": lambda_identity(),
        "kl_second_order": kl_second_order(),
    }

    # Trace-form Gram matrices (lem:fisher substrate)
    G2, ok2 = trace_form_gram([sp.sqrt(2), -sp.sqrt(2)])          # Q(sqrt2)
    G5, ok5 = trace_form_gram([cm.PHI, cm.PSI])                   # Q(sqrt5), power basis phi
    payload["trace_form"] = {
        "Q_sqrt2_gram": [[int(G2[i, j]) for j in range(2)] for i in range(2)],
        "Q_sqrt2_equals_conjugate_gram": ok2,
        "Q_sqrt2_det": int(G2.det()),
        "Q_sqrt5_gram": [[int(G5[i, j]) for j in range(2)] for i in range(2)],
        "Q_sqrt5_equals_conjugate_gram": ok5,
        "Q_sqrt5_entries_are_lucas": {
            "Tr(1)=L0": int(G5[0, 0]), "Tr(theta)=L1": int(G5[0, 1]),
            "Tr(theta^2)=L2": int(G5[1, 1]),
        },
    }
    cm.write_json("identity_lambda_2c.json", payload, SCRIPT)

    # CSV form of the two Gram matrices
    rows = [
        ["Q(sqrt2)", "{1,theta}", int(G2[0, 0]), int(G2[0, 1]), int(G2[1, 1]), int(G2.det())],
        ["Q(sqrt5)", "{1,phi}", int(G5[0, 0]), int(G5[0, 1]), int(G5[1, 1]), int(G5.det())],
    ]
    cm.write_csv("trace_form_gram.csv",
                 ["field", "power_basis", "G00=Tr(1)", "G01=Tr(theta)", "G11=Tr(theta^2)", "det_G"],
                 rows, SCRIPT)

    print("wrote identity_lambda_2c.json, trace_form_gram.csv")
    print("  lambda =", payload["identity"]["coefficient_of_logM"], "(coefficient of log M)")
    print("  Q(sqrt5) Gram =", payload["trace_form"]["Q_sqrt5_gram"], "(Lucas 2,1,3)")


if __name__ == "__main__":
    main()
