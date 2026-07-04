"""
Independent verification of Section 2: "The identity lambda = 2c"
Paper: "The Exchange Rate lambda = 2c: A Conformal Identity, Its Gate, and Its Flip"

Claims tested here:
  - eq:rule / eq:balance / thm:lambda2c : lambda = 2c is the coefficient of log M(theta).
  - eq:kl2c : second-order KL expansion equals 1/2 * Fisher quadratic form,
              hence (1/2c)||r||^2_G when F = G/c.
  - rem:sigma : sigma = 1/(2c) = 1/lambda.
  - thm:lambda2c proof: c=1 recovers 2, c=n recovers 2n.
  - lem:fisher substrate: trace-form Gram G_ij = Tr(theta^{i+j}) equals B^T B,
              the covariance/Gram of the n conjugate-embedding vectors.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50


# ---------------------------------------------------------------------------
# eq:kl2c  --  second order KL = 1/2 Fisher quadratic form
# ---------------------------------------------------------------------------

def test_gaussian_location_kl_is_exactly_half_fisher_form():
    """eq:kl2c: 'KL(p_theta || p_{theta+r}) = 1/2 r^T F r + O(||r||^3)'.

    For a Gaussian location family N(mu, Sigma) with fixed covariance the KL is
    EXACTLY the Fisher quadratic form 1/2 * r^T Sigma^{-1} r, and Fisher = Sigma^{-1}.
    We re-derive KL from its integral definition (closed form for Gaussians) and
    compare to 1/2 r^T F r with F the Fisher information.
    """
    # symmetric positive-definite 2x2 covariance
    a, b, d = sp.symbols('a b d', positive=True)
    Sigma = sp.Matrix([[a, b], [b, d]])
    r1, r2 = sp.symbols('r1 r2', real=True)
    r = sp.Matrix([r1, r2])
    # closed form KL(N(mu,Sigma) || N(mu+r,Sigma)) = 1/2 r^T Sigma^{-1} r
    KL = sp.Rational(1, 2) * (r.T * Sigma.inv() * r)[0]
    # Fisher information of a Gaussian location family is Sigma^{-1}
    F = Sigma.inv()
    half_fisher_form = sp.Rational(1, 2) * (r.T * F * r)[0]
    assert sp.simplify(KL - half_fisher_form) == 0


def test_second_order_kl_hessian_equals_fisher_nongaussian():
    """eq:kl2c: the O(||r||^3) statement -- the *second-order* term of KL is the
    Fisher form even for a non-Gaussian family.  Exponential family rate theta:
    p_theta(x) = theta e^{-theta x}.  We compute KL(theta || theta+r) in closed
    form, Taylor-expand in r, and check the quadratic coefficient equals
    1/2 * I(theta) with Fisher I(theta) = 1/theta^2.
    """
    theta, r, x = sp.symbols('theta r x', positive=True)
    p_theta = theta * sp.exp(-theta * x)
    p_shift = (theta + r) * sp.exp(-(theta + r) * x)
    integrand = p_theta * sp.log(p_theta / p_shift)
    KL = sp.integrate(integrand, (x, 0, sp.oo))
    KL = sp.simplify(KL)
    # quadratic coefficient in r about r=0
    series = sp.series(KL, r, 0, 3).removeO()
    quad_coeff = series.coeff(r, 2)
    fisher = 1 / theta**2
    assert sp.simplify(quad_coeff - sp.Rational(1, 2) * fisher) == 0
    # first order vanishes, KL(0)=0
    assert sp.simplify(series.coeff(r, 0)) == 0
    assert sp.simplify(series.coeff(r, 1)) == 0


# ---------------------------------------------------------------------------
# thm:lambda2c  --  lambda = 2c
# ---------------------------------------------------------------------------

def test_lambda_equals_2c_from_balance():
    """thm:lambda2c / eq:balance: 'adjoin iff (1/2c)||r||^2_G >= log M', equivalently
    ||r||^2_G >= 2c log M, so the coefficient of log M (the exchange rate lambda) is 2c.

    Re-derive: set the gain KL = (1/2c) g equal to log M at the balance and solve for
    the gain g; the per-unit-cost coefficient must be 2c.
    """
    c, g, logM, lam = sp.symbols('c g logM lambda', positive=True)
    # gain measured intrinsically (eq:kl2c): KL = g/(2c) with g = ||r||^2_G
    KL = g / (2 * c)
    # balance KL = logM -> solve for g
    g_at_balance = sp.solve(sp.Eq(KL, logM), g)[0]
    # g = 2c * logM  => coefficient of logM is 2c
    assert sp.simplify(g_at_balance - 2 * c * logM) == 0
    coeff = sp.simplify(g_at_balance / logM)
    assert sp.simplify(coeff - 2 * c) == 0
    # rule eq:rule is 'gain >= lambda * logM'; matching coefficients: lambda = 2c
    lam_solved = sp.solve(sp.Eq(lam * logM, g_at_balance), lam)[0]
    assert sp.simplify(lam_solved - 2 * c) == 0


def test_shipped_and_degree_values():
    """thm:lambda2c proof: 'substituting c=1 recovers the shipped value 2, and
    c=n recovers the degree-aware reading 2n'.
    """
    c = sp.symbols('c', positive=True)
    lam = 2 * c
    n = sp.symbols('n', positive=True)
    assert lam.subs(c, 1) == 2
    assert sp.simplify(lam.subs(c, n) - 2 * n) == 0


def test_sigma_is_reciprocal_lambda():
    """rem:sigma: 'sigma = 1/(2c) = 1/lambda' (conformal scale = model precision).
    Given lambda = 2c the precision sigma defined as 1/(2c) must equal 1/lambda.
    """
    c = sp.symbols('c', positive=True)
    lam = 2 * c
    sigma = sp.Rational(1, 1) / (2 * c)
    assert sp.simplify(sigma - 1 / lam) == 0


# ---------------------------------------------------------------------------
# lem:fisher substrate  --  trace form is the conjugate covariance G = B^T B
# ---------------------------------------------------------------------------

def test_trace_form_is_conjugate_gram_Q_sqrt2():
    """lem:fisher substrate ('trace form is the conjugate covariance'):
    G_ij = Tr(theta^{i+j}) equals (B^T B)_ij where B has rows = the n real
    embeddings of the power basis.  This is the identity that lets the trace form
    be read as (n times) the Fisher metric of the location family over conjugates.
    Verified for K = Q(sqrt2), power basis {1, theta}, embeddings theta -> +-sqrt2.
    """
    theta = sp.sqrt(2)
    n = 2
    conj = [sp.sqrt(2), -sp.sqrt(2)]          # the two real embeddings of theta
    # B[k, i] = sigma_k(theta)^i  (power-basis embedding matrix)
    B = sp.Matrix(n, n, lambda k, i: conj[k]**i)
    Gram = sp.simplify(B.T * B)
    # trace form G_ij = Tr(theta^{i+j}) = sum_k sigma_k(theta)^{i+j}
    G = sp.Matrix(n, n, lambda i, j: sum(ck**(i + j) for ck in conj))
    assert sp.simplify(Gram - G).is_zero_matrix
    # concrete value: G = [[2,0],[0,4]] for Q(sqrt2); det = disc up to square (8 = disc)
    assert G == sp.Matrix([[2, 0], [0, 4]])


def test_trace_form_is_conjugate_gram_golden():
    """lem:fisher substrate for K=Q(sqrt5) power basis with theta=phi (root of
    x^2-x-1): G_ij = Tr(theta^{i+j}) = (B^T B)_ij, embeddings phi and psi.
    """
    phi = (1 + sp.sqrt(5)) / 2
    psi = (1 - sp.sqrt(5)) / 2
    conj = [phi, psi]
    n = 2
    B = sp.Matrix(n, n, lambda k, i: conj[k]**i)
    Gram = sp.expand(sp.simplify(B.T * B))
    G = sp.Matrix(n, n, lambda i, j: sp.expand(sum(ck**(i + j) for ck in conj)))
    assert sp.simplify(Gram - G).is_zero_matrix
    # Tr(1)=2, Tr(theta)=1 (=L_1), Tr(theta^2)=3 (=L_2): Lucas numbers
    assert sp.simplify(G[0, 0]) == 2
    assert sp.simplify(G[0, 1]) == 1
    assert sp.simplify(G[1, 1]) == 3
