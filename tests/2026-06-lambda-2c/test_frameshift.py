"""
Independent verification of Sections 4 & 6: the three canonicalizations of c and
the frame-shift value.  Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - canon table : c in {1, n, sqrt(1+4C)/(2C)} (Jeffreys / degree / frame-shift).
  - def:frameshift / eq:frameshift : solving 2cC = sqrt(1+4C) gives c = sqrt(1+4C)/(2C);
                    at C=1: c = sqrt5/2, lambda = 2c = sqrt5 = phi - psi.
"""
import sympy as sp

C = sp.symbols('C', positive=True)


def test_frameshift_solves_gate_balance():
    """eq:frameshift: 'the frame-shift value of c is the solution of 2cC = sqrt(1+4C)',
    namely c = sqrt(1+4C)/(2C).  Re-derive by solving the gate balance for c.
    """
    c = sp.symbols('c', positive=True)
    sol = sp.solve(sp.Eq(2 * c * C, sp.sqrt(1 + 4 * C)), c)
    assert len(sol) == 1
    assert sp.simplify(sol[0] - sp.sqrt(1 + 4 * C) / (2 * C)) == 0


def test_frameshift_at_golden_gate():
    """def:frameshift: 'at the golden gate C=1: c = sqrt5/2 and lambda = 2c = sqrt5'.
    """
    c_frame = sp.sqrt(1 + 4 * C) / (2 * C)
    c_at_1 = sp.simplify(c_frame.subs(C, 1))
    assert sp.simplify(c_at_1 - sp.sqrt(5) / 2) == 0
    lam = sp.simplify(2 * c_at_1)
    assert sp.simplify(lam - sp.sqrt(5)) == 0


def test_lambda_is_spectral_gap_phi_minus_psi():
    """def:frameshift: 'lambda = sqrt5 = phi - psi' (the spectral gap of the self-action).
    """
    phi = (1 + sp.sqrt(5)) / 2
    psi = (1 - sp.sqrt(5)) / 2
    assert sp.simplify((phi - psi) - sp.sqrt(5)) == 0


def test_three_canonicalizations_values():
    """canon table: Jeffreys c=1, degree-invariant c=n, frame-shift
    c = sqrt(1+4C)/(2C).  Check the three selectors and that at C=1 the
    frame-shift value coincides with sqrt5/2 (cost floor lambda = sqrt5).
    """
    n = sp.symbols('n', positive=True)
    c_jeffreys = sp.Integer(1)
    c_degree = n
    c_frame = sp.sqrt(1 + 4 * C) / (2 * C)
    assert c_jeffreys == 1
    assert sp.simplify(c_degree - n) == 0
    assert sp.simplify(c_frame.subs(C, 1) - sp.sqrt(5) / 2) == 0
    # cost floors 2 log muS and 2n log muS are just lambda*log muS at c=1, c=n
    muS = sp.symbols('muS', positive=True)
    assert sp.simplify(2 * c_jeffreys * sp.log(muS) - 2 * sp.log(muS)) == 0
    assert sp.simplify(2 * c_degree * sp.log(muS) - 2 * n * sp.log(muS)) == 0
