"""
Independent verification of Section 'entry' (the superselection sector):
the angle charge A of the catalog generators (Prop), the tensor products
phi (x) phi and phi^4 (x) phi^4, the lattice closure of the operators, and
the off-lattice lift for beta_4.
"""
import sympy as sp
import mpmath as mp
import cmath

x = sp.symbols('x')


def _args_over_pi(poly):
    """arguments (in units of pi) of all roots of a polynomial."""
    return [cmath.phase(complex(r)) / cmath.pi for r in sp.nroots(poly, n=30)]


def _on_half_integer_lattice(a_over_pi, tol=1e-9):
    """True iff arg/pi is an integer multiple of 1/2 (i.e. angle in (pi/2)Z)."""
    return abs(2 * a_over_pi - round(2 * a_over_pi)) < tol


def test_K_generator_angles_on_lattice():
    """Prop: K = x^4 + 5x^2 - 5 has all roots at arguments in (pi/2)Z:
    a real pair at {0,pi} and an imaginary pair +-i*beta at +-pi/2."""
    K = x**4 + 5*x**2 - 5
    args = _args_over_pi(K)
    assert all(_on_half_integer_lattice(a) for a in args)
    # exactly one imaginary pair at +-1/2 and one real pair at {0,1}
    rounded = sorted(round(2 * a) for a in args)   # in half-pi units
    assert rounded == [-1, 0, 1, 2]                # {-pi/2, 0, pi/2, pi}


def test_real_seed_angles_on_lattice():
    """Prop: the real seed phi (root of x^2-x-1) has roots at args {0,pi}."""
    seed = x**2 - x - 1
    args = _args_over_pi(seed)
    assert all(_on_half_integer_lattice(a) for a in args)
    assert sorted(round(a) for a in args) == [0, 1]   # 0 and pi


def test_tensor_phi_phi():
    """Sec entry: phi (x) phi (pairwise products of the roots of x^2-x-1)
    = (x+1)^2 (x^2-3x+1); its on-circle root is -1, a 2nd root of unity."""
    seed = x**2 - x - 1
    roots = list(sp.roots(seed).keys())
    prods = [a * b for a in roots for b in roots]
    poly = sp.expand(sp.prod([x - p for p in prods]))
    assert sp.factor(poly) == (x + 1)**2 * (x**2 - 3*x + 1)
    # on-circle root is -1 (a 2nd root of unity, order 2): read it directly
    # from the exact pairwise products (avoids numerics on the double root)
    neg_ones = [p for p in prods if sp.simplify(p + 1) == 0]
    assert len(neg_ones) == 2                 # phi*psi and psi*phi
    assert sp.Abs(sp.simplify(neg_ones[0])) == 1   # modulus 1 => on the unit circle


def test_tensor_phi4_phi4():
    """Sec entry: phi^4 (x) phi^4 = (x-1)^2 (x^2-47x+1); its on-circle root
    is +1, a 1st root of unity (order 1). phi^4's minimal polynomial is
    derived independently."""
    phi = (1 + sp.sqrt(5)) / 2
    seed4 = sp.minimal_polynomial(phi**4, x)     # = x^2 - 7x + 1
    assert sp.expand(seed4 - (x**2 - 7*x + 1)) == 0
    roots = list(sp.roots(seed4).keys())
    prods = [sp.simplify(a * b) for a in roots for b in roots]
    poly = sp.expand(sp.prod([x - p for p in prods]))
    assert sp.factor(poly) == (x - 1)**2 * (x**2 - 47*x + 1)
    # on-circle root is +1 (order 1): read it directly from the exact products
    ones = [p for p in prods if sp.simplify(p - 1) == 0]
    assert len(ones) == 2                     # phi^4*psi^4 and psi^4*phi^4
    assert sp.Abs(sp.simplify(ones[0])) == 1


def test_operator_lattice_closure():
    """Thm (superselection): (pi/2)Z is closed under squaring's angle
    doubling (2 * (pi/2)Z ⊆ (pi/2)Z) and tensor's angle addition
    ((pi/2)Z + (pi/2)Z ⊆ (pi/2)Z). In units of pi/2 these are the
    statements 2k in Z and j+k in Z."""
    for k in range(-6, 7):
        assert (2 * k) == int(2 * k)                       # doubling stays integer
    for j in range(-4, 5):
        for k in range(-4, 5):
            assert (j + k) == int(j + k)                   # addition stays integer
    # doubled fourth-root-of-unity lattice is contained in itself (mod 4)
    lattice = {0, 1, 2, 3}                                 # units of pi/2 mod 2pi
    assert {(2 * k) % 4 for k in lattice} <= lattice
    assert {(j + k) % 4 for j in lattice for k in lattice} <= lattice


def test_beta4_offlattice_trace_root_lifts_to_salem_conjugate():
    """Sec entry: for beta_4 the trace root (1-sqrt13)/2 = -1.3028 is an
    OFF-lattice captured position; its lift x^2 - t x + 1 gives an
    on-circle conjugate at argument 0.726 pi (not in (pi/2)Z)."""
    mp.mp.dps = 40
    tr = (1 - mp.sqrt(13)) / 2
    assert abs(tr - mp.mpf('-1.3027756')) < 5e-7
    root = (tr + mp.sqrt(tr**2 - 4)) / 2                   # complex lift
    assert abs(abs(root) - 1) < mp.mpf(10)**(-30)          # on the unit circle
    arg_over_pi = mp.arg(root) / mp.pi
    assert abs(arg_over_pi - mp.mpf('0.725813')) < 1e-5
    # off the (pi/2)Z lattice: 2*arg/pi is not an integer
    assert abs(2 * arg_over_pi - round(float(2 * arg_over_pi))) > 1e-3
