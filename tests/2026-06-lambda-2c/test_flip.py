"""
Independent verification of Sections 9 & 13: "The flip: D = 1+4C" and the second flip.
Paper: "The Exchange Rate lambda = 2c ...".

Claims:
  - thm:flip table : for C in {1,0,-1/4,-1,-2}, D=1+4C in {5,1,0,-3,-7}; roots of
                     x^2+x-C; gap +-sqrt D real (D>0) or imaginary (D<0); double root
                     -1/2 at C=-1/4; e^{+-2pi i/3} at C=-1.
  - prop:detG      : in the gap basis {1, sqrt D}, G = diag(2, 2D), det G = 4D;
                     PD iff D>0, Lorentzian iff D<0, degenerate at D=0.
                     Q(i): basis {1,i} -> diag(2,-2), det=-4=4D with D=-1.
  - prop:meet      : at C=-1, x^2+x+1 has roots e^{+-2pi i/3} on the unit circle
                     (Mah=1), discriminant D=-3.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
x = sp.symbols('x')


def test_flip_table_discriminants_and_roots():
    """thm:flip: the flip table.  For each C compute D=1+4C and the roots of
    x^2+x-C; check the stated D and the real/complex character of the gap.
    """
    cases = {
        sp.Integer(1): 5,
        sp.Integer(0): 1,
        sp.Rational(-1, 4): 0,
        sp.Integer(-1): -3,
        sp.Integer(-2): -7,
    }
    for Cval, Dexp in cases.items():
        D = 1 + 4 * Cval
        assert D == Dexp
        roots = sp.solve(x**2 + x - Cval, x)
        if Dexp > 0:
            assert all(sp.im(sp.nsimplify(r)) == 0 for r in roots)     # real
        elif Dexp == 0:
            assert roots == [sp.Rational(-1, 2)]                        # double root
        else:
            assert all(sp.im(r) != 0 for r in roots)                    # complex


def test_flip_C0_roots():
    """thm:flip: 'C=0, D=1, real roots {-1, 0}, gap +-1'."""
    roots = set(sp.solve(x**2 + x - 0, x))
    assert roots == {0, -1}


def test_flip_double_root_at_minus_quarter():
    """thm:flip: 'C=-1/4, D=0, double root -1/2, gap 0'."""
    Cval = sp.Rational(-1, 4)
    p = sp.Poly(x**2 + x - Cval, x)
    assert sp.discriminant(p) == 0
    roots = sp.roots(p)                       # dict root->multiplicity
    assert roots == {sp.Rational(-1, 2): 2}


def test_flip_cube_roots_of_unity_at_minus_one():
    """thm:flip / prop:meet: 'C=-1, x^2+x+1 has roots e^{+-2pi i/3}, gap +-i sqrt3,
    D=-3'.
    """
    Cval = sp.Integer(-1)
    poly = x**2 + x - Cval        # = x^2 + x + 1
    assert sp.expand(poly) == x**2 + x + 1
    roots = sp.solve(poly, x)
    # e^{+-2pi i/3} in rectangular form (expand_complex ties exp to the algebraic form)
    z1 = sp.expand_complex(sp.exp(2 * sp.pi * sp.I / 3))
    z2 = sp.expand_complex(sp.exp(-2 * sp.pi * sp.I / 3))
    assert sp.simplify(z1 - (sp.Rational(-1, 2) + sp.sqrt(3) / 2 * sp.I)) == 0
    got = {sp.simplify(r) for r in roots}
    assert any(sp.simplify(r - z1) == 0 for r in got)
    assert any(sp.simplify(r - z2) == 0 for r in got)
    # each is a primitive cube root of unity: r^3 = 1, r != 1
    for r in roots:
        assert sp.simplify(r**3 - 1) == 0
        assert sp.simplify(r - 1) != 0
        assert sp.simplify(sp.Abs(r) - 1) == 0        # on the unit circle
    # gap +- i sqrt3
    assert sp.discriminant(sp.Poly(poly, x)) == -3


def test_gap_imaginary_at_C_minus_two():
    """thm:flip: 'C=-2, D=-7, complex roots, gap +- i sqrt7'."""
    Cval = sp.Integer(-2)
    p = sp.Poly(x**2 + x - Cval, x)
    assert sp.discriminant(p) == -7
    gap = sp.sqrt(sp.discriminant(p))
    assert sp.simplify(gap - sp.I * sp.sqrt(7)) == 0


# ---------------------------------------------------------------------------
# prop:detG  --  metric signature flip
# ---------------------------------------------------------------------------

def test_detG_equals_4D_gap_basis():
    """prop:detG / eq:detG: in the gap basis {1, sqrt D} with sqrt D = 2 theta + 1,
    theta a root of x^2+x-C, the trace-form Gram is diag(2, 2D) with det G = 4D.
    Re-derive the traces: Tr(1)=2, Tr(sqrt D)=2 Tr(theta)+Tr(1)=2(-1)+2=0, Tr(D)=2D.
    """
    C = sp.symbols('C')
    D = 1 + 4 * C
    # theta root of x^2+x-C: sum of the two conjugate thetas = -1 (Vieta)
    tr_theta = sp.Integer(-1)
    tr_one = sp.Integer(2)                       # degree-2 field: Tr(1)=2
    tr_sqrtD = 2 * tr_theta + tr_one             # sqrt D = 2 theta + 1
    assert tr_sqrtD == 0
    tr_D = 2 * D                                 # (sqrt D)^2 = D, Tr = 2D
    G = sp.Matrix([[tr_one, tr_sqrtD], [tr_sqrtD, tr_D]])
    assert G == sp.Matrix([[2, 0], [0, 2 * D]])
    assert sp.expand(G.det() - 4 * D) == 0


def test_signature_flip_signs():
    """prop:detG: 'PD (Riemannian) iff D>0, indefinite (1,1) Lorentzian iff D<0,
    degenerate at D=0'.  Check eigenvalue signs of G=diag(2,2D).
    """
    for Dval, expect in [(5, 'pd'), (1, 'pd'), (0, 'deg'), (-3, 'lor'), (-7, 'lor')]:
        G = sp.Matrix([[2, 0], [0, 2 * Dval]])
        eigs = list(G.eigenvals().keys())
        assert G.det() == 4 * Dval
        if expect == 'pd':
            assert all(e > 0 for e in eigs)
        elif expect == 'deg':
            assert 0 in eigs
        else:
            assert any(e > 0 for e in eigs) and any(e < 0 for e in eigs)


def test_Qi_trace_form_is_diag_2_minus2():
    """prop:detG canonical instance: 'Q(i): basis {1,i} gives diag(2,-2), det=-4=4D
    with D=-1 (since i^2=-1)'.  Tr_{Q(i)/Q}(1)=2, Tr(i)=0, Tr(i^2)=Tr(-1)=-2.
    """
    conj = [sp.I, -sp.I]                          # the two embeddings of i
    tr = lambda a: sum(a_k for a_k in a)
    # basis {1, i}: Gram_{jk} = Tr(b_j b_k) with b in {1, i}
    b = [[sp.Integer(1), sp.Integer(1)], [conj[0], conj[1]]]   # embeddings of 1 and i
    G = sp.Matrix(2, 2, lambda j, k: sp.simplify(sum(b[j][e] * b[k][e] for e in (0, 1))))
    assert G == sp.Matrix([[2, 0], [0, -2]])
    assert G.det() == -4
    assert G.det() == 4 * (-1)                    # D = -1


def test_two_flips_meet_mahler_one():
    """prop:meet: 'at C=-1 the additive flip's rotation regime (D<0) meets the
    multiplicative flip's marginal vertex |lambda|=1: roots on the unit circle,
    Mah=1'.
    """
    roots = mp.polyroots([1, 1, 1])              # x^2+x+1
    M = mp.mpf(1)
    for r in roots:
        if abs(r) > 1:
            M *= abs(r)
    assert abs(M - 1) < mp.mpf(10) ** (-30)      # cyclotomic, Mahler 1
    for r in roots:
        assert abs(abs(r) - 1) < mp.mpf(10) ** (-30)
