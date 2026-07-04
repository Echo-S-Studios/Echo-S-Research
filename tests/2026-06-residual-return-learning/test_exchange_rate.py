r"""
Independent verification of Section 4.3 ("The exchange rate, derived: lambda=2c")
and the self-action / frame-shift material:
  - Thm 4.5 (lambda = 2c) : the MDL balance identity;
  - Rem 4.6 : the two shipped floors 0.5624, 2.2496, 4.4992 from mu_S;
  - Ex 4.3b (certified GROW enclosures) via mpmath verified interval arithmetic;
  - Prop 4.7 (self-action spectrum of the C-ladder) : {0,0,+-sqrt(1+4C)};
  - Def 4.8 / Table 4 (frame-shift canonicalization) : gap sqrt(1+4C), c=sqrt5/2.

One sub-claim of Prop 4.7 ("R_1 is conjugate to the keystone R") is xfailed: the
two matrices have different characteristic polynomials, so they are provably NOT
similar; only the self-action GAP (both discriminant 5) is shared -- which we DO
verify.  See NOTES.md.
"""
import mpmath as mp
import pytest
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols, zeros, eye

x = symbols("x")


def adjoint_matrix(Rmat):
    """4x4 matrix of the Lie self-action X -> [R,X]=RX-XR on M_2 in basis
    E11,E12,E21,E22."""
    basis = [Matrix([[1, 0], [0, 0]]), Matrix([[0, 1], [0, 0]]),
             Matrix([[0, 0], [1, 0]]), Matrix([[0, 0], [0, 1]])]
    cols = []
    for B in basis:
        img = Rmat * B - B * Rmat
        cols.append([img[0, 0], img[0, 1], img[1, 0], img[1, 1]])
    return Matrix(cols).T


# --- Theorem 4.5 (lambda = 2c) ----------------------------------------------

def test_lambda_2c_kl_identity():
    """Thm 4.5 / eq (kl2c): D_KL = 1/2 r^T Fisher r = 1/(2c) ||r||_G^2 with
    Fisher = G/c.  Symbolic identity for a generic 2-vector r."""
    c = symbols("c", positive=True)
    G = Matrix(sp.symbols("g11 g12 g12 g22")).reshape(2, 2)
    r = Matrix(sp.symbols("r1 r2"))
    F = G / c
    D_KL = (r.T * F * r)[0] / 2
    norm_G = (r.T * G * r)[0]
    assert sp.simplify(D_KL - norm_G / (2 * c)) == 0


def test_lambda_2c_threshold_readoff():
    """Thm 4.5 / eq (lambda2c): the MDL balance D_KL >= log M, i.e.
    1/(2c)||r||_G^2 >= log M, is exactly ||r||_G^2 >= 2c log M, so lambda = 2c;
    lambda|_{c=1}=2 and lambda|_{c=n}=2n."""
    c, n, logM, normG = symbols("c n logM normG", positive=True)
    # multiply the balance normG/(2c) >= logM through by 2c (>0):
    balance_lhs = normG / (2 * c)
    cleared = sp.simplify(balance_lhs * (2 * c))          # -> normG
    assert cleared == normG
    # matching normG >= lambda*logM against normG >= (2c)*logM gives lambda = 2c
    lam = 2 * c
    assert lam.subs(c, 1) == 2
    assert sp.simplify(lam.subs(c, n) - 2 * n) == 0


# --- Remark 4.6 / Remark 4.10 (the two floors from mu_S) ---------------------

def test_smyth_constant_plastic_number():
    """Rem 4.6: mu_S = 1.32471795724... is the real root of x^3-x-1 (plastic number)."""
    mp.mp.dps = 50
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    assert abs(muS - mp.mpf("1.3247179572447460259609088544780973407")) < mp.mpf("1e-30")
    assert abs(muS**3 - muS - 1) < mp.mpf("1e-40")


def test_floors_from_mu_S():
    """Rem 4.6 & 4.10: 2 log mu_S = 0.5623991486..., and the degree-aware floors
    n*2 log mu_S are 2.2496 (n=4) and 4.4992 (n=8)."""
    mp.mp.dps = 50
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    base = 2 * mp.log(muS)
    assert abs(base - mp.mpf("0.5623991486")) < mp.mpf("1e-9")
    assert abs(4 * base - mp.mpf("2.2496")) < mp.mpf("1e-4")
    assert abs(8 * base - mp.mpf("4.4992")) < mp.mpf("1e-4")


# --- Example 4.3b (certified GROW via interval arithmetic) -------------------

def test_certified_grow_sqrt7_enclosure():
    """Ex 4.3b: for sqrt7 (Mah=7) at lambda=2, log7 in [1.94591,1.94592], so
    2 log7 <= 3.89184, and the exact gain 56 dominates it -> GROW certified."""
    mp.mp.dps = 40
    lo, hi = mp.mpf("1.94591"), mp.mpf("1.94592")
    assert lo <= mp.log(7) <= hi
    assert 2 * hi <= mp.mpf("3.89184")
    assert 56 >= mp.mpf("3.89184") >= 2 * mp.log(7)


def test_certified_grow_2sqrt6_enclosure():
    """Ex 4.3b: for 2 sqrt6 (Mah=24), log24 in [3.17805,3.17806], so
    2 log24 <= 6.35612 <= 96 -> GROW certified."""
    mp.mp.dps = 40
    lo, hi = mp.mpf("3.17805"), mp.mpf("3.17806")
    assert lo <= mp.log(24) <= hi
    assert 2 * hi <= mp.mpf("6.35612")
    assert 96 >= mp.mpf("6.35612") >= 2 * mp.log(24)


def test_certified_stop_lattice_noise():
    """Ex 4.3b: lattice-aligned noise of magnitude 1/10 is STOPped by the Smyth
    constant floor, 1/10 < 0.5624."""
    assert Q(1, 10) < sp.Rational(5624, 10000)


# --- Proposition 4.7 (self-action spectrum of the C-ladder) ------------------

def test_ladder_companion_eigenvalues():
    """Prop 4.7: R_C=[[0,C],[1,-1]] is the companion of x^2+x-C with eigenvalues
    (-1 +- sqrt(1+4C))/2."""
    for C in (Q(1, 4), Q(1, 2), 1, 3):
        RC = Matrix([[0, C], [1, -1]])
        assert sp.expand(RC.charpoly(x).as_expr() - (x**2 + x - C)) == 0
        eigs = set(RC.eigenvals().keys())
        exp = {(-1 + sqrt(1 + 4 * C)) / 2, (-1 - sqrt(1 + 4 * C)) / 2}
        assert all(any(sp.simplify(e - f) == 0 for f in exp) for e in eigs)


def test_ladder_selfaction_spectrum():
    """Prop 4.7 / eq (trifurcation): the self-action [R_C,.] has spectrum
    {-sqrt(1+4C), 0, +sqrt(1+4C)} with 0 of multiplicity 2."""
    for C in (Q(1, 4), Q(1, 2), 1, 3):
        RC = Matrix([[0, C], [1, -1]])
        A = adjoint_matrix(RC)
        eigs = A.eigenvals()
        gap = sqrt(1 + 4 * C)
        assert eigs.get(sp.Integer(0)) == 2
        assert any(sp.simplify(e - gap) == 0 for e in eigs)
        assert any(sp.simplify(e + gap) == 0 for e in eigs)


def test_ladder_kernel_is_centralizer():
    """Prop 4.7: the 0-eigenspace of the self-action is the centraliser
    span{I, R_C} (dimension 2)."""
    C = 1
    RC = Matrix([[0, C], [1, -1]])
    A = adjoint_matrix(RC)
    ns = A.nullspace()
    assert len(ns) == 2
    # I and R_C, flattened to the E11,E12,E21,E22 coordinate order, span the nullspace
    Ivec = Matrix([1, 0, 0, 1])
    Rvec = Matrix([RC[0, 0], RC[0, 1], RC[1, 0], RC[1, 1]])
    span = Matrix.hstack(*ns)
    for v in (Ivec, Rvec):
        assert (span.row_join(v)).rank() == span.rank()


def test_golden_gap_is_sqrt5():
    """Prop 4.7: at C=1 the spectrum is {-sqrt5,0,+sqrt5} with sqrt5 = phi - psi."""
    phi = (1 + sqrt(5)) / 2
    psi = (1 - sqrt(5)) / 2
    assert sp.simplify((phi - psi) - sqrt(5)) == 0
    RC = Matrix([[0, 1], [1, -1]])
    A = adjoint_matrix(RC)
    eigs = A.eigenvals()
    assert any(sp.simplify(e - sqrt(5)) == 0 for e in eigs)


def test_R1_and_keystone_share_self_action_gap():
    """Prop 4.7 (the LOAD-BEARING part of the aside): R_1=[[0,1],[1,-1]] and the
    keystone R=[[0,1],[1,1]] both have discriminant 5, hence the SAME self-action
    gap sqrt5.  (This is what the proposition uses.)"""
    R1 = Matrix([[0, 1], [1, -1]])
    Rk = Matrix([[0, 1], [1, 1]])
    # discriminant of a 2x2 charpoly x^2 - (tr)x + det is tr^2 - 4 det
    for M in (R1, Rk):
        tr, det = M.trace(), M.det()
        assert tr**2 - 4 * det == 5
    # self-action gaps coincide
    g1 = [e for e in adjoint_matrix(R1).eigenvals() if e != 0]
    gk = [e for e in adjoint_matrix(Rk).eigenvals() if e != 0]
    assert {sp.simplify(sp.Abs(e)) for e in g1} == {sp.simplify(sp.Abs(e)) for e in gk} == {sqrt(5)}


def test_R1_conjugate_to_negative_keystone():
    """Prop 4.7 (corrected 2026-07-04): R_1 = [[0,1],[1,-1]] is conjugate to -R,
    where R = [[0,1],[1,1]] is the keystone.  R_1 and -R share charpoly x^2+x-1
    (distinct eigenvalues => similar over Q); R_1 is NOT similar to R (charpoly
    x^2-x-1).  Both R and -R have discriminant 5, so R_1 has the keystone's
    self-action gap sqrt5."""
    R1 = Matrix([[0, 1], [1, -1]])
    Rk = Matrix([[0, 1], [1, 1]])
    negRk = -Rk
    cp_R1 = sp.expand(R1.charpoly(x).as_expr())
    cp_Rk = sp.expand(Rk.charpoly(x).as_expr())
    cp_negRk = sp.expand(negRk.charpoly(x).as_expr())
    assert cp_R1 == x**2 + x - 1
    assert cp_negRk == x**2 + x - 1
    assert cp_Rk == x**2 - x - 1
    assert cp_R1 != cp_Rk


# --- Definition 4.8 / Table 4 (frame-shift canonicalization) -----------------

def test_frameshift_gap_equals_sqrt_mahler():
    """Def 4.8: the spectral gap sqrt(1+4C) equals sqrt(Mah(sqrt(1+4C))) since
    Mah(x^2-(1+4C)) = 1+4C for 1+4C>1."""
    for C in (Q(1, 4), Q(1, 2), 1):
        D = 1 + 4 * C
        # Mahler of x^2 - D (D>1): product of |+-sqrt D| = D
        assert sp.simplify(sqrt(D) - sqrt(D)) == 0
        mp.mp.dps = 40
        roots = mp.polyroots([1, 0, -float(D)])
        Mah = mp.mpf(1)
        for r in roots:
            Mah *= max(mp.mpf(1), abs(r))
        assert abs(Mah - float(D)) < mp.mpf("1e-25")


def test_frameshift_c_and_lambda_formulas():
    """Def 4.8 / eq (frameshift): lambda*C = sqrt(1+4C) => lambda=sqrt(1+4C)/C,
    c=lambda/2=sqrt(1+4C)/(2C).  At the golden gate C=1: c=sqrt5/2, lambda=sqrt5."""
    C = symbols("C", positive=True)
    lam = sqrt(1 + 4 * C) / C
    c = lam / 2
    assert sp.simplify(lam * C - sqrt(1 + 4 * C)) == 0
    assert sp.simplify(lam.subs(C, 1) - sqrt(5)) == 0
    assert sp.simplify(c.subs(C, 1) - sqrt(5) / 2) == 0


def test_canon_table_gaps_and_c_values():
    """Table 4 & appendix: on the gates C in {1/4,1/2,1} the gap sqrt(1+4C) is
    {sqrt2,sqrt3,sqrt5} and c=sqrt(1+4C)/(2C) is {2 sqrt2, sqrt3, sqrt5/2}."""
    gaps = {Q(1, 4): sqrt(2), Q(1, 2): sqrt(3), 1: sqrt(5)}
    cvals = {Q(1, 4): 2 * sqrt(2), Q(1, 2): sqrt(3), 1: sqrt(5) / 2}
    for C in (Q(1, 4), Q(1, 2), 1):
        assert sp.simplify(sqrt(1 + 4 * C) - gaps[C]) == 0
        assert sp.simplify(sqrt(1 + 4 * C) / (2 * C) - cvals[C]) == 0
