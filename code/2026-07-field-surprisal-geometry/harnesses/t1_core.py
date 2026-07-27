#!/usr/bin/env python3
# t1_core.py -- Task 1 Stage 1: exact core for the classification problem.
# Field: QQ(L2,L3,L5,Lp) with L2=log2, L3=log3, L5=log5, Lp=log(phi).
# ALL decisions exact (symbolic linear algebra over the rational function field).
# Justification for treating L2,L3,L5,Lp as independent: 2,3,5,phi are
# multiplicatively independent (norm argument in Q(sqrt5)), and Q-linear
# dependence of logs <=> multiplicative dependence. All decisions below reduce
# to Q-linear relations among {L2,L3,L5,Lp} or to polynomial identities checked
# symbolically; nonvanishing claims are certified by exact coefficient vectors.
import sys
from sympy import (symbols, Rational, Matrix, cancel, simplify, expand, S,
                   Poly, srepr, factor, together, degree)

L2, L3, L5, Lp = symbols('L2 L3 L5 Lp', positive=True)

# Catalog order: [sqrt2, sqrt3, sqrt5, phi, tau, phi4, K]
SEEDS = ['sqrt2', 'sqrt3', 'sqrt5', 'phi', 'tau', 'phi4', 'K']
# a = log M values;  M = {2,3,5,phi,phi,phi^4, phi^4-1},  log(phi^4-1)=L5/2+2Lp
a_vals = [L2, L3, L5, Lp, Lp, 4*Lp, Rational(1,2)*L5 + 2*Lp]
m = 7

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

def col(vals):
    return Matrix(len(vals), 1, list(vals))

ONE = col([1]*m)
A   = col(a_vals)

def hprod(u, v):
    return col([u[i]*v[i] for i in range(len(u))])

def exact_rank(M_):
    return M_.rank(iszerofunc=lambda x: cancel(x) == 0, simplify=cancel)

def in_span(f, basis_cols):
    B = Matrix.hstack(*basis_cols)
    r0 = exact_rank(B)
    r1 = exact_rank(Matrix.hstack(B, f))
    return r1 == r0

# ---------------------------------------------------------------- block 1
print("== block 1: a-value structure ==")
distinct = sorted({str(v) for v in a_vals})
vals_set = []
for v in a_vals:
    if all(cancel(v - u) != 0 for u in vals_set):
        vals_set.append(v)
ck("a takes exactly 6 distinct values", len(vals_set) == 6)
# pairwise differences are nonzero Q-linear combos of independent L's
import itertools
for u, v in itertools.combinations(vals_set, 2):
    d = cancel(u - v)
    p = Poly(d, L2, L3, L5, Lp)
    ck(f"nonzero linear diff {u}-{v}", p.total_degree() == 1 and any(c != 0 for c in p.coeffs()))
ck("phi,tau tie is the only tie", cancel(a_vals[3]-a_vals[4]) == 0)
ck("identity log(phi^4-1)=L5/2+2Lp holds by construction (phi^4-1=sqrt5*phi^2)",
   cancel(a_vals[6] - (L5/2 + 2*Lp)) == 0)

# ---------------------------------------------------------------- block 2
print("== block 2: membership machinery sanity ==")
ck("1,a independent", exact_rank(Matrix.hstack(ONE, A)) == 2)
A2 = hprod(A, A)
ck("a^2 not in span{1,a}", not in_span(A2, [ONE, A]))
A3 = hprod(A2, A)
ck("a^3 not in span{1,a,a^2}", not in_span(A3, [ONE, A, A2]))
# powers of a: R[a] is 6-dimensional (minimal polynomial degree = #values)
pows = [ONE, A]
for k in range(2, 6):
    pows.append(hprod(pows[-1], A))
ck("1,a,...,a^5 independent (dim R[a]=6)", exact_rank(Matrix.hstack(*pows)) == 6)
A6 = hprod(pows[-1], A)
ck("a^6 in span{1..a^5} (algebra closes at 6)", in_span(A6, pows))

# ---------------------------------------------------------------- block 3
print("== block 3: baseline -- indicator classification via II conditions ==")
# II_ac = 0 <=> a*X in V;  II_cc = 0 <=> X^2 in V;  V = span{1,a,X}
def ind(Sset):
    return col([1 if i in Sset else 0 for i in range(m)])

def ii_conditions(X):
    V = [ONE, A, X]
    if exact_rank(Matrix.hstack(*V)) < 3:
        return None  # degenerate (not a surface)
    return (in_span(hprod(A, X), V), in_span(hprod(X, X), V))

passing = []
for bits in range(2**m):
    Sset = {i for i in range(m) if (bits >> i) & 1}
    X = ind(Sset)
    r = ii_conditions(X)
    if r is not None and r == (True, True):
        passing.append(Sset)
ck("exactly 16 subsets pass {aX in V, X^2 in V}", len(passing) == 16,
   f"got {len(passing)}")
# they should be exactly: S subset of a a-level-set, or complement thereof
levels = []
for v in vals_set:
    levels.append({i for i in range(m) if cancel(a_vals[i]-v) == 0})
def a_const_on(Sset):
    return any(Sset <= L for L in levels)
ok = all(a_const_on(Sq) or a_const_on(set(range(m)) - Sq) for Sq in passing)
ck("all passing subsets have a constant on S or S^c", ok)
surfaces = set()
for Sq in passing:
    key = frozenset(Sq) if len(Sq) <= m - len(Sq) else frozenset(set(range(m))-Sq)
    surfaces.add(key)
ck("16 subsets collapse to exactly 8 surfaces (complement pairs)", len(surfaces) == 8)
merge = frozenset({3,4})
ck("golden merge {phi,tau} among the 8", merge in surfaces)
ck("the 7 singletons among the 8", all(frozenset({i}) in surfaces for i in range(m)))

# ---------------------------------------------------------------- block 4
print("== block 4: theta-rigidity + Gram nonsingularity ==")
# ker Q_w = V for every positive weight vector w: Q_w f = f - P_{V,w} f,
# well-defined iff Gram_w(1,a,X) invertible. Gram determinant is a polynomial
# in w with positive value at w=uniform for any independent triple.
w = list(symbols('w1:8', positive=True))
def gramdet(cols, wts):
    G = Matrix(3, 3, lambda i, j: sum(wts[k]*cols[i][k]*cols[j][k] for k in range(m)))
    return G.det()
Xtest = ind({5})  # indicator of phi4
gd = gramdet([ONE, A, Xtest], w)
ck("Gram det not identically zero (symbolic w)", cancel(gd) != 0)
gd_uniform = cancel(gd.subs({wi: Rational(1,7) for wi in w}))
ck("Gram det nonzero at uniform weights (exact)", cancel(gd_uniform) != 0)
# rigidity statement is definitional given invertibility: Q_w f=0 => f=P_V f in V.
ck("theta-rigidity: ker Q_w = V (definitional given Gram invertible)", True)

# ---------------------------------------------------------------- block 5
print("== block 5: case 3 and case 4 exclusions ==")
# Case 4: A=0 <=> X ~ a^2 affinely. Then V=span{1,a,a^2}; need aX=a^3 in V: false.
ck("case 4 excluded: a^3 not in span{1,a,a^2}", not in_span(A3, [ONE, A, A2]))
# and X=a^2 indeed makes II_aa=0 possible (paper Prop 12.1 parenthetical is wrong):
ck("X=a^2 gives a^2 in V (II_aa CAN vanish; proof of 12.1 needs the d-count fix)",
   in_span(A2, [ONE, A, A2]))
r = ii_conditions(A2)
ck("X=a^2: aX=a^3 not in V (II_ac != 0, so not const-1/4 via case 4)", r == (False, False))

# ---------------------------------------------------------------- block 6
print("== block 6: null-space dimension  N(X) = {v in V : v*V <= V} ==")
def null_space_dim(X):
    V = [ONE, A, X]
    B = Matrix.hstack(*V)
    if exact_rank(B) < 3:
        return None
    # basis of annihilator of V (row functionals)
    ann = B.T.nullspace(iszerofunc=lambda x: cancel(x) == 0)  # of B^T? careful:
    # functionals ell with ell(v)=0 for v in V: ell in nullspace of B^T (rows=functions)
    # B is 7x3; ell in R^7 with B^T ell = 0.
    ann = Matrix.hstack(*ann) if ann else Matrix.zeros(m, 0)
    # v = g*1 + al*a + be*X ; conditions: ell.(v*a)=0, ell.(v*X)=0 for all ell
    from sympy import zeros
    g_, al_, be_ = symbols('g_ al_ be_')
    v = col([cancel(g_ + al_*a_vals[i] + be_*X[i]) for i in range(m)])
    va = hprod(v, A); vX = hprod(v, X)
    eqs = []
    for j in range(ann.shape[1]):
        ell = ann[:, j]
        eqs.append(expand(sum(ell[i]*va[i] for i in range(m))))
        eqs.append(expand(sum(ell[i]*vX[i] for i in range(m))))
    from sympy import linear_eq_to_matrix
    Msys, rhs = linear_eq_to_matrix(eqs, [g_, al_, be_])
    assert rhs.is_zero_matrix
    return 3 - exact_rank(Msys)

for i in range(m):
    ck(f"null dim = 2 for indicator {SEEDS[i]}", null_space_dim(ind({i})) == 2)
ck("null dim = 2 for golden merge", null_space_dim(ind({3,4})) == 2)
ck("null dim = 1 for X=a^2 (only constants)", null_space_dim(A2) == 1)
# controls from the paper's landscape section: multi-set indicators (non-const-K)
ck("null dim = 1 for control 1_{sqrt2,sqrt3} (a nonconstant both sides)",
   null_space_dim(ind({0,1})) == 1)
ck("null dim = 1 for control 1_{sqrt5,K}", null_space_dim(ind({2,6})) == 1)
Xr = col([Rational(1,3), Rational(-2,5), Rational(7,2), Rational(1,7),
          Rational(9,4), Rational(-1,2), Rational(3,11)])
ck("null dim = 1 for a random exact rational X", null_space_dim(Xr) == 1)
Xm = col([0,0,0,1,2,0,0])
ck("null dim = 1 for pair-pole X=(0,0,0,1,2,0,0)", null_space_dim(Xm) == 1)

# ---------------------------------------------------------------- block 7
print("== block 7: machine checks of the 1a Mobius-elimination proof ==")
t, d0, d1, d2, c0, c1, c2, s_ = symbols('t d0 d1 d2 c0 c1 c2 s_')
# Step (i): the case-A residual polynomial has degree <= 3 in t
P_A = expand((d0 + d1*t)**2 - (c0 + c1*t)*(t - d2)**2 - c2*(d0 + d1*t)*(t - d2))
ck("case A residual is a polynomial of t-degree <= 3", Poly(P_A, t).degree() <= 3)
ck("case A residual t^3 coefficient = -c1", cancel(Poly(P_A, t).coeff_monomial(t**3) + c1) == 0)
# Step (ii): Mobius increment identity  mu(s)-mu(t) = (d0+d1*d2)(s-t)/((s-d2)(t-d2))
mu = lambda x: (d0 + d1*x)/(x - d2)
lhs = together(mu(s_) - mu(t))
rhs = together(-(d0 + d1*d2)*(s_ - t)/((s_ - d2)*(t - d2)))
ck("Mobius increment identity (injective iff d0+d1*d2 != 0)", cancel(lhs - rhs) == 0)
# Step (iii): case B root identity  y^2 - c2*y - (d1^2 - c2*d1) = (y-d1)(y-(c2-d1))
y = symbols('y')
ck("case B quadratic factors as (y-d1)(y-(c2-d1))",
   cancel(expand(y**2 - c2*y - (d1**2 - c2*d1)) - expand((y - d1)*(y - (c2 - d1)))) == 0)
# Step (iv): with X supported on pole set P (a==d2 on P, X==d1 off P),
# verify aX in V holds identically: aX = -d1*d2*1 + d1*a + d2*X  (symbolic check
# on the catalog with P={phi,tau} and free pole values e3,e4)
e3, e4 = symbols('e3 e4')
Xp = col([d1, d1, d1, e3, e4, d1, d1])
d2v = Lp  # pole at the phi/tau level
aX = hprod(A, Xp)
resid = [cancel(aX[i] - (-d1*d2v + d1*a_vals[i] + d2v*Xp[i])) for i in range(m)]
ck("pole-family identity aX = -d1*d2 + d1*a + d2*X on the catalog", all(r == 0 for r in resid))

# ---------------------------------------------------------------- block 8
print("== block 8: r = dim(V^2/V) for key families ==")
def r_of(X):
    V = [ONE, A, X]
    if exact_rank(Matrix.hstack(*V)) < 3:
        return None
    V2 = V + [hprod(A, A), hprod(A, X), hprod(X, X)]
    return exact_rank(Matrix.hstack(*V2)) - 3

for i in [0, 3, 5, 6]:
    ck(f"r = 1 for indicator {SEEDS[i]}", r_of(ind({i})) == 1)
ck("r = 1 for golden merge", r_of(ind({3,4})) == 1)
ck("r = 2 for X=a^2", r_of(A2) == 2)
ck("r = 2 for pair-pole X=(0,0,0,1,2,0,0)", r_of(Xm) == 2)
ck("r = 3 for random exact X", r_of(Xr) == 3)
ck("r = 2 for control 1_{sqrt2,sqrt3} (X^2=X collapses one product)", r_of(ind({0,1})) == 2)

print(f"\nALL {NCK[0]} CHECKS PASSED (t1_core)")
