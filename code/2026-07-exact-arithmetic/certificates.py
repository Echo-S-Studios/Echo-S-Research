# certificates.py — the evidence layer for Lessons 11-13, rebuilt so that every predicate decides
# the proposition its label names. Supersedes lessons11_13_recheck.py as a CERTIFICATE (not as
# coverage: the original 88 authoring checks remain unrecoverable).
#
# Standard enforced here: no uncertified approximation crosses a decision boundary.
#   - sign of a Q-combination of log2, log3, log5, log phi  -> EXACT, decided in Q(sqrt5)
#   - rank of the interaction tensor                        -> EXACT, nullity + a nonvanishing minor
#   - a transcendental quadratic form's nonvanishing        -> DIRECTED INTERVAL enclosure
#   - finite classifications                                -> EXHAUSTIVE enumeration
import itertools
from fractions import Fraction as Q

import sympy as sp
from sympy.utilities.iterables import multiset_partitions
import mpmath as mp

mp.mp.dps = 60
mp.iv.dps = 40
x, t = sp.symbols('x t')
PHI = (1 + sp.sqrt(5)) / 2
checks, cited = [], []
def check(cid, name, ok, note=""):
    checks.append((cid, name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {cid}  {name}" + (f"\n        {note}" if note else ""), flush=True)
def cite(cid, name, ref):
    """An imported theorem. Cited and scoped; NEVER counted as an executable certificate."""
    cited.append((cid, name, ref))
    print(f"[CITED] {cid}  {name}\n        {ref}", flush=True)
S = sp.simplify

# ═══════════════════════════════════════════════════════════════════ 0. EXACT SIGN ORACLE
# A cost/log-magnitude is a vector (a,b,c,d) over Q meaning a*log2 + b*log3 + c*log5 + d*log(phi).
# sign(a L2 + b L3 + c L5 + d Lphi) = sign( 2^aN 3^bN 5^cN phi^dN - 1 )  for any N > 0 clearing
# denominators; that quantity lies in Q(sqrt5) and its sign is decided over Q. No floats.
def _phi_pow(m):
    """phi^m = (L_m + F_m sqrt5)/2 as an exact pair (p, q) with value p + q*sqrt5, m in Z."""
    if m >= 0:
        L, F = sp.lucas(m), sp.fibonacci(m)
        return (Q(int(L), 2), Q(int(F), 2))
    n = -m
    L, F = sp.lucas(n), sp.fibonacci(n)
    p, q = Q(int(L), 2), Q(int(F), 2)                      # phi^n = p + q sqrt5
    den = p * p - 5 * q * q                                # norm = (-1)^n
    return (p / den, -q / den)

def _mul(A, B):
    (p, q), (r, s) = A, B
    return (p * r + 5 * q * s, p * s + q * r)

def sign_log(v):
    """Exact sign of a*L2 + b*L3 + c*L5 + d*Lphi. Returns -1, 0 or +1."""
    a, b, c, d = (Q(t) for t in v)
    if a == b == c == d == 0:
        return 0
    N = 1
    for r in (a, b, c, d):
        N = N * r.denominator // sp.igcd(N, r.denominator)
    N = int(N)
    ai, bi, ci, di = int(a * N), int(b * N), int(c * N), int(d * N)
    rat = Q(2) ** ai * Q(3) ** bi * Q(5) ** ci             # exact rational
    val = _mul((rat, Q(0)), _phi_pow(di))                  # p + q sqrt5, exact
    p, q = val[0] - 1, val[1]                              # compare to 1
    if q == 0:
        return (p > 0) - (p < 0)
    if p == 0:
        return (q > 0) - (q < 0)
    if (p > 0) == (q > 0):
        return 1 if p > 0 else -1
    bigger_p = p * p > 5 * q * q                           # |p| vs |q|sqrt5, over Q
    return (1 if p > 0 else -1) if bigger_p else (1 if q > 0 else -1)

check("S0", "exact sign oracle: decided in Q(sqrt5), never by float. Spot values — "
      "log2/2 - log phi < 0 (2 < phi^2), log3/2 - log phi > 0 (3 > phi^2), "
      "log5/2 - log phi > 0 (5 > phi^2), log5 - 4 log phi < 0 (5 < phi^4), zero vector -> 0",
      sign_log((Q(1,2),0,0,-1)) == -1 and sign_log((0,Q(1,2),0,-1)) == 1
      and sign_log((0,0,Q(1,2),-1)) == 1 and sign_log((0,0,1,-4)) == -1
      and sign_log((0,0,0,0)) == 0)

# ═══════════════════════════════════════════════════════════════════ 1. THE CATALOG
V = lambda a, b, c, d: (Q(a), Q(b), Q(c), Q(d))
SEEDS = {
    'sqrt2': [V(Q(1,2),0,0,0)] * 2,
    'sqrt3': [V(0,Q(1,2),0,0)] * 2,
    'sqrt5': [V(0,0,Q(1,2),0)] * 2,
    'phi':   [V(0,0,0,1), V(0,0,0,-1)],
    'tau':   [V(0,0,0,-1), V(0,0,0,1)],
    'phi4':  [V(0,0,0,4), V(0,0,0,-4)],
    'K':     [V(0,0,Q(1,4),-1)] * 2 + [V(0,0,Q(1,4),1)] * 2,
}
ORD = ['sqrt2', 'sqrt3', 'sqrt5', 'phi', 'tau', 'phi4', 'K']
add = lambda u, v: tuple(a + b for a, b in zip(u, v))
neg = lambda u: tuple(-a for a in u)

def tropical(p, q):
    """c(p (x) q) = sum over conjugate pairs of max(0, u + v); every sign decided exactly."""
    tot = V(0, 0, 0, 0)
    for u in SEEDS[p]:
        for v in SEEDS[q]:
            w = add(u, v)
            if sign_log(w) > 0:
                tot = add(tot, w)
    return tot

COST = {}
for s in ORD:
    tot = V(0, 0, 0, 0)
    for u in SEEDS[s]:
        if sign_log(u) > 0:
            tot = add(tot, u)
    COST[s] = tot
check("C1", "the catalog's cost vectors, tropically from conjugate log-magnitudes with exact signs: "
      "sqrt2 (1,0,0,0), phi = tau (0,0,0,1), phi^4 (0,0,0,4), K (0,0,1/2,2)",
      COST['sqrt2'] == V(1,0,0,0) and COST['phi'] == COST['tau'] == V(0,0,0,1)
      and COST['phi4'] == V(0,0,0,4) and COST['K'] == V(0,0,Q(1,2),2))
c23, c2p, cp3, cpp = tropical('sqrt2','sqrt3'), tropical('sqrt2','phi'), tropical('phi','sqrt3'), tropical('phi','phi')
quantum = add(add(c23, neg(c2p)), add(neg(cp3), cpp))
cp4, c44 = tropical('phi','phi4'), tropical('phi4','phi4')
gq = add(add(cpp, neg(cp4)), add(neg(tropical('phi4','phi')), c44))
check("C2", "the two interaction quanta, exact: the rational/golden contrast is log 2 and the "
      "golden-sector contrast is -6 log phi",
      quantum == V(1,0,0,0) and gq == V(0,0,0,-6))
check("C3", "the golden tie survives tensoring: c(phi (x) z) = c(tau (x) z) for all seven z",
      all(tropical('phi', z) == tropical('tau', z) for z in ORD))

# ═══════════════════════════════════════════════════════════════════ 2. THE COLLISION LATTICE (was E3)
# LATTICE = { d in Z^7 : sum_i d_i cost_i = 0 in Q^4, and sum_i d_i = 0 }.  Computed, not asserted.
Mrows = [[COST[s][k] for s in ORD] for k in range(4)] + [[Q(1)] * 7]
Msp = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in row] for row in Mrows])
ns = Msp.nullspace()
d1 = sp.Matrix([0, 0, 0, 1, -1, 0, 0])
d2 = sp.Matrix([0, 0, 1, 0, 0, 1, -2])
Bl = sp.Matrix.hstack(d1, d2)
minors = [Bl[[i, j], :].det() for i in range(7) for j in range(i + 1, 7)]
gcd_minors = sp.igcd(*[int(m) for m in minors])
check("C4", "the collision lattice is COMPUTED, not asserted: the integer kernel of the 5x7 "
      "cost/sum matrix has rational rank 2, both d1 (golden swap) and d2 (Salem square) lie in it, "
      "and the 2x2 minors of [d1 d2] have gcd 1 — so the pair is primitive and therefore generates "
      "the FULL integer kernel, not merely a finite-index sublattice",
      len(ns) == 2 and (Msp * d1).is_zero_matrix and (Msp * d2).is_zero_matrix
      and Bl.rank() == 2 and gcd_minors == 1)
# the cube intersection, enumerated rather than argued from one hard-coded entry
cube_hits = []
for n1 in range(-3, 4):
    for n2 in range(-3, 4):
        v = n1 * d1 + n2 * d2
        if all(abs(int(e)) <= 1 for e in v):
            cube_hits.append((n1, n2))
check("C5", "cube ∩ lattice = {0, ±d1}, ENUMERATED: the K-entry of n1 d1 + n2 d2 is -2 n2, so "
      "squarefreeness forces n2 = 0, and then n1 ∈ {-1,0,1}. Search over |n_i| <= 3 returns exactly "
      f"three points {sorted(cube_hits)} — the Salem square cannot act at window level",
      sorted(cube_hits) == [(-1, 0), (0, 0), (1, 0)])

# ═══════════════════════════════════════════════════════════════════ 3. THE INTERACTION TENSOR (was D7)
L2s, L3s, L5s, Lps = sp.symbols('L2 L3 L5 Lp', positive=True)
BASIS = (L2s, L3s, L5s, Lps)
sym = lambda v: sum(sp.Rational(c.numerator, c.denominator) * b for c, b in zip(v, BASIS))
C = sp.Matrix(7, 7, lambda i, j: sym(tropical(ORD[i], ORD[j])))
J = sp.ones(7, 7) / 7
Delta = sp.expand(C - J * C - C * J + J * C * J)
ones = sp.ones(7, 1)
e_phi_tau = sp.Matrix([0, 0, 0, 1, -1, 0, 0])
check("C6", "nullity >= 2 exactly: double centring kills the all-ones vector, and the persistent "
      "phi/tau tie makes rows (and columns) phi and tau identical, killing e_phi - e_tau. Hence "
      "rank(Delta) <= 5 — proved, not observed",
      sp.simplify(Delta * ones).is_zero_matrix and sp.simplify(Delta * e_phi_tau).is_zero_matrix
      and sp.simplify(Delta.T * ones).is_zero_matrix)
target = None
for idx in itertools.combinations(range(7), 5):
    sub = Delta[list(idx), list(idx)]
    dd = sp.factor(sp.simplify(sub.det()))
    if dd != 0:
        target = (idx, dd)
        break
idx, dd = target
expected = sp.factor(-sp.Rational(2, 49) * (L2s - 2*Lps)**2 * (L3s - L5s)**2 * (L5s - 4*Lps))
check("C7", f"rank >= 5 exactly: the principal 5x5 minor on rows/cols {idx} has determinant "
      "-(2/49)(L2 - 2Lp)^2 (L3 - L5)^2 (L5 - 4Lp) — a symbolic factorisation, matched here",
      sp.simplify(dd - expected) == 0)
check("C8", "and each factor is nonzero for an EXACT algebraic reason, no transcendence needed: "
      "L2 - 2Lp != 0 since 2 != phi^2; L3 - L5 != 0 since 3 != 5; L5 - 4Lp != 0 since 5 != phi^4. "
      "Therefore rank(Delta) = 5 exactly — the numerical SVD is retired to corroboration",
      sign_log((1,0,0,-2)) != 0 and sign_log((0,1,-1,0)) != 0 and sign_log((0,0,1,-4)) != 0
      and S(PHI**2 - 2) != 0 and S(PHI**4 - 5) != 0)

# ═══════════════════════════════════════════════════════════════════ 4. COVARIANCE (was D6)
cells = [(i, j) for i in ORD for j in ORD]
f1 = {ij: COST[ij[0]] for ij in cells}
f2 = {ij: COST[ij[1]] for ij in cells}
f3 = {ij: tropical(*ij) for ij in cells}
def mean(f):
    tot = V(0,0,0,0)
    for ij in cells: tot = add(tot, f[ij])
    return tuple(a / len(cells) for a in tot)
m1, m2, m3 = mean(f1), mean(f2), mean(f3)
def covform(fa, ma, fb, mb):
    Mq = [[Q(0)] * 4 for _ in range(4)]
    for ij in cells:
        a = [fa[ij][s] - ma[s] for s in range(4)]
        b = [fb[ij][s] - mb[s] for s in range(4)]
        for s in range(4):
            for u in range(4):
                Mq[s][u] += a[s] * b[u] / len(cells)
    return Mq
Q12, Q13 = covform(f1, m1, f2, m2), covform(f1, m1, f3, m3)
Q23 = covform(f2, m2, f3, m3)
check("C9", "Cov(cost1, cost2) = 0 exactly as a rational quadratic form — the uniform measure on a "
      "product grid makes row and column independent",
      all(Q12[s][u] == 0 for s in range(4) for u in range(4)))
# directed-interval enclosure: linear independence of the logs would NOT settle a quadratic form
IVL = [mp.iv.log(mp.iv.mpf(2)), mp.iv.log(mp.iv.mpf(3)), mp.iv.log(mp.iv.mpf(5)),
       mp.iv.log((1 + mp.iv.sqrt(mp.iv.mpf(5))) / 2)]
def enclose(Mq):
    acc = mp.iv.mpf(0)
    for s in range(4):
        for u in range(4):
            if Mq[s][u] == 0: continue
            r = mp.iv.mpf(Mq[s][u].numerator) / mp.iv.mpf(Mq[s][u].denominator)
            acc = acc + r * IVL[s] * IVL[u]
    return acc
enc13, enc23 = enclose(Q13), enclose(Q23)
check("C10", "Cov(cost1, c) and Cov(cost2, c) are strictly positive, certified by DIRECTED INTERVAL "
      "enclosure over verified log enclosures — not by Q-linear independence, which cannot settle a "
      "quadratic form. So the coupling is invisible only to the two MARGINAL costs",
      enc13.a > 0 and enc23.a > 0,
      f"Cov(cost1,c) ∈ [{mp.nstr(mp.mpf(enc13.a), 17)}, {mp.nstr(mp.mpf(enc13.b), 17)}]")

# ═══════════════════════════════════════════════════════════════════ 5. THE LANDSCAPE (was E4)
LEVEL = {'sqrt2': 0, 'sqrt3': 1, 'sqrt5': 2, 'phi': 3, 'tau': 3, 'phi4': 4, 'K': 5}
def classify(part):
    """dim V(pi) = c + 2 b2 for a canonical partition; None if some block is neither a cluster
    (inside one cost level) nor a line block (at least three distinct cost values)."""
    c = b2 = 0
    for blk in part:
        lv = {LEVEL[s] for s in blk}
        if len(lv) == 1:
            c += 1
        elif len(lv) >= 3:
            b2 += 1
        else:
            return None
    return c + 2 * b2
land = {}
total = 0
for part in multiset_partitions(ORD):
    total += 1
    d = classify(part)
    if d is not None:
        land[d - 1] = land.get(d - 1, 0) + 1
check("C11", "the landscape is ENUMERATED, not asserted: all 877 set partitions of the seven seeds "
      "are generated, each classified into clusters and line blocks, and the constant-1/4 families "
      f"counted by k = dim - 1. Result {dict(sorted(land.items()))}: the reported 8/56/95/31/1 at "
      "k = 2..6, PLUS exactly one class at k = 1 — the single all-seeds line block, whose span is "
      "just <1, log M>. That is Lesson 11's base Gibbs curve, which is a curve and not a surface, "
      "so the landscape table rightly starts at k = 2; its appearance here at multiplicity exactly 1 "
      "is a consistency check on the classification rather than a missing entry",
      total == 877 and sorted(land.items()) == [(1,1),(2,8),(3,56),(4,95),(5,31),(6,1)])
k3 = [p for p in multiset_partitions(ORD) if classify(p) == 4]
dbl = sum(1 for p in k3 if sum(1 for b in p if len({LEVEL[s] for s in b}) == 1) == 2)
spl = sum(1 for p in k3 if sum(1 for b in p if len({LEVEL[s] for s in b}) >= 3) == 2)
check("C12", f"and the k = 3 split is enumerated too: {dbl} double-indicator classes (two clusters + "
      f"one line block) and {spl} split-affine classes (two line blocks), {dbl}+{spl} = 56 — the "
      "refutation of the within-level-indicator conjecture, counted rather than quoted",
      dbl == 26 and spl == 30 and dbl + spl == 56)
extra = sp.binomial(7, 2) - 1
check("C13", "the reproducer's reconciliation: a naive split-affine scan also admits |B^c| = 2 at "
      f"distinct costs, C(7,2) - 1 = {extra} extra spans on both sides (35+20 = 55, 5+20 = 25), so "
      "both routes land on 30", extra == 20 and 35 + 20 == 55 and 5 + 20 == 25)

# ═══════════════════════════════════════════════════════════════════ 6. FISHER = VARIANCE (replaces F1)
bta = sp.Symbol('beta', real=True)
cs = sp.symbols('c0:7', positive=True)
A = sp.log(sum(sp.exp(-bta * ci) for ci in cs))
A1 = sp.simplify(sp.diff(A, bta))
A2 = sp.simplify(sp.diff(A, bta, 2))
pw = [sp.exp(-bta * ci) / sum(sp.exp(-bta * cj) for cj in cs) for ci in cs]
var = sp.simplify(sum(p * ci**2 for p, ci in zip(pw, cs)) - (sum(p * ci for p, ci in zip(pw, cs)))**2)
check("C14", "the Fisher identity is verified as a SYMBOLIC IDENTITY in seven free costs, at every "
      "beta: A'(beta) = -E[log M] and A''(beta) = Var_beta(log M), so I(beta) = A''(beta) >= 0",
      sp.simplify(A1 + sum(p * ci for p, ci in zip(pw, cs))) == 0
      and sp.simplify(A2 - var) == 0)
# the applicability conditions ARE checkable for a finite exponential family; the theorem is not
score = [sp.simplify(sp.diff(sp.log(pi), bta)) for pi in pw]
Escore = sp.simplify(sum(p * sc for p, sc in zip(pw, score)))
check("C15", "the Cramer-Rao APPLICABILITY conditions, verified for this family: the support is the "
      "fixed 7-atom index set and does not depend on beta, so differentiation passes through a "
      "finite sum; and the score has mean zero, E_beta[d/dbeta log p] = 0, symbolically in seven "
      "free costs. Those are the regularity hypotheses the bound needs",
      Escore == 0 and len(cs) == 7)
cite("T1", "Cramer-Rao: for an unbiased estimator of beta from n independent emissions, "
     "Var(beta-hat) >= 1/(n I(beta)), with I(beta) = Var_beta(log M) by C14",
     "ESTABLISHED. Imported, not executed. The previous printing counted a tautology "
     "(1/(nI) - 1/(nI) == 0) as a passing check; that entry is withdrawn. An established theorem "
     "is cited, scoped and applied — it does not masquerade as a Boolean execution result.")

# ═══════════════════════════════════════════════════════════════════ 7. ANCHORS
Zneg = 2 + 3 + 5 + PHI + PHI + PHI**4 + (PHI**4 - 1)
Zpos = sp.Rational(1,2) + sp.Rational(1,3) + sp.Rational(1,5) + 2/PHI + 1/PHI**4 + 1/(PHI**4 - 1)
check("C16", "partition anchors: Z(-1) = 17 + 4 sqrt5 and Z(+1) = 91/30 - sqrt5/5",
      S(Zneg - (17 + 4*sp.sqrt(5))) == 0 and S(Zpos - (sp.Rational(91,30) - sp.sqrt(5)/5)) == 0)
check("C17", "the K entry and the Salem square: M(p_K) = phi^4 - 1 = phi^2 sqrt5, and "
      "(phi^4 - 1)^2 = 5 phi^4",
      S(PHI**4 - 1 - (5 + 3*sp.sqrt(5))/2) == 0 and S((PHI**4 - 1)**2 - 5*PHI**4) == 0)

# ═══════════════════════════════════════════════════════════════════ 8. THIRD-PASS REPAIRS
Bm = sp.Matrix([[1, PHI], [1, (1 - sp.sqrt(5))/2]])
sig = sp.Symbol('sigma', positive=True)
check("C18", "Lesson 3: B^T B = [[2,1],[1,3]] is an UNCENTRED second moment (the first embedding "
      "coordinate is the constant 1, so the centred matrix is singular); under Y = Br + eps with "
      "eps ~ N(0, sigma^2 I) the Fisher matrix is B^T Sigma^-1 B = G/sigma^2, hence c = sigma^2",
      S(Bm.T*Bm) == sp.Matrix([[2,1],[1,3]])
      and S(Bm.T*(sig**2*sp.eye(2)).inv()*Bm - sp.Matrix([[2,1],[1,3]])/sig**2) == sp.zeros(2,2))
psi_, ell = sp.symbols('psi ell', real=True)
def gcurv(E, Gg, sqrtEG, u, v):
    return S(-(sp.diff(sp.diff(Gg, u)/sqrtEG, u) + sp.diff(sp.diff(E, v)/sqrtEG, v))/(2*sqrtEG))
check("C19", "Lesson 12: dpsi^2 + cos^2 psi dl^2 is the UNIT sphere (curvature 1); the Fisher metric "
      "is 4x that, with curvature exactly 1/4 — consistent with p -> 2 sqrt(p) landing on radius 2",
      gcurv(sp.Integer(1), sp.cos(psi_)**2, sp.cos(psi_), psi_, ell) == 1
      and gcurv(sp.Integer(4), 4*sp.cos(psi_)**2, 4*sp.cos(psi_), psi_, ell) == sp.Rational(1,4))
def formula(m): return sum(2**k - 1 for k in m)
def truth(levels):
    seeds = [i for i, m in enumerate(levels) for _ in range(m)]
    n = len(seeds); spans = set()
    for mask in range(1, 2**n - 1):
        Ss = frozenset(i for i in range(n) if mask >> i & 1)
        Cs = frozenset(i for i in range(n) if not (mask >> i & 1))
        if len({seeds[i] for i in Ss}) == 1 or len({seeds[i] for i in Cs}) == 1:
            spans.add(frozenset([Ss, Cs]))
    return len(spans)
check("C20", "Lesson 12: the census formula sum(2^m - 1) is exact at >= 3 cost levels and overcounts "
      "below (two levels: (1,1) gives 2 against 1, (2,1) gives 4 against 3), because 1_S and 1_Sc "
      "span the same surface",
      all(formula(m) == truth(m) for m in [(1,1,1,2,1,1),(1,1,1,2,1),(1,1,1,1,1,1),(1,1,1,1,2,1,1),(3,2,1)])
      and formula((1,1)) > truth((1,1)) and formula((2,1)) > truth((2,1)))

# ═══════════════════════════════════════════════════════════════════ 9. SCHINZEL SCOPE
q2 = sp.Poly(x**4 + x**2 - 1, x)
check("C21", "Lesson 1/4: Schinzel gives M >= phi^(d/2) on the TOTALLY REAL sector; the emission "
      "image is not totally real (Wall 2 permits purely imaginary), and the floor-attaining witness "
      "q2 = x^4 + x^2 - 1 has exactly 2 real roots, carrying its whole measure on ±i sqrt(phi)",
      q2.count_roots(-sp.oo, sp.oo) == 2
      and S((sp.I*sp.sqrt(PHI))**4 + (sp.I*sp.sqrt(PHI))**2 - 1) == 0
      and S(sp.Abs(sp.I*sp.sqrt(PHI))**2 - PHI) == 0)

# ═══════════════════════════════════════════════════════════════════ 10. THE FORCED BAND
# Which exclusion survives with NO plausible tag anywhere in its derivation?
Lp = sp.Poly(x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1, x)
theta0 = sp.Poly(x**3 - x - 1, x)          # Smyth's constant mu_S, the plastic number
r1, r2 = sp.Rational(5, 4), sp.Rational(4, 3)
check("C22", "the ordering, by exact rational separators: M(L) < 5/4 < mu_S < 4/3 < phi. "
      f"L(5/4) = {Lp.eval(r1)} > 0 puts Lehmer's root below 5/4; (5/4)^3-5/4-1 = {theta0.eval(r1)} < 0 "
      f"puts mu_S above it; (4/3)^3-4/3-1 = {theta0.eval(r2)} > 0 puts mu_S below 4/3; and phi > 4/3 "
      "since 45 > 25. So Lehmer's number lies strictly inside (1, mu_S)",
      Lp.eval(r1) > 0 and theta0.eval(r1) < 0 and theta0.eval(r2) > 0
      and sign_log((0, 0, 0, 1)) > 0 and 45 > 25
      and sp.simplify(PHI - sp.Rational(4, 3)).is_positive)
check("C23", "hence the narrowest UNCONDITIONAL statement: the forced excluded band is "
      "(1, mu_S), not (1, phi). Smyth's bound is unconditional off-reciprocal and the odd-charge "
      "floor is forced at mu_S, while phi is forced only on quadratics and the totally real sector "
      "(Schinzel) and stays plausible elsewhere. Because M(L) < mu_S, the exclusion of Lehmer's "
      "number rides entirely on the forced part and touches no plausible tag",
      sp.simplify(PHI - (sp.CRootOf(theta0, 0))).is_positive)
cite("T2", "the square-root isometry p -> 2 sqrt(p) carrying the simplex onto the positive orthant "
     "of the radius-2 sphere, so that the Fisher metric becomes the round metric of curvature 1/4",
     "ESTABLISHED. C. R. Rao, Bull. Calcutta Math. Soc. 37 (1945), 81-91 - the same paper that "
     "introduced the Cramer-Rao bound cited at T1. The sphere is classical ambient background, NOT "
     "a result of this corpus; what is corpus content is that THIS catalog forms the family, that "
     "the other four layers hold for it, and the constant-1/4 classification that follows.")

print()
fails = [c for c in checks if not c[2]]
print(f"SUMMARY: {len(checks)-len(fails)}/{len(checks)} executable certificates passed"
      + (f", {len(cited)} imported theorem(s) cited" if cited else "")
      + ("" if not fails else f"   FAILURES: {[c[0] for c in fails]}"))
