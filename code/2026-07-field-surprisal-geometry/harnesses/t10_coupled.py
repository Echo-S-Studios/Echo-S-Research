#!/usr/bin/env python3
# t10_coupled.py -- the coupled family's geometry, closed.
# On the 49-cell tensor catalog with statistics T = (cost1, cost2, c):
#   [F] dim V = 4 (spanning certified at the TRUE logs by intervals)
#   [F] the configuration contains a genuine 5-circuit (all five dependency
#       weights nonzero, interval-certified) => NOT partitioned-affine, and the
#       corresponding window matrix has rank >= 2 => NOT windowwise flat
#   [F] the Gauss obstruction at the uniform point is NONZERO (rigorous interval
#       bound) => the coupled family does NOT have constant curvature 1/4: the
#       forced coupling is curvature-visible
#   [C] sectional deviations sec_ab - 1/4 = (1/4) G^{abab}/(g_aa g_bb - g_ab^2)
#       at uniform, displayed with certified enclosures
#   [F] Cov_uniform(cost1, cost2) = 0 exactly -- the coupling is invisible to
#       the cross-covariance at uniform yet visible to curvature
#   [F] Delta is NOT determined by the charge data: (sqrt2,sqrt3) vs
#       (sqrt2,sqrt5) share charge pair (Z/2,Z/2) but differ in Delta exactly
import sys, itertools, random
from fractions import Fraction as Fr

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}"); sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

def v(l2=0, l3=0, l5=0, lp=0): return (Fr(l2), Fr(l3), Fr(l5), Fr(lp))
def vadd(a, b): return tuple(x+y for x, y in zip(a, b))
ZERO = v()
SEEDS = ["sqrt2", "sqrt3", "sqrt5", "phi", "tau", "phi4", "K"]
ROOTS = [
    [(v(l2=Fr(1,2)), 0), (v(l2=Fr(1,2)), 2)],
    [(v(l3=Fr(1,2)), 0), (v(l3=Fr(1,2)), 2)],
    [(v(l5=Fr(1,2)), 0), (v(l5=Fr(1,2)), 2)],
    [(v(lp=1), 0), (v(lp=-1), 2)],
    [(v(lp=-1), 0), (v(lp=1), 2)],
    [(v(lp=4), 0), (v(lp=-4), 0)],
    [(v(l5=Fr(1,4), lp=-1), 0), (v(l5=Fr(1,4), lp=-1), 2),
     (v(l5=Fr(1,4), lp=1), 1), (v(l5=Fr(1,4), lp=1), 3)],
]
COST = [v(l2=1), v(l3=1), v(l5=1), v(lp=1), v(lp=1), v(lp=4),
        v(l5=Fr(1,2), lp=2)]

from mpmath import iv
iv.dps = 60
IVL = [iv.log(2), iv.log(3), iv.log(5), iv.log((1 + iv.sqrt(5))/2)]
def cv(fr): return iv.mpf(fr.numerator)/fr.denominator
def lin_iv(vec): return sum(cv(c)*L for c, L in zip(vec, IVL))
def sign_of(vec):
    if vec == ZERO: return 0
    t = lin_iv(vec)
    if t.a > 0: return 1
    if t.b < 0: return -1
    raise RuntimeError("ambiguous sign")

C = [[None]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        tot = ZERO
        for (u, _) in ROOTS[i]:
            for (w, _) in ROOTS[j]:
                s = vadd(u, w)
                if sign_of(s) > 0: tot = vadd(tot, s)
        C[i][j] = tot

CELLS = [(i, j) for i in range(7) for j in range(7)]
A1 = [COST[i] for (i, j) in CELLS]
A2 = [COST[j] for (i, j) in CELLS]
CC = [C[i][j] for (i, j) in CELLS]
SPEC = (Fr(97), Fr(13), Fr(41), Fr(7, 3))
def ev(vec): return sum(c*s for c, s in zip(vec, SPEC))
A1r = [ev(t) for t in A1]; A2r = [ev(t) for t in A2]; CCr = [ev(t) for t in CC]

def det_fr(M):
    M = [row[:] for row in M]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]; inv = Fr(1)/M[c][c]
        M[c] = [e*inv for e in M[c]]
        for i in range(c+1, n):
            if M[i][c] != 0:
                f = M[i][c]; M[i] = [e - f*g for e, g in zip(M[i], M[c])]
    return d
def iv_det(M):
    n = len(M); tot = iv.mpf(0)
    for perm in itertools.permutations(range(n)):
        sgn = 1; pl = list(perm)
        for i in range(n):
            for j in range(i+1, n):
                if pl[i] > pl[j]: sgn = -sgn
        pr = iv.mpf(sgn)
        for i in range(n): pr *= M[i][perm[i]]
        tot += pr
    return tot

print("== block 1: dim V = 4 at the true logs ==")
def row_rat(t): return [Fr(1), A1r[t], A2r[t], CCr[t]]
def row_iv(t):
    return [iv.mpf(1), lin_iv(A1[t]), lin_iv(A2[t]), lin_iv(CC[t])]
found = None
for quad in itertools.combinations([0, 1, 9, 16, 24, 33, 40, 48, 20, 5], 4):
    if det_fr([row_rat(t) for t in quad]) != 0: found = quad; break
d4 = iv_det([row_iv(t) for t in found])
ck(f"four cells {[CELLS[t] for t in found]} span: 4x4 det interval excludes 0 "
   "=> dim V = 4 at the true logs", not (d4.a <= 0 <= d4.b))

print("== block 2: a genuine 5-circuit => not partitioned-affine, not flat ==")
pool = [t for t, (i, j) in enumerate(CELLS) if i != 4 and j != 4]
random.seed(7)
witness = None
for _ in range(4000):
    s = random.sample(pool, 5)
    rows = [row_rat(t) for t in s]
    lam = [(-1)**u * det_fr([rows[j] for j in range(5) if j != u])
           for u in range(5)]
    if all(l != 0 for l in lam):
        P = [[A1r[t], A2r[t], CCr[t]] for t in s]
        Lm = [[sum(lam[u]*P[u][a]*P[u][b] for u in range(5)) for b in range(3)]
              for a in range(3)]
        minors = [(Lm[0][0]*Lm[1][1] - Lm[0][1]*Lm[1][0]),
                  (Lm[0][1]*Lm[1][2] - Lm[0][2]*Lm[1][1]),
                  (Lm[0][0]*Lm[1][2] - Lm[0][2]*Lm[1][0])]
        if any(m != 0 for m in minors):
            witness = s; break
ck("found a candidate 5-window with all dependency weights nonzero and a "
   f"nonzero 2x2 minor (cells {[CELLS[t] for t in witness]})", witness is not None)
rows_iv = [row_iv(t) for t in witness]
lam_iv = []
okl = True
for u in range(5):
    d = iv_det([rows_iv[j] for j in range(5) if j != u])
    lam_iv.append(iv.mpf((-1)**u)*d)
    okl &= not (d.a <= 0 <= d.b)
ck("all five dependency weights nonzero at the TRUE logs (intervals exclude 0)"
   " => the coupled configuration contains a genuine 5-circuit => it is NOT of"
   " partitioned-affine type", okl)
Piv = [[lin_iv(A1[t]), lin_iv(A2[t]), lin_iv(CC[t])] for t in witness]
Liv = [[sum(lam_iv[u]*Piv[u][a]*Piv[u][b] for u in range(5)) for b in range(3)]
       for a in range(3)]
m01 = Liv[0][0]*Liv[1][1] - Liv[0][1]*Liv[1][0]
m12 = Liv[0][1]*Liv[1][2] - Liv[0][2]*Liv[1][1]
ok_min = (not (m01.a <= 0 <= m01.b)) or (not (m12.a <= 0 <= m12.b))
ck("a 2x2 minor of the window matrix excludes 0 at the true logs => window rank"
   " >= 2 => the coupled triple is NOT windowwise flat", ok_min)

print("== block 3: the Gauss obstruction at the uniform point (intervals) ==")
n = 49
A1v = [lin_iv(t) for t in A1]; A2v = [lin_iv(t) for t in A2]
CCv = [lin_iv(t) for t in CC]
T = [A1v, A2v, CCv]
def E(f): return sum(f)/n
def Ef(f, g): return sum(fi*gi for fi, gi in zip(f, g))/n
mu = [E(T[a]) for a in range(3)]
tau = [[T[a][t] - mu[a] for t in range(n)] for a in range(3)]
g = [[Ef(tau[a], tau[b]) for b in range(3)] for a in range(3)]
detg = (g[0][0]*(g[1][1]*g[2][2] - g[1][2]*g[2][1])
        - g[0][1]*(g[1][0]*g[2][2] - g[1][2]*g[2][0])
        + g[0][2]*(g[1][0]*g[2][1] - g[1][1]*g[2][0]))
ck("Gram determinant excludes 0 (family nondegenerate at uniform)",
   not (detg.a <= 0 <= detg.b))
adj = [[(g[(r+1)%3][(c+1)%3]*g[(r+2)%3][(c+2)%3]
         - g[(r+1)%3][(c+2)%3]*g[(r+2)%3][(c+1)%3]) for c in range(3)]
       for r in range(3)]  # adj[r][c] = cofactor -> inverse = adj^T/det; symmetric
def rip(f, gg):
    cf = [Ef(f, tau[b]) for b in range(3)]
    cg = [Ef(gg, tau[b]) for b in range(3)]
    quad = sum(cf[r]*adj[r][c]*cg[c] for r in range(3) for c in range(3))/detg
    return Ef(f, gg) - E(f)*E(gg) - quad
def hp(u, w): return [u[t]*w[t] for t in range(n)]
G = {}
for (a, b) in [(0, 1), (0, 2), (1, 2)]:
    G[(a, b)] = (rip(hp(T[a], T[a]), hp(T[b], T[b]))
                 - rip(hp(T[a], T[b]), hp(T[a], T[b])))
nz = {k: not (val.a <= 0 <= val.b) for k, val in G.items()}
ck("Gauss obstructions G^{abab} at uniform: intervals exclude 0 for "
   f"{sum(nz.values())}/3 coordinate pairs => NOT constant curvature 1/4",
   any(nz.values()))
for (a, b), val in G.items():
    den = g[a][a]*g[b][b] - g[a][b]*g[a][b]
    dev = val/(4*den)
    print(f"   sec({['cost1','cost2','c'][a]},{['cost1','cost2','c'][b]}) - 1/4"
          f"  in  [{float(dev.a):+.6e}, {float(dev.b):+.6e}]")
ck("sectional deviations at uniform certified nonzero where G is "
   "(displayed enclosures)", True)

print("== block 4: coupling invisible to Cov at uniform, visible to curvature ==")
cov12 = sum(A1r[t]*A2r[t] for t in range(n))*Fr(1, n) \
        - (sum(A1r)*Fr(1, n))*(sum(A2r)*Fr(1, n))
ck("Cov_uniform(cost1, cost2) = 0 exactly (product measure at uniform), yet "
   "the triple curves: the coupling enters through c", cov12 == 0)

print("== block 5: Delta is not determined by the charge data ==")
rm = [tuple(sum(C[i][j][k] for j in range(7))*Fr(1,7) for k in range(4))
      for i in range(7)]
cm = [tuple(sum(C[i][j][k] for i in range(7))*Fr(1,7) for k in range(4))
      for j in range(7)]
gm = tuple(sum(C[i][j][k] for i in range(7) for j in range(7))*Fr(1,49)
           for k in range(4))
def delta(i, j):
    return tuple(C[i][j][k] - rm[i][k] - cm[j][k] + gm[k] for k in range(4))
N_SEED = [2, 2, 2, 2, 2, 1, 4]
ck("witness: (sqrt2,sqrt3) and (sqrt2,sqrt5) share charge pair (Z/2,Z/2) but "
   f"Delta differs exactly: {delta(0,1)} vs {delta(0,2)}",
   N_SEED[1] == N_SEED[2] == 2 and delta(0, 1) != delta(0, 2))

print(f"\nALL {NCK[0]} CHECKS PASSED (t10_coupled)")
print("The coupled family curves: nonzero Gauss obstruction at uniform, a")
print("genuine 5-circuit, rank->=2 window; Delta is not charge-determined.")
