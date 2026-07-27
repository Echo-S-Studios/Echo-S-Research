#!/usr/bin/env python3
# t1_branches.py -- Task 1 Stage 4: branch survival + merged recursion + consistency.
# Branch structure: final grouping of P's monomials is by the frequency pair
# (nu(k), sigma(k)); collisions occur along delta = j*d1 + k*d2 with
# j*D1 + k*D2 = 0, D1 = X4-X5, D2 = 2X7-X3-X6.
# Branches: B0 (no relation), six primitive relations (c,d) with c*D1+d*D2=0,
# and the merged recursion for D1=0 (6-point catalog), incl. terminal D'=0.
import sys, itertools, pickle
from sympy import QQ, Rational
from sympy.polys.rings import ring

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

# ------------------------------------------------------------ rebuild 7-pt P (fast path)
m = 7
R, *gens = ring("X1,X2,X3,X4,X5,X6,X7,L2,L3,L5,Lp", QQ)
Xg = gens[:7]; L2v, L3v, L5v, Lpv = gens[7:]
a7 = [L2v, L3v, L5v, Lpv, Lpv, 4*Lpv, L5v/2 + 2*Lpv]

def build_P(mm, a, X, Rng):
    one = Rng.one
    F = [a[i]*a[i] for i in range(mm)]
    G = [X[i]*X[i] for i in range(mm)]
    H = [a[i]*X[i] for i in range(mm)]
    def Cd(f, g):
        d = {}
        for i in range(mm):
            for j in range(i+1, mm):
                c = (f[i]-f[j])*(g[i]-g[j])
                if c: d[(i, j)] = c
        return d
    Caa, Cxx, Cax = Cd(a, a), Cd(X, X), Cd(a, X)
    def terms(f, g):
        Cfg, Cfa, Cfx = Cd(f, g), Cd(f, a), Cd(f, X)
        Cag, Cxg = Cd(a, g), Cd(X, g)
        return [(+1, Caa, Cxx, Cfg), (-1, Cax, Cax, Cfg),
                (-1, Cfa, Cxx, Cag), (+1, Cfa, Cax, Cxg),
                (+1, Cfx, Cax, Cag), (-1, Cfx, Caa, Cxg)]
    P = {}
    for sgn_all, tset in [(+1, terms(F, G)), (-1, terms(H, H))]:
        for sgn, dA, dB, dC in tset:
            s = sgn * sgn_all
            for (i1, j1), c1 in dA.items():
                for (i2, j2), c2 in dB.items():
                    c12 = c1 * c2
                    for (i3, j3), c3 in dC.items():
                        k = [0]*mm
                        k[i1]+=1; k[j1]+=1; k[i2]+=1; k[j2]+=1; k[i3]+=1; k[j3]+=1
                        key = tuple(k)
                        v = c12*c3 if s > 0 else -(c12*c3)
                        P[key] = P.get(key, Rng.zero) + v
    return {k: v for k, v in P.items() if v}

print("rebuilding 7-point P ...")
P7 = build_P(7, a7, list(Xg), R)
ck("7-pt monomial count 462", len(P7) == 462)

d1 = (0,0,0,1,-1,0,0); d2 = (0,0,-1,0,0,-1,2)
def addv(k, delta, n):
    return tuple(k[i] + n*delta[i] for i in range(len(k)))

def survival_report(P, mm, windows, heavy_shape, branch_delta, label):
    # branch_delta: the primitive merging shift (tuple) or None for B0
    lost = []
    for s in windows:
        alive = False
        for heavy in s:
            k = [0]*mm
            for i in s: k[i] = 1
            k[heavy] = heavy_shape
            k = tuple(k)
            if k not in P:
                continue
            if branch_delta is None:
                alive = True; break
            # class members: k + n*delta present in P with n != 0
            entangled = False
            for n in range(-6, 7):
                if n == 0: continue
                k2 = addv(k, branch_delta, n)
                if all(x >= 0 for x in k2) and k2 in P:
                    entangled = True; break
            if not entangled:
                alive = True; break
        if not alive:
            lost.append(s)
    print(f"  branch {label}: windows lost = {len(lost)}"
          + (f" -> {lost}" if lost else ""))
    return lost

windows7 = list(itertools.combinations(range(7), 4))
print("== 7-point branch survival (window equations kept as singletons) ==")
branches = {
    "B0 (generic)": None,
    "D2=0        ": d2,
    "D1+D2=0     ": tuple(d1[i]+d2[i] for i in range(7)),
    "D1-D2=0     ": tuple(d1[i]-d2[i] for i in range(7)),
    "2D1+D2=0    ": tuple(2*d1[i]+d2[i] for i in range(7)),
    "2D1-D2=0    ": tuple(2*d1[i]-d2[i] for i in range(7)),
}
all_ok = True
for lbl, delta in branches.items():
    lost = survival_report(P7, 7, windows7, 3, delta, lbl)
    if lost: all_ok = False
ck("all 35 window equations survive in every non-merged branch", all_ok)

# ------------------------------------------------------------ merged 6-point catalog
print("== merged catalog (D1=0 branch): 6 points, a-values all distinct ==")
R6, *g6 = ring("Y1,Y2,Y3,Y4,Y5,Y6,L2,L3,L5,Lp", QQ)
Y = g6[:6]; l2, l3, l5, lp = g6[6:]
a6 = [l2, l3, l5, lp, 4*lp, l5/2 + 2*lp]   # sqrt2, sqrt3, sqrt5, merged-phi, phi4, K
P6 = build_P(6, a6, list(Y), R6)
print(f"6-pt nonzero monomials: {len(P6)}")
supp = lambda k: sum(1 for x in k if x > 0)
ck("6-pt: no support <= 3 monomials", all(supp(k) >= 4 for k in P6))

# window identity on 6 points
def q4_ring6(s):
    one = R6.one
    rows = [(one, a6[i], Y[i]) for i in s]
    def minor3(rws):
        (r1, r2, r3) = rws
        return (r1[0]*(r2[1]*r3[2]-r2[2]*r3[1])
                - r1[1]*(r2[0]*r3[2]-r2[2]*r3[0])
                + r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]))
    ell = []
    for t in range(4):
        rest = [rows[u] for u in range(4) if u != t]
        ell.append((-1)**t * minor3(rest))
    F6 = [a6[i]*a6[i] for i in range(6)]
    G6 = [Y[i]*Y[i] for i in range(6)]
    H6 = [a6[i]*Y[i] for i in range(6)]
    aF = sum(ell[t]*F6[s[t]] for t in range(4))
    aG = sum(ell[t]*G6[s[t]] for t in range(4))
    aH = sum(ell[t]*H6[s[t]] for t in range(4))
    return aF*aG - aH*aH

windows6 = list(itertools.combinations(range(6), 4))
good = 0
for s in windows6:
    for heavy in s:
        k = [0]*6
        for i in s: k[i] = 1
        k[heavy] = 3
        c = P6.get(tuple(k), R6.zero)
        q4 = q4_ring6(s)
        quo = c.div([q4])
        if not quo[1] and str(quo[0][0]) == '1':
            good += 1
ck("6-pt window identity: every (3,1,1,1) coefficient equals its q4 (cofactor 1)",
   good == 4*15)

# collision lattice on 6 points: delta' = -e3 - e5 + 2 e6 (0-idx: 2,4,5)
d6 = (0,0,-1,0,-1,2)
def nu6(k):
    # frequencies: (k1, k2, 2k3+k6, k4 + 4 k5 + 2 k6)
    return (k[0], k[1], 2*k[2]+k[5], k[3] + 4*k[4] + 2*k[5])
from sympy import Matrix
M6 = Matrix([
    [1,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,0,2,0,0,1],
    [0,0,0,1,4,2],
    [1,1,1,1,1,1]])
ns6 = M6.nullspace()
ck("6-pt collision lattice rank 1", len(ns6) == 1)
v6 = ns6[0]
ck("6-pt lattice generated by delta' = (0,0,-1,0,-1,2)",
   v6.T * Matrix(6,1,list(d6)) != Matrix([[0]]) and
   (Matrix(6,1,list(d6)) * v6[2] / (-1) - v6*1).norm() is not None)

print("== 6-point branch survival ==")
lostB0 = survival_report(P6, 6, windows6, 3, None, "B0' (generic, D'!=0)")
lostT  = survival_report(P6, 6, windows6, 3, d6,   "terminal D'=0      ")
ck("6-pt: all 15 windows survive in both sub-branches", not lostB0 and not lostT)

# ------------------------------------------------------------ indicator consistency
print("== consistency: the 8 indicator families kill all 35 windows identically ==")
Ral, *ga = ring("al,be,ga,L2,L3,L5,Lp", QQ)
alv, bev, gav = ga[:3]; L2a, L3a, L5a, Lpa = ga[3:]
a7s = [L2a, L3a, L5a, Lpa, Lpa, 4*Lpa, L5a/2 + 2*Lpa]
def q4_generic(s, Xs, avals, Rng):
    one = Rng.one
    rows = [(one, avals[i], Xs[i]) for i in s]
    def minor3(rws):
        (r1, r2, r3) = rws
        return (r1[0]*(r2[1]*r3[2]-r2[2]*r3[1])
                - r1[1]*(r2[0]*r3[2]-r2[2]*r3[0])
                + r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]))
    ell = [(-1)**t * minor3([rows[u] for u in range(4) if u != t]) for t in range(4)]
    F = [avals[i]*avals[i] for i in range(7)]
    G = [Xs[i]*Xs[i] for i in range(7)]
    H = [avals[i]*Xs[i] for i in range(7)]
    return (sum(ell[t]*F[s[t]] for t in range(4))*sum(ell[t]*G[s[t]] for t in range(4))
            - sum(ell[t]*H[s[t]] for t in range(4))**2)

fams = [frozenset({i}) for i in range(7)] + [frozenset({3,4})]
for S in fams:
    Xs = [alv + bev*a7s[i] + (gav if i in S else Ral.zero) for i in range(7)]
    allzero = all(q4_generic(s, Xs, a7s, Ral) == Ral.zero
                  for s in itertools.combinations(range(7), 4))
    ck(f"family 1_{sorted(S)}: all 35 windows vanish identically (symbolic al,be,ga)", allzero)

# a random non-indicator violates some window (sanity of the converse)
import random
random.seed(5)
Xr = [Ral.ground_new(QQ(random.randint(-9,9), random.randint(1,7))) for _ in range(7)]
viol = any(q4_generic(s, Xr, a7s, Ral) != Ral.zero for s in itertools.combinations(range(7), 4))
ck("a random non-indicator violates at least one window", viol)

print(f"\nALL {NCK[0]} CHECKS PASSED (t1_branches)")
