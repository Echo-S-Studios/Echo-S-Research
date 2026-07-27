#!/usr/bin/env python3
# t1_reduction.py -- Task 1 Stage 2: reduction machinery, exact.
#   q(theta) = <Q a^2, Q X^2> - |Q aX|^2  (vanishing <=> constant curvature 1/4)
#   Z^2 D q = P(w) := R(F,G) - R(H,H),  F=a^2, G=X^2, H=aX
#   R(f,g) = D C(f,g) - C(f,.)^T adj(C) C(.,g),  C(f,g)=sum_{i<j} w_i w_j (f_i-f_j)(g_i-g_j)
#   q = Tr(Gamma(theta) N) with N = sym([a^2][X^2]^T) - [aX][aX]^T in Sym_r
# All decisions exact over QQ (rational instances) or QQ(L)[X] (engine).
import sys, itertools
from sympy import (symbols, Rational, Matrix, cancel, expand, S, QQ, zeros,
                   linear_eq_to_matrix, Poly)
from sympy.polys.rings import ring

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

m = 7

# ---------------------------------------------------------------- helpers (exact, Matrix based)
def covC(w, f, g):
    Z = sum(w)
    return cancel(Z*sum(w[i]*f[i]*g[i] for i in range(len(w)))
                  - sum(w[i]*f[i] for i in range(len(w)))*sum(w[i]*g[i] for i in range(len(w))))

def covC_pairs(w, f, g):
    return cancel(sum(w[i]*w[j]*(f[i]-f[j])*(g[i]-g[j])
                      for i in range(len(w)) for j in range(i+1, len(w))))

def q_moment(w, a, X):
    # residual-moment (Q-projection) route, exact
    n = len(w); Z = sum(w); p = [cancel(wi/Z) for wi in w]
    def E(h): return cancel(sum(p[i]*h[i] for i in range(n)))
    def cen(h):
        mh = E(h); return [cancel(h[i]-mh) for h in [h]][0] if False else [cancel(hi-mh) for hi in h]
    ta, tX = cen(a), cen(X)
    G = Matrix(2,2, lambda r,c: E([ [ta,tX][r][i]*[ta,tX][c][i] for i in range(n)]))
    def resid_ip(f, g):
        cf = [E([f[i]*ta[i] for i in range(n)]), E([f[i]*tX[i] for i in range(n)])]
        cg = [E([g[i]*ta[i] for i in range(n)]), E([g[i]*tX[i] for i in range(n)])]
        base = cancel(E([f[i]*g[i] for i in range(n)]) - E(f)*E(g))
        corr = (Matrix(1,2,cf) * G.inv() * Matrix(2,1,cg))[0,0]
        return cancel(base - corr)
    F = [a[i]*a[i] for i in range(n)]
    Gf = [X[i]*X[i] for i in range(n)]
    H = [a[i]*X[i] for i in range(n)]
    return cancel(resid_ip(F, Gf) - resid_ip(H, H))

def P_numer(w, a, X):
    n = len(w)
    F = [a[i]*a[i] for i in range(n)]
    Gf = [X[i]*X[i] for i in range(n)]
    H = [a[i]*X[i] for i in range(n)]
    Caa = covC(w, a, a); Cxx = covC(w, X, X); Cax = covC(w, a, X)
    D = cancel(Caa*Cxx - Cax**2)
    def R(f, g):
        Cfg = covC(w, f, g)
        Cfa = covC(w, f, a); Cfx = covC(w, f, X)
        Cag = covC(w, a, g); Cxg = covC(w, X, g)
        return cancel(D*Cfg - (Cfa*Cxx*Cag - Cfa*Cax*Cxg - Cfx*Cax*Cag + Cfx*Caa*Cxg))
    return cancel(R(F, Gf) - R(H, H)), D

# ---------------------------------------------------------------- block 1: identities at exact instances
print("== block 1: exact identity checks (rational instances) ==")
import random
random.seed(20260722)
def rand_q(lo=-9, hi=9):
    while True:
        v = Rational(random.randint(lo, hi), random.randint(1, 9))
        if v != 0: return v

wsym = [Rational(x) for x in [2, 3, 5, 7, 11, 13, 17]]
for trial in range(3):
    a = [rand_q() for _ in range(m)]
    # ensure a distinct
    while len({str(x) for x in a}) < m: a = [rand_q() for _ in range(m)]
    X = [rand_q() for _ in range(m)]
    w = [Rational(random.randint(1, 20)) for _ in range(m)]
    Z = sum(w)
    qm = q_moment(w, a, X)
    Pn, D = P_numer(w, a, X)
    ck(f"trial {trial}: Z^2 D q == P (exact)", cancel(Z**2 * D * qm - Pn) == 0)
    ck(f"trial {trial}: pair-sum identity for C", 
       cancel(covC(w, a, X) - covC_pairs(w, a, X)) == 0)

# ---------------------------------------------------------------- block 2: Tr(Gamma N) reduction at instances
print("== block 2: q = Tr(Gamma N) at exact instances ==")
def exact_rank(M_):
    return M_.rank(iszerofunc=lambda x: cancel(x) == 0, simplify=cancel)

def trGammaN_check(w, a, X):
    n = len(w)
    ONE = [S(1)]*n
    F = [a[i]*a[i] for i in range(n)]
    Gf = [X[i]*X[i] for i in range(n)]
    H = [a[i]*X[i] for i in range(n)]
    Vc = [ONE, a, X]
    B = Matrix([[Vc[j][i] for j in range(3)] for i in range(n)])
    # complement basis of V inside V^2: pick from {F,H,G} the classes that extend rank
    cand = [F, H, Gf]
    comp = []
    cur = B
    for u in cand:
        u_col = Matrix(n, 1, u)
        if exact_rank(Matrix.hstack(cur, u_col)) > exact_rank(cur):
            comp.append(u); cur = Matrix.hstack(cur, u_col)
    r = len(comp)
    # coordinates of [F],[H],[G] in the complement basis (mod V): solve  F = V-part + sum f_al u_al
    def coords(target):
        A_ = Matrix.hstack(B, Matrix([[c[i] for c in comp] for i in range(n)]))
        sol = A_.solve(Matrix(n, 1, target))
        return Matrix(r, 1, sol[3:])
    fF, fH, fG = coords(F), coords(H), coords(Gf)
    N = Rational(1,2)*(fF*fG.T + fG*fF.T) - fH*fH.T
    # Gamma: residual inner products of the complement basis under p = w/Z
    Z = sum(w); p = [cancel(wi/Z) for wi in w]
    def E(h): return cancel(sum(p[i]*h[i] for i in range(n)))
    ta = [cancel(a[i]-E(a)) for i in range(n)]
    tX = [cancel(X[i]-E(X)) for i in range(n)]
    Gm = Matrix(2,2, lambda rr,cc: E([[ta,tX][rr][i]*[ta,tX][cc][i] for i in range(n)]))
    def resid_ip(f, g):
        cf = Matrix(2,1,[E([f[i]*ta[i] for i in range(n)]), E([f[i]*tX[i] for i in range(n)])])
        cg = Matrix(2,1,[E([g[i]*ta[i] for i in range(n)]), E([g[i]*tX[i] for i in range(n)])])
        return cancel(E([f[i]*g[i] for i in range(n)]) - E(f)*E(g) - (cf.T*Gm.inv()*cg)[0,0])
    Gam = Matrix(r, r, lambda rr, cc: resid_ip(comp[rr], comp[cc]))
    lhs = q_moment(w, a, X)
    rhs = cancel((Gam*N).trace())
    return cancel(lhs - rhs) == 0, r

for trial in range(2):
    a = [rand_q() for _ in range(m)]
    while len({str(x) for x in a}) < m: a = [rand_q() for _ in range(m)]
    X = [rand_q() for _ in range(m)]
    w = [Rational(random.randint(1, 15)) for _ in range(m)]
    ok, r = trGammaN_check(w, a, X)
    ck(f"trial {trial}: q = Tr(Gamma N) exact (r={r})", ok and r == 3)

# r = 2 instance: pair-pole X on the true catalog values needs L symbols; use a rational surrogate
a_sur = [Rational(7,10), Rational(11,10), Rational(16,10), Rational(48,100),
         Rational(48,100), Rational(192,100), Rational(176,100)]
Xpp = [S(0), S(0), S(0), S(1), S(2), S(0), S(0)]
w = [Rational(random.randint(1, 15)) for _ in range(m)]
ok, r = trGammaN_check(w, a_sur, Xpp)
ck(f"pair-pole surrogate: q = Tr(Gamma N) exact (r={r})", ok and r == 2)

# ---------------------------------------------------------------- block 3: support <= 3 vanishing + 4-window
print("== block 3: face structure ==")
a = [rand_q() for _ in range(m)]
while len({str(x) for x in a}) < m: a = [rand_q() for _ in range(m)]
X = [rand_q() for _ in range(m)]
w3 = [Rational(3), Rational(5), Rational(7), 0, 0, 0, 0]
Pn3, D3 = P_numer(w3, a, X)
ck("P vanishes on 3-point supports (exact instance)", cancel(Pn3) == 0)

def ann4(s, a, X):
    B = Matrix([[1, a[i], X[i]] for i in s])   # 4x3
    ns = B.T.nullspace(iszerofunc=lambda x: cancel(x)==0)
    # nullspace of B^T? we need ell in R^4 with ell^T B = 0 -> B^T ell = 0
    ns = Matrix(4,1,[0]*4) if not ns else ns[0]
    return ns

s = (0, 2, 4, 6)
w4 = [0]*m
for idx, i in enumerate(s): w4[i] = Rational([3,5,7,11][idx])
Pn4, D4 = P_numer(w4, a, X)
ell = ann4(s, a, X)
def annval(ell, s, h):
    return cancel(sum(ell[k]*h[s[k]] for k in range(4)))
F = [a[i]*a[i] for i in range(m)]; Gf = [X[i]*X[i] for i in range(m)]; H = [a[i]*X[i] for i in range(m)]
q4 = cancel(annval(ell,s,F)*annval(ell,s,Gf) - annval(ell,s,H)**2)
qval = cancel(Pn4/( (sum(w4))**2 * D4 ))
# on a 4-point support residuals are rank-1: q = kappa * q4 with kappa > 0
# kappa = <Q u, Q u>/ann(u)^2 for any u with ann(u) != 0
uF = F
Zw = sum(w4)
p4 = [cancel(wi/Zw) for wi in w4]
def E4(h): return cancel(sum(p4[i]*h[i] for i in range(m)))
ta = [cancel(a[i]-E4(a)) for i in range(m)]; tX = [cancel(X[i]-E4(X)) for i in range(m)]
G2 = Matrix(2,2, lambda rr,cc: E4([[ta,tX][rr][i]*[ta,tX][cc][i] for i in range(m)]))
def rip(f,g):
    cf = Matrix(2,1,[E4([f[i]*ta[i] for i in range(m)]), E4([f[i]*tX[i] for i in range(m)])])
    cg = Matrix(2,1,[E4([g[i]*ta[i] for i in range(m)]), E4([g[i]*tX[i] for i in range(m)])])
    return cancel(E4([f[i]*g[i] for i in range(m)]) - E4(f)*E4(g) - (cf.T*G2.inv()*cg)[0,0])
kappa = cancel(rip(uF,uF)/annval(ell,s,F)**2)
ck("4-window: q = kappa * q4 with kappa = <Qu,Qu>/ann(u)^2 (exact)",
   cancel(qval - kappa*q4) == 0)
ck("4-window: kappa > 0 at the instance", kappa > 0)
ck("4-window: rank-1 structure <Qf,Qg> = kappa ann(f)ann(g) (cross check F,H)",
   cancel(rip(F,H) - kappa*annval(ell,s,F)*annval(ell,s,H)) == 0)

# ---------------------------------------------------------------- block 4: collision lattice of the catalog
print("== block 4: theta_1 frequency collision lattice (7-point catalog) ==")
# a-frequency of monomial k: nu(k) = (k1, k2, k3 + k7/2, k4 + k5 + 4 k6 + 2 k7)
# scaled integer: (k1, k2, 2 k3 + k7, k4 + k5 + 4 k6 + 2 k7)
def nu(k):
    return (k[0], k[1], 2*k[2] + k[6], k[3] + k[4] + 4*k[5] + 2*k[6])
# lattice of differences delta with nu(delta)=0, sum(delta)=0 -- verify rank 2 with
# generators d1 = e4 - e5 (phi/tau swap), d2 = -e3 - e6 + 2 e7 (Salem square)
d1 = (0,0,0,1,-1,0,0); d2 = (0,0,-1,0,0,-1,2)
def nu_lin(d):
    return (d[0], d[1], 2*d[2]+d[6], d[3]+d[4]+4*d[5]+2*d[6])
ck("generator d1 in lattice", nu_lin(d1) == (0,0,0,0) and sum(d1) == 0)
ck("generator d2 in lattice", nu_lin(d2) == (0,0,0,0) and sum(d2) == 0)
# rank of the constraint system: deltas satisfy 5 linear conditions on 7 ints ->
# solution rank should be 2; verify by solving the homogeneous system over QQ
Msys = Matrix([
    [1,0,0,0,0,0,0],
    [0,1,0,0,0,0,0],
    [0,0,2,0,0,0,1],
    [0,0,0,1,1,4,2],
    [1,1,1,1,1,1,1]])
ns = Msys.nullspace()
ck("collision lattice has rank exactly 2", len(ns) == 2)
sp = Matrix.hstack(*ns)
ck("d1, d2 span the lattice", exact_rank(Matrix.hstack(sp, Matrix(7,1,d1), Matrix(7,1,d2))) == 2)

print(f"\nALL {NCK[0]} CHECKS PASSED (t1_reduction)")
