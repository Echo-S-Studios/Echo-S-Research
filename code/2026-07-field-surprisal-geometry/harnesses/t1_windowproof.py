#!/usr/bin/env python3
# t1_windowproof.py -- OP-1 closure: the window identity is a determinant identity.
#
# THEOREM (master window identity).  With C(f,g) = sum_{i<j} w_i w_j (f_i-f_j)(g_i-g_j),
# D = C(a,a)C(X,X)-C(a,X)^2, R(f,g) = D C(f,g) - C(f,.)^T adj(C) C(.,g), F=a^2, G=X^2, H=aX:
#
#   P(w) := R(F,G) - R(H,H)  =  Z^2 * sum_{|s|=4} q4(s) * w_s ,     Z = sum_i w_i,  w_s = prod_{i in s} w_i,
#   q4(s) = l_s(F) l_s(G) - l_s(H)^2,   l_s(f) := det[(1,a,X,f)|_s]  (4x4).
#
# Proof chain, each step machine-verified below:
#   (S1) R(f,g) = det Chat(f,g), the bordered 3x3 determinant of C-covariances of (a,X,f) vs (a,X,g).
#   (S2) Sylvester's determinant identity, pivot Z = B_11:  det Chat(f,g) = Z^{4-2} det B(f,g),
#        B(f,g) = mixed Gram of (1,a,X,f) against (1,a,X,g) under weights w.
#   (S3) mixed Cauchy-Binet:  det B(f,g) = sum_{|s|=4} w_s l_s(f) l_s(g).
#   => P = Z^2 sum_s [l_s(F)l_s(G)-l_s(H)^2] w_s = Z^2 sum_s q4(s) w_s.          [F]
#
# COROLLARIES (all verified):
#   (C1) every monomial has degree 6 and support >= 4; the shape census 140/210/105/7 = 462
#        is forced: (3,1,1,1):35*4, (2,2,1,1):35*6, (2,1,1,1,1):21*5, (1^6):7.
#   (C2) heavy coefficient: [w_h^3 w_j w_k w_l] P = q4({h,j,k,l}), cofactor exactly 1
#        (unique factorization w_h^3 w_j w_k w_l = (w_h * w_h) * w_s, s = {h,j,k,l}).
#   (C3) full coefficient law: [w^k]P = sum over 4-subsets s of supp(k) with k-1_s >= 0
#        of mult(k-1_s) q4(s), mult(2e_a)=1, mult(e_a+e_b)=2.
#   (C4) D = Z * sum_{|t|=3} w_t Delta_t^2, Delta_t = det[(1,a,X)|_t]  (Sylvester+Cauchy-Binet on A);
#        hence q = [sum_s w_s q4(s)] / [Z sum_t w_t Delta_t^2]: the alpha=0 curvature obstruction is
#        a ratio of consecutive Cauchy-Binet compound sums.
#   (C5) on a 4-point face w = w|_s:  q = kappa q4(s) with the EXPLICIT kappa = w_s / (Z sum_{t in s,|t|=3} w_t Delta_t^2) > 0.
#   (C6) factorization: q4(s) = - prod_{t in s} Delta_{s \ {t}}  (product of the four triple minors);
#        the four-point trichotomy (q4 = 0 iff some three of the four points are affinely dependent)
#        is an immediate corollary.
# Exact arithmetic throughout; floats never touch a decision.
import sys, itertools, pickle
from sympy import QQ, Rational, Matrix, cancel, S
from sympy.polys.rings import ring

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

# ---------------------------------------------------------------- generic ring machinery
def det_ring(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = None
    for j in range(n):
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        term = M[0][j] * det_ring(minor)
        if j % 2: term = -term
        tot = term if tot is None else tot + term
    return tot

def build_all(Rng, w, a, X):
    m = len(w)
    one = Rng.one
    Z = sum(w)
    F = [a[i]*a[i] for i in range(m)]
    G = [X[i]*X[i] for i in range(m)]
    H = [a[i]*X[i] for i in range(m)]
    def Sf(f): return sum(w[i]*f[i] for i in range(m))
    def Sfg(f, g): return sum(w[i]*f[i]*g[i] for i in range(m))
    def C(f, g): return Z*Sfg(f, g) - Sf(f)*Sf(g)
    return m, one, Z, F, G, H, Sf, Sfg, C

def P_via_R(Rng, w, a, X):
    m, one, Z, F, G, H, Sf, Sfg, C = build_all(Rng, w, a, X)
    Caa, Cxx, Cax = C(a, a), C(X, X), C(a, X)
    D = Caa*Cxx - Cax*Cax
    def Rr(f, g):
        return (D*C(f, g) - (C(f, a)*Cxx*C(a, g) - C(f, a)*Cax*C(X, g)
                             - C(f, X)*Cax*C(a, g) + C(f, X)*Caa*C(X, g)))
    return Rr(F, G) - Rr(H, H), D, Rr, C

def Chat_det(Rng, w, a, X, f, g):
    m, one, Z, F, G, H, Sf, Sfg, C = build_all(Rng, w, a, X)
    M = [[C(a, a), C(a, X), C(a, g)],
         [C(X, a), C(X, X), C(X, g)],
         [C(f, a), C(f, X), C(f, g)]]
    return det_ring(M)

def B_det(Rng, w, a, X, f, g):
    m, one, Z, F, G, H, Sf, Sfg, C = build_all(Rng, w, a, X)
    u = [[one]*m, a, X, f]; v = [[one]*m, a, X, g]
    M = [[Sfg(u[r], v[c]) for c in range(4)] for r in range(4)]
    return det_ring(M)

def ell(Rng, one, a, X, f, s):
    M = [[one, a[i], X[i], f[i]] for i in s]
    return det_ring(M)

def Delta(Rng, one, a, X, t):
    M = [[one, a[i], X[i]] for i in t]
    return det_ring(M)

def RHS_master(Rng, w, a, X):
    m, one, Z, F, G, H, Sf, Sfg, C = build_all(Rng, w, a, X)
    tot = Rng.zero
    for s in itertools.combinations(range(m), 4):
        lF = ell(Rng, one, a, X, F, s)
        lG = ell(Rng, one, a, X, G, s)
        lH = ell(Rng, one, a, X, H, s)
        ws = w[s[0]]*w[s[1]]*w[s[2]]*w[s[3]]
        tot = tot + (lF*lG - lH*lH)*ws
    return Z*Z*tot

# ================================================================ block 1: fully symbolic proof steps
print("== block 1: fully symbolic verification of the proof chain (generic m=4 and m=5) ==")
for m in (4, 5):
    names = (",".join(f"w{i+1}" for i in range(m)) + ","
             + ",".join(f"a{i+1}" for i in range(m)) + ","
             + ",".join(f"x{i+1}" for i in range(m)))
    Rng, *gens = ring(names, QQ)
    w = gens[:m]; a = gens[m:2*m]; X = gens[2*m:]
    _, one, Z, F, G, H, Sf, Sfg, C = build_all(Rng, w, a, X)
    # (S1) R = det Chat, for both needed pairs
    P_lhs, D, Rr, Cfun = P_via_R(Rng, w, a, X)
    ck(f"m={m} (S1) R(F,G) = det Chat(F,G) symbolically",
       Rr(F, G) - Chat_det(Rng, w, a, X, F, G) == Rng.zero)
    ck(f"m={m} (S1) R(H,H) = det Chat(H,H) symbolically",
       Rr(H, H) - Chat_det(Rng, w, a, X, H, H) == Rng.zero)
    # (S2) Sylvester pivot Z: det Chat = Z^2 det B
    ck(f"m={m} (S2) det Chat(F,G) = Z^2 det B(F,G) symbolically",
       Chat_det(Rng, w, a, X, F, G) - Z*Z*B_det(Rng, w, a, X, F, G) == Rng.zero)
    ck(f"m={m} (S2) det Chat(H,H) = Z^2 det B(H,H) symbolically",
       Chat_det(Rng, w, a, X, H, H) - Z*Z*B_det(Rng, w, a, X, H, H) == Rng.zero)
    # (S3) Cauchy-Binet: det B = sum_s w_s l_s(f) l_s(g)
    for (f, g, nm) in ((F, G, "F,G"), (H, H, "H,H")):
        cb = Rng.zero
        for s in itertools.combinations(range(m), 4):
            ws = w[s[0]]*w[s[1]]*w[s[2]]*w[s[3]]
            cb = cb + ws*ell(Rng, one, a, X, f, s)*ell(Rng, one, a, X, g, s)
        ck(f"m={m} (S3) det B({nm}) = sum_s w_s l_s l_s symbolically",
           B_det(Rng, w, a, X, f, g) - cb == Rng.zero)
    # MASTER
    ck(f"m={m} MASTER  P = Z^2 sum_s q4(s) w_s symbolically",
       P_lhs - RHS_master(Rng, w, a, X) == Rng.zero)
    # (C4) D = Z sum_t w_t Delta_t^2
    dd = Rng.zero
    for t in itertools.combinations(range(m), 3):
        wt = w[t[0]]*w[t[1]]*w[t[2]]
        dt = Delta(Rng, one, a, X, t)
        dd = dd + wt*dt*dt
    ck(f"m={m} (C4) D = Z sum_t w_t Delta_t^2 symbolically", D - Z*dd == Rng.zero)

# ================================================================ block 2: the factorization q4 = -prod Delta
print("== block 2: q4(s) = - prod_{t in s} Delta_(s minus t), fully symbolic (4 generic points) ==")
Rng4, *g4 = ring("a1,a2,a3,a4,x1,x2,x3,x4", QQ)
a4 = g4[:4]; X4 = g4[4:]
one4 = Rng4.one
F4 = [a4[i]*a4[i] for i in range(4)]
G4 = [X4[i]*X4[i] for i in range(4)]
H4 = [a4[i]*X4[i] for i in range(4)]
s = (0, 1, 2, 3)
lF = ell(Rng4, one4, a4, X4, F4, s)
lG = ell(Rng4, one4, a4, X4, G4, s)
lH = ell(Rng4, one4, a4, X4, H4, s)
q4sym = lF*lG - lH*lH
prodD = Rng4.one
for t in itertools.combinations(range(4), 3):
    prodD = prodD * Delta(Rng4, one4, a4, X4, t)
ck("q4(s) + prod of the four triple minors = 0 (exact symbolic)", q4sym + prodD == Rng4.zero)
# trichotomy corollary instances
def q4_inst(pts):
    aa = [Rational(p[0]) for p in pts]; xx = [Rational(p[1]) for p in pts]
    M3 = lambda t: Matrix([[1, aa[i], xx[i]] for i in t]).det()
    l = lambda f: Matrix([[1, aa[i], xx[i], f[i]] for i in range(4)]).det()
    Fv = [aa[i]**2 for i in range(4)]; Gv = [xx[i]**2 for i in range(4)]; Hv = [aa[i]*xx[i] for i in range(4)]
    return cancel(l(Fv)*l(Gv) - l(Hv)**2)
ck("trichotomy: collinear triple => q4 = 0",
   q4_inst([(0, 0), (1, 1), (2, 2), (5, 7)]) == 0)
ck("trichotomy: coincident pair => q4 = 0",
   q4_inst([(3, 4), (3, 4), (1, 0), (0, 1)]) == 0)
ck("trichotomy: generic quadruple => q4 != 0",
   q4_inst([(0, 0), (1, 0), (0, 1), (Rational(1, 3), Rational(1, 5))]) != 0)
ck("trichotomy: vertical pair alone does NOT kill q4 (cost tie only)",
   q4_inst([(1, 0), (1, 5), (0, 1), (2, 3)]) != 0)

# ================================================================ block 3: exact rational instances m=6,7 (all corollaries)
print("== block 3: exact rational instance checks, m = 6 and 7 ==")
import random
random.seed(20260722)
def rq():
    while True:
        v = Rational(random.randint(-9, 9), random.randint(1, 9))
        if v != 0: return v

def q_moment(w, a, X):
    n = len(w); Z = sum(w); p = [cancel(wi/Z) for wi in w]
    def E(h): return cancel(sum(p[i]*h[i] for i in range(n)))
    ta = [cancel(a[i]-E(a)) for i in range(n)]
    tX = [cancel(X[i]-E(X)) for i in range(n)]
    Gm = Matrix(2, 2, lambda r, c: E([[ta, tX][r][i]*[ta, tX][c][i] for i in range(n)]))
    def rip(f, g):
        cf = Matrix(2, 1, [E([f[i]*ta[i] for i in range(n)]), E([f[i]*tX[i] for i in range(n)])])
        cg = Matrix(2, 1, [E([g[i]*ta[i] for i in range(n)]), E([g[i]*tX[i] for i in range(n)])])
        return cancel(E([f[i]*g[i] for i in range(n)]) - E(f)*E(g) - (cf.T*Gm.inv()*cg)[0, 0])
    F = [a[i]*a[i] for i in range(n)]; G = [X[i]*X[i] for i in range(n)]; H = [a[i]*X[i] for i in range(n)]
    return cancel(rip(F, G) - rip(H, H))

def master_rhs_rat(w, a, X):
    n = len(w); Z = sum(w)
    def l(f, s): return Matrix([[1, a[i], X[i], f[i]] for i in s]).det()
    F = [a[i]**2 for i in range(n)]; G = [X[i]**2 for i in range(n)]; H = [a[i]*X[i] for i in range(n)]
    tot = S(0)
    for s in itertools.combinations(range(n), 4):
        ws = w[s[0]]*w[s[1]]*w[s[2]]*w[s[3]]
        tot += ws*cancel(l(F, s)*l(G, s) - l(H, s)**2)
    return cancel(Z**2*tot)

def P_rat(w, a, X):
    n = len(w)
    def covC(f, g):
        Z = sum(w)
        return cancel(Z*sum(w[i]*f[i]*g[i] for i in range(n))
                      - sum(w[i]*f[i] for i in range(n))*sum(w[i]*g[i] for i in range(n)))
    F = [a[i]**2 for i in range(n)]; G = [X[i]**2 for i in range(n)]; H = [a[i]*X[i] for i in range(n)]
    Caa, Cxx, Cax = covC(a, a), covC(X, X), covC(a, X)
    D = cancel(Caa*Cxx - Cax**2)
    def Rr(f, g):
        return cancel(D*covC(f, g) - (covC(f, a)*Cxx*covC(a, g) - covC(f, a)*Cax*covC(X, g)
                                      - covC(f, X)*Cax*covC(a, g) + covC(f, X)*Caa*covC(X, g)))
    return cancel(Rr(F, G) - Rr(H, H)), D

for n in (6, 7):
    for trial in range(2):
        a = [rq() for _ in range(n)]
        while len({str(v) for v in a}) < n:
            a = [rq() for _ in range(n)]
        X = [rq() for _ in range(n)]
        w = [Rational(random.randint(1, 20)) for _ in range(n)]
        Pn, D = P_rat(w, a, X)
        ck(f"m={n} trial {trial}: MASTER P = Z^2 sum q4 w_s (exact instance)",
           cancel(Pn - master_rhs_rat(w, a, X)) == 0)
        Z = sum(w)
        dd = sum(w[t[0]]*w[t[1]]*w[t[2]]*Matrix([[1, a[i], X[i]] for i in t]).det()**2
                 for t in itertools.combinations(range(n), 3))
        ck(f"m={n} trial {trial}: D = Z sum w_t Delta_t^2 (exact instance)",
           cancel(D - Z*dd) == 0)
        qm = q_moment(w, a, X)
        num = sum(w[s[0]]*w[s[1]]*w[s[2]]*w[s[3]]
                  *cancel(Matrix([[1, a[i], X[i], a[i]**2] for i in s]).det()
                          *Matrix([[1, a[i], X[i], X[i]**2] for i in s]).det()
                          - Matrix([[1, a[i], X[i], a[i]*X[i]] for i in s]).det()**2)
                  for s in itertools.combinations(range(n), 4))
        ck(f"m={n} trial {trial}: q = [sum w_s q4]/[Z sum w_t Delta_t^2] (compound ratio, exact)",
           cancel(qm*Z*dd - num) == 0)

# (C5) explicit kappa on a 4-point face
n = 7
a = [rq() for _ in range(n)]
while len({str(v) for v in a}) < n: a = [rq() for _ in range(n)]
X = [rq() for _ in range(n)]
sface = (1, 3, 4, 6)
w = [S(0)]*n
for idx, i in enumerate(sface): w[i] = Rational([3, 5, 7, 11][idx])
Z = sum(w)
qm = q_moment(w, a, X)
q4v = cancel(Matrix([[1, a[i], X[i], a[i]**2] for i in sface]).det()
             *Matrix([[1, a[i], X[i], X[i]**2] for i in sface]).det()
             - Matrix([[1, a[i], X[i], a[i]*X[i]] for i in sface]).det()**2)
den = cancel(Z*sum(w[t[0]]*w[t[1]]*w[t[2]]*Matrix([[1, a[i], X[i]] for i in t]).det()**2
                   for t in itertools.combinations(sface, 3)))
ws = w[sface[0]]*w[sface[1]]*w[sface[2]]*w[sface[3]]
kappa = cancel(ws/den)
ck("(C5) 4-point face: q = kappa q4 with kappa = w_s/(Z sum_{t in s} w_t Delta_t^2), exact",
   cancel(qm - kappa*q4v) == 0)
ck("(C5) kappa > 0 on the face", kappa > 0)

# ================================================================ block 4: the m=7 catalog ring -- dict-level equality
print("== block 4: catalog ring (m=7): rebuild the engine dictionary of P and equate with Z^2 sum q4 w_s ==")
m = 7
Rng7, *gens = ring("X1,X2,X3,X4,X5,X6,X7,L2,L3,L5,Lp", QQ)
Xg = gens[:7]
L2v, L3v, L5v, Lpv = gens[7:]
acat = [L2v, L3v, L5v, Lpv, Lpv, 4*Lpv, L5v/2 + 2*Lpv]
Xs = list(Xg)
one7 = Rng7.one
Fc = [acat[i]*acat[i] for i in range(m)]
Gc = [Xs[i]*Xs[i] for i in range(m)]
Hc = [acat[i]*Xs[i] for i in range(m)]

def Cdict(f, g):
    d = {}
    for i in range(m):
        for j in range(i+1, m):
            c = (f[i]-f[j])*(g[i]-g[j])
            if c: d[(i, j)] = c
    return d
Caa = Cdict(acat, acat); Cxx = Cdict(Xs, Xs); Cax = Cdict(acat, Xs)
def triple_terms(fg):
    f, g = fg
    Cfg = Cdict(f, g); Cfa = Cdict(f, acat); Cfx = Cdict(f, Xs)
    Cag = Cdict(acat, g); Cxg = Cdict(Xs, g)
    return [(+1, Caa, Cxx, Cfg), (-1, Cax, Cax, Cfg),
            (-1, Cfa, Cxx, Cag), (+1, Cfa, Cax, Cxg),
            (+1, Cfx, Cax, Cag), (-1, Cfx, Caa, Cxg)]
def accumulate(terms, sgn_all, out):
    for sgn, dA, dB, dC in terms:
        sg = sgn*sgn_all
        for (i1, j1), c1 in dA.items():
            for (i2, j2), c2 in dB.items():
                c12 = c1*c2
                for (i3, j3), c3 in dC.items():
                    k = [0]*m
                    k[i1] += 1; k[j1] += 1; k[i2] += 1; k[j2] += 1; k[i3] += 1; k[j3] += 1
                    v = c12*c3
                    if sg < 0: v = -v
                    key = tuple(k)
                    out[key] = out.get(key, Rng7.zero) + v
print("rebuilding P (12 triple convolutions)...")
P = {}
accumulate(triple_terms((Fc, Gc)), +1, P)
accumulate(triple_terms((Hc, Hc)), -1, P)
P = {k: v for k, v in P.items() if v}
ck("rebuilt P has 462 monomials", len(P) == 462)

print("building RHS dictionary Z^2 sum_s q4(s) w_s ...")
RHSd = {}
q4cache = {}
for s in itertools.combinations(range(m), 4):
    lF = ell(Rng7, one7, acat, Xs, Fc, s)
    lG = ell(Rng7, one7, acat, Xs, Gc, s)
    lH = ell(Rng7, one7, acat, Xs, Hc, s)
    q4cache[s] = lF*lG - lH*lH
    for al in range(m):
        for be in range(al, m):
            k = [0]*m
            for i in s: k[i] += 1
            k[al] += 1; k[be] += 1
            coeff = q4cache[s] if al == be else 2*q4cache[s]
            key = tuple(k)
            RHSd[key] = RHSd.get(key, Rng7.zero) + coeff
RHSd = {k: v for k, v in RHSd.items() if v}
ck("RHS dictionary also has exactly 462 monomials", len(RHSd) == 462)
ck("MASTER on the catalog ring: P == Z^2 sum_s q4(s) w_s, all 462 coefficients equal",
   set(P.keys()) == set(RHSd.keys()) and all(P[k] - RHSd[k] == Rng7.zero for k in P))

# shape census forced by the identity
from collections import Counter
shape_count = Counter(tuple(sorted((x for x in k if x), reverse=True)) for k in P)
ck("shape census = {(3,1,1,1):140, (2,2,1,1):210, (2,1,1,1,1):105, (1^6):7}",
   shape_count == Counter({(2, 2, 1, 1): 210, (3, 1, 1, 1): 140, (2, 1, 1, 1, 1): 105,
                           (1, 1, 1, 1, 1, 1): 7}))
ck("support >= 4 and degree 6 for every monomial (corollary C1)",
   all(sum(k) == 6 and sum(1 for x in k if x) >= 4 for k in P))

# (C2) heavy coefficients: cofactor exactly 1, all 140 placements, directly from the identity
heavy_ok = 0
for s in itertools.combinations(range(m), 4):
    for h in s:
        k = [1 if i in s else 0 for i in range(m)]
        k[h] = 3
        if P.get(tuple(k), Rng7.zero) - q4cache[s] == Rng7.zero:
            heavy_ok += 1
ck("(C2) all 140 heavy coefficients equal q4(s) with cofactor exactly 1", heavy_ok == 140)

# (C3) full coefficient law spot-verified on every shape
def law_coeff(k):
    supp = [i for i in range(m) if k[i] > 0]
    tot = Rng7.zero
    for s in itertools.combinations(supp, 4):
        rem = list(k)
        ok = True
        for i in s:
            rem[i] -= 1
            if rem[i] < 0: ok = False
        if not ok or sum(rem) != 2: continue
        nz = [i for i in range(m) if rem[i] > 0]
        mult = 1 if len(nz) == 1 else 2
        tot = tot + mult*q4cache[tuple(s)]
    return tot
random.seed(7)
sample = random.sample(list(P.keys()), 40)
ck("(C3) coefficient law reproduces 40 random monomial coefficients exactly",
   all(P[k] - law_coeff(k) == Rng7.zero for k in sample))
ck("(C3) coefficient law reproduces the (1^6) and (2,2,1,1) shapes exactly",
   all(P[k] - law_coeff(k) == Rng7.zero for k in P
       if tuple(sorted((x for x in k if x), reverse=True)) in ((1, 1, 1, 1, 1, 1), (2, 2, 1, 1))))

# cross-check against the artifact saved by t1_engine.py, if present
try:
    with open('P_coeffs.pkl', 'rb') as fh:
        saved = pickle.load(fh)
    same_keys = set(saved.keys()) == set(P.keys())
    same_vals = all(str(P[k]) == saved[k] for k in P)
    ck("rebuilt P matches the artifact saved by t1_engine.py (string-exact)", same_keys and same_vals)
except FileNotFoundError:
    print("      (P_coeffs.pkl not present; artifact cross-check skipped, not counted)")

print(f"\nALL {NCK[0]} CHECKS PASSED (t1_windowproof)")
