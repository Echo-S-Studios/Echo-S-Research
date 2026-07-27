#!/usr/bin/env python3
# t3_suspension.py -- Task 3: iterated suspension / join theorem + TG threshold.
# Part A: exact symbolic -- join metric has constant sectional curvature 1;
#         pullback of the round metric under the join map is the join metric.
# Part B: validated numerics -- sectional curvature routines validated on
#         sphere(+1), hyperbolic(-1), flat(0), full simplex(1/4), then applied
#         to the catalog 3-family (log M, 1_K, 1_{phi4})  [C].
# Part C: exact -- totally geodesic threshold: first TG at k=5, uniquely R[log M].
import sys, itertools
from sympy import (symbols, sin, cos, tan, Function, diff, simplify, trigsimp,
                   Matrix, Rational, cancel, S, exp, log, sqrt)

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

# ---------------------------------------------------------------- Part A: symbolic
def riemann_lower(coords, g):
    n = len(coords)
    ginv = g.inv()
    Gam = [[[S(0)]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                s = S(0)
                for m_ in range(n):
                    s += ginv[i, m_]*(diff(g[m_, j], coords[k]) + diff(g[m_, k], coords[j])
                                      - diff(g[j, k], coords[m_]))
                Gam[i][j][k] = simplify(s/2)
    Rl = [[[[S(0)]*n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    expr = diff(Gam[i][l][j], coords[k]) - diff(Gam[i][k][j], coords[l])
                    for m_ in range(n):
                        expr += Gam[i][k][m_]*Gam[m_][l][j] - Gam[i][l][m_]*Gam[m_][k][j]
                    Rl[i][j][k][l] = expr
    # lower first index
    Rlow = [[[[S(0)]*n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    s = S(0)
                    for m_ in range(n):
                        s += g[i, m_]*Rl[m_][j][k][l]
                    Rlow[i][j][k][l] = trigsimp(simplify(s))
    return Rlow

def is_zero_expr(e):
    from sympy import expand_trig, together, exp as _exp
    e1 = simplify(expand_trig(e))
    if e1 == 0:
        return True
    return cancel(together(e.rewrite('exp').expand())) == 0

def is_constant_curvature(coords, g, K):
    R = riemann_lower(coords, g)
    n = len(coords)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    target = K*(g[i, k]*g[j, l] - g[i, l]*g[j, k])
                    if not is_zero_expr(R[i][j][k][l] - target):
                        return False, (i, j, k, l)
    return True, None

print("== Part A: join metrics have constant sectional curvature 1 (exact) ==")
Psi, u1, u2, ell = symbols('Psi u1 u2 ell', real=True)
# k = 3: corner S^1, base curve; g = dPsi^2 + sin^2 Psi du1^2 + cos^2 Psi dell^2
g3 = Matrix([[1, 0, 0], [0, sin(Psi)**2, 0], [0, 0, cos(Psi)**2]])
ok, bad = is_constant_curvature([Psi, u1, ell], g3, 1)
ck("k=3 join metric: R_ijkl = g_ik g_jl - g_il g_jk (sec == 1)", ok, str(bad))
# k = 4: corner S^2 (round), base curve
g4 = Matrix([[1, 0, 0, 0],
             [0, sin(Psi)**2, 0, 0],
             [0, 0, sin(Psi)**2*sin(u1)**2, 0],
             [0, 0, 0, cos(Psi)**2]])
ok, bad = is_constant_curvature([Psi, u1, u2, ell], g4, 1)
ck("k=4 join metric: sec == 1", ok, str(bad))

print("== Part A2: pullback of the round metric under the join map ==")
# Phi(Psi,u,ell) = (sin Psi cos u, sin Psi sin u, cos Psi * x(ell)),
# with x.x = 1, x.x' = 0, x'.x' = 1 (arc length). Verify pullback = g3.
from sympy import Function
n_amb = 5  # abstract; use inner-product symbols instead of components
xx, xxp, xpxp = S(1), S(0), S(1)
# pullback entries via bilinearity:
# dPhi_Psi = (cosPsi cos u, cosPsi sin u, -sinPsi x)
# dPhi_u   = (-sinPsi sin u, sinPsi cos u, 0)
# dPhi_ell = (0, 0, cosPsi x')
E_PsiPsi = simplify(cos(Psi)**2*cos(u1)**2 + cos(Psi)**2*sin(u1)**2 + sin(Psi)**2*xx)
E_uu     = simplify(sin(Psi)**2*sin(u1)**2 + sin(Psi)**2*cos(u1)**2)
E_ll     = simplify(cos(Psi)**2*xpxp)
E_Psiu   = simplify(-cos(Psi)*cos(u1)*sin(Psi)*sin(u1) + cos(Psi)*sin(u1)*sin(Psi)*cos(u1))
E_Psil   = simplify(-sin(Psi)*cos(Psi)*xxp)
E_ul     = S(0)
ck("pullback: E_PsiPsi = 1", simplify(E_PsiPsi - 1) == 0)
ck("pullback: E_uu = sin^2 Psi", simplify(E_uu - sin(Psi)**2) == 0)
ck("pullback: E_ll = cos^2 Psi", simplify(E_ll - cos(Psi)**2) == 0)
ck("pullback: off-diagonals vanish", all(simplify(e) == 0 for e in (E_Psiu, E_Psil, E_ul)))

# ---------------------------------------------------------------- Part B: validated numerics
print("== Part B: validated numeric sectional curvature ==")
import mpmath as mp
mp.mp.dps = 30

def riemann_fd(gfun, x0, n):
    h = mp.mpf(10)**(-7)
    def g(x): return gfun(x)
    def dg(i, x):
        xp = list(x); xm = list(x)
        xp[i] += h; xm[i] -= h
        Gp, Gm = g(xp), g(xm)
        return [[(Gp[a][b]-Gm[a][b])/(2*h) for b in range(n)] for a in range(n)]
    def d2g(i, j, x):
        xp = list(x); xm = list(x)
        xp[j] += h; xm[j] -= h
        Dp, Dm = dg(i, xp), dg(i, xm)
        return [[(Dp[a][b]-Dm[a][b])/(2*h) for b in range(n)] for a in range(n)]
    G0 = g(x0)
    Gm_ = mp.matrix(G0); Ginv = Gm_**-1
    dG = [dg(i, x0) for i in range(n)]
    Gam_l = [[[ (dG[j][m_][k] + dG[k][m_][j] - dG[m_][j][k])/2
               for k in range(n)] for j in range(n)] for m_ in range(n)]
    Gam = [[[ sum(Ginv[i, m_]*Gam_l[m_][j][k] for m_ in range(n))
             for k in range(n)] for j in range(n)] for i in range(n)]
    # dGamma_l via FD of Gamma-lower then convert: easier direct FD of Gam
    def Gam_at(x):
        Gx = g(x)
        dGx = [dg(i, x) for i in range(n)]
        Gmx = mp.matrix(Gx); Gix = Gmx**-1
        Gl = [[[ (dGx[j][m_][k] + dGx[k][m_][j] - dGx[m_][j][k])/2
                for k in range(n)] for j in range(n)] for m_ in range(n)]
        return [[[ sum(Gix[i, m_]*Gl[m_][j][k] for m_ in range(n))
                  for k in range(n)] for j in range(n)] for i in range(n)]
    dGam = []
    for kdir in range(n):
        xp = list(x0); xm = list(x0)
        xp[kdir] += h; xm[kdir] -= h
        Gp, Gm2 = Gam_at(xp), Gam_at(xm)
        dGam.append([[[ (Gp[i][j][k]-Gm2[i][j][k])/(2*h) for k in range(n)]
                     for j in range(n)] for i in range(n)])
    R = [[[[mp.mpf(0)]*n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    expr = dGam[k][i][l][j] - dGam[l][i][k][j]
                    for m_ in range(n):
                        expr += Gam[i][k][m_]*Gam[m_][l][j] - Gam[i][l][m_]*Gam[m_][k][j]
                    R[i][j][k][l] = expr
    Rlow = [[[[sum(G0[i][m_]*R[m_][j][k][l] for m_ in range(n)) for l in range(n)]
              for k in range(n)] for j in range(n)] for i in range(n)]
    return Rlow, G0

def sec_from_R(Rlow, G, u, v, n):
    def q(w1, w2, w3, w4):
        s = mp.mpf(0)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        s += Rlow[i][j][k][l]*w1[i]*w2[j]*w3[k]*w4[l]
        return s
    def ip(w1, w2):
        return sum(G[i][j]*w1[i]*w2[j] for i in range(n) for j in range(n))
    num = q(u, v, u, v)
    den = ip(u, u)*ip(v, v) - ip(u, v)**2
    return num/den

# validation metrics
def g_sphere(x):   # radius 1, coords (th, ph): diag(1, sin^2 th)
    th = x[0]
    return [[mp.mpf(1), 0], [0, mp.sin(th)**2]]
def g_hyper(x):    # upper half plane: diag(1/y^2, 1/y^2)
    y = x[1]
    return [[1/y**2, 0], [0, 1/y**2]]
def g_flat(x):
    return [[mp.mpf(1), 0], [0, mp.mpf(1)]]

for name, gf, x0, expect in [("unit sphere", g_sphere, [mp.mpf('1.1'), mp.mpf('0.7')], 1),
                             ("hyperbolic plane", g_hyper, [mp.mpf('0.4'), mp.mpf('1.3')], -1),
                             ("flat plane", g_flat, [mp.mpf('0.2'), mp.mpf('0.5')], 0)]:
    R, G = riemann_fd(gf, x0, 2)
    s = sec_from_R(R, G, [1, 0], [0, 1], 2)
    ck(f"FD routine on {name}: sec = {expect}", abs(s - expect) < mp.mpf(10)**-4, str(s))

# cumulant closed form for exponential families (Amari):
# R_ijkl = (1/4) g^{mn} ( T_ikm T_jln - T_ilm T_jkn ),  T = third central moments
import math
LOG2, LOG3, LOG5 = mp.log(2), mp.log(3), mp.log(5)
LPHI = mp.log((1+mp.sqrt(5))/2)
a_num = [LOG2, LOG3, LOG5, LPHI, LPHI, 4*LPHI, LOG5/2 + 2*LPHI]

def expfam_R(Tstats, theta):
    kdim = len(Tstats); mm = len(Tstats[0])
    logits = [sum(theta[al]*Tstats[al][i] for al in range(kdim)) for i in range(mm)]
    mx = max(logits)
    ws = [mp.e**(l - mx) for l in logits]
    Z = sum(ws); p = [w/Z for w in ws]
    mu = [sum(p[i]*Tstats[al][i] for i in range(mm)) for al in range(kdim)]
    tau = [[Tstats[al][i]-mu[al] for i in range(mm)] for al in range(kdim)]
    g = [[sum(p[i]*tau[a1][i]*tau[a2][i] for i in range(mm)) for a2 in range(kdim)] for a1 in range(kdim)]
    T3 = [[[sum(p[i]*tau[a1][i]*tau[a2][i]*tau[a3][i] for i in range(mm))
            for a3 in range(kdim)] for a2 in range(kdim)] for a1 in range(kdim)]
    Gm_ = mp.matrix(g); Gi = Gm_**-1
    R = [[[[mp.mpf(0)]*kdim for _ in range(kdim)] for _ in range(kdim)] for _ in range(kdim)]
    for i in range(kdim):
        for j in range(kdim):
            for k in range(kdim):
                for l in range(kdim):
                    s = mp.mpf(0)
                    for m_ in range(kdim):
                        for n_ in range(kdim):
                            s += Gi[m_, n_]*(T3[i][l][m_]*T3[j][k][n_] - T3[i][k][m_]*T3[j][l][n_])
                    R[i][j][k][l] = s/4
    return R, g

def gfun_expfam(Tstats):
    kdim = len(Tstats); mm = len(Tstats[0])
    def gf(x):
        logits = [sum(x[al]*Tstats[al][i] for al in range(kdim)) for i in range(mm)]
        mx = max(logits)
        ws = [mp.e**(l - mx) for l in logits]
        Z = sum(ws); p = [w/Z for w in ws]
        mu = [sum(p[i]*Tstats[al][i] for i in range(mm)) for al in range(kdim)]
        return [[sum(p[i]*(Tstats[a1][i]-mu[a1])*(Tstats[a2][i]-mu[a2]) for i in range(mm))
                 for a2 in range(kdim)] for a1 in range(kdim)]
    return gf

IND = lambda S_, mm=7: [1 if i in S_ else 0 for i in range(mm)]
T_cat3 = [a_num, IND({6}), IND({5})]      # (log M, 1_K, 1_{phi4})
theta0 = [mp.mpf('0.31'), mp.mpf('-0.42'), mp.mpf('0.27')]
R_cf, g_cf = expfam_R(T_cat3, theta0)
R_fd, g_fd = riemann_fd(gfun_expfam(T_cat3), theta0, 3)
maxdev = max(abs(R_cf[i][j][k][l]-R_fd[i][j][k][l])
             for i in range(3) for j in range(3) for k in range(3) for l in range(3))
ck("closed-form Riemann matches FD on the catalog 3-family", maxdev < mp.mpf(10)**-4, str(maxdev))

# full simplex control: statistics = 6 singleton indicators on 7 pts -> sec 1/4
T_simplex = [IND({i}) for i in range(6)]
Rs, gs = expfam_R(T_simplex, [mp.mpf('0.1')*i - mp.mpf('0.2') for i in range(6)])
import random
random.seed(11)
devs = []
for _ in range(6):
    u = [mp.mpf(random.uniform(-1, 1)) for _ in range(6)]
    v = [mp.mpf(random.uniform(-1, 1)) for _ in range(6)]
    devs.append(abs(sec_from_R(Rs, gs, u, v, 6) - mp.mpf(1)/4))
ck("simplex control: random sectional curvatures = 1/4", max(devs) < mp.mpf(10)**-18, str(max(devs)))

# product control: two independent binary factors -> cross-plane sec = 0
T_prod = [[1,1,0,0], [1,0,1,0]]
Rp, gp = expfam_R(T_prod, [mp.mpf('0.3'), mp.mpf('-0.7')])
ck("product control: cross-plane sectional curvature = 0",
   abs(sec_from_R(Rp, gp, [1,0], [0,1], 2)) < mp.mpf(10)**-20)

# k=2 suspension control: (log M, 1_K) -> 1/4  (established)
T_susp = [a_num, IND({6})]
Rk2, gk2 = expfam_R(T_susp, [mp.mpf('0.21'), mp.mpf('0.65')])
ck("k=2 suspension control: sec = 1/4",
   abs(sec_from_R(Rk2, gk2, [1,0], [0,1], 2) - mp.mpf(1)/4) < mp.mpf(10)**-18)

# THE catalog 3-family: random theta, random planes -> 1/4
worst = mp.mpf(0)
for t in range(8):
    th = [mp.mpf(random.uniform(-1.2, 1.2)) for _ in range(3)]
    Rc, gc = expfam_R(T_cat3, th)
    for _ in range(4):
        u = [mp.mpf(random.uniform(-1, 1)) for _ in range(3)]
        v = [mp.mpf(random.uniform(-1, 1)) for _ in range(3)]
        worst = max(worst, abs(sec_from_R(Rc, gc, u, v, 3) - mp.mpf(1)/4))
ck("catalog (log M, 1_K, 1_{phi4}): all sampled sectional curvatures = 1/4",
   worst < mp.mpf(10)**-18, str(worst))

# control: non-indicator third statistic -> sec varies
T_bad = [a_num, IND({6}), [0, 1, 0, 2, 0, 0, 1]]
vals = []
for t in range(4):
    th = [mp.mpf(random.uniform(-1, 1)) for _ in range(3)]
    Rb, gb = expfam_R(T_bad, th)
    vals.append(sec_from_R(Rb, gb, [1,0,0], [0,0,1], 3))
spread = max(vals) - min(vals)
ck("control non-indicator 3-family: sectional curvature is NOT constant",
   spread > mp.mpf(10)**-3, str(spread))

# ---------------------------------------------------------------- Part C: TG threshold (exact)
print("== Part C: totally geodesic threshold ==")
from sympy import symbols as syms
L2s, L3s, L5s, Lps = syms('L2 L3 L5 Lp', positive=True)
a_sym = [L2s, L3s, L5s, Lps, Lps, 4*Lps, Rational(1,2)*L5s + 2*Lps]
def colv(vals): return Matrix(len(vals), 1, list(vals))
def hp(u, v): return colv([u[i]*v[i] for i in range(7)])
def exact_rank(M_): return M_.rank(iszerofunc=lambda x: cancel(x) == 0, simplify=cancel)
ONEc = colv([1]*7); Ac = colv(a_sym)
pows = [ONEc, Ac]
for k in range(2, 6):
    pows.append(hp(pows[-1], Ac))
B = Matrix.hstack(*pows)
ck("R[log M] is 6-dimensional", exact_rank(B) == 6)
for k in range(6, 11):
    nxt = hp(pows[-1], Ac) if k == 6 else hp(nxt, Ac)
    r0 = exact_rank(B)
    ck(f"a^{k} in R[log M] (subalgebra closes)", exact_rank(Matrix.hstack(B, nxt)) == r0)
# R[log M] = functions constant on {phi, tau}
Efun = Matrix.hstack(*[colv([1 if i == j else 0 for i in range(7)]) for j in [0,1,2,5,6]]
                     + [colv([1 if i in (3,4) else 0 for i in range(7)])])
ck("R[log M] = {functions constant on the golden pair} (span equality)",
   exact_rank(Matrix.hstack(B, Efun)) == 6)
# no subalgebra of dim <= 5 contains {1, a}: powers force dim >= 6
ck("no TG for k <= 4: any subalgebra containing log M has dim >= 6 (powers independent)",
   exact_rank(Matrix.hstack(*pows)) == 6)
# and the 3-family (a, 1_K, 1_{phi4}) is NOT totally geodesic: a^2 not in V
V3 = Matrix.hstack(ONEc, Ac, colv(IND({6})), colv(IND({5})))
ck("catalog 3-family not TG: a^2 not in span V",
   exact_rank(Matrix.hstack(V3, hp(Ac, Ac))) == 5)

print(f"\nALL {NCK[0]} CHECKS PASSED (t3_suspension)")
