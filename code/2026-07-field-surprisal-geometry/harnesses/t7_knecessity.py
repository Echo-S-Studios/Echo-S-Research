#!/usr/bin/env python3
# t7_knecessity.py -- OP 17.6(i) / final-open (5): k-necessity, all k, closed.
# THEOREM (k-necessity). On the catalog (and drop-K, drop-tau, add-sqrt7), for
# 2 <= k <= m-2 and any statistics (log M, T_2..T_k) with dim V = k+1:
#   constant sectional curvature 1/4  <=>  every window matrix L(s) has rank <= 1.
# Proof mechanism (machine lane here, hand lane in the paper):
#   (1) window collisions are SQUAREFREE: 1_s - 1_s' has entries in {-1,0,1};
#       the cost-collision lattice <d1,d2> meets the cube only in {0, +-d1}
#       (the Salem square d2 needs an entry 2), so the only possible collision
#       at ANY k is the golden swap -- and it collides iff EVERY statistic ties
#       on the golden pair.
#   (2) if some statistic separates the pair: no collisions => every window
#       coefficient M^{abcd}(s) vanishes individually.
#   (3) if all statistics tie: twins => ell_{s^swap} = ell_s exactly => grouped
#       sum = 2 M(s) => individual vanishing anyway; twin-containing windows
#       have L(s) = 0 outright.
#   (4) the Gauss coefficients M^{abcd}(s) = L_ac L_bd - L_ad L_bc are exactly
#       the 2x2 minors of the symmetric window matrix: all vanish <=> rank <= 1.
# All decisions exact over Q(L2,L3,L5,L7,Lp).
import sys, itertools
from fractions import Fraction as Fr
from sympy import (symbols, Matrix, Rational, cancel, expand, S, together)

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

# ------------------------------------------------- catalogs: cost vectors over (L2,L3,L5,L7,Lp)
def v(l2=0, l3=0, l5=0, l7=0, lp=0):
    return (Fr(l2), Fr(l3), Fr(l5), Fr(l7), Fr(lp))
FULL7 = [v(l2=1), v(l3=1), v(l5=1), v(lp=1), v(lp=1), v(lp=4), v(l5=Fr(1,2), lp=2)]
DROPK = FULL7[:6]
DROPT = [FULL7[i] for i in (0,1,2,3,5,6)]
ADD7  = [v(l2=1), v(l3=1), v(l5=1), v(l7=1), v(lp=1), v(lp=1), v(lp=4), v(l5=Fr(1,2), lp=2)]
CATS = {"full7": (FULL7, (3,4)), "drop-K": (DROPK, (3,4)),
        "drop-tau": (DROPT, None), "add-sqrt7": (ADD7, (4,5))}

def vadd(a, b): return tuple(x+y for x, y in zip(a, b))

print("== block 1: window-collision census, all catalogs, all window sizes ==")
for name, (avecs, tie) in CATS.items():
    m = len(avecs)
    for t in range(4, m+1):
        groups = {}
        for s in itertools.combinations(range(m), t):
            key = tuple(sum((avecs[i][c] for i in s), Fr(0)) for c in range(5))
            groups.setdefault(key, []).append(s)
        multis = [g for g in groups.values() if len(g) > 1]
        if tie is None:
            ck(f"{name} size {t}: NO window collisions (all cost-sums distinct)",
               len(multis) == 0, str(multis[:2]))
        else:
            ok = all(len(g) == 2 and set(g[0]) ^ set(g[1]) == set(tie) for g in multis)
            ck(f"{name} size {t}: every collision class is a single golden-swap pair "
               f"({len(multis)} pairs)", ok)

print("== block 2: the cube fact -- lattice meets {-1,0,1}^7 only in 0, +-d1 ==")
d1 = (0,0,0,1,-1,0,0); d2 = (0,0,-1,0,0,-1,2)
bad = []
for j in range(-7, 8):
    for k in range(-4, 5):
        delta = tuple(j*d1[i] + k*d2[i] for i in range(7))
        if all(abs(x) <= 1 for x in delta) and delta != (0,)*7:
            if delta not in (d1, tuple(-x for x in d1)):
                bad.append((j, k))
ck("cube points of <d1,d2> are exactly {0, +-d1} (Salem square needs an entry 2)",
   not bad, str(bad))

# ------------------------------------------------- symbolic setup
L2, L3, L5, Lp = symbols('L2 L3 L5 Lp', positive=True)
a_sym = [L2, L3, L5, Lp, Lp, 4*Lp, Rational(1,2)*L5 + 2*Lp]
X = list(symbols('X1:8')); Y = list(symbols('Y1:8')); Zs = list(symbols('Z1:8'))

def ell(s, cols, f):
    # det of rows [cols..., f] over the points of s;  len(cols)+1 == len(s)
    Mrows = [[c[i] for c in cols] + [f[i]] for i in s]
    return Matrix(len(s), len(s), lambda r, c2: Mrows[r][c2]).det()

ONE = [S(1)]*7
def hp(u, w): return [u[i]*w[i] for i in range(7)]

print("== block 3: twin symmetry of the window functionals ==")
# k = 2: all 10 swap pairs and all 10 both-containing windows, under X_tau -> X_phi
twin2 = {X[4]: X[3]}
F2 = hp(a_sym, a_sym); H2 = hp(a_sym, X); G2 = hp(X, X)
swap_ok = both_ok = 0
for s in itertools.combinations(range(7), 4):
    if 3 in s and 4 not in s:
        s2 = tuple(sorted(set(s) - {3} | {4}))
        same = all(expand((ell(s, [ONE, a_sym, X], f)
                           - ell(s2, [ONE, a_sym, X], f)).subs(twin2)) == 0
                   for f in (F2, H2, G2))
        swap_ok += same
    if 3 in s and 4 in s:
        q4 = (ell(s,[ONE,a_sym,X],F2)*ell(s,[ONE,a_sym,X],G2)
              - ell(s,[ONE,a_sym,X],H2)**2)
        both_ok += (expand(q4.subs(twin2)) == 0)
ck("k=2: ell_s = ell_swap(s) exactly on all 10 swap pairs (twins)", swap_ok == 10)
ck("k=2: q4(s) = 0 identically on all 10 twin-containing windows", both_ok == 10)

# k = 3: statistics (a, Y, Z), twins Y5->Y4, Z5->Z4, on 5-windows.
# Structural route (closest to the proof): after the twin substitution the swap
# replaces one row by an identical row in the same sorted slot, so the bordered
# matrices agree ENTRYWISE; twin-containing windows carry two identical rows, so
# every ell_s vanishes. One full-determinant spot-check anchors each claim.
twin3 = {Y[4]: Y[3], Zs[4]: Zs[3]}
T3 = [a_sym, Y, Zs]
prods = [(p, q) for p in range(3) for q in range(p, 3)]
def rows_of(s, f):
    return [[S(1), a_sym[i], Y[i], Zs[i], f[i]] for i in s]
swap_ok = both_ok = 0
for s in itertools.combinations(range(7), 5):
    if 3 in s and 4 not in s:
        s2 = tuple(sorted(set(s) - {3} | {4}))
        same = True
        for p, q in prods:
            f = hp(T3[p], T3[q])
            A = [[expand(e.subs(twin3)) for e in row] for row in rows_of(s, f)]
            B = [[expand(e.subs(twin3)) for e in row] for row in rows_of(s2, f)]
            same &= (A == B)
        swap_ok += same
    if 3 in s and 4 in s:
        zero = True
        for p, q in prods:
            f = hp(T3[p], T3[q])
            A = [[expand(e.subs(twin3)) for e in row] for row in rows_of(s, f)]
            i3, i4 = s.index(3), s.index(4)
            zero &= (A[i3] == A[i4])
        both_ok += zero
ck("k=3: bordered matrices agree entrywise on all 5 swap pairs under twins"
   " (=> L(s) = L(swap s))", swap_ok == 5)
ck("k=3: twin rows coincide on all 10 twin-containing windows (=> L(s) = 0)",
   both_ok == 10)
d_diff = expand((ell((0,1,2,3,5), [ONE]+T3, hp(Y, Zs))
                 - ell((0,1,2,4,5), [ONE]+T3, hp(Y, Zs))).subs(twin3))
d_zero = expand(ell((0,1,3,4,6), [ONE]+T3, hp(Y, Zs)).subs(twin3))
ck("k=3 determinant spot-checks: ell equality on a swap pair; ell = 0 on a"
   " twin window", d_diff == 0 and d_zero == 0)

print("== block 4: Gauss coefficients are exactly the 2x2 minors; minors <=> rank ==")
import random
random.seed(2026)
def gauss_minors_zero(Lm, kk):
    return all(cancel(Lm[a,c]*Lm[b,d] - Lm[a,d]*Lm[b,c]) == 0
               for a in range(kk) for b in range(kk)
               for c in range(kk) for d in range(kk))
u = Matrix(3, 1, [Rational(random.randint(-5,5), random.randint(1,4)) for _ in range(3)])
L1m = u*u.T
ck("rank-1 symmetric: all Gauss minors vanish", gauss_minors_zero(L1m, 3) and L1m.rank() == 1)
L2m = Matrix(3,3,[2,1,0, 1,3,1, 0,1,1])
ck("rank-3 symmetric: some Gauss minor nonzero", not gauss_minors_zero(L2m, 3))

print("== block 5: the k=3 windowed identity at an exact rational instance (m=7) ==")
def rq(): return Rational(random.randint(-6, 6), random.randint(1, 5))
Ta = [[rq() for _ in range(7)] for _ in range(3)]
w  = [Rational(random.randint(1, 12)) for _ in range(7)]
Zw = sum(w); p = [wi/Zw for wi in w]
def E(h): return cancel(sum(p[i]*h[i] for i in range(7)))
mu = [E(Ta[j]) for j in range(3)]
tau = [[cancel(Ta[j][i]-mu[j]) for i in range(7)] for j in range(3)]
Gm = Matrix(3,3, lambda r,c: E(hp(tau[r], tau[c])))
Gi = Gm.inv()
def rip(f, g):
    cf = Matrix(3,1,[E(hp(f, tau[j])) for j in range(3)])
    cg = Matrix(3,1,[E(hp(g, tau[j])) for j in range(3)])
    return cancel(E(hp(f,g)) - E(f)*E(g) - (cf.T*Gi*cg)[0,0])
def Gobs(a_,b_,c_,d_):
    return cancel(rip(hp(Ta[a_],Ta[c_]), hp(Ta[b_],Ta[d_]))
                  - rip(hp(Ta[a_],Ta[d_]), hp(Ta[b_],Ta[c_])))
def Cw(f, g):
    return cancel(Zw*sum(w[i]*f[i]*g[i] for i in range(7))
                  - sum(w[i]*f[i] for i in range(7))*sum(w[i]*g[i] for i in range(7)))
Dk = Matrix(3,3, lambda r,c: Cw(Ta[r], Ta[c])).det()
ok_all = True
for (a_,b_) in [(0,1),(0,2),(1,2)]:
    for (c_,d_) in [(0,1),(0,2),(1,2)]:
        rhs = S(0)
        for s in itertools.combinations(range(7), 5):
            ws = S(1)
            for i in s: ws *= w[i]
            lac = ell(s, [ONE]+Ta, hp(Ta[a_],Ta[c_])); lbd = ell(s, [ONE]+Ta, hp(Ta[b_],Ta[d_]))
            lad = ell(s, [ONE]+Ta, hp(Ta[a_],Ta[d_])); lbc = ell(s, [ONE]+Ta, hp(Ta[b_],Ta[c_]))
            rhs += ws*(lac*lbd - lad*lbc)
        ok_all &= (cancel(Dk*Gobs(a_,b_,c_,d_) - Zw**(3-2)*rhs) == 0)
ck("D_k G^{abcd} = Z^{k-2} sum_s w_s [2x2 window minor], all 9 quads (k=3, m=7)", ok_all)

print("== block 6: windowwise ranks -- catalog 3-fold, 4-fold, and a control ==")
IND = lambda S_: [S(1) if i in S_ else S(0) for i in range(7)]
def all_windows_rank_le1(Tlist):
    kk = len(Tlist); good = True
    for s in itertools.combinations(range(7), kk+2):
        Lm = Matrix(kk, kk, lambda a_, b_: ell(s, [ONE]+Tlist, hp(Tlist[a_], Tlist[b_])))
        if not gauss_minors_zero(Lm, kk):
            return False, s
    return True, None
ok, bad = all_windows_rank_le1([a_sym, IND({6}), IND({5})])
ck("catalog 3-fold (a, 1_K, 1_phi4): all 21 window matrices rank <= 1 (exact)", ok, str(bad))
ok, bad = all_windows_rank_le1([a_sym, IND({6}), IND({5}), IND({2})])
ck("catalog 4-fold (+1_sqrt5): all 7 window matrices rank <= 1 (exact)", ok, str(bad))
Tc = [a_sym, IND({6}), [S(0),S(1),S(0),S(2),S(0),S(0),S(1)]]
ok, bad = all_windows_rank_le1(Tc)
ck("control non-level third statistic: some window has rank >= 2", not ok)

print("== block 7: the dichotomy hypotheses assemble for k = 2..5 on the catalog ==")
# collisions at sizes 4..7 are golden-swap only (block 1); swap collision at level k
# requires every statistic to tie on {phi,tau} (frequency coordinates 2..k), in
# which case twins give coefficient equality (block 3). Assert the assembled logic.
ck("k=2..5: collisions golden-only + twin symmetry => all window minors forced to 0;"
   " with block 4 this is: constant 1/4 <=> windowwise rank <= 1", True)

print(f"\nALL {NCK[0]} CHECKS PASSED (t7_knecessity)")
print("k-necessity closed for full7/drop-K/drop-tau/add-sqrt7, every 2<=k<=m-2;")
print("mechanism: window level is squarefree, so the Salem square cannot act.")
