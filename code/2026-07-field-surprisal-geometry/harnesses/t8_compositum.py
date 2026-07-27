#!/usr/bin/env python3
# t8_compositum.py -- Open Problem (4): multi-field coupling over the compositum.
# RESOLUTION (first horn): the compositum's own tensor catalog FORCES coupling.
# Construction: for catalog seeds p, q let p (x) q be the tensor polynomial whose
# roots are ALL conjugate products alpha_i beta_j (gauge-free: no root choice),
# and c(p,q) = log M(p (x) q). Every conjugate log-magnitude is a Q-vector over
# (L2, L3, L5, Lphi); c(p,q) is an exact Q-vector; on-circle products (zero
# vector, e.g. phi*psi = -1) are excluded symbolically; all other sign decisions
# are certified by rigorous interval arithmetic (mpmath.iv, dps 60).
# Findings verified here:
#   [F] interaction contrast on ({sqrt2,phi} x {sqrt3,phi}) = log 2 exactly
#   [F] golden-sector contrast on ({phi,phi4})^2 = -6 log phi exactly
#   [F] rational sector {sqrt2,sqrt3,sqrt5} exactly separable
#   [F] golden tie persists: c(phi,.) = c(tau,.)
#   [F] double-centered interaction tensor Delta != 0; exact rank over Q(L)
#   [C] cross-Fisher Cov_beta(cost1, cost2) != 0 at beta in {1,-1,sqrt5}
#   [F-instance] charge tensor law: n(p (x) q) = lcm(n_p, n_q) on all 28 pairs
import sys, itertools, math
from fractions import Fraction as Fr
from sympy import symbols, Matrix, Rational, S, simplify

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

def v(l2=0, l3=0, l5=0, lp=0):
    return (Fr(l2), Fr(l3), Fr(l5), Fr(lp))
def vadd(a, b): return tuple(x+y for x, y in zip(a, b))
ZERO = v()

# roots: (log-magnitude vector, angle in quarter turns: 0=+R, 1=+iR, 2=-R, 3=-iR)
SEEDS = ["sqrt2", "sqrt3", "sqrt5", "phi", "tau", "phi4", "K"]
ROOTS = [
    [(v(l2=Fr(1,2)), 0), (v(l2=Fr(1,2)), 2)],                       # +-sqrt2
    [(v(l3=Fr(1,2)), 0), (v(l3=Fr(1,2)), 2)],                       # +-sqrt3
    [(v(l5=Fr(1,2)), 0), (v(l5=Fr(1,2)), 2)],                       # +-sqrt5
    [(v(lp=1), 0), (v(lp=-1), 2)],                                  # phi, psi=-1/phi
    [(v(lp=-1), 0), (v(lp=1), 2)],                                  # 1/phi, -phi
    [(v(lp=4), 0), (v(lp=-4), 0)],                                  # phi^4, phi^-4
    [(v(l5=Fr(1,4), lp=-1), 0), (v(l5=Fr(1,4), lp=-1), 2),          # +-r  (real, inside)
     (v(l5=Fr(1,4), lp=1), 1), (v(l5=Fr(1,4), lp=1), 3)],           # +-i beta (outside)
]
N_SEED = [2, 2, 2, 2, 2, 1, 4]

from mpmath import iv, mp
iv.dps = 60
LIV = [iv.log(2), iv.log(3), iv.log(5), iv.log((1 + iv.sqrt(5))/2)]
min_margin = [None]
def sign_of(vec):
    if vec == ZERO:
        return 0
    tot = iv.mpf(0)
    for f, Li in zip(vec, LIV):
        tot += (iv.mpf(f.numerator)/f.denominator)*Li
    if tot.a > 0:
        m = float(tot.a)
    elif tot.b < 0:
        m = float(-tot.b)
    else:
        raise RuntimeError(f"ambiguous sign for {vec}")
    min_margin[0] = m if min_margin[0] is None else min(min_margin[0], m)
    return 1 if tot.a > 0 else -1

oncircle = {}
C = [[None]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        tot = ZERO; zc = 0
        for (u, _) in ROOTS[i]:
            for (w, _) in ROOTS[j]:
                s = vadd(u, w)
                sg = sign_of(s)
                if sg == 0: zc += 1
                elif sg > 0: tot = vadd(tot, s)
        C[i][j] = tot
        if zc: oncircle[(i, j)] = zc

print("== block 1: the exact tensor-cost matrix and its anchors ==")
ck("symmetry c(i,j) = c(j,i) on all 49 pairs",
   all(C[i][j] == C[j][i] for i in range(7) for j in range(7)))
anchors = [((0,1), v(l2=2,l3=2)), ((0,3), v(l2=1,lp=2)), ((3,1), v(l3=2)),
           ((3,3), v(lp=2)), ((3,5), v(lp=8)), ((5,5), v(lp=8)),
           ((0,6), v(l2=4,l5=2)), ((3,6), v(l5=Fr(3,2),lp=4)), ((0,0), v(l2=4))]
ck("hand-derived anchors match: c(sqrt2,sqrt3)=2log6, c(sqrt2,phi)=log(2phi^2), "
   "c(phi,sqrt3)=log9, c(phi,phi)=2Lphi, c(phi,phi4)=c(phi4,phi4)=8Lphi, "
   "c(sqrt2,K)=4L2+2L5, c(phi,K)=(3/2)L5+4Lphi, c(sqrt2,sqrt2)=4L2",
   all(C[a][b] == val for (a, b), val in anchors))
ck("on-circle exclusions are exactly the golden/quartic unit products: "
   "10 zero-vector root pairs on the 5 expected seed pairs",
   sum(oncircle.values()) == 10 and
   set(oncircle) == {(3,3),(3,4),(4,3),(4,4),(5,5)})
ck("all certified sign margins comfortably positive (rigorous intervals)",
   min_margin[0] is not None and min_margin[0] > Fr(1, 100))

def contrast(i, i2, j, j2):
    return tuple(C[i][j][c] + C[i2][j2][c] - C[i][j2][c] - C[i2][j][c]
                 for c in range(4))

print("== block 2: forced coupling -- exact nonzero interaction ==")
ck("WITNESS 1: contrast({sqrt2,phi} x {sqrt3,phi}) = log 2 exactly",
   contrast(0, 3, 1, 3) == v(l2=1), str(contrast(0, 3, 1, 3)))
ck("WITNESS 2: golden-sector contrast({phi,phi4}^2) = -6 log phi exactly",
   contrast(3, 5, 3, 5) == v(lp=-6), str(contrast(3, 5, 3, 5)))
ck("rational sector {sqrt2,sqrt3,sqrt5} exactly separable: all 9 contrasts zero",
   all(contrast(i, i2, j, j2) == ZERO
       for i, i2 in itertools.combinations(range(3), 2)
       for j, j2 in itertools.combinations(range(3), 2)))
ck("golden tie persists in the compositum: c(phi, .) = c(tau, .) on all columns",
   all(C[3][j] == C[4][j] for j in range(7)))

print("== block 3: the interaction tensor Delta (double-centered), exact rank ==")
L2s, L3s, L5s, Lps = symbols('L2 L3 L5 Lp', positive=True)
LS = (L2s, L3s, L5s, Lps)
def sym_of(vec): return sum(Rational(f)*Li for f, Li in zip(vec, LS))
Cm = Matrix(7, 7, lambda i, j: sym_of(C[i][j]))
rm = [sum(Cm[i, j] for j in range(7))/7 for i in range(7)]
cm = [sum(Cm[i, j] for i in range(7))/7 for j in range(7)]
gm = sum(Cm[i, j] for i in range(7) for j in range(7))/49
D = Matrix(7, 7, lambda i, j: simplify(Cm[i, j] - rm[i] - cm[j] + gm))
ck("Delta is symmetric and annihilates constants (all row sums zero)",
   D == D.T and all(simplify(sum(D[i, j] for j in range(7))) == 0 for i in range(7)))
ck("Delta != 0: the compositum cost is NOT additively separable",
   any(D[i, j] != 0 for i in range(7) for j in range(7)))
rk = D.rank()
ck(f"exact rank of Delta over Q(L2,L3,L5,Lphi): rank = {rk} (>= 1 coupling modes)",
   1 <= rk <= 7)

print("== block 4: cross-Fisher covariance on the compositum Gibbs law [C] ==")
mp.dps = 50
Lnum = [mp.log(2), mp.log(3), mp.log(5), mp.log((1 + mp.sqrt(5))/2)]
cnum = [[sum(float(C[i][j][k])*Lnum[k] for k in range(4)) for j in range(7)]
        for i in range(7)]
a1n = [mp.log(2), mp.log(3), mp.log(5), Lnum[3], Lnum[3], 4*Lnum[3],
       mp.log(5)/2 + 2*Lnum[3]]
covs = {}
for name, beta in [("1", mp.mpf(1)), ("-1", mp.mpf(-1)), ("sqrt5", mp.sqrt(5))]:
    wts = [[mp.e**(-beta*cnum[i][j]) for j in range(7)] for i in range(7)]
    Zt = sum(sum(r) for r in wts)
    E1 = sum(a1n[i]*wts[i][j] for i in range(7) for j in range(7))/Zt
    E2 = sum(a1n[j]*wts[i][j] for i in range(7) for j in range(7))/Zt
    E12 = sum(a1n[i]*a1n[j]*wts[i][j] for i in range(7) for j in range(7))/Zt
    covs[name] = E12 - E1*E2
ck("Cov_beta(cost1, cost2) nonzero at beta = 1, -1, sqrt5 "
   f"(values {[float(covs[k]) for k in covs]})",
   all(abs(c) > mp.mpf('1e-9') for c in covs.values()))

print("== block 5: charge tensor law on all 28 unordered pairs [F-instance] ==")
def n_of_angles(angs):
    if all(a == 0 for a in angs): return 1
    if all(a in (0, 2) for a in angs): return 2
    return 4
ok = True
for i in range(7):
    for j in range(i, 7):
        angs = [(u[1] + w[1]) % 4 for u in ROOTS[i] for w in ROOTS[j]]
        if n_of_angles(angs) != math.lcm(N_SEED[i], N_SEED[j]):
            ok = False; print("   mismatch at", SEEDS[i], SEEDS[j])
ck("n(p (x) q) = lcm(n_p, n_q) on all 28 unordered pairs "
   "(every product conjugate is real or purely imaginary)", ok)

print("\ncost matrix c(i,j) as (L2,L3,L5,Lphi)-vectors:")
for i in range(7):
    print("  ", SEEDS[i], [tuple(str(x) for x in C[i][j]) for j in range(7)])
print(f"\nALL {NCK[0]} CHECKS PASSED (t8_compositum)")
print("OP (4) resolves on the first horn: the compositum forces the coupling;")
print("interaction quantum log 2 on the golden/rational window, -6 log phi in the")
print("golden sector; the rational sector alone stays separable.")
