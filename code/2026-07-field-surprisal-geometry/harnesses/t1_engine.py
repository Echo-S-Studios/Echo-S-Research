#!/usr/bin/env python3
# t1_engine.py -- Task 1 Stage 3: symbolic coefficient engine for P(w).
# P(w) = R(F,G) - R(H,H), coefficients in QQ[X1..X7, L2,L3,L5,Lp].
# Exact sparse-ring arithmetic throughout.
import sys, pickle, itertools
from sympy import QQ, Rational
from sympy.polys.rings import ring

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

m = 7
R, *gens = ring("X1,X2,X3,X4,X5,X6,X7,L2,L3,L5,Lp", QQ)
Xg = gens[:7]
L2v, L3v, L5v, Lpv = gens[7:]
a = [L2v, L3v, L5v, Lpv, Lpv, 4*Lpv, L5v/2 + 2*Lpv]
X = list(Xg)
one = R.one
F = [a[i]*a[i] for i in range(m)]
G = [X[i]*X[i] for i in range(m)]
H = [a[i]*X[i] for i in range(m)]

def Cdict(f, g):
    d = {}
    for i in range(m):
        for j in range(i+1, m):
            c = (f[i]-f[j])*(g[i]-g[j])
            if c:
                d[(i, j)] = c
    return d

Caa = Cdict(a, a); Cxx = Cdict(X, X); Cax = Cdict(a, X)

def triple_terms(fg_pair):
    f, g = fg_pair
    Cfg = Cdict(f, g)
    Cfa = Cdict(f, a); Cfx = Cdict(f, X)
    Cag = Cdict(a, g); Cxg = Cdict(X, g)
    # R(f,g) = Caa*Cxx*Cfg - Cax*Cax*Cfg - Cfa*Cxx*Cag + Cfa*Cax*Cxg + Cfx*Cax*Cag - Cfx*Caa*Cxg
    return [(+1, Caa, Cxx, Cfg), (-1, Cax, Cax, Cfg),
            (-1, Cfa, Cxx, Cag), (+1, Cfa, Cax, Cxg),
            (+1, Cfx, Cax, Cag), (-1, Cfx, Caa, Cxg)]

def accumulate(terms, sign_all, out):
    for sgn, dA, dB, dC in terms:
        s = sgn * sign_all
        for (i1, j1), c1 in dA.items():
            for (i2, j2), c2 in dB.items():
                c12 = c1 * c2
                for (i3, j3), c3 in dC.items():
                    k = [0]*m
                    k[i1] += 1; k[j1] += 1; k[i2] += 1; k[j2] += 1; k[i3] += 1; k[j3] += 1
                    key = tuple(k)
                    v = c12 * c3
                    if s < 0:
                        v = -v
                    if key in out:
                        out[key] = out[key] + v
                    else:
                        out[key] = v

print("building P(w) coefficient dictionary (12 triple convolutions)...")
P = {}
accumulate(triple_terms((F, G)), +1, P)
accumulate(triple_terms((H, H)), -1, P)
P = {k: v for k, v in P.items() if v}
print(f"nonzero monomials: {len(P)}")

# ---------------------------------------------------------------- structure
supp = lambda k: sum(1 for x in k if x > 0)
ck("all monomials have degree 6", all(sum(k) == 6 for k in P))
ck("no monomial with support <= 3 survives", all(supp(k) >= 4 for k in P))
from collections import Counter
shape_count = Counter(tuple(sorted((x for x in k if x), reverse=True)) for k in P)
print("shapes:", dict(shape_count))

# ---------------------------------------------------------------- instance cross-check
print("== instance cross-check against direct rational computation ==")
import random
random.seed(77)
subs_vals = {}
for gvar, val in zip(gens, [Rational(random.randint(-9, 9), random.randint(1, 7)) for _ in range(7)]
                     + [Rational(6931,10000), Rational(10986,10000), Rational(16094,10000), Rational(4812,10000)]):
    subs_vals[gvar] = val
wr = [Rational(random.randint(1, 12)) for _ in range(m)]

def eval_ring(pol, sv):
    tot = Rational(0)
    for monom, coeff in pol.terms():
        t = Rational(coeff.numerator, coeff.denominator) if hasattr(coeff, 'numerator') else Rational(str(coeff))
        for gi, e in zip(gens, monom):
            if e:
                t *= sv[gi]**e
        tot += t
    return tot

Pval = Rational(0)
for k, c in P.items():
    mono = Rational(1)
    for i, e in enumerate(k):
        if e:
            mono *= wr[i]**e
    Pval += mono * eval_ring(c, subs_vals)

# direct route (Matrix based, reuse formulas)
from sympy import cancel, Matrix, S
av = [subs_vals[g_] for g_ in gens[7:]]
a_num = [av[0], av[1], av[2], av[3], av[3], 4*av[3], av[2]/2 + 2*av[3]]
X_num = [subs_vals[g_] for g_ in gens[:7]]
def covC(w, f, g):
    Z = sum(w)
    return cancel(Z*sum(w[i]*f[i]*g[i] for i in range(len(w)))
                  - sum(w[i]*f[i] for i in range(len(w)))*sum(w[i]*g[i] for i in range(len(w))))
Fn = [a_num[i]**2 for i in range(m)]; Gn = [X_num[i]**2 for i in range(m)]; Hn = [a_num[i]*X_num[i] for i in range(m)]
Caa_n = covC(wr, a_num, a_num); Cxx_n = covC(wr, X_num, X_num); Cax_n = covC(wr, a_num, X_num)
Dn = cancel(Caa_n*Cxx_n - Cax_n**2)
def Rnum(f, g):
    return cancel(Dn*covC(wr, f, g)
                  - (covC(wr, f, a_num)*Cxx_n*covC(wr, a_num, g)
                     - covC(wr, f, a_num)*Cax_n*covC(wr, X_num, g)
                     - covC(wr, f, X_num)*Cax_n*covC(wr, a_num, g)
                     + covC(wr, f, X_num)*Caa_n*covC(wr, X_num, g)))
Pdirect = cancel(Rnum(Fn, Gn) - Rnum(Hn, Hn))
ck("engine P-dict evaluates identically to direct formula", cancel(Pval - Pdirect) == 0)

# ---------------------------------------------------------------- extreme monomial vs window invariant
print("== extreme (3,1,1,1) coefficients vs 4-window invariant q4 ==")
def minor3(rows):
    (r1, r2, r3) = rows
    return (r1[0]*(r2[1]*r3[2]-r2[2]*r3[1])
            - r1[1]*(r2[0]*r3[2]-r2[2]*r3[0])
            + r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]))

def q4_ring(s):
    rows = [(one, a[i], X[i]) for i in s]
    ell = []
    for t in range(4):
        rest = [rows[u] for u in range(4) if u != t]
        ell.append((-1)**t * minor3(rest))
    aF = sum(ell[t]*F[s[t]] for t in range(4))
    aG = sum(ell[t]*G[s[t]] for t in range(4))
    aH = sum(ell[t]*H[s[t]] for t in range(4))
    return aF*aG - aH*aH

tested = 0
ratios = set()
for s in itertools.combinations(range(m), 4):
    for heavy in s:
        k = [0]*m
        for i in s:
            k[i] = 1
        k[heavy] = 3
        key = tuple(k)
        c = P.get(key, R.zero)
        q4 = q4_ring(s)
        # hypothesis: c = -(w-independent) * q4 ... discover exact cofactor by division
        if q4:
            quo = c.div([q4])
            if not quo[1]:  # exact division
                ratios.add(str(quo[0][0]))
                tested += 1
if tested:
    print(f"exact divisions found: {tested}; distinct cofactors: {len(ratios)}")
    for rr in sorted(ratios):
        print("  cofactor:", rr)
ck("every (3,1,1,1) coefficient is an exact multiple of its window q4", tested == 4*35)

# ---------------------------------------------------------------- nu-grouping census
print("== nu-group census ==")
def nu(k):
    return (k[0], k[1], 2*k[2] + k[6], k[3] + k[4] + 4*k[5] + 2*k[6])
groups = {}
for k in P:
    groups.setdefault(nu(k), []).append(k)
sizes = Counter(len(v) for v in groups.values())
print("group-size histogram:", dict(sizes))
multi = {kk: v for kk, v in groups.items() if len(v) > 1}
print(f"groups with >1 monomial: {len(multi)} of {len(groups)}")
# show which delta separations occur inside multi-groups
d1 = (0,0,0,1,-1,0,0); d2 = (0,0,-1,0,0,-1,2)
def as_comb(delta):
    # solve delta = j*d1 + k*d2 over integers
    kk = delta[6] // 2
    jj = delta[3]
    if tuple(jj*d1[i] + kk*d2[i] for i in range(m)) == delta:
        return (jj, kk)
    return None
seps = Counter()
for kk, v in multi.items():
    for A, B in itertools.combinations(v, 2):
        delta = tuple(B[i]-A[i] for i in range(m))
        comb = as_comb(delta) or as_comb(tuple(-x for x in delta))
        seps[comb] += 1
print("pair separations (j,k) counts:", dict(seps))

with open('P_dict.pkl', 'wb') as fh:
    pickle.dump({'monoms': list(P.keys())}, fh)
# store coefficients as strings for reload safety
with open('P_coeffs.pkl', 'wb') as fh:
    pickle.dump({k: str(v) for k, v in P.items()}, fh)
print(f"\nALL {NCK[0]} CHECKS PASSED (t1_engine); P saved")
