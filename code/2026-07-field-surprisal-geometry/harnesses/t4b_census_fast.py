#!/usr/bin/env python3
# t4b_census_fast.py -- LANE 2 of the k=3 double-indicator census (OP-2).
#
# Independent re-derivation and fast exact implementation of block 3 of t4_kwindows.py.
# The decision "windowwise-flat" for a pair {S,T} reduces (witness lemma, dim V = 4) to
# three span memberships over K = QQ(L2,L3,L5,Lp):
#     a*1_S , a*1_T , 1_{S cap T}   in   V = span_K(1, a, 1_S, 1_T).
# LANE-2 DERIVATION (this file): write the membership system
#     c0 + c1*a_i + sigma_i*c2 + tau_i*c3 = rhs_i ,  (sigma,tau) = (chS,chT)(i),
# collapse rows to (block, level) classes, and observe:
#   Step 1 (c1-forcing): any block containing two distinct cost levels forces c1 by a
#     row difference:  c1*(alpha-alpha') = rhs-difference, i.e. c1 = sigma_B (target a*1_S),
#     c1 = tau_B (target a*1_T), c1 = 0 (target 1_{S cap T}).  Conflicting forcings => NOT in span.
#   Step 2, Case A (c1 forced to v in {0,1}): the remaining system A*(c0,c2,c3) = b has a
#     RATIONAL coefficient matrix A (rows [1,sigma,tau]) and b with entries in the QQ^5 space
#     spanned by (1,L2,L3,L5,Lp).  Solvable over K  <=>  every rational left-null vector y of A
#     annihilates b:  y.b = 0 as a QQ^5 vector.  Exact, no symbolic algebra.
#   Step 2, Case B (no multi-level block; c1 unknown): for each left-null y of A put
#     u_y = y.alpha, v_y = y.rhs (QQ^5 linear forms).  Solvability over K  <=>
#     (u_y = 0 => v_y = 0) for all y, and the cross-ratios u_y*v_z - u_z*v_y vanish
#     identically (degree-2 forms in the L-symbols; exact 15-dim rational zero test).
# All decisions are therefore exact rational tests; floats never appear.
#
# The census law (derived by hand, both lanes): a nondegenerate pair is windowwise-flat
# iff V = <1, a, 1_S', 1_T'> with S',T' DISJOINT subsets each inside a single cost level;
# the number of distinct such spans for the catalog levels {r2}{r3}{r5}{phi,tau}{phi4}{K}
# is C(8,2) - 2 = 26 (the two non-disjoint golden pairs are excluded; the lone collapse
# 1_phi + 1_tau = 1_{phi,tau} lives inside the single all-golden class).
#
# Self-validation: sampled cross-check of the fast membership routine against the lane-1
# symbolic-rank membership (sympy, iszerofunc = cancel o expand) on mixed-verdict pairs.

import sys, itertools, random
from fractions import Fraction as Fr

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    tag = f"[{NCK[0]:03d}]"
    if not cond:
        print(f"{tag} FAIL {label} {detail}")
        sys.exit(1)
    print(f"{tag} PASS {label}")

# ---------------- QQ^5 linear forms over basis (1, L2, L3, L5, Lp) ----------------
B = 5
def lin(*c): return tuple(Fr(x) for x in c)
ZERO = lin(0, 0, 0, 0, 0)
def ladd(u, v): return tuple(a + b for a, b in zip(u, v))
def lsub(u, v): return tuple(a - b for a, b in zip(u, v))
def lscale(q, u): return tuple(q * a for a in u)
def lzero(u): return all(x == 0 for x in u)

def lmul(u, v):
    # product of two linear forms: dict {(i,j) i<=j : Fr}, exact
    d = {}
    for i, a in enumerate(u):
        if a == 0: continue
        for j, b in enumerate(v):
            if b == 0: continue
            k = (i, j) if i <= j else (j, i)
            d[k] = d.get(k, Fr(0)) + a * b
    return d
def qsub(d1, d2):
    d = dict(d1)
    for k, x in d2.items(): d[k] = d.get(k, Fr(0)) - x
    return d
def qzero(d): return all(x == 0 for x in d.values())

# ---------------- catalog ----------------
IDX = {"r2": 0, "r3": 1, "r5": 2, "phi": 3, "tau": 4, "phi4": 5, "K": 6}
m = 7
AVAL = [lin(0, 1, 0, 0, 0), lin(0, 0, 1, 0, 0), lin(0, 0, 0, 1, 0),
        lin(0, 0, 0, 0, 1), lin(0, 0, 0, 0, 1), lin(0, 0, 0, 0, 4),
        (Fr(0), Fr(0), Fr(0), Fr(1, 2), Fr(2))]
LEVEL = [0, 1, 2, 3, 3, 4, 5]
LEVELS = [(0,), (1,), (2,), (3, 4), (5,), (6,)]

ck("catalog sanity: a_phi = a_tau and the six level values are pairwise distinct",
   AVAL[3] == AVAL[4] and len({AVAL[i] for i in (0, 1, 2, 3, 5, 6)}) == 6)

# ---------------- small exact linear algebra over QQ ----------------
def rref_rows(rows):
    M = [list(r) for r in rows]
    if not M: return M
    R, C = len(M), len(M[0]); pr = 0
    for c in range(C):
        piv = next((r for r in range(pr, R) if M[r][c] != 0), None)
        if piv is None: continue
        M[pr], M[piv] = M[piv], M[pr]
        inv = M[pr][c]
        M[pr] = [x / inv for x in M[pr]]
        for r in range(R):
            if r != pr and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[pr])]
        pr += 1
        if pr == R: break
    return M
def qrank(rows):
    return sum(1 for r in rref_rows(rows) if any(x != 0 for x in r))
def nullspace(eqs, n):
    # basis of {x in QQ^n : E x = 0}
    M = rref_rows(eqs) if eqs else []
    piv = {}
    for r in M:
        c = next((j for j, x in enumerate(r) if x != 0), None)
        if c is not None: piv[c] = r
    free = [j for j in range(n) if j not in piv]
    basis = []
    for f in free:
        x = [Fr(0)] * n; x[f] = Fr(1)
        for c, r in piv.items():
            x[c] = -r[f]
        basis.append(x)
    return basis

# ---------------- the membership routine (lane 2) ----------------
def membership(chS, chT, target):
    """target in {'aS','aT','ST'}: is the product in span_K(1, a, 1_S, 1_T)?"""
    # collapse to (block, level) row classes
    classes = {}
    for i in range(m):
        key = (chS[i], chT[i], LEVEL[i])
        classes[key] = classes.get(key, 0) + 1
    # rows: (sigma, tau, alpha, rhs)
    rows = []
    blocks = {}
    for (s, t, lv) in classes:
        alpha = AVAL[[i for i in range(m) if LEVEL[i] == lv][0]]
        if target == 'aS':   rhs = lscale(Fr(s), alpha)
        elif target == 'aT': rhs = lscale(Fr(t), alpha)
        else:                rhs = lin(s * t, 0, 0, 0, 0)
        rows.append((s, t, alpha, rhs))
        blocks.setdefault((s, t), []).append(lv)
    # Step 1: c1-forcing from multi-level blocks
    forced = set()
    for (s, t), lvs in blocks.items():
        if len(lvs) >= 2:
            if target == 'aS':   forced.add(s)
            elif target == 'aT': forced.add(t)
            else:                forced.add(0)
    if len(forced) >= 2:
        return False
    A = [[Fr(1), Fr(s), Fr(t)] for (s, t, _, _) in rows]
    AT = [[A[r][c] for r in range(len(A))] for c in range(3)]
    ynull = nullspace(AT, len(A))
    if len(forced) == 1:
        v = forced.pop()
        bvec = [lsub(rhs, lscale(Fr(v), alpha)) for (_, _, alpha, rhs) in rows]
        for y in ynull:
            acc = ZERO
            for yr, br in zip(y, bvec):
                if yr != 0: acc = ladd(acc, lscale(yr, br))
            if not lzero(acc): return False
        return True
    # Case B: c1 unknown in K
    uv = []
    for y in ynull:
        u, w = ZERO, ZERO
        for yr, (s, t, alpha, rhs) in zip(y, rows):
            if yr != 0:
                u = ladd(u, lscale(yr, alpha))
                w = ladd(w, lscale(yr, rhs))
        if lzero(u):
            if not lzero(w): return False
        else:
            uv.append((u, w))
    for i in range(len(uv)):
        for j in range(i + 1, len(uv)):
            if not qzero(qsub(lmul(uv[i][0], uv[j][1]), lmul(uv[j][0], uv[i][1]))):
                return False
    return True

# ---------------- the census ----------------
subsets = [tuple(i for i in range(m) if msk >> i & 1) for msk in range(1, 2 ** m - 1)]
ck("126 nonempty proper subsets; 7875 unordered pairs iterated",
   len(subsets) == 126 and len(subsets) * (len(subsets) - 1) // 2 == 7875)

def is_dim4(chS, chT):
    if qrank([[Fr(1)] * m, [Fr(x) for x in chS], [Fr(x) for x in chT]]) != 3:
        return False
    blocks = {}
    for i in range(m):
        blocks.setdefault((chS[i], chT[i]), set()).add(LEVEL[i])
    return any(len(v) >= 2 for v in blocks.values())

passing, dim4_pairs = [], []
for iS in range(len(subsets)):
    for iT in range(iS + 1, len(subsets)):
        Sset, Tset = subsets[iS], subsets[iT]
        chS = [1 if i in Sset else 0 for i in range(m)]
        chT = [1 if i in Tset else 0 for i in range(m)]
        if not is_dim4(chS, chT): continue
        dim4_pairs.append((Sset, Tset))
        if (membership(chS, chT, 'aS') and membership(chS, chT, 'aT')
                and membership(chS, chT, 'ST')):
            passing.append((Sset, Tset))
print(f"pairs with dim V = 4: {len(dim4_pairs)}; windowwise-flat pairs: {len(passing)}")

def vclass_key(Sset, Tset):
    chS = [Fr(1) if i in Sset else Fr(0) for i in range(m)]
    chT = [Fr(1) if i in Tset else Fr(0) for i in range(m)]
    M = rref_rows([[Fr(1)] * m, chS, chT])
    return tuple(tuple(r) for r in M)

classes = {}
for (Sset, Tset) in passing:
    classes.setdefault(vclass_key(Sset, Tset), []).append((Sset, Tset))
print(f"distinct V-classes among passing pairs: {len(classes)}")

within = []
for lv in LEVELS:
    for r in range(1, len(lv) + 1):
        for comb in itertools.combinations(lv, r):
            within.append(tuple(sorted(comb)))
ck("within-level subsets number 8 = sum over levels of (2^mult - 1)", len(within) == 8)
pred = set()
for i1 in range(len(within)):
    for i2 in range(i1 + 1, len(within)):
        A_, B_ = within[i1], within[i2]
        if set(A_) & set(B_): continue
        if len(A_) + len(B_) == m: continue
        pred.add(vclass_key(A_, B_))
ck("predicted class set has size 26 = C(8,2) - 2 (golden overlaps excluded)", len(pred) == 26)
ck("structural law: passing V-classes == { <1,a,1_S',1_T'> : S',T' disjoint within-level }",
   set(classes.keys()) == pred)
ck("the k=3 windowwise-flat count is exactly 26 V-classes", len(classes) == 26)
ck("k=2 embedding: the pair ({phi},{tau}) is windowwise-flat",
   any({tuple(sorted(p[0])), tuple(sorted(p[1]))} == {(3,), (4,)} for p in passing))
ck("control: ({r2,r3},{r5}) spans two levels in one block and is NOT windowwise-flat",
   ((0, 1), (2,)) not in passing and ((2,), (0, 1)) not in passing)

# ---------------- lane-1 cross-validation on samples (sympy symbolic ranks) ----------------
import sympy as sp
L2, L3, L5, Lp = sp.symbols("L2 L3 L5 Lp")
SYM = [L2, L3, L5, Lp]
def to_expr(u): return u[0] + sum(Fr(u[k + 1]) * SYM[k] for k in range(4))
acat_expr = [to_expr(a) for a in AVAL]
def in_K_span_sym(chS, chT, target):
    cols = [[sp.Integer(1)] * m, acat_expr,
            [sp.Integer(x) for x in chS], [sp.Integer(x) for x in chT]]
    if target == 'aS':   h = [acat_expr[i] * chS[i] for i in range(m)]
    elif target == 'aT': h = [acat_expr[i] * chT[i] for i in range(m)]
    else:                h = [sp.Integer(chS[i] * chT[i]) for i in range(m)]
    Mfull = sp.Matrix(m, 4, lambda i, j: cols[j][i])
    Maug = Mfull.row_join(sp.Matrix(m, 1, lambda i, j: h[i]))
    zf = lambda x: sp.cancel(sp.expand(x)) == 0
    return Maug.rank(iszerofunc=zf) == Mfull.rank(iszerofunc=zf)

random.seed(11)
fails = [p for p in dim4_pairs if p not in passing]
sample = random.sample(passing, 12) + random.sample(fails, 12)
agree = True
for (Sset, Tset) in sample:
    chS = [1 if i in Sset else 0 for i in range(m)]
    chT = [1 if i in Tset else 0 for i in range(m)]
    for tg in ('aS', 'aT', 'ST'):
        if membership(chS, chT, tg) != in_K_span_sym(chS, chT, tg):
            agree = False
ck("two-lane agreement: fast membership == symbolic-rank membership on 24 sampled pairs x 3 targets",
   agree)

print(f"\nALL {NCK[0]} CHECKS PASSED (t4b_census_fast)")
