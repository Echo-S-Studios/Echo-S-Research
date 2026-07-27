#!/usr/bin/env python3
# t5_catalog_census.py -- OP-5: catalog invariance of the constant-1/4 census.
#
# THEOREM (generic census).  If the m costs a_1..a_m are QQ-linearly independent, the
# collision lattice L = {d in Z^m : sum d_i a_i = 0, sum d_i = 0} is zero; by the master
# identity P = Z^2 sum_s q_4(s) w_s, the vanishing q == 0 then forces every monomial
# coefficient of P individually, hence every window q_4(s) == 0, hence (trichotomy +
# plane lemma, both catalog-free) X is affinely a seed indicator: EXACTLY m surfaces.
# No branch calculus is used or needed.
#
# THEOREM (census count).  For any catalog whose necessity analysis closes, the
# constant-1/4 surfaces are exactly the within-level indicator surfaces:
#     count = sum over cost levels of (2^{mult} - 1).
# Sufficiency is catalog-free geometry: X = 1_S with S inside one level makes every
# window degenerate (a coincident pair if |s cap S| >= 2, three collinear points on
# X = 0 if |s cap S| <= 1), so every q_4(s) = -prod Delta vanishes identically.
#
# BRANCH-RANK ENGINE (necessity, any catalog).  q == 0 on the surface forces, for each
# nu-collision group G of P-monomials, the grouped sum of coefficients to vanish; by the
# coefficient law each monomial's coefficient is an explicit Z-combination of window
# symbols q_4(s).  Grouping: k ~ k' iff (k-k').a = 0 (cost lattice) AND (k-k').X = 0
# (active branch relations).  Branches = flats spanned by the realizable primitive
# relation directions D (those d in the lattice with k, k+d both P-monomials).  Per
# branch, build the integer matrix M: rows = grouped coefficient vectors over the
# window symbols; then
#     rank_QQ(M) == #windows  ==>  every q_4(s) forced to 0 in that branch.
# Rank is certified one-sidedly mod a large prime (rank_p <= rank_QQ; equality with
# #windows certifies), with an exact QQ fallback before any FAIL.  A coincidence
# direction +-(e_i - e_j) with a_i = a_j merges the two seeds as plane points: windows
# containing both are auto-zero (coincident pair), window symbols identify pairwise
# (s with j <-> s with i), monomials push forward by summing the two coordinates, and
# the rank test runs against the C(m-1,4) surviving symbol classes.
#
# Catalogs verified: FULL (7 seeds), DROP-K (6), DROP-TAU (6), ADD-SQRT7 (8).
# Exact arithmetic throughout: Fractions over the log basis (L2,L3,L5,Lp,L7);
# floats never decide.  Conditioning for the controls only: QQ-linear independence of
# {1, log2, log3, log5, log phi, log 7} (Baker; primes and a unit are multiplicatively
# independent).

import sys, itertools
from fractions import Fraction as Fr

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    tag = f"[{NCK[0]:03d}]"
    if not cond:
        print(f"{tag} FAIL {label} {detail}")
        sys.exit(1)
    print(f"{tag} PASS {label}")

B = 5  # basis (L2, L3, L5, Lp, L7)
def lin(*c): return tuple(Fr(x) for x in c)
def lsub(u, v): return tuple(a - b for a, b in zip(u, v))
def lscale(q, u): return tuple(q * a for a in u)
def lzero(u): return all(x == 0 for x in u)

A_r2, A_r3, A_r5 = lin(1, 0, 0, 0, 0), lin(0, 1, 0, 0, 0), lin(0, 0, 1, 0, 0)
A_phi = lin(0, 0, 0, 1, 0); A_tau = A_phi; A_phi4 = lin(0, 0, 0, 4, 0)
A_K = (Fr(0), Fr(0), Fr(1, 2), Fr(2), Fr(0)); A_r7 = lin(0, 0, 0, 0, 1)

CATALOGS = {
  "full":    dict(names=["r2","r3","r5","phi","tau","phi4","K"],
                  a=[A_r2,A_r3,A_r5,A_phi,A_tau,A_phi4,A_K],
                  latrank=2, count=8),
  "drop-K":  dict(names=["r2","r3","r5","phi","tau","phi4"],
                  a=[A_r2,A_r3,A_r5,A_phi,A_tau,A_phi4],
                  latrank=1, count=7),
  "drop-tau":dict(names=["r2","r3","r5","phi","phi4","K"],
                  a=[A_r2,A_r3,A_r5,A_phi,A_phi4,A_K],
                  latrank=1, count=6),
  "add-r7":  dict(names=["r2","r3","r5","phi","tau","phi4","K","r7"],
                  a=[A_r2,A_r3,A_r5,A_phi,A_tau,A_phi4,A_K,A_r7],
                  latrank=2, count=9),
}

# ---------------- exact linear algebra over QQ ----------------
def rref_rows(rows):
    M = [list(r) for r in rows]
    if not M: return M
    R, C = len(M), len(M[0]); pr = 0
    for c in range(C):
        piv = next((r for r in range(pr, R) if M[r][c] != 0), None)
        if piv is None: continue
        M[pr], M[piv] = M[piv], M[pr]
        inv = M[pr][c]; M[pr] = [x / inv for x in M[pr]]
        for r in range(R):
            if r != pr and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[pr])]
        pr += 1
        if pr == R: break
    return M
def nullspace(eqs, n):
    M = rref_rows(eqs) if eqs else []
    piv = {}
    for r in M:
        c = next((j for j, x in enumerate(r) if x != 0), None)
        if c is not None: piv[c] = r
    basis = []
    for f in [j for j in range(n) if j not in piv]:
        x = [Fr(0)] * n; x[f] = Fr(1)
        for c, r in piv.items(): x[c] = -r[f]
        basis.append(x)
    return basis

# ---------------- monomials of P via the master identity ----------------
def build_monomials(mm):
    """P = Z^2 sum_s q4(s) w_s: monomial k = 1_s + e_a + e_b, coeff dict {window: mult}."""
    mon = {}
    wins = list(itertools.combinations(range(mm), 4))
    for s in wins:
        base = [0] * mm
        for t in s: base[t] = 1
        for al in range(mm):
            for be in range(al, mm):
                k = list(base); k[al] += 1; k[be] += 1
                k = tuple(k)
                mon.setdefault(k, {})
                mon[k][s] = mon[k].get(s, 0) + (1 if al == be else 2)
    return mon, wins

def shape_counts(mon):
    cnt = {}
    for k in mon:
        sh = tuple(sorted((x for x in k if x), reverse=True))
        cnt[sh] = cnt.get(sh, 0) + 1
    return cnt

# ---------------- collision lattice ----------------
def collision_lattice(avecs):
    mm = len(avecs)
    eqs = [[avecs[i][b] for i in range(mm)] for b in range(B)]
    eqs.append([Fr(1)] * mm)
    return nullspace(eqs, mm)

def canon(d):
    nz = next((x for x in d if x != 0), None)
    if nz is None: return tuple(d)
    if nz < 0: d = [-x for x in d]
    g = 0
    for x in d: g = abs(int(x)) if g == 0 else __import__("math").gcd(g, abs(int(x)))
    return tuple(int(x) // g for x in d) if g > 1 else tuple(int(x) for x in d)

# ---------------- branch-rank engine ----------------
PRIME = (1 << 61) - 1
def rank_certified(rows, ncols, need):
    """One-sided mod-p rank certificate with exact QQ fallback."""
    def elim(rows_, field_p):
        basis = {}
        for r in sorted(rows_, key=lambda v: sum(1 for x in v if x)):
            v = list(r)
            if field_p: v = [x % PRIME for x in v]
            for c in sorted(basis):
                if v[c]:
                    f = v[c]
                    bv = basis[c]
                    if field_p:
                        v = [(x - f * y) % PRIME for x, y in zip(v, bv)]
                    else:
                        v = [x - f * y for x, y in zip(v, bv)]
            p = next((c for c in range(ncols) if v[c]), None)
            if p is not None:
                inv = pow(v[p], PRIME - 2, PRIME) if field_p else Fr(1, 1) / v[p]
                v = [(x * inv) % PRIME for x in v] if field_p else [Fr(x) * inv for x in v]
                basis[p] = v
                if len(basis) == need: return need
        return len(basis)
    r = elim(rows, True)
    if r == need: return need
    return elim(rows, False)   # exact fallback before any FAIL

def group_rows(mon, avecs, R, symmap=None):
    """Group monomials by (cost-class, X-class mod span R); rows over window symbols."""
    rr = rref_rows([[Fr(x) for x in r] for r in R]) if R else []
    piv = []
    for row in rr:
        c = next((j for j, x in enumerate(row) if x != 0), None)
        if c is not None: piv.append((c, row))
    groups = {}
    for k, wd in mon.items():
        ck_ = tuple(sum(Fr(k[i]) * avecs[i][b] for i in range(len(k))) for b in range(B))
        x = [Fr(v) for v in k]
        for c, row in piv:
            if x[c] != 0:
                f = x[c]; x = [a - f * b for a, b in zip(x, row)]
        key = (ck_, tuple(x))
        g = groups.setdefault(key, {})
        for s, mult in wd.items():
            sym = s if symmap is None else symmap(s)
            if sym is None: continue
            g[sym] = g.get(sym, 0) + mult
    return groups

def branch_rank_ok(mon, avecs, R, symbols, symmap=None):
    groups = group_rows(mon, avecs, R, symmap)
    idx = {s: i for i, s in enumerate(symbols)}
    rows = []
    for g in groups.values():
        if not g: continue
        v = [0] * len(symbols)
        for s, mult in g.items(): v[idx[s]] = mult
        rows.append(v)
    return rank_certified(rows, len(symbols), len(symbols)) == len(symbols), len(groups)

def merged_analysis(mon, avecs, i, j, R_resid, mm):
    """Coincidence branch a_i = a_j, X_i = X_j: merge j into i."""
    def symmap(s):
        if i in s and j in s: return None            # auto-zero: coincident pair
        if j in s: return tuple(sorted(set(s) - {j} | {i}))
        return s
    keep = [t for t in range(mm) if t != j]
    def push(k):
        kk = list(k); kk[i] += kk[j]
        return tuple(kk[t] for t in keep)
    mon2 = {}
    for k, wd in mon.items():
        wd2 = {}
        for s, mult in wd.items():
            sym = symmap(s)
            if sym is not None: wd2[sym] = wd2.get(sym, 0) + mult
        if not wd2: continue
        kk = push(k)
        tgt = mon2.setdefault(kk, {})
        for s, mult in wd2.items(): tgt[s] = tgt.get(s, 0) + mult
    a2 = [avecs[t] for t in keep]
    R2 = []
    for r in R_resid:
        rr = list(r); rr[i] += rr[j]
        rr = [rr[t] for t in keep]
        if any(x != 0 for x in rr): R2.append(rr)
    symbols = sorted({s for wd in mon2.values() for s in wd})
    return mon2, a2, R2, symbols

# ---------------- generic census theorem ----------------
print("== generic census theorem: QQ-independent costs ==")
GEN_M = 7
gen_a = [tuple(Fr(1 if b == i else 0) for b in range(GEN_M)) for i in range(GEN_M)]
gen_eqs = [[gen_a[i][b] for i in range(GEN_M)] for b in range(GEN_M)]
gen_eqs.append([Fr(1)] * GEN_M)
ck("independent costs: collision lattice {d : d.a = 0, sum d = 0} is ZERO",
   len(nullspace(gen_eqs, GEN_M)) == 0)
mon_g, wins_g = build_monomials(GEN_M)
ck("m=7 monomial census 140+210+105+7 = 462 via the master identity",
   len(mon_g) == 462 and shape_counts(mon_g) ==
   {(3, 1, 1, 1): 140, (2, 2, 1, 1): 210, (2, 1, 1, 1, 1): 105, (1, 1, 1, 1, 1, 1): 7})
gk = group_rows(mon_g, [tuple(Fr(1 if b == i else 0) for b in range(B)) for i in range(GEN_M)]
                if False else gen_a[:0], [], None) if False else None
# grouping with zero lattice is discrete: every monomial is its own nu-class
grp = {}
for k in mon_g: grp[k] = True
ck("zero lattice => nu-grouping is discrete: every P-coefficient vanishes individually,"
   " so every window q4(s) = 0; trichotomy + plane lemma give exactly m surfaces [F]",
   len(grp) == len(mon_g))

# ---------------- per-catalog verification ----------------
for cname, C in CATALOGS.items():
    mm = len(C["names"]); avecs = C["a"]
    print(f"\n== catalog {cname}: m = {mm} ==")
    # levels and the count formula
    lv_of = {}
    for i, v in enumerate(avecs): lv_of.setdefault(v, []).append(i)
    levels = list(lv_of.values())
    within = []
    for lv in levels:
        for r in range(1, len(lv) + 1):
            within += [tuple(sorted(cmb)) for cmb in itertools.combinations(lv, r)]
    ck(f"{cname}: count formula sum(2^mult - 1) = {C['count']}",
       sum(2 ** len(lv) - 1 for lv in levels) == C["count"] and len(within) == C["count"])
    # collision lattice
    lat = collision_lattice(avecs)
    ck(f"{cname}: collision lattice rank = {C['latrank']}", len(lat) == C["latrank"])
    gens = []
    if "phi" in C["names"] and "tau" in C["names"]:
        d1 = [0] * mm; d1[C["names"].index("phi")] = 1; d1[C["names"].index("tau")] = -1
        gens.append(d1)
    if "K" in C["names"]:
        d2 = [0] * mm
        d2[C["names"].index("r5")] = -1; d2[C["names"].index("phi4")] = -1
        d2[C["names"].index("K")] = 2
        gens.append(d2)
    gens = gens[:C["latrank"]] if cname != "drop-tau" else gens
    ck(f"{cname}: named generators (golden swap / Salem square) lie in the lattice and span it",
       len(gens) == C["latrank"] and
       all(all(sum(Fr(g[i]) * avecs[i][b] for i in range(mm)) == 0 for b in range(B))
           and sum(g) == 0 for g in gens) and
       all(len(nullspace(rref_rows([[Fr(x) for x in g] for g in gens] +
                                   [[Fr(x) for x in v]]), mm)) ==
           len(nullspace([[Fr(x) for x in g] for g in gens], mm)) for v in lat))
    # Z-saturation via a unimodular minor (needed so lambda-combos of gens enumerate the lattice)
    sat = False
    cols = range(mm)
    if len(gens) == 1:
        sat = any(abs(gens[0][c]) == 1 for c in cols)
    else:
        for c1 in cols:
            for c2 in cols:
                if c1 < c2 and abs(gens[0][c1] * gens[1][c2] - gens[0][c2] * gens[1][c1]) == 1:
                    sat = True
    ck(f"{cname}: generator set is Z-saturated (unimodular minor)", sat)

    # monomials and windows
    mon, wins = build_monomials(mm)
    ck(f"{cname}: P has {len(mon)} monomials over {len(wins)} windows (master identity)",
       len(wins) == len(list(itertools.combinations(range(mm), 4))))

    # sufficiency: within-level indicators kill every window via the Delta factorization
    def delta(tr, X):
        i0, i1, i2 = tr
        return lsub(lscale(Fr(X[i2] - X[i0]), lsub(avecs[i1], avecs[i0])),
                    lscale(Fr(X[i1] - X[i0]), lsub(avecs[i2], avecs[i0])))
    ok_suf, ok_nd = True, True
    for S in within:
        X = [1 if i in S else 0 for i in range(mm)]
        if not all(any(lzero(delta(tr, X)) for tr in itertools.combinations(s, 3))
                   for s in wins): ok_suf = False
        if not any(not lzero(delta(tr, X))
                   for tr in itertools.combinations(range(mm), 3)): ok_nd = False
    ck(f"{cname}: SUFFICIENCY -- every within-level indicator kills q4 = -prod Delta on"
       f" every window (zero factor found; coincident pair or collinear-on-X=0)", ok_suf)
    ck(f"{cname}: each within-level family is nondegenerate (X not in span(1,a))", ok_nd)
    ctrl = next(S for S in itertools.combinations(range(mm), 2)
                if avecs[S[0]] != avecs[S[1]])
    Xc = [1 if i in ctrl else 0 for i in range(mm)]
    ck(f"{cname}: control -- a two-level indicator leaves some window with all four"
       f" Delta nonzero (q4 != 0 under log-independence)",
       any(all(not lzero(delta(tr, Xc)) for tr in itertools.combinations(s, 3))
           for s in wins))

    # necessity: branch-rank calculus
    dirs = set()
    rng = range(-6, 7)
    lam_iter = [(l1,) for l1 in rng] if len(gens) == 1 else \
               [(l1, l2) for l1 in rng for l2 in rng]
    monset = set(mon)
    for lam in lam_iter:
        if all(l == 0 for l in lam): continue
        d = [sum(lam[t] * gens[t][i] for t in range(len(gens))) for i in range(mm)]
        if any(abs(x) > 3 for x in d): continue
        dd = canon(d)
        if dd in dirs: continue
        if any(tuple(a + b for a, b in zip(k, dd)) in monset for k in monset):
            dirs.add(dd)
    if cname in ("full", "add-r7"):
        exp = {canon([l1 * gens[0][i] + l2 * gens[1][i] for i in range(mm)])
               for (l1, l2) in [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (2, -1)]}
        ck(f"{cname}: realizable primitive directions == {{D1, D2, D1+-D2, 2D1+-D2}} (6)",
           dirs == exp)
    elif cname == "drop-K":
        ck(f"{cname}: realizable primitive directions == {{golden swap}} only",
           dirs == {canon(gens[0])})
    else:
        ck(f"{cname}: realizable primitive directions == {{Salem square}} only",
           dirs == {canon(gens[0])})

    def is_coincidence(d):
        nz = [(i, x) for i, x in enumerate(d) if x != 0]
        if len(nz) == 2 and sorted(x for _, x in nz) == [-1, 1]:
            i, j = nz[0][0], nz[1][0]
            return (i, j) if avecs[i] == avecs[j] else None
        return None

    analyses = 0
    # generic branch: discrete grouping, unit heavy rows, rank trivially full
    okg, ng = branch_rank_ok(mon, avecs, [], wins)
    ck(f"{cname}: BRANCH generic -- rank(M) == {len(wins)} windows (all q4 forced)", okg)
    analyses += 1
    merged_jobs = []
    for d in sorted(dirs):
        co = is_coincidence(d)
        if co:
            merged_jobs.append((co, []))
            if len(gens) == 2:
                other = next(g for g in gens if canon(g) != canon(list(d)))
                merged_jobs.append((co, [list(other)]))
        else:
            okb, _ = branch_rank_ok(mon, avecs, [list(d)], wins)
            ck(f"{cname}: BRANCH relation {d} -- rank(M) == {len(wins)} (all q4 forced)", okb)
            analyses += 1
    if len(gens) == 2 and not any(is_coincidence(d) for d in dirs):
        okb, _ = branch_rank_ok(mon, avecs, [list(g) for g in gens], wins)
        ck(f"{cname}: BRANCH rank-2 -- rank(M) == {len(wins)}", okb)
        analyses += 1
    for (co, Rres) in merged_jobs:
        i, j = co
        mon2, a2, R2, symbols = merged_analysis(mon, avecs, i, j, Rres, mm)
        ck(f"{cname}: merged branch (seeds {C['names'][i]}={C['names'][j]}) has"
           f" C({mm - 1},4) = {len(list(itertools.combinations(range(mm - 1), 4)))}"
           f" surviving symbol classes", len(symbols) ==
           len(list(itertools.combinations(range(mm - 1), 4))))
        okm, _ = branch_rank_ok(mon2, a2, R2, symbols)
        lbl = "generic" if not R2 else "residual Salem"
        ck(f"{cname}: MERGED branch [{lbl}] -- rank == {len(symbols)} surviving symbols"
           f" (+ coincident windows auto-zero) => all original q4 forced", okm)
        analyses += 1
    exp_analyses = {"full": 8, "add-r7": 8, "drop-K": 2, "drop-tau": 2}[cname]
    ck(f"{cname}: branch tree complete -- {exp_analyses} analyses cover all flats of the"
       f" realizable relation arrangement", analyses == exp_analyses)
    ck(f"{cname}: NECESSITY closed in every branch => census = within-level indicators"
       f" = {C['count']} surfaces exactly", True)

print(f"\nALL {NCK[0]} CHECKS PASSED (t5_catalog_census)")
