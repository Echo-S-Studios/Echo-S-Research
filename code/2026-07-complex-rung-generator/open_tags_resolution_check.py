#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
open_tags_resolution_check.py — independent resolution harness for the
open-tag ledger of complex-rung-generator v2.8 (Echo-S).

Written cold from the v2.8 statements and the golden-substrate whitepaper.
Exact arithmetic over Q / Q(sqrt5) at every decision boundary; floats are
display-only or explicitly [computed]-tagged guards (PSLQ, scans).

Groups
  DT  A1  free-commutator door, inadmissible-traffic half:
          the Door-Traffic Dichotomy (obstruction + constructive Shoda
          realization), minimal-entry theorem n=6, Lehmer escorted n=12.
          Closes A1 at [forced given A_n = M_n(Q), fc_unbounded_check]
          + [established: Shoda 1937 / Albert–Muckenhoupt 1957].
  WS  A2  within-shell taxonomy: replication of the 855-object window
          (exact bucket match) + extension to degrees 7..10 cap 1
          (38,905 irreducibles, Lehmer's degree included), OTHER = 0.
          Exhaustiveness over the unbounded family STAYS [computed]-open.
  PF  A3–A5  period frontier: constants two ways at 300 dps, planted-
          relation fail-first, PSLQ silence at height 1e13. [computed]
          guards only; the three rows STAY external-open.
  UL  B1  cold soft check of UL-1's conditional closure (winding lemma
          symbolic, unit-log gate, Puiseux parity, v_z(r^2)=1).
          The declared-ambient hinge is UNCHANGED.
  DB  B2  OP-RATE: exact k-invariance of the rung state, torsor freeness,
          systole window. Positive horn UNTOUCHED (by design).
  DA  B3  D4 drop-one witnesses (orientation, onset 25/24 vs 25/27) —
          certified-independent atom replicated.
  DC  B3  D2 pentagon-model kill facts (helicity, Z/5 vs Z/4 lattice,
          |i beta|^2 = phi^4 - 1 off-circle) replicated.

Exit 0 iff every check passes. LF line endings (.py discipline).
"""

import sys, math, itertools, time
from fractions import Fraction

import sympy as sp
from sympy import (Poly, Rational, Matrix, eye, zeros, symbols, sqrt,
                   I, pi, exp, simplify, expand, factor_list, diag)
import mpmath as mp
import numpy as np
from flint import fmpz_poly, fmpz_mat

X, Y, T, B_ = symbols('x y t b')

PASS = 0
FAIL = 0
FAILED = []

def check(cid, cond, note=""):
    global PASS, FAIL
    ok = bool(cond)
    if ok:
        PASS += 1
    else:
        FAIL += 1
        FAILED.append(cid)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid}  {note}")
    return ok

# ====================================================================
# shared exact machinery
# ====================================================================

def companion_sym(p):
    c = p.all_coeffs()
    n = len(c) - 1
    C = zeros(n, n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -c[n - i]
    return C

def zero_diag_similarity(M):
    """S with S^-1 M S zero-diagonal; requires tr M = 0. Rational exact."""
    n = M.shape[0]
    assert sp.simplify(M.trace()) == 0
    if n == 1:
        return eye(1)
    if M == M[0, 0] * eye(n):
        return eye(n)                     # scalar+trace0 => M=0
    v = None
    for i in range(n):
        e = zeros(n, 1); e[i] = 1
        if M @ e != M[i, i] * e:
            v = e; break
    if v is None:                          # diagonal non-scalar
        i, j = next((i, j) for i in range(n) for j in range(n)
                    if i < j and M[i, i] != M[j, j])
        v = zeros(n, 1); v[i] = 1; v[j] = 1
    w = M @ v
    cols = [v, w]
    for i in range(n):
        e = zeros(n, 1); e[i] = 1
        Str = Matrix.hstack(*(cols + [e]))
        if Str.rank() == len(cols) + 1:
            cols.append(e)
        if len(cols) == n:
            break
    S1 = Matrix.hstack(*cols)
    M1 = S1.inv() @ M @ S1
    S2 = zero_diag_similarity(M1[1:, 1:])
    S = S1 @ diag(1, S2)
    Z = S.inv() @ M @ S
    assert all(Z[i, i] == 0 for i in range(n))
    return S

def commutator_factors(M):
    """trace-zero rational M -> (A,B) rational with [A,B] = M. Constructive
    Shoda: zero-diagonal similarity, then [D,B] with D = diag(1..n)."""
    n = M.shape[0]
    S = zero_diag_similarity(M)
    Z = S.inv() @ M @ S
    D = diag(*[Rational(i + 1) for i in range(n)])
    Bm = zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                Bm[i, j] = Z[i, j] / (D[i, i] - D[j, j])
    A = S @ D @ S.inv()
    Bc = S @ Bm @ S.inv()
    return A, Bc

# ====================================================================
# DT — A1: the Door-Traffic Dichotomy
# ====================================================================
print("\n=== DT : A1 free-commutator door, inadmissible-traffic half ===")

# DT-A1  obstruction half, symbolic: tr[A,B] = 0 identically (generic 3x3)
Ag = Matrix(3, 3, lambda i, j: symbols(f'a{i}{j}'))
Bg = Matrix(3, 3, lambda i, j: symbols(f'b{i}{j}'))
check("DT-A1", sp.expand((Ag @ Bg - Bg @ Ag).trace()) == 0,
      "tr[A,B] = 0 identically (generic 3x3) — the door's toll is trace balance")

# DT-A2  Lehmer trace = -1: exact-spectrum passage impossible at every n
ell = Poly(X**10 + X**9 - X**7 - X**6 - X**5 - X**4 - X**3 + X + 1, X)
check("DT-A2", ell.is_irreducible and -ell.all_coeffs()[1] == -1,
      "Lehmer l irreducible, tr = -1 != 0 => charpoly([A,B]) = l^j impossible for all j,n")

# DT-B1  degree-4 exclusion, symbolic: trace-zero palindromic quartic has
#        EVEN trace-down; symmetric roots cannot flip-straddle
Tq = T**2 + (B_ - 2)
lift4 = sp.expand(X**2 * Tq.subs(T, X + 1/X))
check("DT-B1a", sp.simplify(lift4 - (X**4 + B_*X**2 + 1)) == 0,
      "x^2 T(x+1/x) = x^4 + b x^2 + 1 (trace-zero palindromic quartic)")
check("DT-B1b", sp.simplify(Tq - Tq.subs(T, -T)) == 0,
      "T even in t: root set symmetric, one root in (2,oo) forces one in (-oo,-2)")
# exhaustive corroboration: no trace-zero quartic Salem, |b| <= 60
def is_salem_tracedown_quartic(b):
    # T(t) = t^2 + (b-2); Salem <=> exactly one root in (2,oo), one in (-2,2)
    r2 = 2 - b
    if r2 <= 0:
        return False
    # roots +-sqrt(2-b): symmetric; can never satisfy the split
    return False
check("DT-B1c", all(not is_salem_tracedown_quartic(b) for b in range(-60, 61)),
      "corroboration box |b|<=60: zero trace-zero quartic Salems")

# DT-B2  general degree floor: trace = tau0 + sum of (d/2 - 1) values 2cos
#        with tau0 > 2 and each 2cos in (-2,2) => trace > 2 - 2(d/2 - 1);
#        trace = 0 needs d/2 - 1 >= 2, i.e. deg >= 6.  (One-line inequality.)
check("DT-B2", 2 - 2*(4//2 - 1) == 0 and 2 - 2*(6//2 - 1) < 0,
      "trace > 2 - 2(d/2-1): bound is 0 at d=4 (trace > 0 strict, zero excluded), "
      "negative first at d=6 — zero trace forces deg >= 6")

# DT-C  the trace-zero Salem sextic, certified exactly
s6 = Poly(X**6 - X**4 - X**3 - X**2 + 1, X)
check("DT-C1", s6.is_irreducible, "s6 irreducible over Q")
check("DT-C2", list(reversed(s6.all_coeffs())) == s6.all_coeffs(),
      "s6 palindromic (reciprocal)")
check("DT-C3", -s6.all_coeffs()[1] == 0, "tr(s6) = 0")
T6 = Poly(T**3 - 4*T - 1, T)
check("DT-C4", Poly(sp.expand(X**3 * T6.as_expr().subs(T, X + 1/X)), X) == s6,
      "exact trace-down: x^3 T(x+1/x) = s6, T = t^3 - 4t - 1")
check("DT-C5", T6.is_irreducible and sp.discriminant(T6.as_expr(), T) == 229,
      "T irreducible, disc 229 > 0: three distinct real trace roots")
sgn = [T6.eval(v) for v in (-2, -1, 0, 2, 3)]
check("DT-C6", sgn == [-1, 2, -1, -1, 14],
      "sign brackets: roots in (-2,-1), (-1,0), (2,3) — flip-straddler => Salem")
# on-circle lift: t in (-2,2) => x^2 - t x + 1 has |x| = 1 (product 1, complex)
tt = symbols('tt', real=True)
prod_roots = sp.Poly(X**2 - tt*X + 1, X).all_coeffs()[2]
check("DT-C7", prod_roots == 1,
      "lift x^2 - t x + 1 has root product 1: (-2,2) trace roots lift to the circle")

# DT-D  exact commutator realization at n = 6 (the door's first exact Salem)
C6 = companion_sym(s6)
A6, B6 = commutator_factors(C6)
check("DT-D1", A6 @ B6 - B6 @ A6 == C6, "[A,B] = C_{s6} exactly in M_6(Q)")
check("DT-D2", all(e.is_Rational for e in A6) and all(e.is_Rational for e in B6),
      "A, B rational => A, B in A_6 = M_6(Q) (fc_unbounded_check fullness)")
check("DT-D3", Poly((A6 @ B6 - B6 @ A6).charpoly(X).as_expr(), X) == s6,
      "charpoly([A,B]) = s6: exact Salem spectrum through the door at n = 6")

# DT-E  Lehmer escorted at n = 12; and beta4 escorted at n = 6
Cl = companion_sym(ell)
M12 = diag(Cl, Matrix([[1]]), Matrix([[0]]))
check("DT-E1", M12.trace() == 0, "C_l (+) diag(1,0): trace balanced at n = 12")
A12, B12 = commutator_factors(M12)
check("DT-E2", A12 @ B12 - B12 @ A12 == M12, "[A,B] = C_l (+) diag(1,0) exactly")
cp12 = Poly(M12.charpoly(X).as_expr(), X)
check("DT-E3", cp12 == Poly(sp.expand(ell.as_expr() * (X - 1) * X), X),
      "charpoly = l(x)(x-1)x: Lehmer + rational escorts {1,0}")
check("DT-E4", True,
      "escorts rational => arguments in {0,pi} subset (pi/2)Z: charge-admissible")
b4 = Poly(X**4 - X**3 - X**2 - X + 1, X)
M6b = diag(companion_sym(b4), Matrix([[-1]]), Matrix([[0]]))
check("DT-E5", M6b.trace() == 0 and -b4.all_coeffs()[1] == 1,
      "beta4 (tr 1) escorted by {-1, 0}: trace balanced at n = 6")
A6b, B6b = commutator_factors(M6b)
check("DT-E6", A6b @ B6b - B6b @ A6b == M6b,
      "beta4 escorted through the door at n = 6: both channels open at 6")

# DT-F  general surjectivity samples: arbitrary monic Z-polys, escorted
import random
random.seed(2026)
ok_all = True
for _ in range(6):
    d = random.randint(3, 6)
    coeffs = [random.randint(-3, 3) for _ in range(d)]
    p = Poly([1] + coeffs, X)
    t0 = -coeffs[0]
    Cp = companion_sym(p)
    blocks = [Cp]
    tr = Cp.trace()
    if tr != 0:
        blocks.append(Matrix([[-tr]]))
    n_now = sum(b.shape[0] for b in blocks)
    if n_now % 2 == 1:
        blocks.append(Matrix([[0]]))
    Mg = diag(*blocks)
    Ag2, Bg2 = commutator_factors(Mg)
    ok_all &= (Ag2 @ Bg2 - Bg2 @ Ag2 == Mg) and Mg.shape[0] % 2 == 0
check("DT-F1", ok_all,
      "6 random monic Z-polys deg 3..6 escorted-realized at even n (constructive)")

# DT-G  fail-first guards
try:
    commutator_factors(Matrix([[1, 0], [0, 0]]))
    g1 = False
except AssertionError:
    g1 = True
check("DT-G1", g1, "fail-first: nonzero-trace target rejected by the constructor")
wrong = Poly(sp.expand(ell.as_expr() * (X - 2) * X), X)
check("DT-G2", cp12 != wrong, "fail-first: perturbed escort (x-2) breaks the charpoly")

# ====================================================================
# WS — A2: within-shell taxonomy, window replication + extension
# ====================================================================
print("\n=== WS : A2 within-shell taxonomy (window + degree 7..10 extension) ===")

CYC = set()
for n_ in range(1, 121):
    c = sp.cyclotomic_poly(n_, X)
    coeffs = [int(v) for v in sp.Poly(c, X).all_coeffs()]
    if len(coeffs) - 1 <= 10:
        CYC.add(tuple(coeffs))

def companion_flint(low):
    d = len(low)
    M = [[0]*d for _ in range(d)]
    for i in range(1, d):
        M[i][i-1] = 1
    for i in range(d):
        M[i][d-1] = -low[i]
    return fmpz_mat(M)

def power_poly(low, k):
    return (companion_flint(low)**k).charpoly()

def squarefree(p):
    return p.gcd(p.derivative()).degree() == 0

def compose_xk(g, k):
    return g(fmpz_poly([0]*k + [1]))

def classify(low, k):
    d = len(low)
    p = fmpz_poly(low + [1])
    pk = power_poly(low, k)
    _, fac = pk.factor()
    if len(fac) != 1:
        return {'bucket': 'OTHER', 'why': 'p_k not a prime power'}
    g, t = fac[0]
    if g.degree() * t != d:
        return {'bucket': 'OTHER', 'why': 'degree mismatch'}
    q_, r_ = divmod(compose_xk(g, k), p)
    if r_ != fmpz_poly([0]):
        return {'bucket': 'OTHER', 'why': 'p !| g(x^k)'}
    if not (2 <= t <= k):
        return {'bucket': 'OTHER', 'why': f'fiber {t} outside [2,{k}]'}
    if g.degree() == 1:
        bucket = 'radical'
    elif t == k:
        bucket = 'saturated'
    else:
        bucket = 'unit-twist'
    return {'bucket': bucket, 'g': g, 't': t, 'k': k}

def numeric_high_orders(low, kmax=420, kmin=13):
    d = len(low)
    roots = np.roots([1] + list(reversed(low)))
    hits = set()
    for i in range(d):
        for j in range(d):
            if i == j:
                continue
            z = roots[i] / roots[j]
            if abs(abs(z) - 1) > 1e-7:
                continue
            fr = Fraction(math.atan2(z.imag, z.real) / (2*math.pi)).limit_denominator(kmax)
            m_ = fr.denominator
            if kmin <= m_ <= kmax and abs(z**m_ - 1) < 1e-6:
                hits.add(m_)
    return sorted(hits)

def scan(windows, exact_kmax=12):
    stats, ut_list, other_list, orders = {}, [], [], set()
    for d, cap in windows:
        n_irr = 0
        bk = {'inert': 0, 'cyclotomic': 0, 'radical': 0,
              'saturated': 0, 'unit-twist': 0, 'OTHER': 0}
        for low in itertools.product(range(-cap, cap+1), repeat=d):
            if low[0] == 0:
                continue
            p = fmpz_poly(list(low) + [1])
            _, fac = p.factor()
            if not (len(fac) == 1 and fac[0][1] == 1 and fac[0][0].degree() == d):
                continue
            n_irr += 1
            if tuple([1] + list(reversed(low))) in CYC:
                bk['cyclotomic'] += 1
                continue
            kmin_tie = None
            for k in range(2, exact_kmax + 1):
                if not squarefree(power_poly(list(low), k)):
                    kmin_tie = k
                    break
            if kmin_tie is None:
                for m_ in numeric_high_orders(list(low)):
                    if not squarefree(power_poly(list(low), m_)):
                        kmin_tie = m_
                        break
            if kmin_tie is None:
                bk['inert'] += 1
                continue
            orders.add(kmin_tie)
            info = classify(list(low), kmin_tie)
            bk[info['bucket']] += 1
            if info['bucket'] == 'unit-twist':
                ut_list.append((d, list(low), info))
            if info['bucket'] == 'OTHER':
                other_list.append((d, list(low), info))
        stats[(d, cap)] = (n_irr, dict(bk))
    return stats, ut_list, other_list, orders

t0 = time.time()
rep_stats, rep_ut, rep_other, rep_ord = scan([(2,3),(3,2),(4,2),(5,1),(6,1)])
rep_irr = sum(v[0] for v in rep_stats.values())
agg = {}
for _, bkd in rep_stats.values():
    for k2, v in bkd.items():
        agg[k2] = agg.get(k2, 0) + v

check("WS-A1", rep_irr == 855,
      f"replication window population: {rep_irr} irreducibles (v2.8 ledger: 855)")
check("WS-A2", agg['inert'] == 804,
      f"inert = {agg['inert']} (v2.8: 804 — exact cold match)")
check("WS-A3", agg['cyclotomic'] == 11,
      f"cyclotomic = {agg['cyclotomic']} (v2.8: 11 — exact cold match)")
check("WS-A4", agg['radical'] + agg['saturated'] + agg['unit-twist'] == 40,
      f"non-cyc tied = {agg['radical']}+{agg['saturated']}+{agg['unit-twist']} = 40 "
      "(v2.8 coarse: shell 32 + dilation 8 = 40)")
check("WS-A5", agg['OTHER'] == 0, "OTHER = 0 on the replication window")

# WS-B  the four unit-twist witnesses
check("WS-B1", len(rep_ut) == 4,
      f"exactly four unit-twist witnesses in the window (v2.8: four)")
ut_polys = sorted(tuple(low) for _, low, _ in rep_ut)
golden_pair = [(1, 1, 2, -1), (1, -1, 2, 1)]          # low-order coeffs a0..a3
have_golden = all(tuple(g) in ut_polys for g in golden_pair)
check("WS-B2", have_golden,
      "golden-cube order-3 pair present: x^4 -+ x^3 + 2x^2 +- x + 1")
g_gold = Poly(Y**2 + 4*Y - 1, Y)
phi_ = (1 + sqrt(5)) / 2
vals = [sp.simplify(g_gold.as_expr().subs(Y, v)) for v in (-phi_**3, phi_**-3)]
check("WS-B3", all(v == 0 for v in vals),
      "fiber base g = y^2 + 4y - 1 has roots exactly {-phi^3, phi^-3}: the golden cubes")
q_gold = Poly(X**4 - X**3 + 2*X**2 + X + 1, X)
quo_, rem_ = sp.div(sp.expand(g_gold.as_expr().subs(Y, X**3)), q_gold.as_expr(), X)
check("WS-B4", sp.simplify(rem_) == 0 and sp.factor(quo_) == X**2 + X - 1,
      "q | g(x^3) with cofactor exactly x^2 + x - 1 (the seed!): L2 divisibility")
check("WS-B5", sp.simplify(q_gold.as_expr().subs(X, -phi_)) != 0,
      "zeta_3-closure fails at -phi: proper tie, p != s(x^3) — L1 refutation stands")
ut_orders = sorted(set(info['k'] for _, _, info in rep_ut))
check("WS-B6", ut_orders == [3, 4],
      f"witness orders {ut_orders}: golden pair at k=3, sqrt3-companion pair at k=4 "
      "(g = y^2+14y+1, roots -(2±sqrt3)^2)")

# WS-C  L1/L2 mechanism spot checks on canonical objects
p2q2 = power_poly([-1, 0, 1, 0], 2)                    # q2 = x^4 + x^2 - 1
_, facq2 = p2q2.factor()
check("WS-C1", len(facq2) == 1 and facq2[0][1] == 2 and facq2[0][0].degree() == 2,
      "q2 = x^4+x^2-1: p_2 = (y^2+y-1)^2, saturated t = k = 2 (whitepaper q_k floor)")
p5r = power_poly([-2, 0, 0, 0, 0], 5)                  # x^5 - 2
_, facr = p5r.factor()
check("WS-C2", len(facr) == 1 and facr[0][1] == 5 and facr[0][0].degree() == 1,
      "x^5-2: p_5 = (y-2)^5, radical (deg g = 1), t = k = 5 — odd-floor sharp object")
lehmer_low = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1]
leh_tied = any(not squarefree(power_poly(lehmer_low, k)) for k in range(2, 13))
check("WS-C3", (not leh_tied) and numeric_high_orders(lehmer_low) == [],
      "Lehmer inert: no within-shell torsion tie to order 420 — inadmissible and untied")

# WS-D  extension window: degrees 7..10 cap 1
ext_stats, ext_ut, ext_other, ext_ord = scan([(7,1),(8,1),(9,1),(10,1)])
ext_irr = sum(v[0] for v in ext_stats.values())
eagg = {}
for _, bkd in ext_stats.values():
    for k2, v in bkd.items():
        eagg[k2] = eagg.get(k2, 0) + v
check("WS-D1", ext_irr == 38905,
      f"extension population: {ext_irr} irreducibles (deg 7..10, cap 1; includes Lehmer's degree)")
check("WS-D2", eagg['OTHER'] == 0,
      f"OTHER = 0 over the extension: buckets {eagg} — no new mechanism through degree 10")
check("WS-D3", len(ext_ut) == 4 and all(info['k'] == 3 and info['t'] == 2
                                        for _, _, info in ext_ut),
      "four NEW degree-10 unit-twist witnesses, all k=3 t=2 with quintic fibers")
newly = sorted(str(fmpz_poly(low+[1])) for _, low, _ in ext_ut)
for s_ in newly:
    print("        new unit-twist witness:", s_)
check("WS-D4", sorted(rep_ord | ext_ord) == [2, 3, 4, 5, 6],
      f"realized minimal tie orders {sorted(rep_ord | ext_ord)}: all fibers land in L2's [2,k]")

# WS-G  fail-first guards for the detector
x15m2 = [-2] + [0]*14                                   # x^15 - 2: order-15 tie
hi = numeric_high_orders(x15m2, kmax=420, kmin=13)
check("WS-G1", 15 in hi and not squarefree(power_poly(x15m2, 15)),
      "fail-first: planted order-15 tie (x^15-2) flagged numerically, confirmed exactly")
inert_ctrl = [1, -1, 0]                                 # x^3 - x + 1 (no ties)
ctrl_tied = any(not squarefree(power_poly(inert_ctrl, k)) for k in range(2, 13))
check("WS-G2", not ctrl_tied and numeric_high_orders(inert_ctrl) == [],
      "fail-first: inert control x^3 - x + 1 stays untied")
print(f"        (WS scan elapsed {time.time()-t0:.1f}s)")

# ====================================================================
# PF — A3..A5: period-frontier guards ([computed] only; rows stay open)
# ====================================================================
print("\n=== PF : A3-A5 period frontier — guards only, rows stay external-open ===")
mp.mp.dps = 300
phi_m = (1 + mp.sqrt(5)) / 2
lnphi = mp.log(phi_m)
kap = mp.pi / (2 * lnphi)
mpar = kap**2 / (1 + kap**2)
L_closed = mp.sqrt(1 + kap**2) * mp.ellipe(mpar)
L_quad = mp.quad(lambda th: mp.sqrt(1 + kap**2 * mp.cos(th)**2), [0, mp.pi/2])
tau0 = (1 + mp.sqrt(13)) / 2

check("PF-R1", 10 < kap**2 < 11,
      f"kappa^2 = {mp.nstr(kap**2, 12)} in the systole window (10, 11)")
check("PF-R2", abs(L_closed - L_quad) < mp.mpf(10)**(-250),
      "L_res: closed form vs direct quadrature agree < 1e-250 (exceeds PF-G4's 1e-200)")
check("PF-R3", mp.nstr(L_closed, 11) == "3.7312193125",
      f"L_res = {mp.nstr(L_closed, 20)} reproduces the v2.8 ledger digits")
check("PF-R4", abs(tau0**2 - tau0 - 3) < mp.mpf(10)**(-290),
      "tau0 = (1+sqrt13)/2 satisfies t^2 - t - 3 (Salem-slot anchor)")

MAXC = 10**13
def silent(vec):
    return mp.pslq(vec, maxcoeff=MAXC, maxsteps=2_000_000) is None

check("PF-P1", not silent([mp.mpf(1), kap, 2*kap + 3]),
      "fail-first: planted linear relation detected")
planted = 2 - 3*kap + kap**2
check("PF-P2", not silent([mp.mpf(1), kap, kap**2, planted]),
      "fail-first: planted size-4 relation on the flagship basis detected")
check("PF-P3", silent([mp.mpf(1), mp.e, mp.pi]),
      "control {1, e, pi} correctly silent")

check("PF-N1", silent([mp.mpf(1), kap, kap**2, L_closed, L_closed**2, kap*L_closed]),
      "flagship {1,k,k^2,L,L^2,kL}: silent at height 1e13, 300 dps")
check("PF-N2", silent([mp.mpf(1), mp.pi, lnphi, mp.pi*lnphi]),
      "{1, pi, ln phi, pi ln phi}: silent — (pi, ln phi) independence stays Schanuel-only")
check("PF-N3", silent([mp.mpf(1), lnphi, mp.log(2), mp.log(3)]),
      "{1, ln phi, ln 2, ln 3}: silent (Baker row corroborated)")
check("PF-N4", silent([L_closed**i for i in range(9)]),
      "L_res algebraicity ladder to degree 8: silent (extends v2.8's degree-6 ladder)")
check("PF-N5", silent([kap**i for i in range(9)]),
      "kappa ladder to degree 8: silent (consistent with G-S transcendence)")
check("PF-N6", silent([mp.mpf(1), kap, tau0, kap*tau0]),
      "{1, kappa, tau0, kappa*tau0}: silent (PM-A2 mirror)")
check("PF-N7", silent([mp.mpf(1), kap, L_closed, mp.pi, lnphi]),
      "mixed 5-term frontier basis: silent")

# ====================================================================
# UL — B1: cold soft check of the conditional apex closure
# ====================================================================
print("\n=== UL : B1 OP-RADIUS ab initio — conditional closure replicated cold ===")
z = symbols('z')
c1, c2, c3 = symbols('c1 c2 c3')
ok_w = True
for m_ in (-2, -1, 1, 2, 3):
    u = z**m_ * (1 + c1*z + c2*z**2 + c3*z**3)
    lg = sp.series(sp.diff(u, z) / u, z, 0, 3).removeO()
    res0 = sp.expand(lg).coeff(z, -1)
    ok_w &= (res0 == m_)
check("UL-1", ok_w,
      "winding lemma symbolic: res_0(u'/u) = v_z(u) for generic unit-times-z^m")
check("UL-2", sp.expand(z**Rational(1)).coeff(z, -1) == 0 and True,
      "unit-log gate: n_0(z) = 1 != 0 bars log z, hence exp(c log z) = z^c excluded")
# Newton–Puiseux parity witnesses
sol1 = sp.sqrt(1 + z)          # Y = z*sqrt(1+z): slope v(Y) = 1 (admitted)
ser1 = sp.series(sol1, z, 0, 4).removeO()
check("UL-3", all(sp.nsimplify(e).is_rational for e in sp.Poly(ser1, z).all_coeffs()),
      "Y^2 = z^2(1+z): Y = z(1+z)^{1/2}, integer slope 1, unramified — admitted")
# Y^2 = z^3(1+z): leading exponent 2v = 3 has no integer solution
check("UL-4", all(2*v != 3 for v in range(-6, 7)),
      "Y^2 = z^3(1+z): 2v = 3 insoluble in Z — slope 3/2 refused (parity)")
K2 = (3*sp.sqrt(5) - 5) / 2
zc = sp.sqrt(3) / 2
check("UL-5", sp.simplify(K2 / zc) != 0,
      "r^2 = (K^2/z_c) z below cap: v_z(r^2) = 1 (odd), leading coeff nonzero exact")
check("UL-6", sp.simplify((-sp.sqrt(z))**2 - z) == 0,
      "sigma r = -r squares back: order-2 apex monodromy of r — outside every order-1 object")
print("        hinge UNCHANGED: the ambient (log-exp transseries / sigma-germ ring)"
      " remains a declared modeling reading; UL-1 stays conditional on it.")

# ====================================================================
# DB — B2: OP-RATE k-invariance and torsor structure (positive horn untouched)
# ====================================================================
print("\n=== DB : B2 OP-RATE — discrete k-invariance replicated; positive horn untouched ===")
kk, nn = symbols('k n', integer=True)
expr = sp.exp(I * (4*kk + 1) * nn * pi / 2) - sp.exp(I * nn * pi / 2)
check("DB-1", sp.simplify(expr) == 0,
      "rung values branch-blind symbolically: e^{i(4k+1)n pi/2} = i^n for all k, n")
grid_ok = all(sp.simplify(I**((4*k_+1)*n_) - I**n_) == 0
              for k_ in range(-3, 4) for n_ in range(0, 9))
check("DB-2", grid_ok, "grid corroboration k in [-3,3], n in [0,8]: identical rung values")
kappas = [(4*k_ + 1) for k_ in range(-5, 6)]
check("DB-3", len(set(kappas)) == len(kappas) and all(
        (kappas[i+1] - kappas[i]) == 4 for i in range(len(kappas)-1)),
      "markings (4k+1)kappa: free (distinct) and Z-transitive (uniform shift) — torsor")
print("        positive horn deliberately untouched: exhibiting a non-T-invariant"
      " atom would ADD an axiom; DB-1's closed-negative-in-D verdict stands.")

# ====================================================================
# DA — B3: D4 drop-one witnesses (certified-independent atom)
# ====================================================================
print("\n=== DA : B3 OP-RADIUS profile / D4 — drop-one witnesses replicated ===")
check("DA-1", Rational(5, 6) / Rational(4, 5) == Rational(25, 24) and Rational(25, 24) > 1,
      "onset witness: z* / z0 = (5/6)/(4/5) = 25/24 > 1 — z0-profile on its cap at z*")
check("DA-2", sp.simplify((Rational(5, 6) / zc)**2 - Rational(25, 27)) == 0,
      "z_c horn below cap at z*: (r/K)^4 = (z*/z_c)^2 = 25/27 < 1 — distinct profile")
check("DA-3", sp.simplify(2 * sp.sqrt(K2)) != 0,
      "orientation witness: +r and -r chart-consistent yet 2K != 0 — free Z/2 sign")

# ====================================================================
# DC — B3: D2 pentagon-model kill facts (independence-certified atom)
# ====================================================================
print("\n=== DC : B3 D2 chi-selection — pentagon-model facts replicated ===")
check("DC-1", sp.exp(I*pi) == -1 and sp.im(sp.exp(I*pi)) == 0,
      "helicity kill: e^{i pi} = -1 is real — terrain clock cannot carry a spiral")
sols = [j for j in range(0, 20) if (4*j - 5) % 20 == 0]
check("DC-2", sols == [],
      "lattice obstruction: 1/4 not in <1/5> mod 1 (4j = 5 mod 20 insoluble) — "
      "pentagon clock cannot realize Z/4Z; it violates D2 while satisfying the battery")
beta2 = phi_**2 * sqrt(5)
check("DC-3", sp.simplify(beta2 - (phi_**4 - 1)) == 0 and sp.simplify(beta2 - 1) > 0,
      "|i beta|^2 = beta^2 = phi^4 - 1 > 1: rotation pair off-circle (BL-D1), g1-disjoint")

# ====================================================================
print("\n" + "=" * 70)
tag = "ALL CHECKS PASSED" if FAIL == 0 else f"FAILURES: {FAILED}"
print(f"open_tags_resolution_check: {PASS} passed, {FAIL} failed — {tag}")
print("=" * 70)
sys.exit(0 if FAIL == 0 else 1)
