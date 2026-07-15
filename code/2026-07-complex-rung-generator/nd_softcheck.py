#!/usr/bin/env python3
"""
nd_softcheck.py -- ODD-2 by relative-norm descent.  Written cold, 2026-07-15.
Register prefix ND.  Supersedes the prime-m reduction proposed in the session of
2026-07-15 (printed in odd2_box_softcheck.py [E]), whose hypothesis (R1) is refuted below.

PROVENANCE.  v2.0 of the paper never asserted (R1); it filed ODD-2 as an open
named lemma (blk:odd2) and is vindicated on that point.  (R1) was the PRIOR
SESSION'S proposed reduction, not the paper's.  What v2.0 does get wrong is the
framing that the residual gap [mu_S, 2) *is* the SSS/Lehmer core -- see [S].

WHAT BROKE
  That reduction split on  k := [Q(alpha):Q(alpha^m)]  dividing m, so that
  prime m forced k in {1, m}.  reduction_audit.py exhibits alpha = 2*zeta_3
  (min poly x^2+2x+4, on R_3): m = 3, d = 2, d' = 1, k = 2, and 2 does NOT divide 3.
  Degrees of irreducible factors of x^m - beta need not divide m.  The split is void,
  prime m is NOT closed by it, and the fallback M >= phi^(d/m) leaves a window
  d <= floor(1.4404 m) that GROWS with m -- an infinite family of finite searches,
  not one unfinished computation.

THE REPAIR -- never mentions k | m
  alpha an algebraic integer, not a root of unity, every conjugate on
  R_m = {z != 0 : arg z in (2pi/m)Z}, m odd.  Put beta = alpha^m, F = Q(beta),
  K = Q(alpha), d = [K:Q], d' = [F:Q], k = [K:F], so d = k d'.

    (1) beta is totally positive, hence F is totally real.
    (2) N := N_{K/F}(alpha) is an algebraic integer of F, and is TOTALLY POSITIVE:
        each Q-conjugate of N is a product of k conjugates of alpha, so its argument
        lies in (2pi/m)Z; it is also real; and pi is not in (2pi/m)Z mod 2pi for m odd.
    (3) The F-conjugates of alpha are exactly the alpha_j with alpha_j^m = beta, so
        each is alpha*zeta^{j_i} and  N = alpha^k * zeta^J.  Hence N = 1 would make
        alpha^{km} = 1.  So N != 1.
    (4) TP-2: a totally positive algebraic integer N != 1 has M(N) >= 2^r >= 2,
        where r = #{conjugates of N exceeding 1}; equality iff N = 2.
        ELEMENTARY -- no Schinzel, no external input at all.  Let g = minpoly(N),
        deg e.  N != 1 and g irreducible => g(1) != 0 => |g(1)| >= 1 (rational integer).
        Write A = {N_i > 1} (r = |A|), B = {N_i < 1}; no N_i = 1.  Then
            |g(1)| = prod_B (1 - N_i) * prod_A (N_i - 1),   prod_B (1 - N_i) <= 1,
        since every factor lies in (0,1).  Hence  prod_A (N_i - 1) >= |g(1)| >= 1.
        A is nonempty: |Norm N| = |g(0)| >= 1 forces some N_i >= 1, and none is 1.
        t |-> log(1 + e^t) is convex (second derivative 1/(4 cosh^2(t/2)) > 0), so
        Jensen gives the Mahler product inequality
            prod_A (1 + u_j) >= (1 + (prod_A u_j)^(1/r))^r    (u_j := N_j - 1 > 0),
        whence  M(N) = prod_A N_j = prod_A (1 + u_j) >= (1 + 1)^r = 2^r >= 2.  []
        [Schinzel's phi^(d/2) is retained in group [E] as CORROBORATION ONLY.]
    (5) h(N) <= sum_i h(sigma_i alpha) = k h(alpha), and h(N) = log M(N)/e, so
            M(alpha) = exp(d h(alpha)) >= M(N)^(d/(k e)) = M(N)^(d'/e).
    (6) Q(N) <= F so e | d', so d'/e is a positive integer and M(alpha) >= 2.  QED

  No divisibility hypothesis, no window, no external theorem, every odd m at once.
  The whole chain is elementary and internal: ODD-2 closes unconditionally.  Sharp: x^m - 2 has
  N = 2, e = d' = 1, bound 2^(1/1) = 2 = M.  Also sharp at x^6-3x^3+1 (N = phi^2,
  e = d' = 2, refined bound M(N)^(d'/e) = phi^2 = M).

Exit 0 on PASS.
"""

import json, sys, time, itertools, os
import numpy as np
import mpmath as mp
from sympy import (Poly, symbols, resultant, factor_list, CRootOf, ZZ,
                   minimal_polynomial, Rational, nsimplify)
from sympy import symbols as sym_symbols, diff, log, exp, cosh, simplify
from math import comb

mp.mp.dps = 140
x, y = symbols('x y')
PHI = (1 + mp.sqrt(5)) / 2
TOL = mp.mpf(10) ** (-40)

PASS, FAIL = [], []


def chk(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}: {msg}")
    return bool(cond)


def roots_hp(coeffs):
    return mp.polyroots([mp.mpf(int(c)) for c in coeffs], maxsteps=600, extraprec=600)


def mahler(coeffs):
    return mp.fprod([max(mp.mpf(1), abs(r)) for r in roots_hp(coeffs)])


def on_ray(rts, m):
    for r in rts:
        a = mp.arg(r) * m / (2 * mp.pi)
        if abs(a - mp.nint(a)) > mp.mpf('1e-30'):
            return False
    return True


def certify(coeffs, m):
    """Full norm-descent certificate for one on-ray object. Returns dict or raises AssertionError."""
    f = Poly(coeffs, x)
    d = f.degree()
    rts = roots_hp(coeffs)
    assert on_ray(rts, m), "not on R_m at high precision"

    # --- beta = alpha^m : exact.  Res_x(f(x), y - x^m) = prod_i (y - alpha_i^m) = minpoly_beta^k
    g_full = Poly(resultant(f.as_expr().subs(x, x), y - x**m, x), y)
    fl = factor_list(g_full.as_expr(), y)
    facs = [(Poly(p, y), e) for p, e in fl[1]]
    assert len(facs) == 1, f"Res factors into {len(facs)} distinct irreducibles, expected a pure power"
    gb, k = facs[0]
    dprime = gb.degree()
    assert d == k * dprime, f"tower law fails: d={d} k={k} d'={dprime}"

    # beta = alpha^m for the distinguished root alpha := rts[0]
    alpha = rts[0]
    beta = alpha ** m
    assert abs(mp.im(beta)) < TOL and mp.re(beta) > 0, "beta not positive real"
    beta = mp.re(beta)

    # totally positive: every root of gb is a positive real
    gbc = [int(c) for c in gb.all_coeffs()]
    brts = roots_hp(gbc)
    tot_pos = all(abs(mp.im(b)) < TOL and mp.re(b) > TOL for b in brts)
    assert tot_pos, "beta not totally positive"

    # --- N = prod of the F-conjugates of alpha = { alpha_j : alpha_j^m = beta }
    S = [r for r in rts if abs(r**m - beta) < mp.mpf('1e-30')]
    assert len(S) == k, f"|{{alpha_j : alpha_j^m = beta}}| = {len(S)}, expected k = {k}"
    N = mp.fprod(S)
    assert abs(mp.im(N)) < TOL, "N not real"
    N = mp.re(N)
    assert N > 0, "N not positive"

    # identity  N = alpha^k zeta^J
    w = N / alpha**k
    assert abs(w**m - 1) < TOL, "N / alpha^k is not an m-th root of unity"

    # --- min poly of N.  However found, it is CERTIFIED below: monic, irreducible
    #     over Z, vanishing at N far beyond root separation, degree dividing d'.
    #     Those four facts force it to be the minimal polynomial.
    if k == 1:
        gN = f                                   # N = alpha itself
    elif dprime == 1:
        Ni = mp.nint(N)
        assert abs(N - Ni) < TOL, "N not a rational integer though d'=1"
        gN = Poly([1, -int(Ni)], x)
    else:
        gN = None
        for e0 in [t for t in range(1, dprime + 1) if dprime % t == 0]:
            rel = mp.pslq([N**j for j in range(e0, -1, -1)],
                          tol=mp.mpf(10) ** (-(mp.mp.dps - 45)), maxcoeff=10**22, maxsteps=10**6)
            if not rel or abs(rel[0]) != 1:
                continue
            sgn = 1 if rel[0] > 0 else -1
            cand = Poly([int(c) * sgn for c in rel], x)
            resid = abs(mp.polyval([mp.mpf(int(c)) for c in cand.all_coeffs()], N))
            if cand.is_irreducible and resid < mp.mpf('1e-45'):
                gN = cand                      # monic + irreducible + vanishes => IS the min poly
                break
        assert gN is not None, "no monic irreducible relation for N of degree dividing d'"

    assert gN.LC() == 1, "N not an algebraic integer (min poly not monic)"
    assert gN.is_irreducible, "candidate min poly of N is reducible"
    assert abs(mp.polyval([mp.mpf(int(c)) for c in gN.all_coeffs()], N)) < mp.mpf('1e-30'), "min poly of N does not vanish at N"
    e = gN.degree()
    assert dprime % e == 0, f"e={e} does not divide d'={dprime}"

    Nrts = roots_hp([int(c) for c in gN.all_coeffs()])
    assert all(abs(mp.im(t)) < TOL and mp.re(t) > TOL for t in Nrts), "N not totally positive"
    assert not (e == 1 and abs(N - 1) < TOL), "N = 1 (alpha would be a root of unity)"

    M_a, M_N = mahler(coeffs), mahler([int(c) for c in gN.all_coeffs()])
    assert M_N >= 2 - TOL, f"TP-2 violated: M(N) = {M_N}"

    # --- TP-2 by the ELEMENTARY chain, certified on this object (no Schinzel) ---
    gN1 = abs(mp.polyval([mp.mpf(int(c)) for c in gN.all_coeffs()], mp.mpf(1)))
    assert gN1 >= 1 - TOL, f"|g_N(1)| = {gN1} < 1 -- N = 1?"
    Are = [mp.re(t) for t in Nrts if mp.re(t) > 1]
    Bre = [mp.re(t) for t in Nrts if mp.re(t) < 1]
    rA = len(Are)
    assert rA >= 1, "no conjugate of N exceeds 1 (contradicts |Norm N| >= 1)"
    prodB = mp.fprod([1 - b for b in Bre]) if Bre else mp.mpf(1)
    prodA = mp.fprod([a - 1 for a in Are])
    assert prodB <= 1 + TOL, f"prod_B (1-N_i) = {prodB} > 1"
    assert prodA >= 1 - TOL, f"prod_A (N_i-1) = {prodA} < 1"
    assert M_N >= 2**rA - mp.mpf('1e-20'), f"elementary floor 2^r violated: M(N)={M_N} r={rA}"

    bound_new = M_N ** (mp.mpf(dprime) / e)
    bound_old = PHI ** (mp.mpf(d) / m)
    assert M_a >= bound_new - mp.mpf('1e-25'), f"norm-descent bound violated: M={M_a} < {bound_new}"
    assert M_a >= 2 - mp.mpf('1e-25'), f"ODD-2 violated: M = {M_a}"
    # tower identity R2 (survives the R1 refutation)
    assert abs(M_a**m - mahler(gbc)**k) < mp.mpf('1e-20') * M_a**m, "R2 tower identity failed"

    return dict(m=m, d=d, dp=dprime, k=k, e=e, M=M_a, MN=M_N, r=rA, g1=gN1, pa=prodA,
                new=bound_new, old=bound_old, f=f.as_expr(), gN=gN.as_expr())


print("=" * 78)
print("ODD-2 BY RELATIVE-NORM DESCENT  |  closes the paper's named-lemma blocker (blk:odd2)")
print("=" * 78)
t0 = time.time()

# ---------------------------------------------------------------- [R] refutation lock
print("\n[R] refutation lock -- what the repair must NOT depend on")
r = certify([1, 2, 4], 3)
chk("ND-R1", r['k'] == 2 and 3 % r['k'] != 0,
    f"alpha = 2*zeta_3 on R_3: d={r['d']} d'={r['dp']} k={r['k']}, k does NOT divide m=3  [R1 dead]")
chk("ND-R2", True, "R2 (M(alpha)^m = M(beta)^k) is the tower law d = k d' only -- survives, asserted per object below")
gaps = [(m, int(np.ceil(float(mp.log(2)/mp.log(PHI)) * m)) - 1) for m in (5, 7, 9, 11, 13, 15)]
chk("ND-R3", all(g > 5 for _, g in gaps),
    f"fallback phi^(d/m) window grows without bound: {', '.join(f'm={m}:d<={g}' for m, g in gaps)} -- box route cannot close ODD-2")

# ---------------------------------------------------------------- [T] TP-2, complete boxes
print("\n[T] TP-2  (M(N) >= 2 for totally positive N != 1) -- complete forced boxes")
# every violator of M >= 2 has |a_0| < 2 and |a_j| <= C(d,j)*2 (Landau/Mahler), so the box is forced
for dd in (1, 2, 3, 4):
    axes = [range(-2 * comb(dd, j), 2 * comb(dd, j) + 1) for j in range(dd)]
    viol, mins, munits, ntot = [], [], [], 0
    for t in itertools.product(*axes):
        if t[0] == 0:
            continue
        ntot += 1
        c = [1] + list(reversed(t))
        rr = np.roots(c)
        if np.max(np.abs(rr.imag)) > 1e-9 or np.min(rr.real) <= 1e-9:
            continue
        pp = Poly(c, x)
        if not pp.is_irreducible:
            continue
        M = mahler(c)
        if M <= 1 + mp.mpf('1e-30'):
            continue
        mins.append((float(M), pp.as_expr()))
        if abs(c[-1]) == 1:
            munits.append((float(M), pp.as_expr()))
        if M < 2 - mp.mpf('1e-20'):
            viol.append((pp.as_expr(), float(M)))
    mins.sort(key=lambda t: t[0]); munits.sort(key=lambda t: t[0])
    ok = chk(f"ND-T{dd}", len(viol) == 0,
             f"d'={dd}: box {ntot:>7,} | violators of M>=2: {len(viol)} | "
             f"min M = {mins[0][0]:.6f} at {mins[0][1]}" if mins else
             f"d'={dd}: box {ntot:>7,} | violators of M>=2: {len(viol)} | no non-cyclotomic witness in box")
    if munits:
        print(f"         units only: min M = {munits[0][0]:.6f} at {munits[0][1]}")

# ---------------------------------------------------------------- [E] TP-2 is ELEMENTARY
print("\n[E] TP-2 without Schinzel -- the external dependency is removed, not merely cited")

# ND-E1: the only inequality the elementary proof uses, certified symbolically.
tt = sym_symbols('tt', real=True)
d2 = simplify(diff(log(1 + exp(tt)), tt, 2))
conv = simplify(d2 - 1 / (4 * cosh(tt / 2) ** 2)) == 0
rng = np.random.default_rng(20260715)
jen = True
for _ in range(20000):
    rr_ = int(rng.integers(1, 7))
    u = np.exp(rng.normal(0, 2.2, rr_))                      # u_j > 0, wide dynamic range
    lhs = float(np.prod(1.0 + u))
    rhs = float((1.0 + np.prod(u) ** (1.0 / rr_)) ** rr_)
    if lhs < rhs * (1 - 1e-12):
        jen = False
        break
chk("ND-E1", bool(conv) and jen,
    f"Jensen step: d2/dt2 log(1+e^t) = {d2} = 1/(4cosh^2(t/2)) > 0 [convex] "
    f"=> prod(1+u_j) >= (1+(prod u_j)^(1/r))^r; 20,000 randomized instances, 0 violations")

# ND-E2: the elementary chain, per object, over a WIDE totally-positive box.
# (Wider than the M<2 violator box of [T]: it exercises r >= 2, which the [T] box cannot.)
tot_e, byr, minM_e, tight = 0, {}, mp.inf, []
for dd, cap in ((1, 12), (2, 10), (3, 8), (4, 6)):
    for t in itertools.product(*([range(-cap, cap + 1)] * dd)):
        if t[0] == 0:
            continue
        c = [1] + list(reversed(t))
        rts_ = np.roots(c)
        if np.max(np.abs(rts_.imag)) > 1e-9 or np.min(rts_.real) <= 1e-9:
            continue
        pp = Poly(c, x)
        if not pp.is_irreducible:
            continue
        M = mahler(c)
        if M <= 1 + mp.mpf('1e-30'):
            continue                                          # N = 1 / cyclotomic
        rr_ = np.sort(rts_.real)
        g1 = abs(int(pp.eval(1)))
        A = [v for v in rr_ if v > 1]
        B = [v for v in rr_ if v < 1]
        pa = float(np.prod([v - 1 for v in A]))
        pb = float(np.prod([1 - v for v in B])) if B else 1.0
        ok = (g1 >= 1 and len(A) >= 1 and pb <= 1 + 1e-9 and pa >= 1 - 1e-9
              and M >= 2 ** len(A) - mp.mpf('1e-9'))
        assert ok, f"elementary chain broke at {pp.as_expr()}"
        tot_e += 1
        byr[len(A)] = byr.get(len(A), 0) + 1
        if M < minM_e:
            minM_e, argmin_e = M, pp.as_expr()
        if M < 2 ** len(A) * mp.mpf('1.3'):
            tight.append((float(M), len(A), pp.as_expr()))
chk("ND-E2", tot_e > 0,
    f"elementary chain (|g(1)|>=1, A nonempty, prod_B<=1, prod_A>=1, M>=2^r) holds on "
    f"{tot_e} totally-positive irreducibles over the wide box; r-census {dict(sorted(byr.items()))}")
chk("ND-E3", abs(minM_e - 2) < TOL and str(argmin_e) == 'x - 2',
    f"equality analysis: min M = {float(minM_e):.6f} at {argmin_e}; the elementary floor is attained "
    f"only at N = 2 (r=1, prod_A = 1)")
for v in sorted(tight)[:3]:
    print(f"         tight vs 2^r: M = {v[0]:.4f}  r = {v[1]}   {v[2]}")

# ND-E5 (Schinzel, DEMOTED to corroboration): totally-real complete boxes.
# Statement as published -- Smyth, "The Mahler measure of algebraic numbers: a survey",
# arXiv:math/0701397, eq. (E-S), citing Schinzel, Acta Arith. 24 (1973) 385-399,
# Addendum ibid. 26 (1974/75) 329-331:  alpha totally real, alpha != 0,+-1, deg d
# => M(alpha) >= phi^(d/2).  NOT load-bearing anywhere in this file.
viol_s, nbox = [], 0
for dd in (2, 3, 4):
    Bd = float(PHI ** (mp.mpf(dd) / 2))
    axes = [range(-int(np.floor(comb(dd, j) * Bd)), int(np.floor(comb(dd, j) * Bd)) + 1) for j in range(dd)]
    for t in itertools.product(*axes):
        if t[0] == 0:
            continue
        c = [1] + list(reversed(t))
        nbox += 1
        rts_ = np.roots(c)
        if np.max(np.abs(rts_.imag)) > 1e-9:
            continue                                          # not totally real
        # cheap float screen: only a candidate violator earns the exact work
        Mf = float(np.prod(np.maximum(1.0, np.abs(rts_))))
        if Mf > Bd + 1e-6:
            continue
        pp = Poly(c, x)
        if not pp.is_irreducible:
            continue
        M = mahler(c)
        if M <= 1 + mp.mpf('1e-30'):
            continue                                          # M=1: cyclotomic/+-1 (Kronecker)
        if M < PHI ** (mp.mpf(dd) / 2) - mp.mpf('1e-20'):
            viol_s.append((pp.as_expr(), float(M)))
chk("ND-E5", len(viol_s) == 0,
    f"Schinzel M >= phi^(d/2) corroborated on complete totally-real boxes d<=4 ({nbox:,} polys, "
    f"{len(viol_s)} violators) -- CORROBORATION ONLY; no step above depends on it")

# ---------------------------------------------------------------- [N] certificate over the forced population
print("\n[N] norm-descent certificate over the forced on-ray population")
_POP = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pop.json')  # next to this file (CWD-independent)
pop = {int(m): [tuple(c) for c in v] for m, v in json.load(open(_POP)).items()}
rows, bad, byshape = [], [], {}
for m in sorted(pop):
    for c in pop[m]:
        try:
            res = certify(list(c), m)
        except AssertionError as ex:
            bad.append((m, c, str(ex)))
            continue
        rows.append(res)
        byshape.setdefault((res['m'], res['k'], res['dp'], res['e']), 0)
        byshape[(res['m'], res['k'], res['dp'], res['e'])] += 1
chk("ND-N1", len(bad) == 0, f"{len(rows)} on-ray objects certified, {len(bad)} certificate failures")
for b in bad[:5]:
    print("      !!!", b)
chk("ND-N2", all(r['M'] >= 2 - mp.mpf('1e-25') for r in rows), f"ODD-2 holds on all {len(rows)}: min M = {float(min(r['M'] for r in rows)):.6f}")
chk("ND-N3", all(r['M'] >= r['new'] - mp.mpf('1e-25') for r in rows), "norm-descent bound M >= M(N)^(d'/e) holds on all")
chk("ND-N4", all(r['dp'] % r['e'] == 0 for r in rows), "e | d' on all")
chk("ND-E4", all(r['MN'] >= 2**r['r'] - mp.mpf('1e-20') and r['pa'] >= 1 - TOL and r['g1'] >= 1 - TOL
                 for r in rows),
    f"every N arising in the {len(rows)} descents clears TP-2 by the ELEMENTARY floor 2^r "
    f"(r-census {dict(sorted({rr: sum(1 for q in rows if q['r'] == rr) for rr in set(q['r'] for q in rows)}.items()))}): "
    f"Schinzel is invoked nowhere in the certificate")
nontriv = [r for r in rows if 1 < r['k'] < r['m']]
chk("ND-N5", len(nontriv) > 0 and all(r['M'] >= 2 - mp.mpf('1e-25') for r in nontriv),
    f"{len(nontriv)} objects with 1 < k < m (the composite gap the old split could not reach) all satisfy ODD-2")

# ---------------------------------------------------------------- [S] strict dominance
print("\n[S] separation: where the old fallback dies and the descent does not")
sep = [r for r in rows if r['old'] < 2]
chk("ND-S1", len(sep) > 0 and all(r['new'] >= 2 - TOL for r in sep),
    f"{len(sep)}/{len(rows)} objects have phi^(d/m) < 2 (old bound useless); descent gives >= 2 on every one")
seen = set()
print("    %3s %3s %3s %3s %3s %10s %13s %10s   min poly" % ("m","d","dp","k","e","phi^(d/m)","M(N)^(dp/e)","M(alpha)"))
_uniq = []
for r in sorted(sep, key=lambda r: (r['m'], r['d'], r['k'], r['e'], float(r['M']))):
    key = (r['m'], r['d'], r['k'], r['e'])
    if key in seen:
        continue
    seen.add(key); _uniq.append(r)
for r in _uniq[:14]:
    print(f"    {r['m']:>3} {r['d']:>3} {r['dp']:>3} {r['k']:>3} {r['e']:>3} {float(r['old']):>10.4f} "
          f"{float(r['new']):>13.4f} {float(r['M']):>10.4f}   {r['f']}")

# ---------------------------------------------------------------- [X] sharpness
print("\n[X] sharpness of the descent")
for m in (3, 5, 7, 9, 11, 15):
    c = [1] + [0] * (m - 1) + [-2]
    r = certify(c, m)
    chk(f"ND-X{m}", abs(r['M'] - 2) < TOL and abs(r['new'] - 2) < TOL,
        f"x^{m}-2: N={r['gN']}, e={r['e']}, d'={r['dp']} -> bound = {float(r['new']):.10f} = M = {float(r['M']):.10f}  [equality]")
r = certify([1, 0, 0, -3, 0, 0, 1], 3)
chk("ND-XP", abs(r['M'] - PHI**2) < TOL and abs(r['new'] - PHI**2) < TOL,
    f"x^6-3x^3+1: N={r['gN']} = phi^2, e=d'={r['e']} -> refined bound = M = phi^2 = {float(r['M']):.10f}  [equality]")

# ---------------------------------------------------------------- [A] anchor normalisation
print("\n[A] anchor normalisation -- from Delta <= Z/m (pairwise differences) to on-ray")
# Delta is generated by pairwise differences of t(alpha_i) = Arg(alpha_i)/2pi (v2.0 l.1299).
# So Delta <= Z/m only forces  t(alpha_i) in theta0 + (1/m)Z  for a common theta0, NOT theta0 = 0.
# Conjugation symmetry of a real min poly gives {t_i} = {-t_i}, hence 2*theta0 in (1/m)Z, hence
# theta0 in (1/2m)Z.  For m ODD the index-2 quotient leaves exactly two cosets:
#   theta0 in (1/m)Z          -> alpha  is on R_m
#   theta0 in 1/2m + (1/m)Z   -> t(-alpha_i) = t_i + 1/2 in (m+1)/2m + (1/m)Z = (1/m)Z, since 2|(m+1)
#                             -> -alpha is on R_m,  and  M(-alpha) = M(alpha).
# So the sign twist is forced AND exhaustive: exactly one of +-alpha is on R_m.
BOXA = {2: 12, 3: 10, 4: 8}
tw_tot, tw_cert, both, neither = 0, 0, 0, 0
minM_tw = mp.inf
for m in (3, 5, 7, 9):
    for dd, cap in BOXA.items():
        for t in itertools.product(*([range(-cap, cap + 1)] * dd)):
            if t[0] == 0:
                continue
            c = [1] + list(reversed(t))
            rr = np.roots(c)
            tt = np.angle(rr) * m / (2 * np.pi)
            diffs = tt[:, None] - tt[None, :]
            if not np.all(np.abs(diffs - np.round(diffs)) < 1e-7):
                continue                                     # Delta not <= Z/m
            pp = Poly(c, x)
            if not pp.is_irreducible or mahler(c) <= 1 + mp.mpf('1e-30'):
                continue
            plus = np.all(np.abs(tt - np.round(tt)) < 1e-7)
            cneg = [v * (-1) ** i for i, v in enumerate(c)]   # monic rep of f(-x): (-1)^d f(-x)
            tn = np.angle(np.roots(cneg)) * m / (2 * np.pi)
            minus = np.all(np.abs(tn - np.round(tn)) < 1e-7)
            if plus and minus:
                both += 1
            if not plus and not minus:
                neither += 1
            if minus and not plus:                            # genuinely twisted anchor
                tw_tot += 1
                r = certify([int(v) for v in cneg], m)
                assert abs(mahler(c) - r['M']) < TOL, "M(-alpha) != M(alpha)"
                minM_tw = min(minM_tw, r['M'])
                tw_cert += 1
chk("ND-A1", both == 0 and neither == 0,
    f"dichotomy exhaustive on the forced box: 0 objects with both +-alpha on-ray, 0 with neither "
    f"(both={both}, neither={neither})")
chk("ND-A2", tw_tot > 0 and tw_cert == tw_tot,
    f"{tw_tot} twisted-anchor objects (only -alpha on-ray); all certified, M(-alpha) = M(alpha) on every one")
chk("ND-A3", minM_tw >= 2 - mp.mpf('1e-25'),
    f"ODD-2 holds on the twisted branch too: min M = {float(minM_tw):.6f}")

# ---------------------------------------------------------------- shape census
print("\n[census] (m, k, d', e) shapes realised in the forced population")
print("    %3s %3s %3s %3s %6s   k|m ?" % ("m","k","dp","e","count"))
for (m, k, dp, ee), n in sorted(byshape.items()):
    print(f"    {m:>3} {k:>3} {dp:>3} {ee:>3} {n:>6}   {'yes' if m % k == 0 else 'NO  <-- R1 fails here'}")

print("\n" + "=" * 78)
print(f"ND: {len(PASS)}/{len(PASS)+len(FAIL)} checks passed | {len(rows)} objects certified | {time.time()-t0:.1f}s")
print("=" * 78)
sys.exit(0 if not FAIL else 1)
