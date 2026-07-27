# lesson6_census.py — independent fourth-operator replication of the quintic census (stage 1 + 2).
import sys
import sympy as sp
from sympy import Rational as R_
import time

QUICK = '--quick' in sys.argv   # smoke mode: the a = -1 slice of the box, ~40 s instead of ~20 min

x, y = sp.symbols('x y')
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)

def ratio_poly(p_expr):
    Rr = sp.resultant(p_expr.subs(x, y), sp.expand(p_expr.subs(x, x*y)), y)
    return sp.Poly(sp.expand(Rr), x).primitive()[1]

def totients_upto(N):
    phi = list(range(N + 1))
    for i in range(2, N + 1):
        if phi[i] == i:
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi

def contacts(P, phi_cache={}):
    d = P.degree()
    N = 2 * d * d
    if N not in phi_cache: phi_cache[N] = totients_upto(N)
    phi = phi_cache[N]
    hits = {}
    for m in range(1, N + 1):
        if phi[m] > d: continue
        cm = sp.Poly(sp.cyclotomic_poly(m, x), x)
        T, mult = P, 0
        while True:
            q, r = sp.div(T, cm)
            if r.is_zero: mult, T = mult + 1, q
            else: break
        if mult: hits[m] = mult
    return hits

def _complex_intervals(P, eps):
    """Return only the complex root rectangles, tolerating both SymPy shapes:
    older builds hand back one flat list, current builds hand back (reals, complexes)."""
    r = P.intervals(all=True, eps=eps)
    if isinstance(r, tuple) and len(r) == 2 and isinstance(r[0], list):
        return r[1]
    out = []
    for item in r:
        iv = item[0]
        if isinstance(iv, tuple) and len(iv) == 2 and not isinstance(iv[0], tuple):
            z1, z2 = iv
            if sp.im(z1) != 0 or sp.im(z2) != 0:
                out.append(item)
    return out

def disk_inside_all_complex(P):
    """Every non-real root strictly inside |z| < 1, certified by exact rational rectangles."""
    eps = R_(1, 10**3)
    for _ in range(7):
        undecided = False
        cplx = _complex_intervals(P, eps)
        for (z1, z2), _m in cplx:
            ax, ay = sp.re(z1), sp.im(z1)
            bx, by = sp.re(z2), sp.im(z2)
            mx = max(ax*ax, bx*bx) + max(ay*ay, by*by)
            lo_x = R_(0) if ax <= 0 <= bx else min(ax*ax, bx*bx)
            lo_y = R_(0) if ay <= 0 <= by else min(ay*ay, by*by)
            mn = lo_x + lo_y
            if mx < 1: continue
            if mn > 1: return False
            undecided = True
        if not undecided: return True
        eps = eps / 100
    raise ValueError("disk certificate did not terminate")

# ---------------- stage 1 ----------------
print("--- stage 1: certification cascade over the [-2,2]^5 box (3125 quintics) ---", flush=True)
t1 = time.time()
tal = dict(e0=0, recip=0, red=0, pattern=0, disk=0, pisot=0)
patterns = dict(real5=0, mixed=0, twopair=0)
pisots = []
rng = range(-2, 3)
arange = [-1] if QUICK else list(rng)
if QUICK:
    print("    [--quick] a = -1 slice only; stage-1 tallies below are NOT the published census",
          flush=True)
for a in arange:
 for b in rng:
  for c in rng:
   for d in rng:
    for e in rng:
        if e == 0: tal['e0'] += 1; continue
        v = [1, a, b, c, d, e]
        if v == v[::-1] or v == [-t for t in v[::-1]]:
            tal['recip'] += 1; continue
        P = sp.Poly([1, a, b, c, d, e], x)
        if not P.is_irreducible: tal['red'] += 1; continue
        if not (P.count_roots(1, sp.oo) == 1 and P.count_roots(-sp.oo, -1) == 0):
            tal['pattern'] += 1; continue
        nreal = P.count_roots(-sp.oo, sp.oo)
        try:
            inside = disk_inside_all_complex(P)
        except ValueError:
            inside = False
        if not inside: tal['disk'] += 1; continue
        tal['pisot'] += 1
        pisots.append((a, b, c, d, e, nreal))
        if nreal == 5: patterns['real5'] += 1
        elif nreal == 3: patterns['mixed'] += 1
        else: patterns['twopair'] += 1
print(f"    tallies: {tal}   patterns: {patterns}   ({time.time()-t1:.0f}s)", flush=True)
check("D1  stage-1 cascade matches the paper EXACTLY: 625/50/638/1318/411 -> 83 Pisot; "
      "patterns 0/16/67   [independent operator]",
      True if QUICK else
      ((tal['e0'], tal['recip'], tal['red'], tal['pattern'], tal['disk'], tal['pisot'])
       == (625, 50, 638, 1318, 411, 83)
       and (patterns['real5'], patterns['mixed'], patterns['twopair']) == (0, 16, 67)),
      "[--quick] cascade executed on a slice; full-box tallies not asserted" if QUICK else "")

# ---------------- stage 2 ----------------
print("--- stage 2: level-1 scans; Rat^o facts; shell detector ---", flush=True)
t2 = time.time()
lvl1_ok = rat0_ok = detector_ok = sqfree_ok = 0
reducible = []
c2_sample = []
tt = sp.symbols('tt')
zk = {0: sp.Integer(2), 1: tt}
for k in range(2, 11): zk[k] = sp.expand(tt*zk[k-1] - zk[k-2])
for (a, b, c, d, e, nreal) in pisots:
    p = x**5 + a*x**4 + b*x**3 + c*x**2 + d*x + e
    Rp = ratio_poly(p)
    if contacts(Rp) == {1: 5}: lvl1_ok += 1
    if nreal == 1:
        Rat0 = sp.Poly(sp.quo(Rp, sp.Poly((x - 1)**5, x)), x)
        sqfree = sp.gcd(Rat0, Rat0.diff()).degree() == 0
        irred = Rat0.is_irreducible
        if sqfree and Rat0.degree() == 20: sqfree_ok += 1
        if sqfree and irred and Rat0.degree() == 20: rat0_ok += 1
        if sqfree and not irred: reducible.append(sp.sstr(p))
        cf = Rat0.all_coeffs()
        T = sp.expand(sum(cf[i]*zk[10 - i] for i in range(10)) + cf[10])
        if 2 * sp.Poly(T, tt).count_roots(-2, 2) == 4: detector_ok += 1
        if (a, b, c, d, e) == (-2, -2, -2, -2, -2) or len(c2_sample) < 2 and lvl1_ok % 41 == 7:
            c2_sample.append((p, Rat0))
print(f"    level-1 {lvl1_ok}/83; Rat^o sqfree {sqfree_ok}/67, irreducible {rat0_ok}/67; "
      f"detector=4: {detector_ok}/67 "
      f"({time.time()-t2:.0f}s)", flush=True)
check("D2  all 83 level-1 scans return {Phi_1^5} (Thm 3.1's forced prediction, P3)",
      lvl1_ok == len(pisots) if QUICK else lvl1_ok == 83,
      "[--quick] slice only: level-1 holds on every Pisot found in the slice" if QUICK else "")
check("D3  ** THE ERRATUM, ENCODED ** Rat_p^o is squarefree of degree 20 on all 67 two-pair "
      "instances but irreducible on only 66. The single exception is the D5 instance "
      "x^5 - x^3 - 2x^2 - 2x - 1, whose Rat^o splits 10+10 into irreducible self-reciprocal factors. "
      "This check asserts the CORRECTED counts; the note's '67/67' is the defect it records",
      True if QUICK else (sqfree_ok == 67 and rat0_ok == 66
                          and reducible == ['x**5 - x**3 - 2*x**2 - 2*x - 1']),
      "[--quick] full-box counts not asserted" if QUICK
      else f"squarefree {sqfree_ok}/67, irreducible {rat0_ok}/67, exception {reducible}")
check("D4  shell detector reads 4 on all 67 (distinct shells; P8's same-shell case absent from the box)",
      True if QUICK else detector_ok == 67,
      "[--quick] full-box count not asserted" if QUICK else "")

print("--- C2 composed-square scans on a 2-instance sample (incl. the smallest, x^5-2x^4-...-2) ---",
      flush=True)
t3 = time.time()
c2_ok = 0
phi400 = totients_upto(1 if QUICK else 320000)
cands = [] if QUICK else [m for m in range(1, 320001) if phi400[m] <= 400]
print(f"    candidate cyclotomics with phi(m) <= 400: {len(cands)}"
      + (" [--quick] skipped" if QUICK else f" (largest m = {max(cands)})"), flush=True)
for p, S in ([] if QUICK else c2_sample[:2]):
    dS = S.degree()
    C2 = sp.Poly(sp.expand(sp.resultant(S.as_expr().subs(x, y),
                                        sp.expand(y**dS * S.as_expr().subs(x, x/y)), y)), x).primitive()[1]
    dd = C2.degree()
    hits = {}
    for m in cands:
        cm = sp.Poly(sp.cyclotomic_poly(m, x), x)
        T, mult = C2, 0
        while True:
            q, r = sp.div(T, cm)
            if r.is_zero: mult, T = mult + 1, q
            else: break
        if mult: hits[m] = mult
    print(f"      {sp.sstr(p)}: deg C2 = {dd}, contacts = {hits}   ({time.time()-t3:.0f}s)", flush=True)
    if dd == 400 and hits == {1: 20}: c2_ok += 1
check("D5  C2 negative certificate on the sample: degree 400, contacts exactly {Phi_1^20} — "
      "no mirrored cross-shell class (paper: 67/67)", True if QUICK else c2_ok == 2,
      "[--quick] C2 composed-square scans skipped" if QUICK else "")
check("D6  Burnside forced arithmetic: (729+27)/2 = 378 and (15625+125)/2 = 7875",
      (729 + 27)//2 == 378 and (15625 + 125)//2 == 7875)

print()
fails = [c for c in checks if not c[1]]
print(("[--quick] SMOKE RUN — the published census requires --full\n" if QUICK else "")
      + f"CENSUS SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
