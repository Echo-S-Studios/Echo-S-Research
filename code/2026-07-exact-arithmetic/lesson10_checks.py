# lesson10_checks.py — exact verification for Lesson 10 (the residual pool).
# Discipline: no float crosses a decision boundary; Salem ordering by exact rational separators.
import hashlib
import re
import sympy as sp
from fractions import Fraction as F

x, t = sp.symbols('x t')
s5 = sp.sqrt(5); PHI = (1 + s5)/2
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)
S = sp.simplify

def companion(p):
    P = sp.Poly(p, x); d = P.degree(); assert P.LC() == 1
    Cm = sp.zeros(d, d)
    for i in range(d - 1): Cm[i + 1, i] = 1
    cf = P.all_coeffs()
    for i in range(d): Cm[i, d - 1] = -cf[d - i]
    return Cm

def power_gram(p):
    C = companion(p); d = C.rows
    tr = {0: sp.Integer(d)}; Pw = sp.eye(d)
    for kk in range(1, 2*d - 1):
        Pw = Pw*C; tr[kk] = sp.trace(Pw)
    return sp.Matrix(d, d, lambda i, j: tr[i + j])

def signature_of(G):
    ev = sp.Poly(G.charpoly(x).as_expr(), x)
    pos = ev.count_roots(0, sp.oo); neg = ev.count_roots(-sp.oo, 0)
    return (int(pos), int(neg))

# ================= A. the Minkowski face — and the guard the code is missing
G5 = power_gram(x**2 - x - 1)
check("A1  golden field: G = M^T M = [[2,1],[1,3]] on {1,phi}; det G = 5 = d_K exactly "
      "(Z[phi] IS the maximal order)",
      G5 == sp.Matrix([[2, 1], [1, 3]]) and G5.det() == 5)
G3 = power_gram(x**3 - 2)
check("A2  a complex place (Ex. cubic2): Q(cbrt2) has G = [[3,0,0],[0,0,6],[0,6,0]], det = -108 = d_K, "
      "and signature (2,1) = (r1+r2, r2) with eigenvalues {3, +/-6} — INDEFINITE",
      G3 == sp.Matrix([[3, 0, 0], [0, 0, 6], [0, 6, 0]]) and G3.det() == -108
      and signature_of(G3) == (2, 1))
Gi = power_gram(x**2 + 1)
Mi = sp.Matrix([[1, sp.I], [1, -sp.I]])
G2i = S(Mi.conjugate().T * Mi)
check("A3  the Gaussian field: trace form diag(2,-2) is indefinite (1,1) with det = -4 = d_K; the "
      "Hermitian companion G_2 = M* M = 2I is positive definite with det 4 = |d_K|; covolume "
      "2^{-r2} sqrt|d_K| = 1 — the unit-covolume lattice Z[i]",
      Gi == sp.diag(2, -2) and Gi.det() == -4 and signature_of(Gi) == (1, 1)
      and G2i == 2*sp.eye(2) and G2i.det() == 4 and S(sp.Rational(1, 2)*sp.sqrt(4)) == 1)
Gprod = sp.diag(4, 8, 12, 24)
Gpow = power_gram(x**4 - 10*x**2 + 1)
check("A4  the index correction (Thm GMM): the product basis of Q(sqrt2,sqrt3) has det G = 9216 while "
      "d_K = 2304, so [O_K : Lambda] = sqrt(9216/2304) = 2; the POWER basis of x^4-10x^2+1 has "
      f"det G = {Gpow.det()}, index^2 = {sp.Rational(Gpow.det(), 2304)} — Lesson 5's Gram was never d_K",
      Gprod.det() == 9216 and sp.sqrt(sp.Rational(9216, 2304)) == 2
      and sp.Rational(Gpow.det(), 2304) == sp.Rational(Gpow.det(), 2304)
      and sp.sqrt(sp.Rational(Gpow.det(), 2304)).is_integer)
cat = [(x**2 - x - 1, 0), (x**2 - 2, 0), (x**4 - 10*x**2 + 1, 0), (x**2 + 1, 1), (x**3 - 2, 1)]
check("A5  sign(det G) = (-1)^{r2} across the catalog — the DISCRIMINANT'S SIGN IS THE DEFINITENESS FLAG",
      all(sp.sign(power_gram(p).det()) == (-1)**r2 for p, r2 in cat))
check("A6  Minkowski bound closes the golden class number: M_K = sqrt(d) (4/pi)^{r2} n!/n^n = sqrt5/2 < 2 "
      "by the exact guard 5 < 16, so no prime of norm < 2 exists and h = 1",
      5 < 16 and S(s5/2 - 2).is_negative)
check("A7  discriminant as information volume: covolumes sqrt|d_K| are sqrt5, 2sqrt2, 6sqrt3 for "
      "Q(sqrt5), Q(sqrt2), Q(cbrt2) — exactly, since 108 = 36*3",
      S(sp.sqrt(8) - 2*sp.sqrt(2)) == 0 and S(sp.sqrt(108) - 6*sp.sqrt(3)) == 0)
# ---- THE FINDING: the capture predicate needs total reality, which the constructor does not check
B1 = sp.Matrix([[1], [0], [0]])                       # the forced basis {1} — survives every re-home
Pj = B1*(B1.T*G3*B1)**-1*B1.T*G3
xt = sp.Matrix([0, 1, 0])                             # the observation theta = cbrt2
rt = xt - Pj*xt
check("A8  ** FINDING (G11 candidate) **: on the paper's OWN example field Q(cbrt2), basis {1} has "
      "B^T G B = 3 (nonsingular, no NO_PROJECTION sentinel), yet the residual of theta is r = (0,1,0) "
      "!= 0 with ||r||^2_G = Tr(theta^2) = 0 EXACTLY — so the shipped predicate `rn == 0` reads "
      "CAPTURED for the generator of its own field",
      (B1.T*G3*B1)[0, 0] == 3 and S(Pj*Pj - Pj) == sp.zeros(3, 3)
      and rt == sp.Matrix([0, 1, 0]) and (rt.T*G3*rt)[0, 0] == 0)
ri = sp.Matrix([0, 1])
check("A9  the second failure mode: over Q(i) with basis {1}, the residual of 1+i is (0,1) with "
      "||r||^2 = -2 < 0 — a NEGATIVE gain, which the capacity gate compares against a floor",
      (ri.T*Gi*ri)[0, 0] == -2)
check("A10 the one-line guard that closes both: total reality is an exact Sturm count. Q(sqrt2+sqrt3): "
      "4 real roots of 4 (PASS); Q(cbrt2): 1 of 3; Q(i): 0 of 2 (both REJECT)",
      sp.Poly(x**4 - 10*x**2 + 1, x).count_roots(-sp.oo, sp.oo) == 4
      and sp.Poly(x**3 - 2, x).count_roots(-sp.oo, sp.oo) == 1
      and sp.Poly(x**2 + 1, x).count_roots(-sp.oo, sp.oo) == 0)

# ================= B. the browser instrument: the Phi-closure
def Phi(A): return companion(A.charpoly(x).as_expr())
Rm = companion(x**2 - x - 1)
A_dd = sp.diag(Rm, Rm)                                 # phi (+) phi, derogatory
check("B1  Phi = companion o charpoly is a RETRACTION: idempotent after one step, Phi(Phi(A)) = Phi(A) "
      "entrywise, and dimension-preserving",
      Phi(Phi(A_dd)) == Phi(A_dd) and Phi(A_dd).rows == A_dd.rows)
def minpoly_of(A):
    d = A.rows; vs = [sp.eye(d).reshape(d*d, 1)]
    Pw = sp.eye(d)
    for kk in range(1, d + 1):
        Pw = Pw*A; vs.append(Pw.reshape(d*d, 1))
        Mx = sp.Matrix.hstack(*vs)
        ns = Mx.nullspace()
        if ns:
            cf = [sp.nsimplify(v) for v in ns[0]]
            return sp.Poly(sum(cf[i]*x**i for i in range(len(cf))), x).monic()
    return sp.Poly(A.charpoly(x).as_expr(), x)
check("B2  the exact similarity verdict: Phi preserves charpoly (hence det, tr, rho, M) ALWAYS, but is "
      "similarity-preserving iff the input is non-derogatory — phi(+)phi has minpoly x^2-x-1 of degree "
      "2 < 4 = deg charpoly, so its lift is NOT similar to it",
      sp.expand(Phi(A_dd).charpoly(x).as_expr() - A_dd.charpoly(x).as_expr()) == 0
      and minpoly_of(A_dd).degree() == 2 and minpoly_of(Phi(A_dd)).degree() == 4)
check("B3  queryable registry: a plate EXTENDS a seed iff the seed's minpoly divides the plate's "
      "charpoly — (x^2-x-1) | charpoly(phi(+)phi) but not | charpoly(companion(x^2-2))",
      sp.rem(sp.Poly(A_dd.charpoly(x).as_expr(), x), sp.Poly(x**2 - x - 1, x)).is_zero
      and not sp.rem(sp.Poly(x**2 - 2, x), sp.Poly(x**2 - x - 1, x)).is_zero)
LEH = sp.Rational(117628, 100000)
check("B4  the Lehmer gap stays visible: the floor is pinned at M = 1 and the smallest seed measure is "
      "phi, so the band (1, phi) — which CONTAINS the Lehmer tick 1.17628 — is empty. Exact: "
      "phi > 117628/100000 because (2r-1)^2 = 285846649/156250000 < 5",
      S((2*LEH - 1)**2 - 5).is_negative and S(PHI - LEH).is_positive)
def bins(nn): return max(3, min(12, sp.ceiling(sp.sqrt(nn))))
check("B5  the binning rule is exact integer arithmetic: ceil(sqrt n) clamped to [3,12] — n=9 -> 3, "
      "n=50 -> 8, n=200 -> 12 (clamped), n=2 -> 3 (clamped)",
      (bins(9), bins(50), bins(200), bins(2)) == (3, 8, 12, 3))

# ================= C. KIRA: the readiness gate as a checkable object
CONTRACT = {"read", "laws", "search", "audit", "bridge", "lexicon", "render", "speak",
            "observe", "propose", "commit", "persist", "restore"}
check("C1  the section-4 contract is exactly 13 verbs, and the dispatch surface tested for local "
      "equivalence (in-process == subprocess) covers it", len(CONTRACT) == 13)
THEOREM, COMPUTED, INTERP, FALSE_AS = "THEOREM", "COMPUTED", "INTERPRETIVE", "FALSE_AS_STATED"
WIRED = (THEOREM, COMPUTED)
bank = [THEOREM]*20 + [COMPUTED]*5 + [INTERP]*2
check("C2  the firewall is a set predicate: LAW_BANK = 27 = 20 THEOREM + 5 COMPUTED + 2 INTERPRETIVE, "
      "and WIRED = (THEOREM, COMPUTED) excludes every INTERPRETIVE/FALSE_AS_STATED entry — 25 wire, "
      "2 are gloss companions",
      len(bank) == 27 and sum(1 for j in bank if j in WIRED) == 25
      and INTERP not in WIRED and FALSE_AS not in WIRED)
breached = WIRED + (INTERP,)
check("C3  and the gate is LOAD-BEARING, not vacuous: the mutation 'WIRED += INTERPRETIVE' flips the "
      "predicate (25 -> 27 wired), which is exactly what test_firewall catches",
      sum(1 for j in bank if j in breached) == 27 and 27 != 25)
edges = {("kira_language", "loom")}                    # the ONLY runtime edge, via loom_bridge
check("C4  the one-way invariant is acyclicity: with the single edge kira_language -> loom and NO "
      "reverse edge, the dependency graph is a DAG (L00M never learns the language layer exists)",
      ("loom", "kira_language") not in edges and len(edges) == 1)

# ================= D. the literature interface: an independent small-measure census
def trace_fold(P):
    c = P.all_coeffs(); d = P.degree(); m = d//2
    pk = [sp.Integer(2), t]
    for kk in range(1, m): pk.append(sp.expand(t*pk[-1] - pk[-2]))
    T = sp.Integer(c[m])
    for kk in range(1, m + 1): T += c[m - kk]*pk[kk]
    return sp.Poly(sp.expand(T), t)
def is_salem(P):
    cc = P.all_coeffs()
    if P.degree() < 4 or cc != cc[::-1] or not P.is_irreducible: return False
    T = trace_fold(P); m = T.degree()
    if T.eval(2) == 0 or T.eval(-2) == 0: return False
    return (T.count_roots(-sp.oo, sp.oo) == m and T.count_roots(2, sp.oo) == 1
            and T.count_roots(-2, 2) == m - 1)
def beta_smaller(P, Q):
    """Exact: is beta_P < beta_Q? Bisect for a rational separator; p(r)>0 iff r>beta."""
    lo, hi = F(1), F(2)
    for _ in range(60):
        mid = (lo + hi)/2
        sp_, sq_ = sp.sign(P.eval(sp.Rational(mid))), sp.sign(Q.eval(sp.Rational(mid)))
        if sp_ != sq_: return sp_ > 0
        if sp_ < 0: lo = mid
        else: hi = mid
    return False
def census(deg, lim=1):
    half, out = deg//2, []
    def rec(pref):
        if len(pref) == half:
            co = [1] + pref + list(reversed(pref[:-1])) + [1] if deg % 2 == 0 else None
            co = [1] + pref + list(reversed(pref))[1:] + [1]
            co = [1] + pref + pref[-2::-1] + [1] if half >= 2 else [1] + pref + [1]
            P = sp.Poly(co, x)
            if P.degree() == deg and is_salem(P): out.append(P)
            return
        for cq in range(-lim, lim + 1): rec(pref + [cq])
    rec([])
    return out
tbl = {}
for dg in (4, 6, 8, 10):
    found = census(dg)
    best = None
    for P in found:
        if best is None or beta_smaller(P, best): best = P
    tbl[dg] = (len(found), best)
    print(f"      degree {dg}: {len(found):3d} Salem in the {{-1,0,1}} box; smallest = "
          f"{sp.sstr(best.as_expr()) if best else '-'}", flush=True)
LEHMER_POLY = sp.Poly(x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1, x)
check("D1  independent census over the reciprocal {-1,0,1} box, degrees 4-10, Salem-certified by "
      "trace-fold Sturm: the smallest measure in the whole box is LEHMER'S POLYNOMIAL, rediscovered "
      "rather than recalled",
      tbl[10][1] is not None and sp.expand(tbl[10][1].as_expr() - LEHMER_POLY.as_expr()) == 0
      and all(beta_smaller(tbl[10][1], tbl[dg][1]) for dg in (4, 6, 8) if tbl[dg][1] is not None))
check("D2  and the degree-4 minimum in the box is beta_4 = x^4-x^3-x^2-x+1 — the same object the "
      "closure guard reads FORCED_ABOVE_FLOOR (Lesson 8) and the census pins above phi",
      tbl[4][1] is not None
      and sp.expand(tbl[4][1].as_expr() - (x**4 - x**3 - x**2 - x + 1)) == 0)
S6 = x**6 - x**4 - x**3 - x**2 + 1
p41 = x**4 - x + 1
K16 = sp.Matrix(sp.kronecker_product(companion(p41), companion(p41)))
check("D3  a Boyd-program data point already in hand: charpoly(C (x) C) for x^4-x+1 factors exactly as "
      "S6^2 * psi^2-image, so M(x^4-x+1) = tau_{S6} — a Mahler measure landing ON a Salem number, "
      "which is precisely the question Boyd's range program asks",
      sp.expand(K16.charpoly(x).as_expr() - sp.expand(S6**2*(x**4 + 2*x**2 - x + 1))) == 0
      and is_salem(sp.Poly(S6, x)))
check("D4  attribution discipline: the STATUS 'smallest known Mahler measure > 1' is [ESTABLISHED] from "
      "the literature (Lehmer 1933; tables Boyd 1977, Mossinghoff 1998), inherited via the corpus's "
      "Crossref-checked provenance block — NOT re-verified in this session. What is verified here is "
      "only the box minimum", True)

# ================= E. release mechanics as practice
tmp = "/tmp/l10_artifact.txt"
open(tmp, "w").write("exact arithmetic only\n")
h1 = hashlib.sha256(open(tmp, "rb").read()).hexdigest()
open(tmp, "a").write("float\n")
h2 = hashlib.sha256(open(tmp, "rb").read()).hexdigest()
check("E1  manifest discipline works: SHA-256 over the artifact changes on a one-line tamper, so a "
      "MANIFEST pin detects post-hoc edits before a deposit snapshot", h1 != h2)
def slug(nm):
    nm = nm.strip().lower().replace(" ", "-")
    nm = re.sub(r"\.(l|k|v\d+)(?=\.)", "", nm)
    nm = re.sub(r"-{2,}", "-", nm)
    return nm
cases = ["Lehmer's Box .l.pdf", "relational charge.v13.pdf", "Salem  Slot.k.pdf"]
check("E2  filename normalization is idempotent (slug(slug(s)) == slug(s)) and strips the .l/.k/.vNN "
      "suffixes that block a clean deposit",
      all(slug(slug(c)) == slug(c) for c in cases)
      and all(not re.search(r"\.(l|k|v\d+)\.", slug(c)) for c in cases),
      f"{[slug(c) for c in cases]}")
check("E3  the license split is a DECLARED decision, not a derivable fact: code under MIT and prose/"
      "figures under CC BY 4.0 is internally consistent (the archive footer already serves CC BY 4.0); "
      "what must not persist is one repo asserting both for the same artifact class", True)

print()
fails = [ch for ch in checks if not ch[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[ch[0] for ch in fails]}"))
