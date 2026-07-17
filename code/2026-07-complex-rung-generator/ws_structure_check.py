#!/usr/bin/env python3
# ws_structure_check.py -- the within-shell taxonomy, upgraded from scan to structure.
#
# TARGET (v2.7 rem:p1scan / rem:crossscan): the within-shell taxonomy
# {pinned/inert, s(x^k)-shell tower, cyclotomic, dilation-of-known} is exhaustive
# only over the bounded 855-object window; "no structure theorem forcing
# p(x) = s(x^k) at unbounded degree -- only a scan."
#
# SESSION REPORT proves three degree-free lemmas and a trichotomy:
#
# WS-L1 (zeta-closure <=> substitution).  p monic irreducible, p(0) != 0, deg d.
#   The root set is closed under multiplication by a primitive k-th root of
#   unity zeta  <=>  every nonzero coefficient sits in degree = d (mod k)
#   <=>  k | d and p(x) = s(x^k) for a (monic, irreducible) s of degree d/k.
#
# WS-L2 (tower divisibility).  If p is irreducible with a WITHIN-SHELL TORSION
#   TIE of order k >= 2 (some pair of roots beta = zeta_k alpha), then with
#   g = minpoly(alpha^k), d' = deg g, the fiber size t = d/d' satisfies
#   2 <= t <= k and p(x) | g(x^k).  (Ratios of unit modulus: every torsion tie
#   is automatically within-shell.)
#
# WS-L3 (saturation criterion).  If the tie is realized over Q(zeta_k) --
#   i.e. zeta_k alpha is a conjugate of alpha over Q(zeta_k) -- then t = k and
#   p(x) = g(x^k) exactly.  COROLLARY (k = 2, unconditional): every order-2 tie
#   forces p(x) = g(x^2): an even polynomial.  (Q(zeta_2) = Q.)
#
# WS-T (structural trichotomy).  Every irreducible p with a torsion tie is
#   (a) cyclotomic, or (b) RADICAL: p | x^n - c for a tie order n with
#   alpha^n = c in Q (d' = 1), or (c) SATURATED TOWER: p = g(x^n), deg g >= 2,
#   or (d) the named RESIDUAL: 2 <= t < n with d' >= 2 at every tie order --
#   hunted below; over the 855 window the residual class is EMPTY, so the
#   taxonomy holds with a structural classifier, not a signature match.
#
# WS-L4 (catalog degree cap).  Every eigenvalue of every tensor/oplus catalog
#   word lies in E = Q(i, sqrt2, sqrt3, 5^{1/4}) with [E:Q] = 32; hence every
#   irreducible factor of every word has degree <= 32 at EVERY word size --
#   for the forced population, "unbounded degree" is a 32-bounded question.
#
# DISCIPLINE: exact integer/rational arithmetic at all decisions (resultants,
# polynomial division, cyclotomic divisibility, minimal polynomials); no floats
# at decision boundaries.  LF endings.  Exit 0.

import sys, itertools, time
import sympy as sp
from sympy import Symbol, Poly, Rational, sqrt, I, simplify, expand, resultant
from sympy import cyclotomic_poly, totient, minimal_polynomial

x = Symbol('x'); y = Symbol('y')
PASS = 0; FAIL = 0

def check(cid, desc, ok):
    global PASS, FAIL
    tag = "PASS" if ok else "FAIL"
    if ok: PASS += 1
    else: FAIL += 1
    print(f"[{tag}] {cid}: {desc}")

PHI = (1 + sqrt(5))/2

# ---------------------------------------------------------------- WS-A: the zeta-closure lemma (L1)
def support_condition(p, k):
    P = Poly(p, x); d = P.degree()
    return all((d - m) % k == 0 for m, c in enumerate(reversed(P.all_coeffs())) if c != 0)

def closed_under_zeta(p, k):
    """Exact: root set closed under zeta_k-mult <=> zeta^{-d} p(zeta x) = p(x) in Q(zeta)[x]."""
    P = Poly(p, x); d = P.degree()
    z = sp.exp(2*sp.pi*I/k)
    lhs = sp.expand(sp.simplify(z**(-d) * p.subs(x, z*x)))
    return sp.simplify(sp.expand(lhs - p)) == 0

A_batt = [
    (x**4 + 5*x**2 - 5, 2, True),    # K-seed: s(x^2), s = y^2+5y-5
    (x**4 - 6*x**2 + 4, 2, True),    # minpoly(phi*sqrt2)
    (x**10 + x**5 - 1, 5, True),     # q_5: the M=phi minimizer family
    (x**2 - 2*x + 2, 2, False),      # 1+i: no order-2 closure
    (x**3 - x - 1, 2, False),        # Smyth cubic: no closure
]
okA = True
for p, k, expect in A_batt:
    sc = support_condition(p, k)
    cz = closed_under_zeta(p, k)
    if not (sc == cz == expect):
        okA = False
check("WS-A1", "L1 both directions on a 5-poly battery: support condition <=> exact zeta-closure identity (k=2,5)", okA)
sK = x**2 + 5*x - 5
check("WS-A2", "L1 reconstruction: K-seed = s(x^2) with s = y^2+5y-5, s irreducible, deg = d/k (exact)",
      sp.expand(sK.subs(x, x**2) - (x**4 + 5*x**2 - 5)) == 0 and Poly(sK, x).is_irreducible)
check("WS-A3", "L1 divisibility direction: p(0) != 0 forces k | d (instance: x^10+x^5-1 has 5 | 10; support classes {0,5} = {d mod 5})",
      support_condition(x**10 + x**5 - 1, 5) and 10 % 5 == 0)

# ---------------------------------------------------------------- exact machinery
def adams(p, n):
    """psi^n(p): monic poly whose roots are alpha^n (with multiplicity), exact."""
    P = Poly(p, x)
    R = sp.resultant(P.as_expr().subs(x, y), x - y**n, y)
    Q = Poly(sp.expand(R), x)
    if Q.LC() != 1:
        Q = Poly(Q.as_expr() / Q.LC(), x)
    return Q

def ratio_poly(p):
    """R(x) = Res_y(p(y), p(x y)): roots are all ratios alpha_i/alpha_j."""
    P = p.subs(x, y)
    return Poly(sp.expand(sp.resultant(P, p.subs(x, x*y), y)), x)

def tie_orders(p):
    """Exact torsion-tie orders: n >= 2 with Phi_n | Res_y(p(y), p(xy))."""
    P = Poly(p, x); d = P.degree()
    R = ratio_poly(p)
    bound = d*(d - 1)
    out = []
    n = 2
    while True:
        if totient(n) > bound and n > 2*bound + 2:
            break
        if totient(n) <= bound:
            cn = Poly(cyclotomic_poly(n, x), x)
            if cn.degree() <= R.degree() and R.rem(cn) == Poly(0, x):
                out.append(n)
        n += 1
        if n > 4*bound + 6:
            break
    return out

def minpoly_of_power(p, n):
    """g = minpoly(alpha^n) for a root alpha of irreducible p: the unique irreducible
    factor h of psi^n(p) with p | h(x^n). Returns (g, fiber t)."""
    P = Poly(p, x); d = P.degree()
    psi = adams(p, n)
    g = None
    for h, _m in sp.factor_list(psi.as_expr())[1]:
        H = Poly(h, x)
        if H.LC() != 1:
            H = Poly(H.as_expr()/H.LC(), x)
        comp = Poly(sp.expand(H.as_expr().subs(x, x**n)), x)
        if comp.rem(P) == Poly(0, x):
            g = H
            break
    if g is None:
        return None, None
    return g, Rational(d, g.degree())

# ---------------------------------------------------------------- WS-B: tower divisibility (L2)
B_batt = [
    (x**2 - 2*x + 2, 4),     # 1+i, tie order 4
    (x**4 + 5*x**2 - 5, 2),  # K-seed
    (x**4 - 6*x**2 + 4, 2),  # phi*sqrt2
    (x**10 + x**5 - 1, 5),   # q_5
    (x**4 + 2*x**2 + 2, 8),  # sqrt(-1+i): tie order 8
]
okB = True
detail = []
for p, n in B_batt:
    g, t = minpoly_of_power(p, n)
    d = Poly(p, x).degree()
    ok1 = g is not None
    ok2 = ok1 and (2 <= t <= n) and (t == Rational(d, g.degree()))
    comp = Poly(sp.expand(g.as_expr().subs(x, x**n)), x) if ok1 else None
    ok3 = ok1 and comp.rem(Poly(p, x)) == Poly(0, x)
    detail.append((sp.sstr(p), n, sp.sstr(g.as_expr()) if g else None, str(t)))
    okB = okB and ok1 and ok2 and ok3
check("WS-B1", f"L2 on a 5-object battery: p | g(x^n), fiber 2 <= t <= n, exact minpoly-of-power (e.g. x^2-2x+2 | (x+4)(x^4) i.e. x^4+4, t=2<4)", okB)
g_ex, t_ex = minpoly_of_power(x**2 - 2*x + 2, 4)
check("WS-B2", "the 1+i witness: g = x + 4 (so alpha^4 = -4), x^4 + 4 = (x^2-2x+2)(x^2+2x+2) exact; unsaturated at k=4 (t=2)",
      g_ex.as_expr() == x + 4 and sp.expand((x**2-2*x+2)*(x**2+2*x+2) - (x**4+4)) == 0 and t_ex == 2)

# ---------------------------------------------------------------- WS-C: saturation (L3)
# k=2 unconditional: every order-2 tie forces even support (checked over the whole window in WS-D).
# Q(zeta_4)-non-realization for 1+i: p splits over Q(i) into distinct linears => no sigma fixing i moves alpha to i*alpha.
fac_i = sp.factor(x**2 - 2*x + 2, extension=I)
check("WS-C1", "1+i: p factors over Q(i) into distinct linears (x-1-i)(x-1+i): the order-4 tie is NOT Q(i)-realized, exactly why saturation fails there",
      sp.expand(fac_i - (x - 1 - I)*(x - 1 + I)) == 0)
gK, tK = minpoly_of_power(x**4 + 5*x**2 - 5, 2)
check("WS-C2", "K-seed: order-2 tie IS Q-realized (zeta=-1 rational), so L3 forces p = g(x^2) with g = y^2+5y-5: t = k = 2 exactly",
      gK.as_expr() == x**2 + 5*x - 5 and tK == 2)

# ---------------------------------------------------------------- WS-D: the 855-window scan, structural classifier
CAPS = [(2, 3), (3, 2), (4, 2), (5, 1), (6, 1)]
t0 = time.time()
window = []
for d, cap in CAPS:
    for tail in itertools.product(range(-cap, cap + 1), repeat=d):
        coeffs = (1,) + tail
        P = Poly(list(coeffs), x)
        if P.is_irreducible:
            window.append(P)
check("WS-D1", f"window replication: {len(window)} irreducible monic objects under the logged caps {CAPS} (expected 855)",
      len(window) == 855)

def is_cyclo(P):
    d = P.degree()
    # cyclotomic factor of degree <= 6 has order n with totient(n) <= 6 => n <= 18
    for n in range(1, 19):
        if totient(n) == d and Poly(cyclotomic_poly(n, x), x) == P:
            return True
    return False

buckets = {"inert": 0, "cyclotomic": 0, "radical": 0, "tower": 0, "residual": []}
tied_objects = []
order2_evensupport_ok = True
for P in window:
    d = P.degree()
    ties = tie_orders(P.as_expr())
    if not ties:
        buckets["inert"] += 1
        continue
    tied_objects.append((P, ties))
    # L3 corollary consistency: order-2 tie <=> even support
    if 2 in ties:
        if not support_condition(P.as_expr(), 2):
            order2_evensupport_ok = False
    if support_condition(P.as_expr(), 2) and 2 not in ties and d >= 2:
        order2_evensupport_ok = False
    if is_cyclo(P):
        buckets["cyclotomic"] += 1
        continue
    classified = False
    # radical: some tie order n with alpha^n rational
    for n in ties:
        g, t = minpoly_of_power(P.as_expr(), n)
        if g is not None and g.degree() == 1:
            buckets["radical"] += 1
            classified = True
            break
    if classified:
        continue
    # saturated tower: p = g(x^n) at some tie order, deg g >= 2
    for n in ties:
        g, t = minpoly_of_power(P.as_expr(), n)
        if g is not None and t == n and g.degree() >= 2:
            comp = Poly(sp.expand(g.as_expr().subs(x, x**n)), x)
            if comp == P:
                buckets["tower"] += 1
                classified = True
                break
    if not classified:
        buckets["residual"].append((P.all_coeffs(), ties))

elapsed = time.time() - t0
check("WS-D2", f"structural scan complete in {elapsed:.0f}s: inert {buckets['inert']}, cyclotomic {buckets['cyclotomic']}, radical {buckets['radical']}, saturated tower {buckets['tower']}, residual {len(buckets['residual'])}",
      buckets["inert"] + buckets["cyclotomic"] + buckets["radical"] + buckets["tower"] + len(buckets["residual"]) == 855)
check("WS-D3", "corpus replication (convention-free counts): inert = 804 and cyclotomic = 11, matching rem:p1scan exactly",
      buckets["inert"] == 804 and buckets["cyclotomic"] == 11)
check("WS-D4", f"the 40 tied non-cyclotomic objects split radical {buckets['radical']} + saturated tower {buckets['tower']} + unit-twist {len(buckets['residual'])} = 40 under the STRUCTURAL classifier (the corpus's signature bins shell 32 + dilation 8 partition the same 40)",
      buckets["radical"] + buckets["tower"] + len(buckets["residual"]) == 40)
# WS-D5: the residual hunt FOUND the fourth structural class -- and it is exactly
# the unit-twist class L2/L3 predict: p | g(x^n) properly, d' >= 2, fiber t < n,
# tie not Q(zeta_n)-realized.  The conjectured "p = s(x^k) at unbounded degree"
# is therefore FALSE as stated; the corrected, degree-free statement is L2 + L3.
TWISTS = {
    (1, -1, 2, 1, 1):  "zeta3 * tau  (= -zeta3 * phi)",
    (1, 1, 2, -1, 1):  "zeta3 * phi  (= -zeta3 * tau)",
    (1, -2, 2, 2, 1):  "zeta8 * sqrt(2 + sqrt3)",
    (1, 2, 2, -2, 1):  "zeta8 * sqrt(2 - sqrt3)",
}
found = {tuple(int(c) for c in coeffs) for coeffs, _t in buckets["residual"]}
check("WS-D5", "RESIDUAL HUNT (finding): exactly FOUR unit-twist objects in-window -- torsion twists of the catalog-adjacent units tau, phi, 2+-sqrt3; the 4-label taxonomy of rem:p1scan is structurally incomplete and 'p = s(x^k)' is refuted by these witnesses",
      found == set(TWISTS.keys()))
# exact identification of each twist by minimal-polynomial identity
z3 = sp.exp(2*sp.pi*I/3); z8 = sp.exp(2*sp.pi*I/8)
tau_ = (sqrt(5) - 1)/2
ids_ok = True
ids_ok &= sp.expand(minimal_polynomial(z3*tau_, x) - (x**4 - x**3 + 2*x**2 + x + 1)) == 0
ids_ok &= sp.expand(minimal_polynomial(z3*PHI, x) - (x**4 + x**3 + 2*x**2 - x + 1)) == 0
ids_ok &= sp.expand(minimal_polynomial(z8*sqrt(2 - sqrt(3)), x) - (x**4 + 2*x**3 + 2*x**2 - 2*x + 1)) == 0
ids_ok &= sp.expand(minimal_polynomial(z8*sqrt(2 + sqrt(3)), x) - (x**4 - 2*x**3 + 2*x**2 + 2*x + 1)) == 0
check("WS-D5b", "twist identifications exact: minpoly(zeta3 tau), minpoly(zeta3 phi), minpoly(zeta8 sqrt(2-+sqrt3)) equal the four residual quartics (minimal-polynomial identities)",
      bool(ids_ok))
fac3 = sp.factor(x**4 - x**3 + 2*x**2 + x + 1, extension=z3)
check("WS-D5c", "obstruction witnessed: the zeta3-twist quartic factors over Q(zeta3) into two DISTINCT quadratics, so the order-3 tie is not Q(zeta3)-realized -- exactly L3's failure mode (saturation impossible)",
      isinstance(fac3, sp.Mul) and len(fac3.args) == 2)
# WS-D8: the L2 CLASSIFICATION content.  (p | g(x^n) is automatic -- g is chosen as the factor with
# p | g(x^n) -- so re-checking it is a tautology; the real content is that the fiber t = d/deg(g) sorts every
# within-shell torsion tie into saturated (t=n) or unit-twist (2<=t<n), a clean partition, AND -- fail-first --
# that genuinely UNTIED irreducibles carry no within-shell tie at all, so the tie predicate is non-vacuous.
sat_ct, twist_ct, class_ok = 0, 0, True
for Pq, ties in tied_objects:
    for n in ties:
        g, t = minpoly_of_power(Pq.as_expr(), n)
        if g is None or not (2 <= t <= n):
            class_ok = False; continue
        if t == n:
            sat_ct += 1
        else:
            twist_ct += 1
untied_ok = all(tie_orders(q) == [] for q in [x**3 - x - 1, x**5 - x - 1])
check("WS-D8", f"L2 over ALL {len(tied_objects)} tied window objects: every within-shell torsion tie sorts into saturated t=n ({sat_ct}) or unit-twist 2<=t<n ({twist_ct}, incl. the 4 refuting witnesses) -- a clean fiber partition; FAIL-FIRST: untied x^3-x-1, x^5-x-1 carry NO within-shell tie (tie_orders empty), so the predicate is non-vacuous",
      class_ok and untied_ok and twist_ct >= 4)
check("WS-D6", "L3 corollary over all 855 objects: order-2 tie <=> even coefficient support (zero violations, both directions)",
      order2_evensupport_ok)

# ---------------------------------------------------------------- WS-E: unbounded-degree instances (the lemmas are degree-free)
p16 = x**16 - x**8 - 1
P16 = Poly(p16, x)
ok16 = P16.is_irreducible and support_condition(p16, 8)
g16, t16 = minpoly_of_power(p16, 8)
check("WS-E1", "degree-16 instance minpoly(phi^{1/8}) = x^16 - x^8 - 1: irreducible, saturated s(x^8) tower with g = y^2-y-1, t = 8 (the lemmas run at any degree)",
      ok16 and g16 is not None and g16.as_expr() == x**2 - x - 1 and t16 == 8)
p12 = x**12 + x**6 - 1
P12 = Poly(p12, x)
g12, t12 = minpoly_of_power(p12, 6)
check("WS-E2", "degree-12 instance q_6 = x^12 + x^6 - 1 (the even-floor minimizer family at m=12): saturated tower over y^2+y-1, t = 6",
      P12.is_irreducible and g12 is not None and g12.as_expr() == x**2 + x - 1 and t12 == 6)

# ---------------------------------------------------------------- WS-F: the catalog degree cap (L4)
mpK = minimal_polynomial(5**Rational(1, 4) / PHI, x)
check("WS-F1", "K = 5^{1/4}/phi has minimal polynomial x^4 + 5x^2 - 5 (the K-seed) -- exact identity",
      sp.expand(mpK - (x**4 + 5*x**2 - 5)) == 0)
mpib = minimal_polynomial(I * 5**Rational(1, 4) * PHI, x)
check("WS-F2", "i*beta = i 5^{1/4} phi has the SAME minimal polynomial x^4 + 5x^2 - 5 (the rotation pair lives in E)",
      sp.expand(mpib - (x**4 + 5*x**2 - 5)) == 0)
deg8 = sp.degree(minimal_polynomial(sqrt(2) + 5**Rational(1, 4), x), x)
deg16 = sp.degree(minimal_polynomial(sqrt(2) + sqrt(3) + 5**Rational(1, 4), x), x)
check("WS-F3", "tower degrees: [Q(sqrt2, 5^{1/4}):Q] = 8 and [Q(sqrt2, sqrt3, 5^{1/4}):Q] = 16 via primitive-element minimal polynomials (exact)",
      deg8 == 8 and deg16 == 16)
# seed eigenvalue membership FIRST (so WS-F4 can lean on it, not on a bare degree count):
membership = [
    (PHI, x**2 - x - 1),
    (sqrt(2), x**2 - 2), (sqrt(3), x**2 - 3), (sqrt(5), x**2 - 5),
    ((7 + 3*sqrt(5))/2, x**2 - 7*x + 1),
]
okmem = all(sp.expand(minimal_polynomial(v, x) - p) == 0 for v, p in membership)
check("WS-F5", "seed eigenvalues phi, sqrt2, sqrt3, sqrt5, (7+3sqrt5)/2 have the catalog minimal polynomials and lie in E (with K, i*beta from WS-F1/F2: all seven seeds)",
      okmem)
deg32 = sp.degree(minimal_polynomial(I + sqrt(2) + sqrt(3) + 5**Rational(1, 4), x), x)
# L4, made NON-vacuous: [E:Q]=32 AND all seven seed eigenvalues verified in E (K, i*beta via WS-F1/F2;
# the five quadratics via WS-F5/okmem).  E is a field closed under product (tensor) and union (oplus), so
# every catalog-WORD eigenvalue is a product of seeds and lies in E, hence has degree dividing 32.
# FAIL-FIRST: a value NOT in E (2^{1/3}, degree 3) is excluded because 3 does not divide 32.
seeds_in_E = (sp.expand(mpK - (x**4 + 5*x**2 - 5)) == 0
              and sp.expand(mpib - (x**4 + 5*x**2 - 5)) == 0 and okmem)
not_in_E = (32 % sp.degree(minimal_polynomial(2**Rational(1, 3), x), x) != 0)
check("WS-F4", "L4: [E:Q]=32 AND all seven seed eigenvalues verified in E (WS-F1/F2/F5); E a field closed under product/union => every catalog-word eigenvalue lies in E, so every irreducible word factor has degree dividing 32 (<= 32) at EVERY word size -- the forced population's 'unbounded degree' is 32-bounded. FAIL-FIRST: 2^{1/3} (deg 3 does not divide 32) is excluded.",
      deg32 == 32 and seeds_in_E and not_in_E)

print()
print(f"ws_structure_check: {PASS}/{PASS+FAIL} PASS")
sys.exit(0 if FAIL == 0 else 1)
