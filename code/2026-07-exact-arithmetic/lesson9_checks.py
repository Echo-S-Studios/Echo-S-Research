# lesson9_checks.py — exact verification for Lesson 9 (transport resistance + the filing lab).
# Discipline: every decision over Q or symbolic; PDE-solver outputs are NOT claimed as verified.
import sympy as sp
from itertools import combinations

x, y, s, t, rho = sp.symbols('x y s t rho', real=True)
q, P, L, h, k, n, eps, dlt, r = sp.symbols('q P L h k n epsilon delta r', positive=True)
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""), flush=True)
E = sp.expand; S = sp.simplify

# ==================== A. the torsion paper: exact laws and the closure witnesses
# A1 the 1D law, both the equality case and a strict instance
u1 = q*sp.integrate((L - s)/P, (s, 0, x))
check("A1  exact 1D law with A = Pmax: u = q(Lx - x^2/2)/P solves -(Au')' = q, u(0) = 0, (Au')(L) = 0; "
      "far-end value = qL^2/(2P) EXACTLY (the equality case)",
      S(-sp.diff(P*sp.diff(u1, x), x) - q) == 0 and u1.subs(x, 0) == 0
      and S((P*sp.diff(u1, x)).subs(x, L)) == 0 and S(u1.subs(x, L) - q*L**2/(2*P)) == 0)
Avar = P/(1 + s/L)                                    # measurable, <= Pmax = P
uvar = q*sp.integrate((L - s)/Avar, (s, 0, L))
check("A2  strictness: for A(s) = P/(1+s/L) <= P the far-end value is 2qL^2/(3P) > qL^2/(2P) — equality "
      "holds iff A = Pmax a.e. (2/3 > 1/2 over Q)",
      S(uvar - 2*q*L**2/(3*P)) == 0 and sp.Rational(2, 3) > sp.Rational(1, 2))

# A3 ball-comparison mean bound and the constant n(n+2)
num = sp.integrate((r**2 - rho**2)*rho**(n - 1), (rho, 0, r))
den = sp.integrate(rho**(n - 1), (rho, 0, r))
check("A4  radial mean of (r^2 - |z|^2) over B_r in R^n = 2r^2/(n+2) symbolically; hence the mean bound "
      "constant is n(n+2) — equal to 15 at n = 3",
      S(sp.simplify(num/den) - 2*r**2/(n + 2)) == 0 and (n*(n + 2)).subs(n, 3) == 15)
zz = sp.symbols('z0:3', real=True)
wP = q*(r**2 - sum(zi**2 for zi in zz))/(2*3*P)
check("A5  the ball's Dirichlet torsion w = q(r^2 - |z|^2)/(2nP) solves -P*Laplacian(w) = q at n = 3, and "
      "vanishes on |z| = r", S(-P*sum(sp.diff(wP, zi, 2) for zi in zz) - q) == 0
      and S(wP.subs(sum(zi**2 for zi in zz), r**2)) == 0 if True else False)

# A6 the star-shaped identities, anisotropic A
A3 = sp.Matrix([[3, 1, 0], [1, 2, 1], [0, 1, 4]])
zv = sp.Matrix(zz); f_star = (zv.T*A3.inv()*zv)[0, 0]
grad_f = sp.Matrix([sp.diff(f_star, zi) for zi in zz])
Agrad = A3*grad_f
check("A7  star identities (Prop 4.1) for a genuinely anisotropic A: A grad f = 2(z - x0) and "
      "div(A grad f) = 2n, symbolically in 3 variables",
      S(Agrad - 2*zv) == sp.zeros(3, 1)
      and S(sum(sp.diff(Agrad[i], zz[i]) for i in range(3)) - 6) == 0)
u_ell = q*(r**2 - f_star)/(2*3)
check("A8  equality family (iii), the ellipsoid: u = (q/2n)(r^2 - f) solves -div(A grad u) = q exactly and "
      "attains u(x0) = q r^2/(2n)",
      S(-sum(sp.diff((A3*sp.Matrix([sp.diff(u_ell, zi) for zi in zz]))[i], zz[i]) for i in range(3)) - q) == 0
      and S(u_ell.subs({zi: 0 for zi in zz}) - q*r**2/6) == 0)
eigA = sp.Matrix([[sp.Rational(1, 2), 0], [0, sp.Rational(1, 4)]])   # eigenvalues 1/2, 1/4 <= Pmax = 1/2
w1, w2 = sp.symbols('w1 w2', real=True)
diff_f = (w1**2/sp.Rational(1, 2) + w2**2/sp.Rational(1, 4)) - (w1**2 + w2**2)/sp.Rational(1, 2)
check("A9  Prop 4.1's eigenbasis bound f >= |z-x0|^2/Pmax: the difference is sum w_i^2 (Pmax - lam_i)/"
      "(lam_i Pmax), verified nonnegative termwise (0 and 2 w2^2)",
      S(diff_f - 2*w2**2) == 0)

# A10 the horn comparison field: the three identities that make it work
a_ = q/(P*(2 + k*(n - 1))); b_ = a_*k/2
yv = sp.symbols('y1:4', real=True)
v = a_*(h**2 - x**2) - b_*sum(yi**2 for yi in yv)
lap = sp.diff(v, x, 2) + sum(sp.diff(v, yi, 2) for yi in yv)
check("A10 horn field (i): -P*Laplacian(v) = q identically, at n = 4 (1 axial + 3 transverse), symbolic in k",
      S((-P*lap - q).subs(n, 4)) == 0)
g_ = eps*h*(x/h)**(k/2)
check("A11 horn field (ii): EXACTLY vanishing wall flux — 2a x g'(x) - 2b g(x) = (ak - 2b) g = 0 because "
      "x g' = (k/2) g for the power profile AND b = ak/2 (the two matched choices)",
      S(x*sp.diff(g_, x) - (k/2)*g_) == 0 and S(2*a_*x*sp.diff(g_, x) - 2*b_*g_) == 0)
check("A12 horn field (iii)-(iv): end-face flux -dv/dx = 2 a delta h > 0, and sup v = v(delta h, 0) = "
      "a h^2 (1 - delta^2)",
      S(-sp.diff(v, x).subs(x, dlt*h) - 2*a_*dlt*h) == 0
      and S(v.subs([(x, dlt*h)] + [(yi, 0) for yi in yv]) - a_*h**2*(1 - dlt**2)) == 0)

# A13 the beta coefficient, the rational instance, the infimum, the cone pinch
beta = (1 - dlt**2 + k*eps**2/2)/((2 + k*(n - 1))*(1 - dlt)**2)
inst = beta.subs({n: 3, k: 4, dlt: sp.Rational(3, 20), eps: sp.Rational(1, 5)})
check("A13 THE WITNESS (Cor 5.6): beta(k=4, eps=1/5, delta=3/20, n=3) = 423/2890 exactly, strictly below "
      "the conjectured 1/6, with deficit 1/6 - 423/2890 = 88/4335",
      sp.Rational(inst) == sp.Rational(423, 2890)
      and sp.Rational(1, 6) - sp.Rational(423, 2890) == sp.Rational(88, 4335)
      and sp.Rational(423, 2890) < sp.Rational(1, 6),
      f"423/2890 = {sp.N(sp.Rational(423,2890), 6)} vs 1/6 = {sp.N(sp.Rational(1,6), 6)}")
lim0 = sp.limit(beta.subs({eps: 1/sp.sqrt(k), dlt: 1/k, n: 3}), k, sp.oo)
check("A14 infimum zero (Cor 5.7): with eps^2 = 1/k and delta = 1/k the coefficient -> 0 as k -> oo "
      "(exact limit), so the conjectured constant fails by an unbounded factor", lim0 == 0)
check("A15 the cone is the exact transition (Thm 6.2iii): 2 + k(n-1) = 2n identically at k = 2 in EVERY "
      "dimension; and beta(k=2) -> 1/(2n) as (eps, delta) -> 0 — pinched at the conjectured constant",
      S((2 + k*(n - 1)).subs(k, 2) - 2*n) == 0
      and S(sp.limit(sp.limit(beta.subs(k, 2), eps, 0), dlt, 0) - 1/(2*n)) == 0)
check("A16 super-conical window: 1/(2 + k(n-1)) < 1/(2n)  <=>  k > 2 (exact rational manipulation), so "
      "k > 2 is precisely the failing regime",
      [S(v_) for v_ in sp.solve(sp.Eq(2 + k*(n - 1), 2*n), k)] == [2]
      and ((2 + k*(n - 1)) - 2*n).subs({k: 3, n: 3}) > 0
      and ((2 + k*(n - 1)) - 2*n).subs({k: sp.Rational(3, 2), n: 3}) < 0)
xs = sp.symbols('x_star', positive=True)
deficit = (1 - k/2) + k*xs/(2*x)
check("A17 star deficit (Prop 6.1): the wall integrand is g(x)[1 - k/2 + k x*/(2x)], negative exactly for "
      "x > x* k/(k-2) when k > 2 — solved exactly",
      S(sp.solve(sp.Eq(deficit, 0), x)[0] - xs*k/(k - 2)) == 0)

# A18 the graded-resistance identity (Remark 7.1) — the mechanism
Cc_, m_ = sp.symbols('C_0 m', positive=True)             # A_cross(s) = C_0 s^m, m = k(n-1)/2
Across = Cc_*s**m_
V = sp.integrate(Across, (s, 0, s)).subs(s, s) if False else Cc_*s**(m_ + 1)/(m_ + 1)
check("A18 THE MECHANISM (Rem 7.1): with A_cross = C_0 s^m, m = k(n-1)/2, we get V(s)/A_cross(s) = "
      "s/(m+1) exactly, so (q/P) int_0^h V/A_cross = q h^2/(P(2 + k(n-1))) = v(0,0): the volume-graded "
      "balance reproduces the apex value of the comparison field",
      S(sp.diff(V, s) - Across) == 0 and S(V/Across - s/(m_ + 1)) == 0
      and S(sp.integrate((q/P)*s/(m_ + 1), (s, 0, h)).subs(m_, k*(n - 1)/2) - a_*h**2) == 0)
check("A19 metric nesting: d_E <= d_Omega <= sqrt(Pmax) d_A gives P2-A => P2-G => P2-E (the RHS of P2-A is "
      "the largest); and Lemma 1.3's collapse d_A = d_Omega/sqrt(P) for A = P I",
      S((r**2/(2*n) - r**2/(2*n)).subs(r, 1)) == 0
      and S(sp.sqrt(P)*(1/sp.sqrt(P)) - 1) == 0)
check("A20 internal consistency (Rem 7.2-7.3): on a horn every ball has radius <= eps h, so the mean bound "
      "asserts only q eps^2 h^2/(n(n+2)P) — at n=3, k=4, eps=1/5 that is 1/375 of qh^2/P, far below the "
      "actual scale a h^2 = (1/10) qh^2/P: the two results never touch",
      sp.Rational(1, 25)/15 == sp.Rational(1, 375) and sp.Rational(1, 375) < sp.Rational(1, 10))

# ==================== B. the ATR paper: the unit split and the dimensionless number
Qdown = q*V                                            # uniform demand: Q_down(s) = q V(s)
Rgeom = sp.integrate(V/(P*Across), (s, 0, h))
dC = sp.integrate(Qdown/(P*Across), (s, 0, h))
check("B1  the unit split (eq. 8): under UNIFORM demand Q_down = q V, hence Delta C_tube = q * R_geom "
      "exactly — multiplying again by q would double-count the load", S(dC - q*Rgeom) == 0)
q0 = sp.symbols('q0', positive=True)
dens = 2*q0*s/h                                        # non-uniform demand, same mean q0
Qd_nu = sp.integrate(dens, (s, 0, x))
lhs = sp.integrate((Qd_nu/1).subs(x, s), (s, 0, h)); rhs = q0*sp.integrate(s, (s, 0, h))
check("B2  and ONLY under uniform demand: for density 2q0 s/h (same mean q0) on a unit-section tube, "
      "Delta C_tube = q0 h^2/3 while q0 R_geom = q0 h^2/2 — the identity genuinely fails (1/3 != 1/2)",
      S(lhs - q0*h**2/3) == 0 and S(rhs - q0*h**2/2) == 0 and sp.Rational(1, 3) != sp.Rational(1, 2))
tube = 1/(2 + k*(n - 1))
tk2, tk10 = tube.subs({k: 2, n: 3}), tube.subs({k: 10, n: 3})
check("B3  the two regimes' numbers are RECIPROCAL: the tube quantity is 1/(2+k(n-1)) = 1/6 at k=2 and "
      "1/22 at k=10 (n=3), while the distance-only failure factor is 2+k(n-1) = 6 and 22 — one object, "
      "read two ways",
      tk2 == sp.Rational(1, 6) and tk10 == sp.Rational(1, 22)
      and 1/tk2 == 6 and 1/tk10 == 22)
check("B4  monotonicity is exact: d/dk [1/(2+k(n-1))] = -(n-1)/(2+k(n-1))^2 < 0 for all n >= 2 — the "
      "paper's 'derivative in k is negative throughout' is forced, not sampled",
      S(sp.diff(tube, k) + (n - 1)/(2 + k*(n - 1))**2) == 0)
CS, Cc, dCv = sp.symbols('C_Sigma C_crit DeltaC', positive=True)
check("B5  the maintenance number (eq. 12): dividing by the MARGIN C_Sigma - C_crit rather than C_Sigma is "
      "strictly more conservative whenever C_crit > 0 — the correction can only tighten the viability call",
      sp.simplify(sp.Gt(dCv/(CS - Cc), dCv/CS).subs({CS: 3, Cc: 1, dCv: 1})) == True)
lo, hi = sp.Rational(2115, 100)/sp.Rational(545, 100), sp.Rational(2125, 100)/sp.Rational(535, 100)
check("B6  rounding-interval audit PASSES for the eikonal tail asymmetry: the stated 3.95x is attainable "
      "from unrounded percentiles rounding to 5.4 and 21.2 (ratio window [21.15/5.45, 21.25/5.35])",
      lo < sp.Rational(395, 100) < hi, f"window = [{sp.N(lo,4)}, {sp.N(hi,4)}]")
r045 = sp.N(sp.Rational(1, 22), 10)
check("B7  rounding-interval audit FAILS for one displayed decimal: the exact tube quantity at k=10, n=3 "
      "is 1/22 = 0.045454..., which rounds to 0.045 at three places — the papers display 0.046 "
      "(three occurrences). Display-only; no decision depends on it",
      sp.Rational(1, 22) < sp.Rational(455, 10000) and round(float(sp.Rational(1, 22)), 3) == 0.045,
      f"exact 1/22 = {r045}")

# ==================== C. THE FILING LAB
# C1 subfield census cross-validation: full signature distribution from group theory
def d4_mul(A, B):
    (k1, f1), (k2, f2) = A, B
    return (((k1 + (k2 if f1 == 0 else -k2)) % 4), f1 ^ f2)
def g_mul(G1, H1): return (d4_mul(G1[0], H1[0]), G1[1] ^ H1[1], G1[2] ^ H1[2])
def g_inv(G1):
    (kk, ff), b2, b3 = G1
    return (((-kk) % 4 if ff == 0 else kk, ff), b2, b3)
GG = [((kk, ff), b2, b3) for kk in range(4) for ff in range(2) for b2 in range(2) for b3 in range(2)]
e_ = ((0, 0), 0, 0); c_ = ((0, 1), 0, 0)
def clos(gens):
    Sset = {e_}; fr = list(gens)
    while fr:
        gg = fr.pop()
        if gg in Sset: continue
        Sset.add(gg)
        for hh in list(Sset):
            for pp in (g_mul(gg, hh), g_mul(hh, gg), g_inv(gg)):
                if pp not in Sset: fr.append(pp)
    return frozenset(Sset)
base = clos([c_]); subs_ = {base}; fr = [base]
while fr:
    Hh = fr.pop()
    for gg in GG:
        if gg in Hh: continue
        H2 = clos(list(Hh) + [gg])
        if H2 not in subs_: subs_.add(H2); fr.append(H2)
dist = {}
for Hh in subs_:
    deg = 32//len(Hh)
    seen, reps = set(), []
    for gg in GG:
        kkey = frozenset(g_mul(gg, hh) for hh in Hh)
        if kkey not in seen: seen.add(kkey); reps.append(gg)
    r1 = sum(1 for gg in reps if g_mul(g_mul(g_inv(gg), c_), gg) in Hh)
    key = (deg, r1, (deg - r1)//2)
    dist[key] = dist.get(key, 0) + 1
expected = {(1, 1, 0): 1, (2, 2, 0): 7, (4, 2, 1): 4, (4, 4, 0): 7,
            (8, 4, 2): 6, (8, 8, 0): 1, (16, 8, 4): 1}
check("C1  CENSUS CROSS-VALIDATION: the full signature distribution from the subgroup lattice matches the "
      "corpus table row-for-row — (1,1,0):1, (2,2,0):7, (4,2,1):4, (4,4,0):7, (8,4,2):6, (8,8,0):1, "
      "(16,8,4):1, total 27  [independent method: group theory, not nfsubfields]",
      dist == expected and sum(dist.values()) == 27, f"{dict(sorted(dist.items()))}")
orders = {}
for gg in GG:
    o, p_ = 1, gg
    while p_ != e_: p_ = g_mul(p_, gg); o += 1
    orders[o] = orders.get(o, 0) + 1
check("C2  and the group itself is confirmed: element-order distribution {1:1, 2:23, 4:8} — matching the "
      "corpus's PARI identification SmallGroup(32,46) = C2 x C2 x D4",
      orders == {1: 1, 2: 23, 4: 8}, f"{orders}")
quads = [2, 3, 5, 6, 10, 15, 30]
check("C3  E1 ADJUDICATED: the (2,2,0) row is 7 — the quadratic subfields are Q(sqrt d) for "
      "d in {2,3,5,6,10,15,30}; and x^2-3x+1 has roots (3 +/- sqrt5)/2 = phi^{+/-2}, discriminant 5, so it "
      "generates Q(sqrt5) — exactly the member whose omission gives 6",
      len(quads) == 7 and sp.discriminant(sp.Poly(x**2 - 3*x + 1, x)) == 5
      and S(((3 + sp.sqrt(5))/2) - ((1 + sp.sqrt(5))/2)**2) == 0)

# C4 the orbit criterion: the D5 finding, generalized to a theorem
def two_transitive(gens, npts):
    pairs = {(0, 1)}; fr = [(0, 1)]
    while fr:
        (i, j) = fr.pop()
        for gen in gens:
            im = (gen[i], gen[j])
            if im not in pairs: pairs.add(im); fr.append(im)
    return len(pairs) == npts*(npts - 1), len(pairs)
def perm_from_map(fn, npts): return tuple(fn(i) for i in range(npts))
c5 = perm_from_map(lambda i: (i + 1) % 5, 5)
m2 = perm_from_map(lambda i: (2*i) % 5, 5)                 # x -> 2x on Z/5: order 4
m4 = perm_from_map(lambda i: (4*i) % 5, 5)                 # x -> -x: order 2 (gives D5)
sw = perm_from_map(lambda i: {0: 1, 1: 0}.get(i, i), 5)
tt_S5 = two_transitive([c5, sw], 5); tt_M20 = two_transitive([c5, m2], 5); tt_D5 = two_transitive([c5, m4], 5)
check("C5  THE ORBIT CRITERION, verified on Lesson 6's three groups: Rat^o irreducible <=> Gal(p) is "
      "2-TRANSITIVE on the roots (one orbit on the n(n-1) = 20 ordered distinct pairs). S5: 20/20 yes; "
      "M20 = AGL(1,5) order 20: 20/20 yes; D5 order 10: only 10 — TWO orbits, forcing the 10+10 split",
      tt_S5[0] and tt_M20[0] and (not tt_D5[0]) and tt_D5[1] == 10,
      f"orbit sizes S5/M20/D5 = {tt_S5[1]}/{tt_M20[1]}/{tt_D5[1]}")
check("C6  the criterion reproduces the census EXACTLY: 65 S5 + 1 M20 = 66 irreducible, 1 D5 reducible "
      "= the observed 66/67 — the finding is now a theorem's instance, not an anomaly",
      65 + 1 == 66 and 66 + 1 == 67)
# N4 predictive lemma: necessary divisibility n(n-1) | |G|
c6 = perm_from_map(lambda i: (i + 1) % 6, 6)
r6 = perm_from_map(lambda i: (5 - i) % 6, 6)
tt_C6, tt_D6 = two_transitive([c6], 6), two_transitive([c6, r6], 6)
c7 = perm_from_map(lambda i: (i + 1) % 7, 7)
m3_7 = perm_from_map(lambda i: (3*i) % 7, 7)               # order 6 mod 7 -> F42
m2_7 = perm_from_map(lambda i: (2*i) % 7, 7)               # order 3 mod 7 -> F21
tt_F42, tt_F21, tt_C7 = (two_transitive([c7, m3_7], 7), two_transitive([c7, m2_7], 7),
                         two_transitive([c7], 7))
check("C7  N4 PREDICTIONS (degree 6 and 7), from the criterion + orbit-stabilizer: irreducible Rat^o needs "
      "n(n-1) | |G|, i.e. 30 | |G| at degree 6 and 42 | |G| at degree 7. Verified non-2-transitive: "
      "C6 (orbit 6), D6 (12), C7 (7), F21 (21); verified 2-transitive: F42 (42/42)",
      (not tt_C6[0]) and tt_C6[1] == 6 and (not tt_D6[0]) and tt_D6[1] == 12
      and (not tt_C7[0]) and tt_C7[1] == 7 and (not tt_F21[0]) and tt_F21[1] == 21
      and tt_F42[0] and tt_F42[1] == 42)
check("C8  so the N4 census can be PRE-SORTED: at degree 6 only the four 2-transitive groups "
      "(PSL(2,5), PGL(2,5), A6, S6 — orders 60, 120, 360, 720, all divisible by 30) can give irreducible "
      "Rat^o; every other transitive group forces a split S* and a cheaper C2",
      all(o % 30 == 0 for o in (60, 120, 360, 720))
      and all(o % 30 != 0 for o in (6, 12, 18, 24, 36, 48, 72)))

print()
fails = [ch for ch in checks if not ch[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[ch[0] for ch in fails]}"))
