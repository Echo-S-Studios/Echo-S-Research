# lesson1_checks.py — exact verification for the teaching session.
# Discipline: no float crosses a decision boundary. Equality/membership decided
# over Q or Q(sqrt 5) via sympy; mpmath-style floats appear ONLY in display strings.
import sympy as sp

x, t = sp.symbols('x t')
PHI = (1 + sp.sqrt(5)) / 2
checks = []

def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

def mahler_exact(p):
    """Exact Mahler measure of a monic integer polynomial via radical roots.
    Only valid when sympy can express roots in radicals (our small examples)."""
    roots = sp.roots(sp.Poly(p, x))  # dict root -> multiplicity, exact
    m = sp.Integer(1)
    for r, mult in roots.items():
        a = sp.Abs(r)
        # decision: is |root| > 1 ?  -- exact comparison
        gt1 = sp.simplify(a - 1) > 0
        if gt1 == sp.true or gt1 is True:
            m *= a**mult
    return sp.radsimp(sp.simplify(m))

# ---------------------------------------------------------------- C1: worked M values
m_phi  = mahler_exact(x**2 - x - 1)
check("C1a  M(x^2 - x - 1) = phi  (exactly)", sp.simplify(m_phi - PHI) == 0,
      f"M = {m_phi}  ~ {float(m_phi):.6f}")
m_sq2  = mahler_exact(x**2 - 2)
check("C1b  M(x^2 - 2) = 2  (exactly)", sp.simplify(m_sq2 - 2) == 0, f"M = {m_sq2}")
m_cyc  = mahler_exact(x**2 + x + 1)
check("C1c  M(x^2 + x + 1) = 1  (Kronecker floor case)", sp.simplify(m_cyc - 1) == 0)

# ---------------------------------------------------------------- C2: Lehmer polynomial structure
L = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1

# C2a reciprocity is an exact polynomial identity: x^10 * L(1/x) == L(x)
recip = sp.expand(x**10 * L.subs(x, 1/x)) - L
check("C2a  Lehmer L is reciprocal: x^10 L(1/x) = L(x)", sp.simplify(recip) == 0)

# C2b irreducible over Q
check("C2b  L irreducible over Q", sp.Poly(L, x).is_irreducible)

# C2c trace substitution: L(z)/z^5 = Q(z + 1/z), Q of degree 5, integer coeffs.
z = sp.symbols('z')
Lz_over_z5 = sp.expand(L.subs(x, z) / z**5)
# Build Q by symmetric reduction: express in t = z + 1/z
Q = sp.simplify(sp.expand(Lz_over_z5.rewrite(sp.Pow)))
# do it constructively: powers z^k + z^-k are Chebyshev-like polys in t
def zk_plus_zmk(k):
    # p_k(t) with z^k + z^-k = p_k(z + 1/z); recurrence p_k = t*p_{k-1} - p_{k-2}
    p0, p1 = sp.Integer(2), t
    if k == 0: return p0
    if k == 1: return p1
    for _ in range(2, k + 1):
        p0, p1 = p1, sp.expand(t*p1 - p0)
    return p1
# L(z)/z^5 = (z^5+z^-5) + (z^4+z^-4) - (z^2+z^-2) - (z+z^-1) - 1
Qt = sp.expand(zk_plus_zmk(5) + zk_plus_zmk(4) - zk_plus_zmk(2) - zk_plus_zmk(1) - 1)
# verify the identity exactly: substitute t = z + 1/z and compare with L(z)/z^5
diff = sp.simplify(sp.expand(Qt.subs(t, z + 1/z)) - Lz_over_z5)
check("C2c  trace form: L(z)/z^5 = Q(z+1/z) with Q = " + str(Qt), diff == 0)

# C2d Sturm count of real roots of Q: inside (-2,2) -> on-circle pairs of L;
#      outside [-2,2] -> real off-circle pairs (beta, 1/beta).
Qpoly = sp.Poly(Qt, t)
n_inside  = Qpoly.count_roots(-2, 2)        # exact (Sturm)
n_total   = Qpoly.count_roots(-sp.oo, sp.oo)
n_outside = n_total - Qpoly.count_roots(-2, 2)
check("C2d  Q has 5 real roots: 4 in (-2,2), 1 outside", (n_total, n_inside, n_outside) == (5, 4, 1),
      f"total={n_total}, inside={n_inside}, outside={n_outside}")
# consequence: L has 8 on-circle roots + one real pair (beta, 1/beta); exactly ONE root outside |z|=1.

# C2e that outside root of Q is > 2 (not < -2): count in (2, oo)
check("C2e  the off-circle trace root lies in (2, +oo)", Qpoly.count_roots(2, sp.oo) == 1)

# C2f Salem beta display value (float for DISPLAY ONLY): the real root of L in (1,2)
beta_iv = sp.Poly(L, x).count_roots(1, 2)
check("C2f  L has exactly one real root in (1,2)  [= Lehmer's number]", beta_iv == 1)
beta_num = sp.nsolve(L, x, sp.Rational(6, 5))
print(f"      display: Lehmer's number beta ~ {sp.N(beta_num, 15)}")

# ---------------------------------------------------------------- C3: Smyth's plastic number side
Pl = x**3 - x - 1
nonrecip = sp.simplify(sp.expand(x**3 * Pl.subs(x, 1/x)) - Pl) != 0 and \
           sp.simplify(sp.expand(x**3 * Pl.subs(x, 1/x)) + Pl) != 0
check("C3a  x^3 - x - 1 is non-reciprocal (exact identity fails both signs)", nonrecip)
check("C3b  x^3 - x - 1 has exactly one real root, in (1,2)",
      sp.Poly(Pl, x).count_roots(-sp.oo, sp.oo) == 1 and sp.Poly(Pl, x).count_roots(1, 2) == 1)
# its two complex roots: |z|^2 = z*conj(z); product of all three roots = 1 (constant term -(-1)=1)
# so |complex pair|^2 * theta0 = 1  -> pair is INSIDE the circle since theta0 > 1  => Pisot.
theta0 = sp.nsolve(Pl, x, sp.Rational(13, 10))
print(f"      display: plastic number theta0 ~ {sp.N(theta0, 12)}  (Smyth floor, non-reciprocal case)")

# ---------------------------------------------------------------- C4: dsum multiplicativity (elementary law)
m_prod = mahler_exact(sp.expand((x**2 - x - 1) * (x**2 - 2)))
check("C4   M(p*q) = M(p) M(q):  M((x^2-x-1)(x^2-2)) = 2*phi exactly",
      sp.simplify(m_prod - 2*PHI) == 0, f"M = {sp.simplify(m_prod)}")

# ---------------------------------------------------------------- C5: the angle wall vs Lehmer's polynomial
# Emitted on-circle eigenvalues are 4th roots of unity [FORCED, test_p2_02_angle].
# If L's Salem conjugates could appear, some root of L would satisfy x^4 = 1.
g = sp.gcd(sp.Poly(L, x), sp.Poly(x**4 - 1, x))
check("C5a  gcd(L, x^4 - 1) = 1  — no root of L is a 4th root of unity", g == sp.Poly(1, x),
      f"gcd = {g.as_expr()}")
# and integrality of char polys drags the WHOLE conjugate set along:
check("C5b  L irreducible => any emitted spectrum containing beta contains all 10 conjugates",
      sp.Poly(L, x).is_irreducible, "min-poly divides char-poly (Z coeffs)")

# ---------------------------------------------------------------- C6: closure of {1} u [phi, oo) under * and ^n  (elementary half)
# decided symbolically: if a,b in {1} u [phi,oo) then ab in {1} u [phi,oo).
a, b = sp.symbols('a b', positive=True)
# case analysis is the proof; verify the only nontrivial numeric fact: phi^2 >= phi  and  phi*phi > phi
check("C6   phi^2 > phi (exact), so [phi,oo) is closed under products/powers",
      sp.simplify(PHI**2 - PHI) > 0, f"phi^2 - phi = {sp.simplify(PHI**2 - PHI)} = 1... wait")
# phi^2 - phi = 1 exactly (phi^2 = phi + 1): show it
check("C6b  phi^2 = phi + 1 (the defining relation, exact)", sp.simplify(PHI**2 - PHI - 1) == 0)

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed" + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
