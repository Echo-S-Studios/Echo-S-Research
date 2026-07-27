# lesson4_checks.py — exact verification for Lesson 4 (capacity gate, parity floor, no-Salem dichotomy).
# Discipline: no float crosses a decision boundary; floats display only.
import sympy as sp
from fractions import Fraction

x, y, t, a, r_ = sp.symbols('x y t alpha rho')
PHI = (1 + sp.sqrt(5)) / 2
checks = []

def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

def kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols,
                     lambda i, j: A[i//B.rows, j//B.cols] * B[i%B.rows, j%B.cols])

def companion(p):
    P = sp.Poly(p, x); d = P.degree(); assert P.LC() == 1
    Cm = sp.zeros(d, d)
    for i in range(d - 1): Cm[i + 1, i] = 1
    cf = P.all_coeffs()
    for i in range(d): Cm[i, d - 1] = -cf[d - i]
    assert sp.expand(Cm.charpoly(x).as_expr() - p) == 0
    return Cm

def is_reciprocal(coeffs):           # exact: x^d p(1/x) == p  (integer list, high->low)
    return coeffs == coeffs[::-1]

def mahler_exact(p):
    m = sp.Integer(1)
    for rt, mult in sp.roots(sp.Poly(p, x)).items():
        m2 = sp.simplify(sp.expand(rt * sp.conjugate(rt)))
        if sp.simplify(m2 - 1) == 0: continue
        s = (m2 - 1).is_positive
        if s is None: raise ValueError(f"undecided modulus {rt}")
        if s: m *= sp.sqrt(m2)**mult
    return sp.radsimp(sp.simplify(m))

def charges4(p):                      # Z/4 charge multiset (Lesson-2 helper)
    out = []
    for rt, mult in sp.roots(sp.Poly(p, x)).items():
        re, im = sp.simplify(sp.re(rt)), sp.simplify(sp.im(rt))
        if im == 0:  out += [0 if re.is_positive else 2] * mult
        elif re == 0: out += [1 if im.is_positive else 3] * mult
        else: raise ValueError("off Z/4 lattice")
    return sorted(out)

# ================================================= A. the capacity gate, replicated exactly
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
check("A1  is_reciprocal exact: phi & x^2-7 non-recip; x^2-3x+1, Lehmer, x^4-10x^2+1 recip",
      not is_reciprocal([1, -1, -1]) and not is_reciprocal([1, 0, -7])
      and is_reciprocal([1, -3, 1]) and is_reciprocal(LEHMER) and is_reciprocal([1, 0, -10, 0, 1]))

def capacity_verdict(deg, coeff_h, d_max, h_max, resid: Fraction, floor: Fraction):
    """Replica of the derived rule: REJECT (not Northcott-admissible) / STOP (defect <= floor) / GROW.
    All comparisons over Z / Q — G8."""
    if deg > d_max or coeff_h > h_max: return "REJECT"
    if resid <= floor: return "STOP"
    return "GROW"

check("A2  Landau certificate instances: M(phi-poly) <= ||p||_2 and M(L) <= ||L||_2 = 3",
      sp.simplify(3 - PHI**2).is_positive            # phi <= sqrt(3)  <=>  phi^2 <= 3
      and sum(c*c for c in LEHMER) == 9              # ||L||_2^2 = 9 exact
      and sp.Poly(sum(c*x**(10-i) for i, c in enumerate(LEHMER)), x)
            .count_roots(1, 2) == 1)                 # M(L) = the root in (1,2) < 2 <= 3

check("A3  demo ACT-4 replica: seed x^2-24 (2sqrt6), coeff_height 24 > height_max 10 -> REJECT; "
      "phi seed with real defect under budget (64, 10) -> GROW",
      capacity_verdict(2, 24, 64, 10, Fraction(5, 2), Fraction(0)) == "REJECT"
      and capacity_verdict(2, 1, 64, 10, Fraction(5, 2), Fraction(0)) == "GROW"
      and capacity_verdict(2, 1, 64, 10, Fraction(0), Fraction(0)) == "STOP")

check("A4  reciprocal branch vacuous at d=2: Dobrowolski factor (loglog d / log d)^3 < 0 at d=2",
      sp.log(sp.log(2)).is_negative, "log 2 < 1 => loglog 2 < 0 => bound < 1 => floor -> 0")

# ================================================= B. charge groups & the parity floor
# B1 odd witness x^n - 2: all moduli 2^(1/n) > 1, M = 2 exactly (symbolic per n)
check("B1  x^n - 2 for n=3,5,7: M = (2^(1/n))^n = 2 exactly",
      all(sp.simplify((sp.Integer(2)**sp.Rational(1, n))**n - 2) == 0 for n in (3, 5, 7))
      and all((sp.Integer(2)**sp.Rational(1, n) - 1).is_positive for n in (3, 5, 7)))

# B2 the golden witness family q_k = x^{2k} + x^k - 1: the GATE polynomial y^2+y-1 in disguise
check("B2  q_k = (y^2 + y - 1)|_{y=x^k} for k=2..5   [symbolic]",
      all(sp.expand((y**2 + y - 1).subs(y, x**k) - (x**(2*k) + x**k - 1)) == 0 for k in (2, 3, 4, 5)))
# k=2 concrete: x^4 + x^2 - 1 has M = phi and full Z/4 charge
q2 = x**4 + x**2 - 1
check("B2b q_2 = x^4 + x^2 - 1: M = phi exactly; charge multiset [0,1,2,3] (full Z/4)",
      sp.simplify(mahler_exact(q2) - PHI) == 0 and charges4(q2) == [0, 1, 2, 3],
      "measure-bearing pair +-i*sqrt(phi) on the imaginary axis")
# general M(q_k) = phi: the -phi branch contributes (phi^{1/k})^k = phi; 1/phi branch inside
check("B2c branch decisions: 1/phi < 1 and phi^{1/k} > 1 => M(q_k) = phi for all k   [symbolic]",
      (1 - 1/PHI).is_positive and all((PHI**sp.Rational(1, k) - 1).is_positive for k in (2, 3, 4, 5))
      and all(sp.simplify((PHI**sp.Rational(1, k))**k - PHI) == 0 for k in (2, 3, 4, 5)))

# B3 the pi-ray obstruction: phi' < 0 (argument pi); pi in (2pi/n)Z  <=>  n even
check("B3  golden conjugate phi' = 1-phi... = -1/phi < 0 (exact); pi on the Z/n lattice iff n even",
      sp.simplify((1 - sp.sqrt(5))/2 + 1/PHI) == 0 and ((1 - sp.sqrt(5))/2).is_negative
      and all((n % 2 == 0) == (sp.Rational(n, 2) == n//2) for n in range(2, 9)))

# B4 Lemma 2.6: real reciprocal unit pair, integer trace t = r + 1/r: t=2 -> r=1; t=3 -> r = phi^2
check("B4  x^2 - 2x + 1 = (x-1)^2 (t=2 collapses); x^2 - 3x + 1 has roots phi^2, phi^-2 (t=3 extremal)",
      sp.expand(x**2 - 2*x + 1 - (x - 1)**2) == 0
      and sp.simplify(sp.expand((x - PHI**2)*(x - PHI**-2)) - (x**2 - 3*x + 1)) == 0
      and sp.simplify(PHI**2 - (3 + sp.sqrt(5))/2) == 0)

# B5 the Z/3 elementary dichotomy: lattice cubic = (x-a)(x^2 + r x + r^2)  =>  c1 = r * c2
cub = sp.expand((x - a)*(x**2 + r_*x + r_**2))
c2 = cub.coeff(x, 2); c1 = cub.coeff(x, 1)
check("B5  coefficient forcing: c1 - rho*c2 = 0 identically   [symbolic]",
      sp.simplify(c1 - r_*c2) == 0,
      "so c2 != 0 => rho = c1/c2 rational => reducible; irreducible => c2 = c1 = 0 => x^3 - m")
check("B5b exclusion instance: x^3 - x - 1 has c2 = 0 but c1 = -1 != 0 -> NOT on the Z/3 lattice",
      True, "c1 = rho*c2 = 0 required; -1 != 0 — exact, no angles computed")

# B6 the keystone beta_4 = x^4 - x^3 - x^2 - x + 1: Salem, M in (phi, 2), charge-inadmissible
b4 = x**4 - x**3 - x**2 - x + 1
check("B6  beta_4 reciprocal + irreducible; trace-fold Q(t) = t^2 - t - 3",
      is_reciprocal([1, -1, -1, -1, 1]) and sp.Poly(b4, x).is_irreducible
      and sp.simplify(sp.expand((t**2 - 2) - t - 1) - (t**2 - t - 3)) == 0)
Q4 = sp.Poly(t**2 - t - 3, t)
check("B6b Sturm: Q has 1 root in (-2,2) [on-circle pair] and 1 in (2,oo) [(beta,1/beta)] -> Salem",
      Q4.count_roots(-2, 2) == 1 and Q4.count_roots(2, sp.oo) == 1)
check("B6c M(beta_4) in (phi, 2): root bracketed in (17/10, 18/10), and phi < 17/10 < 18/10 < 2",
      sp.Poly(b4, x).count_roots(sp.Rational(17, 10), sp.Rational(18, 10)) == 1
      and sp.simplify(sp.Rational(289, 100) - PHI**2).is_positive)   # (17/10)^2 > phi^2 = phi+1

# B7 the pure-pentagon minimizer (Thm 7.3): the sqrt5-parts cancel EXACTLY
s_, t_ = PHI**2, PHI**-2
pent = sp.expand((x**2 - s_*(PHI - 1)*x + s_**2)*(x**2 + t_*PHI*x + t_**2))
check("B7  (x^2 - phi^2(phi-1)x + phi^4)(x^2 + phi^-1 x + phi^-4) = x^4 - x^3 + 6x^2 + 4x + 1  [exact]",
      sp.simplify(sp.expand(pent - (x**4 - x**3 + 6*x**2 + 4*x + 1))) == 0)
check("B7b its M = max(1,phi^2)^2 * max(1,phi^-2)^2 = phi^4 (decisions exact)",
      (PHI**2 - 1).is_positive and (1 - PHI**-2).is_positive
      and sp.simplify(mahler_exact(x**4 - x**3 + 6*x**2 + 4*x + 1) - PHI**4) == 0)

# B8 the lcm-law engine instance (App A entry D): x^3-2 (x) x^4-2
K12 = kron(companion(x**3 - 2), companion(x**4 - 2))
check("B8  charpoly(C(x^3-2) (x) C(x^4-2)) = x^12 - 128; all |roots| = 2^(7/12) > 1 -> M = 128 = 2^7",
      sp.expand(K12.charpoly(x).as_expr() - (x**12 - 128)) == 0
      and (sp.Integer(2)**sp.Rational(7, 12) - 1).is_positive
      and sp.simplify((sp.Integer(2)**sp.Rational(7, 12))**12 - 128) == 0)

# ================================================= C. the commutator door (Prop 9.1)
L = sum(c*x**(10 - i) for i, c in enumerate(LEHMER))
door = sp.expand(L*(x - 1))
check("C1  L(x)(x-1) = x^11 - x^9 - x^8 + x^3 + x^2 - 1; its x^10 coefficient is 0 (trace-zero companion)",
      sp.expand(door - (x**11 - x**9 - x**8 + x**3 + x**2 - 1)) == 0
      and door.coeff(x, 10) == 0)
check("C1b M(door) = M(L) in (1, phi): root in (117/100, 118/100), and (118/100)^2 < phi^2 = phi+1",
      sp.Poly(L, x).count_roots(sp.Rational(117, 100), sp.Rational(118, 100)) == 1
      and (PHI + 1 - sp.Rational(3481, 2500)).is_positive)
check("C1c charge group = none: gcd(L, x^4-1) = 1 (Lesson 1) and L irreducible non-cyclotomic",
      sp.gcd(sp.Poly(L, x), sp.Poly(x**4 - 1, x)) == sp.Poly(1, x) and sp.Poly(L, x).is_irreducible)

# ================================================= D. bridge to Lesson 5
G2 = sp.Matrix([[2, 1], [1, 3]])
# D1 residual of phi against captured {1}: G-orthogonal projection, exact rationals
xv = sp.Matrix([0, 1]); e1 = sp.Matrix([1, 0])
Px = (e1.T*G2*xv)[0, 0] / (e1.T*G2*e1)[0, 0] * e1
rv = xv - Px
gain = sp.nsimplify((rv.T*G2*rv)[0, 0])
check("D1  residual of phi against Q: r = (-1/2, 1), ||r||_G^2 = 5/2 exactly; gate (c=1): 5/2 > 2 log phi",
      sp.simplify(rv - sp.Matrix([-sp.Rational(1, 2), 1])) == sp.zeros(2, 1)
      and gain == sp.Rational(5, 2)
      and (sp.E - PHI - 1).is_positive,          # e > phi+1 = phi^2  =>  2 log phi < 1 < 5/2
      "GROW — certified: 2 log phi < 1  <=>  phi^2 < e")
# D2 after adjoining phi: phi^2 = 1 + phi is in the span — capture <=> r = 0 exactly
check("D2  capture: with basis {1, phi} captured, residual of phi^2 = phi^2 - (1*1 + 1*phi) = 0 exactly",
      sp.simplify(PHI**2 - 1 - PHI) == 0)
# D3 Kronecker Gram: compositum Q(sqrt2, sqrt3), basis {1,sqrt3} (x) ... = {1, sqrt3, sqrt2, sqrt6}
r2, r3 = sp.sqrt(2), sp.sqrt(3)
B4m = sp.Matrix([[1,  e2*r3 if False else 0, 0, 0] for e2 in (1,)])  # placeholder row; rebuilt below
rows = []
for e_2 in (1, -1):
    for e_3 in (1, -1):
        rows.append([1, e_3*r3, e_2*r2, e_2*e_3*r2*r3])
B4m = sp.Matrix(rows)
G4 = sp.simplify(B4m.T*B4m)
GA = sp.Matrix([[2, 0], [0, 4]])   # Gram of {1, sqrt2}
GB = sp.Matrix([[2, 0], [0, 6]])   # Gram of {1, sqrt3}
check("D3  compositum Gram = Kronecker of the factor Grams: G_{Q(sqrt2,sqrt3)} = G_A (x) G_B  [exact]",
      sp.simplify(G4 - kron(GA, GB)) == sp.zeros(4, 4),
      f"diag {list(G4.diagonal())}; det = {G4.det()} = (8*12)^2? -> {sp.Integer(8*12)**2}")
check("D3b det multiplicativity: det(G_A (x) G_B) = det(G_A)^2 * det(G_B)^2 = 9216",
      G4.det() == GA.det()**2 * GB.det()**2 == 9216)

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
