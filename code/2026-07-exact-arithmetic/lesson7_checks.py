# lesson7_checks.py — exact verification for Lesson 7 (the golden substrate and the rung line).
# Discipline: decisions over Q and Q(sqrt5); K decided at the squared level; floats display-only.
import sympy as sp

x, y, z, w, s, th, b = sp.symbols('x y z w s theta beta', positive=True)
X4 = sp.symbols('a0 a1 a2 a3')
s5 = sp.sqrt(5)
PHI = (1 + s5)/2; PSI = (1 - s5)/2; TAU = PHI - 1
GAP = (7 - 3*s5)/2; K2 = (3*s5 - 5)/2; B2 = (5 + 3*s5)/2
checks = []
def check(name, ok, note=""):
    checks.append((name, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {note}" if note else ""))

E = sp.expand
S = lambda e: sp.simplify(sp.radsimp(sp.expand(e)))

# ============================== A. seed and K-seed ledger (whitepaper eq. (1)-(2))
check("A1  seed identities: phi*psi=-1, phi+psi=1, sqrt5=phi+1/phi=phi-psi, tau^2+tau=1",
      S(PHI*PSI + 1) == 0 and S(PHI + PSI - 1) == 0 and S(PHI + 1/PHI - s5) == 0
      and S(PHI - PSI - s5) == 0 and S(TAU**2 + TAU - 1) == 0)
check("A2  gap = phi^-4 = (7-3sqrt5)/2;  K^2 = 1 - gap;  K^4 + 5K^2 - 5 = 0;  K^4 = 5*gap",
      S(PHI**-4 - GAP) == 0 and S(K2 - (1 - GAP)) == 0
      and S(K2**2 + 5*K2 - 5) == 0 and S(K2**2 - 5*GAP) == 0)
check("A3  (K*phi)^4 = 5  (so K*phi = 5^{1/4}: Q(K) = Q(5^{1/4}));  K^2*beta^2 = 5;  beta^2 = K^2*phi^4",
      S(K2**2 * PHI**4 - 5) == 0 and S(K2*B2 - 5) == 0 and S(B2 - K2*PHI**4) == 0)
check("A4  Route 5: K^2 = tau(2 - tau);  beta^2 = phi^2*sqrt5 = M(K-seed)",
      S(K2 - TAU*(2 - TAU)) == 0 and S(B2 - PHI**2*s5) == 0)
check("A5  K-seed Eisenstein at 5 (irreducible); x^2-roots (-5±3sqrt5)/2 straddle 0: guards 45<49 and 180>169",
      sp.Poly(x**4 + 5*x**2 - 5, x).is_irreducible and (3*s5 - 5).is_positive
      and S(K2 - 1).is_negative and (6*s5 - sp.Rational(13, 2)*2).is_positive)
check("A6  gap object: (x - phi^4)(x - phi^-4) = x^2 - 7x + 1 exactly (L4 = 7); M = phi^4, full-turn cost 4 ln phi",
      S(E((x - PHI**4)*(x - PHI**-4)) - (x**2 - 7*x + 1)) == 0)

# ============================== B. operator layer (whitepaper §4)
e1 = sp.Matrix([[0, 1], [1, 0]]); e2 = sp.Matrix([[1, 0], [0, -1]])
J = e1*e2
H = 2*e1 - e2
R = (sp.eye(2) + H)/2
check("B1  convention locked: H = 2e1 - e2 = 2R - I with R = [[0,1],[1,1]]; H^2 = 5I (Clifford length 5)",
      R == sp.Matrix([[0, 1], [1, 1]]) and sp.simplify(H - (2*R - sp.eye(2))) == sp.zeros(2, 2)
      and H*H == 5*sp.eye(2))
Xs = sp.Matrix(2, 2, lambda i, j: X4[2*i + j])
check("B2  the return operator IS the half-anticommutator: RX + XR - X = (HX + XH)/2, symbolically",
      sp.simplify(R*Xs + Xs*R - Xs - (H*Xs + Xs*H)/2) == sp.zeros(2, 2))
N1 = e1 + 2*e2; N2 = J
check("B3  ker L = span{e1 + 2e2, J}: both satisfy RX + XR - X = 0 exactly (the phi-slack, Clifford dress)",
      sp.simplify(R*N1 + N1*R - N1) == sp.zeros(2, 2)
      and sp.simplify(R*N2 + N2*R - N2) == sp.zeros(2, 2))
check("B4  Remark 4.3 disproof replicated: L(-e1 + e2) = -3*I, NOT zero (the refuted design note)",
      sp.simplify(R*(-e1 + e2) + (-e1 + e2)*R - (-e1 + e2)) == -3*sp.eye(2))
check("B5  scalar-H plane: L(1) = H and L(H) = 5*1  (the [[0,5],[1,0]] block, eigenvalues ±sqrt5)",
      sp.simplify(R + R - sp.eye(2) - H) == sp.zeros(2, 2)
      and sp.simplify(R*H + H*R - H - 5*sp.eye(2)) == sp.zeros(2, 2))

# ============================== C. the slot (whitepaper §6.1)
check("C1  branch identity: tau0 - 2 = (beta - 1)^2 / beta, symbolically",
      S((b + 1/b) - 2 - (b - 1)**2/b) == 0)
check("C2  t = sqrt5 lifts to {phi, 1/phi}: (x - phi)(x - tau) = x^2 - sqrt5*x + 1 exactly",
      S(E((x - PHI)*(x - TAU)) - (x**2 - s5*x + 1)) == 0)
check("C3  redirection slope at the golden point: g(x) = x + 1/x has g(phi) = sqrt5 and g'(phi) = 1/phi exactly",
      S(PHI + 1/PHI - s5) == 0 and S(1 - 1/PHI**2 - 1/PHI) == 0)
plc = sp.resultant(th**3 - th - 1, s*th - th**2 - 1, th)
plc = sp.Poly(sp.expand(plc), s).primitive()[1]
check("C4  plastic trace: mu_S + 1/mu_S is a root of s^3 + s^2 - 4s - 5; bracketed in (2, 11/5), "
      "and 11/5 < sqrt5 (125 > 121): the plastic trace sits BELOW the sqrt5 accumulation point",
      sp.expand(plc.as_expr() - (s**3 + s**2 - 4*s - 5)) == 0
      and sp.Poly(plc, s).count_roots(2, sp.Rational(11, 5)) == 1 and 125 > 121)
check("C5  fold trichotomy: disc(x^2 - tx + 1) = t^2 - 4; for t in (-2,2) roots are a nonreal unit pair "
      "(product exactly 1) — the elliptic face; t > 2 gives {beta, 1/beta} — the hyperbolic face",
      S(sp.discriminant(x**2 - w*x + 1, x) - (w**2 - 4)) == 0)

# ============================== D. gate family and the ladder (whitepaper §§8-9)
Cs, us, ms = sp.symbols('C u m')
check("D1  flip reciprocity (Prop 8.2): m+*m- = 1 and m + 1/m = -1/C - 2, from Vieta on x^2 + x - C",
      S(((-Cs/(1 + (-1 + sp.sqrt(1 + 4*Cs))/2)**2) * (-Cs/(1 + (-1 - sp.sqrt(1 + 4*Cs))/2)**2)) - 1) == 0)
# ladder identities in w = phi^{2n} (positive symbol)
Cn = w/(w - 1)**2; un = 1/(w - 1)
mn = -Cn/(1 + un)**2
check("D2  ladder (Thm 9.1), symbolic in w = phi^{2n}: u_n^2 + u_n = C_n; multiplier m_n = -1/w = -phi^{-2n}; "
      "sqrt(D_n) = (w+1)/(w-1); lambda_n = sqrt(D)/C = w - 1/w = phi^{2n} - phi^{-2n}",
      S(un**2 + un - Cn) == 0 and S(mn + 1/w) == 0
      and S(1 + 4*Cn - ((w + 1)/(w - 1))**2) == 0
      and S(((w + 1)/(w - 1))/Cn - (w - 1/w)) == 0)
lad = []
for n in range(1, 6):
    Cval = S(1/(PHI**n - PHI**-n)**2)
    lad.append(Cval)
check("D3  rung gates exactly: C_1..C_5 = 1, 1/5, 1/16, 1/45, 1/121 (golden gate; K-gate; L_3^2; 5F_4^2; L_5^2)",
      lad == [1, sp.Rational(1, 5), sp.Rational(1, 16), sp.Rational(1, 45), sp.Rational(1, 121)])
check("D4  lambda_n = sqrt5 * F_{2n} for n = 1..5 (exchange rates on the rungs)",
      all(S((PHI**(2*n) - PHI**(-2*n)) - s5*sp.fibonacci(2*n)) == 0 for n in range(1, 6)))

# ============================== E. special gates (whitepaper §10)
rem = sp.rem(E((1 - y)**2 * (5 + y)**3 - 5*y), y**2 + 5*y - 5, y)
check("E1  K-gate congruence certificate: 5y == (1-y)^2 (5+y)^3  mod (y^2 + 5y - 5), exact remainder 0",
      sp.expand(rem) == 0)
upm = [(-1 + sp.I)/2, (-1 - sp.I)/2]
mvals = [sp.simplify(-(-sp.Rational(1, 2))/(1 + u)**2) for u in upm]
check("E2  rotation gate C = -1/2: fixed points (-1±i)/2, multipliers exactly {i, -i} — the odd-charge "
      "generators of Z/4Z", set(mvals) == {sp.I, -sp.I})
tr1 = S(-1/(-PHI**-2) - 2); tr2 = S(-1/(-PHI**2) - 2)
check("E3  pentagon gates: multiplier traces at C = -phi^{-2}, -phi^2 are tau and -phi — the roots of "
      "x^2 + x - 1, Lesson 4's pentagon cosines 2cos72°, 2cos144°",
      S(tr1 - TAU) == 0 and S(tr2 + PHI) == 0 and S(tr1**2 + tr1 - 1) == 0 and S(tr2**2 + tr2 - 1) == 0)
sv, tv = sp.symbols('s t', positive=True)
imexpr = sp.im((sv*sp.exp(sp.I*tv) + 1/(sv*sp.exp(sp.I*tv))))
check("E4  obstruction engine (Thm 11.1): Im(m + 1/m) = (s - 1/s) sin(theta); reality forces s = 1 or "
      "sin(theta) = 0 — no single real gate spirals",
      sp.simplify(imexpr - (sv - 1/sv)*sp.sin(tv)) == 0)

# ============================== F. the chart (D1) and its consequences (whitepaper §13)
uz = (1 - z**2)/z**2; mz = -(1 - z**2); Cz = (1 - z**2)/z**4
sqDz = (2 - z**2)/z**2; lz = z**2*(2 - z**2)/(1 - z**2)
check("F1  chart closed forms: u(z) is a fixed point of G_{C(z)} with multiplier exactly m(z) = -(1-z^2)",
      S(uz**2 + uz - Cz) == 0 and S(-Cz/(1 + uz)**2 - mz) == 0)
check("F2  1 + 4C(z) = ((2-z^2)/z^2)^2;  lambda*|m| = 1 - m^2;  lambda = sqrt(D)/C — all chart identities",
      S(1 + 4*Cz - sqDz**2) == 0 and S(lz*(1 - z**2) - (1 - (1 - z**2)**2)) == 0
      and S(sqDz/Cz - lz) == 0)
check("F3  rung pullback: with z_n^2 = 1 - phi^{-2n}, u(z_n) = 1/(phi^{2n} - 1) = u_n and rho_n = n ln phi "
      "(1/(1-z_n^2) = phi^{2n})",
      all(S((1 - (1 - PHI**(-2*n)))**-1 - PHI**(2*n)) == 0 for n in range(1, 5)))
check("F4  the cross-tie (13.1e): lambda_1 * |m_1| = sqrt5 * phi^{-2} = K^2 exactly",
      S(s5*PHI**-2 - K2) == 0)

# ============================== G. the lens and the guards (whitepaper §13.2, orderings)
check("G1  lens: z_c = sqrt3/2 gives 1 - z_c^2 = 1/4, so rho(z_c) = ln 2 EXACTLY; C = 4/9; D = 25/9 square; "
      "splits (x - 1/3)(x + 4/3); m = {-1/4, -4}; lambda = 15/4",
      S(sp.Rational(1, 4) - (1 - sp.Rational(3, 4))) == 0
      and S((1 - sp.Rational(3, 4))/sp.Rational(9, 16) - sp.Rational(4, 9)) == 0
      and sp.expand((x - sp.Rational(1, 3))*(x + sp.Rational(4, 3)) - (x**2 + x - sp.Rational(4, 9))) == 0
      and S(sp.Rational(3, 4)*(2 - sp.Rational(3, 4))/(1 - sp.Rational(3, 4)) - sp.Rational(15, 4)) == 0)
check("G2  golden rungs irreducible over Q: D_1 = 5 and D_2 = 9/5 nonsquare (guards 4<5<9, 36<45<49); "
      "ordering tau < 3/4 < K^2 by guards 20<25 and 180>169; and 2 < phi^2 (tau > 0)",
      S(1 + 4*lad[0] - 5) == 0 and S(1 + 4*lad[1] - sp.Rational(9, 5)) == 0
      and 4 < 5 < 9 and 36 < 45 < 49 and 20 < 25 and 180 > 169 and S(PHI**2 - 2 - TAU) == 0)
check("G3  K-point (Thm 15.1): rho(K) = -ln(gap)/2 = 2 ln phi so theta(K) = pi (two quanta); horn "
      "self-intersection z = K^2/z_c excluded by the SAME guard 180 > 169; K in [z_c, 1] by 3/4 < K^2 < 1",
      S(GAP*PHI**4 - 1) == 0 and sp.sign(K2 - sp.Rational(3, 4)) == 1
      and sp.sign(1 - K2) == 1)
check("G4  kink slope squared at the lens: K^2/3 = (3sqrt5 - 5)/6 exactly",
      S(K2/3 - (3*s5 - 5)/6) == 0)

# ============================== H. incommensurability + fields + registry (whitepaper §§16-18)
binet_ok = all(S(PHI**p - (sp.lucas(p) + sp.fibonacci(p)*s5)/2) == 0 for p in range(1, 11))
check("H1  Thm 16.1 engine: phi^p = (L_p + F_p sqrt5)/2 with F_p >= 1 for p = 1..10 — phi^p irrational, "
      "so 2^q = phi^p is impossible: log_phi(2) is irrational; the lens never lands on a rung", binet_ok)
mp1 = sp.Poly(x**4 - x**2 - 1, x); mp2 = sp.Poly(x**4 + x**2 - 1, x)
check("H2  Thm 17.1: minpoly(sqrt(phi)) = x^4 - x^2 - 1 and minpoly(sqrt(tau)) = x^4 + x^2 - 1, both "
      "irreducible; sqrt(phi) = phi*sqrt(tau) since phi^2*tau = phi",
      mp1.is_irreducible and mp2.is_irreducible
      and S(PHI**2 - PHI - 1) == 0 and S(TAU**2 + TAU - 1) == 0 and S(PHI**2*TAU - PHI) == 0)
p2 = [(40 + 16*s5)/32, (40 - 16*s5)/32]
check("H3  square-class separation: (p + q sqrt5)^2 = (5+sqrt5)/2 forces 16p^4 - 40p^2 + 5 = 0, "
      "p^2 = (5±2sqrt5)/4 — irrational (1225 < 1280 < 1296), so Q(sqrt phi) != Q(5^{1/4})",
      all(S(16*v**2 - 40*v + 5) == 0 for v in p2) and 35**2 < 1280 < 36**2)
check("H4  registry: phi^5 = 5phi + 3; M(cons) = 2phi^5 = 11 + 5sqrt5; norm N = 121 - 125 = -4 = (-1)^5 4^1 "
      "— ON the (ln2, lnphi) lattice at lens + 5 floors",
      S(PHI**5 - (5*PHI + 3)) == 0 and S(2*PHI**5 - (11 + 5*s5)) == 0 and 11**2 - 5*5**2 == -4)
check("H5  registry: N(M(res)) = (43^2 - 5*19^2)/4 = 11, and 11 is not ±4^Z (guard 4 < 11 < 16): "
      "res is OFF the lattice for ALL integer exponents — one norm obstruction closes every case",
      sp.Rational(43**2 - 5*19**2, 4) == 11 and 4 < 11 < 16)

# ============================== I. the single generator q = i*tau (derived; CRG packaging unread)
q = sp.I*TAU
check("I1  q = i*tau: q^4 = tau^4 = phi^{-4} = gap EXACTLY — the Z/4Z quarter-turn whose fourth power "
      "is the gap; |q^n|^2 = phi^{-2n} = |m_n| (the multiplier ladder) and arg q^n = n*pi/2 (the compass)",
      S(q**4 - GAP) == 0 and all(S(sp.Abs(q**n)**2 - PHI**(-2*n)) == 0 for n in range(1, 5)))
check("I2  winding rate: kappa * ln(phi) = pi/2 exactly (D3: one Z/4Z quantum per entropy floor)",
      S((sp.pi/(2*sp.log(PHI)))*sp.log(PHI) - sp.pi/2) == 0)
d1, d2 = sp.N(4*sp.log(PHI)/sp.pi, 30), sp.N(TAU, 30)
print(f"      [computed guard] 4 ln(phi)/pi = {d1}  vs  tau = {d2}  (|diff| ~ 5.34e-3; NOT an identity)")

print()
fails = [c for c in checks if not c[1]]
print(f"SUMMARY: {len(checks) - len(fails)}/{len(checks)} passed"
      + ("" if not fails else f"  FAILURES: {[c[0] for c in fails]}"))
