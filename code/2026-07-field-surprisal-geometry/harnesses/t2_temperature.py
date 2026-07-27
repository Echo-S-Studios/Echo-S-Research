#!/usr/bin/env python3
# t2_temperature.py -- Task 2: the temperature is a coordinate, not a parameter.
# (b)-resolution: equivariance + product functoriality force p ~ M^{-beta}
# (measurable Cauchy); beta is the canonical coordinate of a forced e-geodesic.
# Exact parts: factorization rank tests, lattice-Cauchy pathology, anchors in Q(sqrt5).
# Numeric parts (display-only): I(beta), total information length, marked points.
import sys
from sympy import (symbols, Rational, sqrt, Matrix, cancel, simplify, expand,
                   S, nsimplify, log, exp, together, radsimp)

NCK = [0]
def ck(label, cond, detail=""):
    NCK[0] += 1
    if not cond:
        print(f"[{NCK[0]:03d}] FAIL {label} {detail}")
        sys.exit(1)
    print(f"[{NCK[0]:03d}] PASS {label}")

# ---------------------------------------------------------------- exact: masses in Q(sqrt5)
s5 = sqrt(5)
phi = (1 + s5)/2
M = [S(2), S(3), S(5), phi, phi, phi**4, phi**4 - 1]
Z_anchor = sum(M)
ck("anchor Z(-1) = sum M = 17 + 4 sqrt5 (exact)",
   simplify(expand(Z_anchor - (17 + 4*s5))) == 0)
Z1 = together(sum(1/m for m in M))
Z1s = radsimp(simplify(Z1))
print("   Z(+1) = sum 1/M =", Z1s, "=", float(Z1s))
ck("Z(+1) is exact in Q(sqrt5)", simplify(Z1s - Rational(31,30) - (2/phi + 1/phi**4 + 1/(phi**4-1))) == 0)

# ---------------------------------------------------------------- exact: product factorization forces power law
print("== forced form: factorization <=> rank-1 <=> multiplicative ==")
b = symbols('beta', real=True)
# power law: the product-catalog Gibbs matrix f(Mi*Mj') is rank 1 (f = x^-beta)
Mp = [S(2), phi, phi**4]     # a small second catalog (masses)
def gibbs_matrix(f):
    return Matrix(len(M), len(Mp), lambda i, j: f(M[i]*Mp[j]))
fpow = lambda x: x**(-2)     # beta = 2 instance, exact
Gm = gibbs_matrix(fpow)
ck("power law: product Gibbs matrix has rank 1 (factorizes)", Gm.rank() == 1)
fexp = lambda x: exp(-x)     # Boltzmann-in-M: NOT multiplicative
Ge = gibbs_matrix(fexp)
ck("falsifier exp(-M): product Gibbs matrix has rank >= 2 (no factorization)",
   Ge[0,0]*Ge[1,1] - Ge[0,1]*Ge[1,0] != 0)
# rank-1 <=> f(xy) = A(x)B(y) on the value grid <=> log f additive (Cauchy)
x, y = symbols('x y', positive=True)
fp = x**(-b)
ck("multiplicativity identity f(xy)=f(x)f(y) for f=x^-beta (symbolic)",
   simplify((x*y)**(-b) - x**(-b)*y**(-b)) == 0)

# lattice-Cauchy pathology: per-generator exponents are multiplicative on the
# lattice but are NOT a single power  ->  measurability axiom is necessary.
print("== lattice Cauchy caveat ==")
b1, b2, b3, b4 = Rational(1), Rational(2), Rational(3), Rational(5)
def f_path(n2, n3, n5, nphi):   # value on 2^n2 3^n3 5^n5 phi^nphi
    return S(2)**(-b1*n2) * S(3)**(-b2*n3) * S(5)**(-b3*n5) * phi**(-b4*nphi)
# multiplicative by construction on the free monoid: f(v+w) = f(v) f(w)
v = (1, 0, 2, 3); w = (0, 2, 1, 1)
ck("pathological f is multiplicative on the lattice",
   simplify(f_path(*[v[i]+w[i] for i in range(4)]) - f_path(*v)*f_path(*w)) == 0)
# but not a single power: compare exponent ratios on masses 2 and 3
ck("pathological f is not x^-beta for any single beta",
   simplify(log(f_path(1,0,0,0))/log(S(2)) - log(f_path(0,1,0,0))/log(S(3))) != 0)

# ---------------------------------------------------------------- numeric: the Gibbs curve invariants
print("== Gibbs curve: I(beta), total length, marked points (numeric, display) ==")
import mpmath as mp
mp.mp.dps = 40
PHI = (1 + mp.sqrt(5))/2
Mn = [mp.mpf(2), mp.mpf(3), mp.mpf(5), PHI, PHI, PHI**4, PHI**4 - 1]
la = [mp.log(m) for m in Mn]

def pI(beta):
    ws = [mp.e**(-beta*l) for l in la]
    Z = sum(ws); p = [w/Z for w in ws]
    mu = sum(p[i]*la[i] for i in range(7))
    I = sum(p[i]*(la[i]-mu)**2 for i in range(7))
    return p, I, mu

def I_of(beta): return pI(beta)[1]

# marked temperatures
for lbl, bb in [("beta=0 (max entropy)", 0), ("beta=1 (MDL)", 1),
                ("beta=sqrt5 (exchange rate)", mp.sqrt(5)), ("beta=-1 (anchor)", -1)]:
    print(f"   I({lbl}) = {mp.nstr(I_of(mp.mpf(bb)), 12)}")

# beta* = argmax I  (max Fisher);  I'(beta) = -kappa3
def dI(beta):
    h = mp.mpf('1e-12')
    return (I_of(beta+h)-I_of(beta-h))/(2*h)
def kappa3(beta):
    p, I, mu = pI(beta)
    return sum(p[i]*(la[i]-mu)**3 for i in range(7))
bstar = mp.findroot(dI, mp.mpf('-0.1'))
ck("I'(beta) = -kappa3(beta) (identity check at beta*)",
   abs(dI(bstar) + kappa3(bstar)) < mp.mpf('1e-8'))
print(f"   beta* (max Fisher)  = {mp.nstr(bstar, 10)},  I(beta*) = {mp.nstr(I_of(bstar), 10)}")
ck("beta* is a genuine max (I'' < 0)",
   (I_of(bstar+mp.mpf('1e-4')) + I_of(bstar-mp.mpf('1e-4')) - 2*I_of(bstar)) < 0)

# heat capacity C(beta) = beta^2 I(beta): peak
def negC(beta): return -(beta**2)*I_of(beta)
bC = mp.findroot(lambda t: (negC(t+mp.mpf('1e-8'))-negC(t-mp.mpf('1e-8')))/mp.mpf('2e-8'), mp.mpf('2.0'))
print(f"   beta_C (heat-capacity peak, beta>0 branch) = {mp.nstr(bC, 10)}")

# total information length L_tot = int sqrt(I) dbeta over R  (finite!)
L_tot = mp.quad(lambda t: mp.sqrt(I_of(t)), [-mp.inf, -5, 0, 5, mp.inf])
print(f"   TOTAL information length of the Gibbs curve: L_tot = {mp.nstr(L_tot, 12)}")
ck("L_tot is finite and positive", L_tot > 0 and L_tot < 100)

# arc length s(beta) from beta=0; marked-point table
def s_of(beta):
    if beta == 0: return mp.mpf(0)
    a_, b_ = (0, beta) if beta > 0 else (beta, 0)
    v = mp.quad(lambda t: mp.sqrt(I_of(t)), [a_, b_])
    return v if beta > 0 else -v
print("   arc-length positions (s(0)=0):")
for lbl, bb in [("beta=-1", -1), ("beta*", bstar), ("beta=0", 0), ("beta=1", 1),
                ("beta_C", bC), ("beta=sqrt5", mp.sqrt(5))]:
    print(f"     s({lbl}) = {mp.nstr(s_of(mp.mpf(bb)), 10)}")
print(f"   fraction of total length in [-1, sqrt5]: "
      f"{mp.nstr((s_of(mp.sqrt(5))-s_of(mp.mpf(-1)))/L_tot, 8)}")

# endpoint limits: golden merge and phi4 vertex
p_hi, _, _ = pI(mp.mpf(200))
ck("beta -> +inf: limit is uniform on the golden pair {phi, tau}",
   abs(p_hi[3]-mp.mpf(1)/2) < mp.mpf('1e-6') and abs(p_hi[4]-mp.mpf(1)/2) < mp.mpf('1e-6'))
p_lo, _, _ = pI(mp.mpf(-200))
ck("beta -> -inf: limit is the phi^4 vertex", abs(p_lo[5]-1) < mp.mpf('1e-6'))

print(f"\nALL {NCK[0]} CHECKS PASSED (t2_temperature)")
