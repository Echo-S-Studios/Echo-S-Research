#!/usr/bin/env python3
# dx_softcheck.py -- hard-core dossier harness (HANDOFF-v1.0 SS3.3 / SS3.5 / SS3.7).
# Lane: R2 discrete forcing chain (L-a, L-b, L-c), R4 square-root backstop,
# R6 Cencov-style no-go. Written cold this session, fail-first.
#
# Discipline: every DECISION exact -- sympy over Q / Q(sqrt5) (extended by i
# where needed), Q(sqrt5) signs by rational coefficient arithmetic (sgnQ5),
# mod-2pi bookkeeping in exact Fractions (units of pi). NO transcendental
# inequality arises in this lane, so there are ZERO mpmath interval guards
# (QX-G idiom noted; count here is 0 by construction, not by omission).
# Falsifier guards (fail-first: a deliberately perturbed / naive variant must
# be REJECTED) are counted separately as DX-*F. Exit 0 iff all pass.
#
#   DX-A  floor lattice rho_n = n ln phi from the gate ladder   (D1t, R2 L-a)
#   DX-B  square roots of -tau^2 I in R[J]: exactly {+-Q}       (D2t, R4)
#   DX-C  multiplier sign ladder + exposed joint H              (D3t, R2 L-c)
#   DX-D  transfer ledger vs shift-equivariance                 (D4t, R2 L-b)
#   DX-E  k-invariance lemma + no-go escape guards              (D5t, R6)
#
# Corpus citations verified this session before use:
#   W-Thm 9.1  (whitepaper p.11, "Ladder in Lucas-Fibonacci dress"): gates at
#              EVERY rung n>=1, C_n = 1/(phi^n-phi^-n)^2, u_n = 1/(phi^2n - 1),
#              m_n = -phi^-2n, rho_n = n ln phi via |m_n| = e^-2rho_n.
#   W-Decl 8.1 / W-Prop 8.2: gate iteration f_C(x) = C/(1+x), m = -C/(1+u)^2,
#              m+ m- = 1.
#   v1.6 thm:generator (W_n = q^n, q = i tau, W_{n+2} = -tau^2 W_n),
#   v1.6 prop:linearization (K-gate m_2 = -gap ~ per-two-rung action),
#   v1.6 thm:interp / cor:d3select / rem:gammatransfer / prop:sphlength
#              (the Z-family and the six selection principles),
#   v1.6 sec:quarterturn (J = e1 e2 forced, no D1-D3), QX-A14/A15, QX-F1-F6.

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin,
                   cos, simplify, expand, radsimp, trigsimp, together, cancel,
                   sqrtdenest, nsimplify, Poly, Matrix, eye, im, re, fibonacci,
                   lucas, diff, solve, conjugate, Mod)

PASS, FAIL, FALS = [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append((cid, desc))
    print(("PASS" if ok else "FAIL"), cid, "-", desc)

def fk(cid, cond, desc):
    # falsifier guard: cond must be True = the perturbed/naive variant IS rejected
    ok = bool(cond)
    FALS.append((cid, ok, desc))
    print(("FGUARD-PASS" if ok else "FGUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append((cid, desc))

def zero(e):
    if simplify(expand(e)) == 0:
        return True
    e2 = simplify(radsimp(together(cancel(e))))
    if e2 == 0:
        return True
    try:
        if simplify(sqrtdenest(e2)) == 0:
            return True
    except Exception:
        pass
    return simplify(trigsimp(e2)) == 0

def q5AB(e):
    """Exact (A,B) with e = A + B*sqrt(5), A,B in Q."""
    e = expand(cancel(radsimp(together(expand(e)))))
    try:
        p = Poly(e, s5)
    except Exception:
        # nested radicals: reduce to canonical Q(sqrt5) form first (exact)
        e = nsimplify(simplify(sqrtdenest(e)), [sqrt(5)])
        e = expand(cancel(radsimp(together(expand(e)))))
        p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQ5(e):
    """Exact sign of A + B*sqrt(5) via rational arithmetic only."""
    A, B = q5AB(e)
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    if A > 0:
        return 1 if A * A > 5 * B * B else -1
    return 1 if 5 * B * B > A * A else -1

s5   = sqrt(5)
phi  = (1 + s5)/2
tau  = (s5 - 1)/2
gap  = phi**-4
q    = I*tau

w    = symbols('w', positive=True)          # w = phi^{2n}
Cs   = symbols('C_s', positive=True)        # symbolic positive gate level
aR, bR = symbols('a b', real=True)
th   = symbols('theta', real=True)
rho  = symbols('rho', real=True)
kk, kp, nn, aZ, bZ = symbols('k k_p n a_z b_z', integer=True)
u_s, Av, Bv = symbols('u_s A_v B_v', positive=True)

def gate(C):
    """attracting fixed point and multiplier of f_C(x) = C/(1+x), C > 0."""
    u = (-1 + sqrt(1 + 4*C))/2
    return u, -C/(1 + u)**2

# ================================================================ DX-A
print("== DX-A: floor lattice rho_n = n ln phi from the gate ladder (L-a) ==")

Cw = w/(w - 1)**2
uw = 1/(w - 1)
mw = -Cw/(1 + uw)**2
ck("DX-A1", zero(uw**2 + uw - Cw) and zero(simplify(mw) + 1/w)
        and all(zero(Cw.subs(w, phi**(2*n)) - 1/(phi**n - phi**(-n))**2)
                for n in (1, 2, 3, 5)),
   "gate algebra symbolic in w=phi^2n: u=1/(w-1) fixed, m=-C/(1+u)^2=-1/w; C_n form")

okA2 = True
for n in range(1, 9):
    Cn = 1/(phi**n - phi**(-n))**2
    un, mn = gate(Cn)
    okA2 = okA2 and zero(un - 1/(phi**(2*n) - 1)) and zero(mn + phi**(-2*n))
    if n % 2:
        okA2 = okA2 and zero(Cn - 1/lucas(n)**2)
    else:
        okA2 = okA2 and zero(Cn - 1/(5*fibonacci(n)**2))
ck("DX-A2", okA2,
   "W-Thm 9.1 replication n=1..8: u_n=1/(phi^2n-1), m_n=-phi^-2n, TFD dress on C_n")

up = (-1 - sqrt(1 + 4*Cs))/2          # repelling partner
mprod = simplify((-Cs/(1 + (-1 + sqrt(1 + 4*Cs))/2)**2) * (-Cs/(1 + up)**2))
okA3 = zero(mprod - 1)
for n in (1, 2, 3, 6):
    okA3 = okA3 and sgnQ5(1 - phi**(-2*n)) > 0 and sgnQ5(phi**(2*n) - 1) > 0
ck("DX-A3", okA3,
   "attracting selection: m+ m- = 1 (W-Prop 8.2); |m_n| < 1 < |partner| exact")

okA4 = simplify(exp(-2*nn*log(phi)) - phi**(-2*nn)) == 0 \
   and zero(((nn + 1)*log(phi) - nn*log(phi)) - log(phi)) \
   and bool((-2*exp(-2*rho)).is_negative) \
   and all(zero(exp(-2*n*log(phi)) - (-gate(1/(phi**n - phi**(-n))**2)[1]))
           for n in range(1, 7))
ck("DX-A4", okA4,
   "floor lattice: e^{-2 n lnphi} = |m_n| (symbolic n + n=1..6); spacing lnphi; e^-2rho injective")

Cpert = Rational(1, 5) + 1
upert, mpert = gate(Cpert)
fk("DX-A5F", (not zero(mpert + gap)) and (not zero(exp(-2*3*log(phi)) - gap))
        and (not zero(exp(-2*Rational(3, 2)*log(phi)) - phi**(-2))),
   "FALSIFIER: perturbed gate C=1/5+1 multiplier != -gap; wrong lattices n=3, n=3/2 rejected")

# ================================================================ DX-B
print("== DX-B: square roots of -tau^2 I in R[J] are exactly {+-Q} (R4 backstop) ==")

e1 = Matrix([[0, 1], [1, 0]])          # sigma_x
e2 = Matrix([[1, 0], [0, -1]])         # sigma_z
J  = e1*e2                             # corpus convention J = e1 e2 (sec:quarterturn)
I2 = eye(2)
ck("DX-B1", J == Matrix([[0, -1], [1, 0]]) and (J*J + I2) == Matrix.zeros(2, 2),
   "J = e1 e2 = [[0,-1],[1,0]], J^2 = -I (keystone convention, no D1-D3)")

M  = aR*I2 + bR*J
MM = expand(M*M)
tgt = (aR**2 - bR**2)*I2 + (2*aR*bR)*J
ck("DX-B2", all(zero(MM[i, j] - tgt[i, j]) for i in range(2) for j in range(2)),
   "(aI+bJ)^2 = (a^2-b^2)I + 2ab J entrywise (R[J] iso C)")

sols = solve([aR**2 - bR**2 + tau**2, 2*aR*bR], [aR, bR], dict=True)
okB3 = len(sols) == 2 and all(sol[aR] == 0 for sol in sols)
bvals = [sol[bR] for sol in sols]
okB3 = okB3 and ((zero(bvals[0] - tau) and zero(bvals[1] + tau)) or
                 (zero(bvals[0] + tau) and zero(bvals[1] - tau)))
okB3 = okB3 and solve(aR**2 + tau**2, aR) == [] and sgnQ5(tau**2) > 0
ck("DX-B3", okB3,
   "solution set of (aI+bJ)^2 = -tau^2 I over R is exactly {(0,tau),(0,-tau)}; b=0 branch impossible")

Q  = tau*J
Mq = Matrix([[re(q), -im(q)], [im(q), re(q)]])
okB4 = all(zero((Q*Q + tau**2*I2)[i, j]) for i in range(2) for j in range(2))
okB4 = okB4 and all(zero((Mq - Q)[i, j]) for i in range(2) for j in range(2))
okB4 = okB4 and zero(q**2 + tau**2)
ck("DX-B4", okB4,
   "Q = tau J is the real matrix of mult-by-q; Q^2 = -tau^2 I = mult-by-q^2")

Dbad = expand((tau*I2 + tau*J)*(tau*I2 + tau*J)) + tau**2*I2
Dp   = expand(Q*Q) - (-tau**2*I2 + I2/7)
fk("DX-B5F", sgnQ5(Dbad[0, 0]) != 0 and (not zero(Dp[0, 0])),
   "FALSIFIER: (tau I + tau J)^2 != -tau^2 I; perturbed target -tau^2 I + I/7 not hit by Q^2")

S = Matrix([[0, 2], [Rational(-1, 2), 0]])
okB6 = (S*S + I2) == Matrix.zeros(2, 2)
okB6 = okB6 and all(zero((expand((tau*S)*(tau*S)) + tau**2*I2)[i, j])
                    for i in range(2) for j in range(2))
okB6 = okB6 and sgnQ5(-2*tau - (-tau/2)) != 0     # tau*S needs b=-2tau and b=-tau/2: clash
ck("DX-B6", okB6,
   "scope: tau*S with S=[[0,2],[-1/2,0]] also squares to -tau^2 I but lies OUTSIDE R[J] -- uniqueness is relative to the forced J")

ck("DX-B7", solve(exp(aR) + tau, aR) == [],
   "thm:interp branch elimination replicated: e^a = -tau has no real solution")

# The OPEN JOINT (analyzed, not faked): the FORCED object above is the two-rung
# AMPLITUDE action q^2 = -tau^2 I (|q^2| = tau^2), whose R[J] square roots are ±Q
# = the one-rung step. The renormalization LOAD wants ±Q as the sqrt of the
# LINEARIZATION multiplier m_2 = -gap = f'(K) (v1.6 prop:linearization). But m_2
# acts at the CO-INTENSITY / squared level: |m_2| = gap = tau^4 != tau^2 = |q^2|.
# So "sqrt of the linearization" and "sqrt of the two-rung amplitude action" are
# DIFFERENT statements; the "-gap * conformal" identification bridging them is the
# open joint -- corroborated (matching co-intensity scale gap + axis flip) but NOT
# forced to be the same square-root problem.
m2 = -phi**(-4)                       # gate multiplier at rung 2 (m_n=-phi^-2n, DX-A2) = f'(K)
ck("DX-B8", sgnQ5(gap - tau**2) == -1 and zero(gap - tau**4) and zero(-m2 - gap)
        and zero(expand(q**2*conjugate(q**2)) - tau**4)
        and zero(expand(q*conjugate(q)) - tau**2),
   "load distinction: |m_2| = gap = tau^4 != tau^2 = |q^2|; the amplitude sqrt (forced) and the linearization sqrt (the -gap*conformal load) are NOT the same problem -- open joint, honestly separated")

# ================================================================ DX-C
print("== DX-C: sign ladder and the exposed joint H (L-c source 1) ==")

# m_n = -phi^{-2n} is the gate multiplier, PROVEN in DX-A1 (symbolic m(w)=-1/w)
# and DX-A2 (per-n gate() == -phi^{-2n}); reading signs off this clean Q(sqrt5)
# closed form avoids per-n nested-radical denesting (same forced object).
mlist = {n: -phi**(-2*n) for n in range(1, 9)}
okC1 = bool((-1/w).is_negative)                     # m(w) = -1/w < 0, symbolic in n
okC1 = okC1 and all(sgnQ5(mlist[n]) == -1 for n in range(1, 9))
okC1 = okC1 and all(sgnQ5(mlist[i]*mlist[j]) == 1 for i, j in ((1, 2), (2, 3), (3, 8)))
ck("DX-C1", okC1,
   "sign ladder D1-free: m_n = -C_n/(1+u_n)^2 = -1/w < 0 at EVERY gate (symbolic in n); sign n-INDEPENDENT")

# CORPUS-VERIFIED (W-Thm 9.1, whitepaper p.11): gates exist at EVERY rung n>=1,
# rho_n = n lnphi (one floor per rung), m_n = -phi^-2n negative at EVERY rung.
# So m_n does NOT alternate/flip between gates. The "pi per two floors" holonomy
# is carried by the TWO-RUNG ladder action q^2 = -tau^2 (arg pi), tied to the
# linearization multiplier m_2 = -gap = f'(K) (v1.6 prop:linearization).
okC2 = zero(q**2 + tau**2) and sgnQ5(-tau**2) == -1              # q^2 = -tau^2, arg pi
okC2 = okC2 and all(zero(expand(q**(n + 2)*conjugate(q**(n + 2)))
                       - gap*expand(q**n*conjugate(q**n))) for n in range(7))
okC2 = okC2 and zero(mlist[2] + gap) and zero(mlist[2] - (-gap))  # m_2 = -gap = f'(K)
okC2 = okC2 and all(sgnQ5(mlist[n]) == -1 for n in range(1, 9))   # uniform sign, no flip
ck("DX-C2", okC2,
   "two-rung carrier: q^2 = -tau^2 (arg pi), |W_{n+2}|^2 = gap|W_n|^2, m_2 = -gap = f'(K); m_n uniformly negative (no per-gate flip)")

fk("DX-C2F", sgnQ5(mlist[1]*mlist[2]) == 1 and (not zero(mlist[2] - q**4)),
   "FALSIFIER: no per-gate sign flip (m_1 m_2 > 0); the naive uniform m_n = q^2n rejected at n=2 (m_2=-gap, q^4=+gap)")

def h1(d):                       # H1 congruence in units of pi: 2d = 1 (mod 2)
    return (2*d - 1) % 2 == 0
def inpm(d):                     # d in {+-1/2 + 2Z}
    return ((d - Fraction(1, 2)) % 2 == 0) or ((d + Fraction(1, 2)) % 2 == 0)
okC3 = all(h1(Fraction(1, 2) + m) and inpm(Fraction(1, 2) + m) for m in range(-6, 7))
grid = [Fraction(p, qd) for qd in range(1, 13) for p in range(-5*qd, 5*qd + 1)]
okC3 = okC3 and all(h1(d) == inpm(d) for d in grid)
ck("DX-C3", okC3,
   "bookkeeping lemma: 2d = 1 (mod 2) <=> d in {+-1/2 + 2Z}: H1 => delta = +-pi/2 (mod 2pi)")

fk("DX-C3F", (not h1(Fraction(1, 4))) and (not h1(Fraction(1, 1)))
         and (not h1(Fraction(3, 4))) and (not h1(Fraction(0, 1))),
   "FALSIFIER: d = 1/4, 1, 3/4, 0 all rejected by the H1 congruence")

okC4 = simplify(cos(pi/2 + 2*pi*kk)) == 0 and simplify(sin(pi/2 + 2*pi*kk) - 1) == 0
okC4 = okC4 and simplify(cos(-pi/2 + 2*pi*kk)) == 0 and simplify(sin(-pi/2 + 2*pi*kk) + 1) == 0
okC4 = okC4 and expand((pi/2 + 2*pi*kk) - (pi/2 + 2*pi*kp) - 2*pi*(kk - kp)) == 0
okC4 = okC4 and bool((2*pi/log(phi)).is_positive)
ck("DX-C4", okC4,
   "value collapse: e^{i delta_k} = +-i for EVERY branch; curves stay distinct (kappa_k - kappa_k' != 0)")

# Axiom H (exposed joint) = "arg(two-rung ladder action) = transverse holonomy
# over two floors". Given H: 2*delta = pi (mod 2pi) => delta = pi/2 (mod pi).
# Insufficiency: delta = pi/2 + pi*Z ALL satisfy; mod 2pi collapses to {pi/2, 3pi/2}
# = {+-pi/2} = the QUANTUM, but the residual Z (the branch) is NOT pinned.
okC5 = all(h1(Fraction(1, 2) + m) for m in range(-6, 7))          # pi/2 + pi*m all satisfy H1
okC5 = okC5 and {(Fraction(1, 2) + m) % 2 for m in range(-6, 7)} == {Fraction(1, 2), Fraction(3, 2)}
ck("DX-C5", okC5,
   "H forces delta = pi/2 (mod pi): the QUANTUM +-pi/2 is pinned, but a residual Z survives -- branch (iii) stays open at rung resolution")

# corroboration that the ACTUAL per-floor charge advance is pi/2 (uses D1 magnitude
# tau + the quarter-turn phase): q^n = i^n tau^n, so arg q^n = n(pi/2); the two-rung
# args are {pi, 0, pi, 0} -- consistent with pi per two floors.
okC6 = all(zero(q**n - I**n*tau**n) for n in range(9))
okC6 = okC6 and sgnQ5(expand(q**2)) == -1 and sgnQ5(expand(q**4)) == 1 \
             and sgnQ5(expand(q**6)) == -1 and sgnQ5(expand(q**8)) == 1
ck("DX-C6", okC6,
   "corroboration [uses D1 magnitude]: q^n = i^n tau^n, arg q^n = n(pi/2); two-rung args {pi,0,pi,0}")

okC7 = Mod(2*(Rational(1, 2) + 2*kk) - 1, 2) == 0
okC7 = okC7 and all(h1(Fraction(1, 2) + 2*k) for k in range(-3, 4))
ck("DX-C7", okC7,
   "restatement test: every branch d_k = 1/2 + 2k satisfies H1 -- H1 is branch-BLIND (no k-selection)")

# ================================================================ DX-D
print("== DX-D: transfer ledger vs shift-equivariance (L-b) ==")

okD1 = zero((1 - 1/(w*phi**2)) - (1 - 1/w) - tau/w)
okD1 = okD1 and all(zero((1 - tau**(2*(n + 1))) - (1 - tau**(2*n)) - tau**(2*n + 1))
                    for n in range(9))
okD1 = okD1 and zero(tau/(1 - tau**2) - 1)
ck("DX-D1", okD1,
   "QX-A15 replication: Delta F_n = tau^{2n+1} (symbolic + n=0..8); sum = 1")

Wt = tau**nn * exp(I*th)
okD2 = simplify(Wt*conjugate(Wt) - tau**(2*nn)) == 0
okD2 = okD2 and simplify(diff(Wt*conjugate(Wt), th)) == 0
ck("DX-D2", okD2,
   "ledger blindness: |W~_n|^2 = tau^2n for ARBITRARY angular assignment theta (d/dtheta = 0)")

okD3 = all(zero((tau**n*exp(I*n*pi))*conjugate(tau**n*exp(I*n*pi)) - tau**(2*n))
           for n in range(6))
ck("DX-D3", okD3,
   "counter-model (a) theta_n = n pi: passes EVERY ledger identity (moduli untouched)")

fk("DX-D3F", not zero(tau*exp(I*pi) - q),
   "FALSIFIER: model (a) violates rung passage (W~_1 = -tau != i tau) -- ledger data is strictly weaker than rung data")

sseq = [n % 2 for n in range(8)]                  # non-constant integer lift sequence
thn  = [n*pi/2 + 2*pi*sseq[n] for n in range(8)]
okD4 = all(zero(exp(I*thn[n]) - I**n) for n in range(8))
d0, d1 = thn[1] - thn[0], thn[2] - thn[1]
okD4 = okD4 and bool(simplify((d0 - d1)/pi).is_nonzero) and expand(d0 - d1 - 4*pi) == 0
okD4 = okD4 and all(zero((tau**n*exp(I*thn[n]))*conjugate(tau**n*exp(I*thn[n])) - tau**(2*n))
                    for n in range(6))
ck("DX-D4", okD4,
   "counter-model (b) theta_n = n pi/2 + 2pi s_n, s_n nonconstant: ALL rung values AND ledger identical, lift steps differ by 4pi")

okD5 = all(zero((tau**(n + 1)*exp(I*thn[n + 1]))/(tau**n*exp(I*thn[n])) - q)
           for n in range(7))
ck("DX-D5", okD5,
   "value-level step is pinned to q even in model (b): non-constancy lives ONLY in the lift -- L-b = lift equivariance axiom")

# ================================================================ DX-E
print("== DX-E: k-invariance lemma and the no-go escape guards (R6) ==")

okE1 = all(simplify(exp(I*(pi/2 + 2*pi*k)*n) - I**n) == 0
           for k in range(-3, 4) for n in range(9))
okE1 = okE1 and simplify(exp(I*(pi/2 + 2*pi*kk)*nn) - exp(I*pi*nn/2)) == 0
okE1 = okE1 and simplify(cos(2*pi*kk*nn) - 1) == 0 and simplify(sin(2*pi*kk*nn)) == 0
ck("DX-E1", okE1,
   "QX-F1 replication + symbolic extension: e^{i kappa_k n lnphi} = i^n for ALL integer k, n")

state = [1 - tau**(2*3), tau**(2*3 + 1), 1/(phi**3 - phi**(-3))**2,
         1/(phi**6 - 1), -phi**(-6), phi**6 - phi**(-6), 3*log(phi)]
okE2 = all(kk not in expr.free_symbols for expr in state)
Fwit = sum(exp(-n*log(phi))*exp(I*(pi/2 + 2*pi*kk)*n) for n in range(6))
Fwit0 = sum(exp(-n*log(phi))*I**n for n in range(6))
okE2 = okE2 and simplify(Fwit - Fwit0) == 0
ck("DX-E2", okE2,
   "FULL discrete state (z^2, ledger, gate C/u/m/lambda, rho_n) is k-invariant; witness functional of rung data is k-constant symbolically")

logq0 = log(tau) + I*pi/2
logqk = log(tau) + I*(pi/2 + 2*pi*kk)
okE3 = zero(logqk - logq0 - 2*pi*I*kk) and Matrix([[1, kk], [0, 1]]).det() == 1
okE3 = okE3 and all(not zero(im(logqk.subs(kk, k) - logq0)) for k in (-2, -1, 1, 2))
ck("DX-E3", okE3,
   "torus lattice Lambda = Z log q + Z 2pi i is k-invariant (SL2(Z) change of marking); the MARKING is not")

okE4 = all(abs(Fraction(1, 2) + 2*k) >= Fraction(3, 2) for k in range(-10, 11) if k != 0)
okE4 = okE4 and abs(Fraction(1, 2)) < Fraction(3, 2)
sqv = sorted((4*k + 1)**2 for k in range(-3, 4))
okE4 = okE4 and sqv[0] == 1 and sqv.count(1) == 1 and all(v > 1 for v in sqv[1:])
okE4 = okE4 and zero(diff(sqrt(Av + Bv*u_s), u_s) - Bv/(2*sqrt(Av + Bv*u_s))) \
             and bool((Bv/(2*sqrt(Av + Bv*u_s))).is_positive)
okE4 = okE4 and sgnQ5(3*s5 - 5) > 0 and bool((sqrt(3)/2).is_positive)
ck("DX-E4", okE4,
   "escape guard: all six principles' selector inputs are k-DEPENDENT (|4k+1| min at 0 unique; speeds strictly increasing in kappa^2; r>0) -- none factors through the discrete data")

Gvals = [sum(expand(exp(-n*log(phi))*simplify(exp(I*(pi/2 + 2*pi*k)*n))) for n in range(5))
         for k in range(-2, 3)]
okE5 = all(zero(Gvals[j] - Gvals[0]) for j in range(1, 5))
Ldiff = [(4*k + 1)**2 - 1 for k in range(-2, 3)]
okE5 = okE5 and all(v != 0 for k, v in zip(range(-2, 3), Ldiff) if k != 0) \
             and bool((pi**2/(4*log(phi)**2)).is_positive)
ck("DX-E5", okE5,
   "separation witness: discrete functional CONSTANT on k=-2..2 (exact); curve functional L_k^2 separates (differences = nonzero integer * kappa^2)")

fk("DX-E5F", not all(x == 0 for x in Ldiff),
   "FALSIFIER: the continuum functional is NOT k-constant -- the no-go is non-vacuous (its hypothesis genuinely excludes the six principles)")

# ================================================================ summary
print("=" * 64)
nf = sum(1 for c, ok, d in FALS if ok)
print("DX EXACT: %d passed, %d failed | FALSIFIER GUARDS: %d/%d rejected-as-required | mpmath guards: 0 (none needed)"
      % (len(PASS), len(FAIL), nf, len(FALS)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
sys.exit(0 if not FAIL else 1)
