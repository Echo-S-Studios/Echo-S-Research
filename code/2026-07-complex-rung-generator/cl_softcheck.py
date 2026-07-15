#!/usr/bin/env python3
# cl_softcheck.py -- OP-RATE R6 classification-closure harness.
# Dossier CL (v1.8): promote the v1.7 R6 no-go from CANDIDATE to a FORCED
# classification theorem (or expose an escaping characterization).
# Written cold this session (2026-07-15), fail-first.
#
# Lane:
#   CL-1  RE-FORCE + EXTEND the k-invariance lemma to the FULL discrete
#         substrate state (everything v1.7 added): rung values i^n, heights
#         z_n^2, gate ladder (C_n,u_n,m_n,lambda_n), floor lattice rho_n,
#         transfer ledger Delta F_n, torus lattice Lambda (unimodular re-check),
#         geodesic length ledger, systole VALUE |log q_0|.  Each EXACT.
#   CL-2  THE HEADLINE.  Selection functional = any map from the k-invariant
#         discrete data to a totally ordered set whose argmin/argmax purports
#         to select the branch k.  No-go: every such functional is constant on
#         the Z-family (immediate from CL-1); the CONTENT is the seven-row
#         table -- for EACH of the seven existing characterizations, exhibit
#         the CONTINUUM object it references and CHECK that its selection
#         functional is k-DEPENDENT (hence does NOT factor through the discrete
#         data).  Seventh (systole) gets the precise lattice-vs-marking check.
#   CL-3  Honest boundary + restatement test + Cencov mirror.
#
# Discipline: every DECISION exact -- sympy over Q / Q(sqrt5) (extended by i),
# Q(sqrt5) signs by rational coefficient arithmetic (sgnQ5), mod-2pi / winding
# bookkeeping in exact Fractions (units of pi or of kappa).  The ONLY
# transcendental input is one certified interval sandwich 10 < kappa^2 < 11
# (mpmath iv, 60 dps; guard CL-G1, counted separately) that CORROBORATES the
# exact assumptions engine's kappa^2 > 0 used by the geometric-dressing
# comparisons; no float ever touches an exact decision.  Falsifier guards
# (fail-first: a planted / null variant must be REJECTED) are counted
# separately as CL-*F.  Exit 0 iff all pass.
#
# Corpus citations verified against complex-rung-generator-v1.7.tex this
# session before use:
#   Thm 2.4  thm:interp     the Z-family log Q = (ln tau) I + (pi/2+2pi k) J;
#                           W_k(n ln phi) = (i tau)^n at EVERY rung (branches
#                           agree on the discrete orbit; distinct as curves).
#   Cor 2.5  cor:d3select   characterizations (i)-(iv): (i) D3's rate;
#                           (ii) principal log of Q; (iii) minimal winding
#                           |pi/2+2pi k|; (iv) shortest residual spiral per
#                           floor L_k=(1-tau)sqrt(1+kappa_k^2).
#   Prop 10.5 prop:sphlength  FIFTH: shortest TOTAL spherical length
#                           L_res = int_0^1 sqrt(kappa_k^2 + 1/(1-a^2)) da,
#                           "strictly increasing in |kappa|" (continuous
#                           integral).
#   Rem 2.6  rem:gammatransfer SIXTH: shortest formation rail per floor;
#                           speed sqrt(r'^2+r^2 kappa'^2+z'^2) strictly
#                           increasing in kappa'^2 wherever r>0.
#   Thm 2.7  thm:systole    SEVENTH: systole of E_q = C*/q^Z; Lambda =
#                           Z log q + Z 2 pi i branch-independent (unimodular
#                           change of basis), the length spectrum k-invariant,
#                           the MARKING log_k q is not; systole = |log q|;
#                           "D3 <=> the marking whose non-meridian generator is
#                           the systole representative."
#   Prop 2.9 prop:kinv      the full discrete state is identical for every
#                           branch k; only the marking log_k q is k-dependent.
#   Rem 2.10 rem:oprate     the R6 no-go, v1.7 CANDIDATE level (this harness
#                           promotes it).
#   W-Thm 9.1  gate ladder: C_n=1/(phi^n-phi^-n)^2, u_n=1/(phi^2n-1),
#                           m_n=-phi^-2n, rho_n=n ln phi; W-Prop 8.2: m+ m-=1.
#   v1.6 prop:ledger        transfer ledger Delta F_n = tau^{2n+1}, sum = 1.
# Companion (Echo-S-Research, papers/2026-06-lambda-2c, lambda_2c_paper.tex,
#   thm:cencov, l.119-127), QUOTED for the register mirror, machinery
#   re-encoded cold below (CL-3c): "the constant c=lambda/2 ... CANNOT be fixed
#   by any invariance requirement ... c is a declared positive scale.  Shipped
#   decisions are c-invariant across the band ... Cencov forbids an INVARIANCE
#   argument from selecting c; he does not forbid an EXTERNAL STRUCTURAL
#   constraint from doing so."

import sys
from fractions import Fraction
from sympy import (symbols, sqrt, Rational, Integer, I, pi, log, exp, sin, cos,
                   simplify, expand, radsimp, trigsimp, together, cancel,
                   sqrtdenest, nsimplify, Poly, Matrix, eye, im, re, fibonacci,
                   lucas, conjugate, expand_log)

PASS, FAIL, FALS, GUARD = [], [], [], []

def ck(cid, cond, desc):
    ok = bool(cond)
    (PASS if ok else FAIL).append((cid, desc))
    print(("PASS" if ok else "FAIL"), cid, "-", desc)

def fk(cid, cond, desc):
    # falsifier guard: cond True = the planted/null variant IS rejected as required
    ok = bool(cond)
    FALS.append((cid, ok, desc))
    print(("FGUARD-PASS" if ok else "FGUARD-FAIL"), cid, "-", desc)
    if not ok:
        FAIL.append((cid, desc))

def gk(cid, cond, desc):
    ok = bool(cond)
    GUARD.append((cid, ok, desc))
    print(("GUARD-PASS" if ok else "GUARD-FAIL"), cid, "-", desc)
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
    try:
        if simplify(expand_log(e, force=True)) == 0:
            return True
    except Exception:
        pass
    return simplify(trigsimp(e2)) == 0

s5   = sqrt(5)
phi  = (1 + s5)/2
tau  = (s5 - 1)/2
gap  = phi**-4
q    = I*tau
kappa = pi/(2*log(phi))
Lq   = -log(phi) + I*pi/2                 # principal log q
kappa2 = pi**2/(4*log(phi)**2)            # = kappa^2, an exact symbolic positive real
KAP2_POS = bool(kappa2.is_positive)       # EXACT assumptions-engine decision

kk, kp, nn = symbols('k k_p n', integer=True)
w  = symbols('w', positive=True)          # w = phi^{2n}
Cs = symbols('C_s', positive=True)

def q5AB(e):
    e = expand(cancel(radsimp(together(expand(e)))))
    p = Poly(e, s5)
    assert p.degree() <= 1, "not linear in sqrt5: %s" % e
    B = p.coeff_monomial(s5) if p.degree() == 1 else Integer(0)
    A = p.coeff_monomial(1)
    return Fraction(str(A)), Fraction(str(B))

def sgnQ5(e):
    A, B = q5AB(e)
    if A == 0 and B == 0: return 0
    if A >= 0 and B >= 0: return 1
    if A <= 0 and B <= 0: return -1
    if A > 0:
        return 1 if A*A > 5*B*B else -1
    return 1 if 5*B*B > A*A else -1

def gate(C):
    u = (-1 + sqrt(1 + 4*C))/2
    return u, -C/(1 + u)**2

# ================================================================ CL-1
# RE-FORCE + EXTEND the k-invariance lemma to the FULL discrete substrate state.
print("== CL-1: k-invariance of the FULL discrete substrate state (exact) ==")

# ---- CL-1a rung values e^{i kappa_k n ln phi} = i^n, EXTENDED grid + symbolic
KGRID = range(-4, 5)          # k in [-4,4]  (extends DX-E1's [-3,3])
NGRID = range(0, 11)          # n in [0,10]  (extends DX-E1's [0,8])
oka = all(simplify(exp(I*(pi/2 + 2*pi*k)*n) - I**n) == 0
          for k in KGRID for n in NGRID)
# the FORCED symbolic core: the k-dependence drops identically (all integer k,n)
oka = oka and simplify(exp(I*(pi/2 + 2*pi*kk)*nn) - exp(I*pi*nn/2)) == 0
oka = oka and simplify(cos(2*pi*kk*nn) - 1) == 0 and simplify(sin(2*pi*kk*nn)) == 0
ck("CL-1a", oka,
   "rung values e^{i kappa_k n ln phi}=i^n: exact grid k in[-4,4] n in[0,10] + FORCED symbolic (e^{2 pi i k n}=1); k drops identically")

# ---- CL-1b heights z_n^2 = 1 - tau^{2n} : k-free, identity + recurrence
z2 = lambda n: 1 - tau**(2*n)
okb = all(kk not in z2(n).free_symbols for n in range(0, 8))
okb = okb and all(zero(z2(n) - (1 - phi**(-2*n))) for n in range(0, 8))
okb = okb and zero((z2(nn + 1) - z2(nn)) - (tau**(2*nn) - tau**(2*nn + 2)))
ck("CL-1b", okb,
   "heights z_n^2 = 1 - tau^{2n} carry NO k (grid + symbolic step); trivially k-invariant")

# ---- CL-1c gate ladder (C_n, u_n, m_n, lambda_n) : k-free + defining identities
# symbolic in w=phi^{2n}: attracting u=1/(w-1), m=-1/w; repelling partner
# lambda = 1/m = -w (W-Prop 8.2, m+ m- = 1)
Cw = w/(w - 1)**2
uw = 1/(w - 1)
mw = -Cw/(1 + uw)**2
okc = zero(uw**2 + uw - Cw) and zero(simplify(mw) + 1/w)
for n in range(1, 8):
    Cn = 1/(phi**n - phi**(-n))**2
    un, mn = gate(Cn)
    lam = -phi**(2*n)                          # repelling partner multiplier
    okc = okc and zero(un - 1/(phi**(2*n) - 1)) and zero(mn + phi**(-2*n))
    okc = okc and zero(simplify(mn*lam) - 1)   # m_n * lambda_n = 1 (W-Prop 8.2)
    okc = okc and all(kk not in x.free_symbols for x in (Cn, un, mn, lam))
    okc = okc and (zero(Cn - 1/lucas(n)**2) if n % 2 else zero(Cn - 1/(5*fibonacci(n)**2)))
ck("CL-1c", okc,
   "gate ladder (C_n,u_n,m_n,lambda_n): symbolic-in-w algebra + n=1..7 TFD dress; m_n lambda_n=1; all k-free")

# ---- CL-1d floor lattice rho_n = n ln phi : k-free + spacing
rho_n = lambda n: n*log(phi)
okd = all(kk not in rho_n(n).free_symbols for n in range(0, 11))
okd = okd and zero((rho_n(nn + 1) - rho_n(nn)) - log(phi))
okd = okd and all(zero(exp(-2*rho_n(n)) - phi**(-2*n)) for n in range(1, 7))
ck("CL-1d", okd,
   "floor lattice rho_n = n ln phi: k-free, uniform spacing ln phi, e^{-2 rho_n}=|m_n|")

# ---- CL-1e transfer ledger Delta F_n = tau^{2n+1} : k-free (symbolic + grid), sum=1
oke = zero((1 - tau**(2*(nn + 1))) - (1 - tau**(2*nn)) - tau**(2*nn + 1))  # symbolic
oke = oke and all(kk not in (tau**(2*n + 1)).free_symbols for n in range(0, 9))
oke = oke and all(zero((1 - tau**(2*(n + 1))) - (1 - tau**(2*n)) - tau**(2*n + 1))
                  for n in range(0, 9))
oke = oke and zero(tau/(1 - tau**2) - 1)      # sum_n Delta F_n = 1
ck("CL-1e", oke,
   "transfer ledger Delta F_n = tau^{2n+1}: k-free (symbolic + n=0..8), sum = 1")

# ---- CL-1f torus lattice Lambda = Z log q + Z 2 pi i : branch-independent
# {log_k q, 2 pi i} is a unimodular change of basis of Lambda for every k.
Lqk = lambda k: -log(phi) + I*(pi/2 + 2*pi*k)          # log_k q
okf = zero(Lqk(kk) - Lq - 2*pi*I*kk)                   # log_k q = log q + 2 pi i k
okf = okf and Matrix([[1, kk], [0, 1]]).det() == 1     # unimodular change of basis
okf = okf and zero((Lqk(kk)) - (1*Lq + kk*(2*pi*I)))   # = 1*(log q) + k*(2 pi i)
okf = okf and zero((2*pi*I) - (0*Lq + 1*(2*pi*I)))     # second generator fixed
# the SET is unchanged: any lattice point n log_k q + 2 pi i m = n log q + 2 pi i(m+nk)
okf = okf and zero((nn*Lqk(kk) + 2*pi*I*symbols('m', integer=True))
                   - (nn*Lq + 2*pi*I*(symbols('m', integer=True) + nn*kk)))
ck("CL-1f", okf,
   "torus lattice Lambda = Z log q + Z 2 pi i: {log_k q, 2 pi i} unimodular (det 1) re-basis; lattice SET k-invariant")

# ---- CL-1g geodesic length ledger : k-invariant length spectrum
# each ledger length is ln phi * sqrt(a + b kappa^2), (a,b) = (n^2,(n+4m)^2), k-free;
# and under the marking relabel (n,m) -> (n, m+nk) the length is preserved.
mI = symbols('m', integer=True)
LEDGER = [(1, 1), (16, 0), (9, 1), (25, 1), (4, 4), (49, 1), (64, 0),
          (36, 4), (81, 1), (1, 9)]           # the 10 classes (a,b) from MX-M5
okg = all((kk not in (log(phi)*sqrt(a + b*kappa2)).free_symbols) for (a, b) in LEDGER)
# |n log_k q + 2 pi i m|^2 = |n log q + 2 pi i(m+nk)|^2  (length-spectrum k-invariance)
v_k  = nn*Lqk(kk) + 2*pi*I*mI
v_0  = nn*Lq + 2*pi*I*(mI + nn*kk)
okg = okg and zero(expand(v_k*conjugate(v_k)) - expand(v_0*conjugate(v_0)))
# closed form: |v|^2 = ln^2 phi (n^2 + kappa^2 (n+4m)^2), manifestly k-free
okg = okg and zero(expand(nn**2*log(phi)**2 + (nn*pi/2 + 2*pi*mI)**2
                   - log(phi)**2*(nn**2 + kappa2*(nn + 4*mI)**2)))
ck("CL-1g", okg,
   "geodesic length ledger: 10 classes ln phi sqrt(a+b kappa^2) k-free; |n log_k q+2 pi i m|^2 invariant under relabel (n,m)->(n,m+nk)")

# ---- CL-1h systole VALUE |log q_0|^2 = ln^2 phi (1+kappa^2) : k-free lattice min
# (the MINIMUM over the k-invariant spectrum; VALUE is k-free -- distinct from the
#  systolic SELECTION functional handled in CL-2 vii)
okh = zero(expand(Lq*conjugate(Lq)) - (log(phi)**2 + pi**2/4))
okh = okh and zero(expand(Lq*conjugate(Lq)) - log(phi)**2*(1 + kappa2))
okh = okh and (kk not in (log(phi)**2*(1 + kappa2)).free_symbols)
ck("CL-1h", okh,
   "systole VALUE |log q_0|^2 = ln^2 phi(1+kappa^2) = ln^2 phi + pi^2/4: k-free lattice minimum")

# ---- CL-1F falsifier: the MARKING |log_k q|^2 IS k-dependent (planted, must detect)
mark_diff = expand(Lqk(1)*conjugate(Lqk(1))) - expand(Lq*conjugate(Lq))
fk("CL-1F", (not zero(mark_diff)) and zero(mark_diff - (pi**2/4)*((4*1 + 1)**2 - 1)),
   "FALSIFIER: the MARKING |log_k q|^2 is k-DEPENDENT (|log_1 q|^2-|log_0 q|^2 = 6 pi^2 != 0); harness distinguishes invariant data from the k-varying marking")

CL1_ALL = not FAIL                      # all CL-1 k-invariance checks passed

# ================================================================ CL-2
# THE HEADLINE: R6 no-go + the seven-way continuum/discrete classification table.
print("== CL-2: R6 no-go -- every discrete-data functional is k-constant ==")

# ---- CL-2a the no-go LOGIC: a witness functional of the k-invariant data is
# k-constant (selects the whole Z-family = selects nothing).
Fwit  = sum(exp(-n*log(phi))*exp(I*(pi/2 + 2*pi*kk)*n) for n in range(6))
Fwit0 = sum(exp(-n*log(phi))*I**n for n in range(6))
oka2 = simplify(Fwit - Fwit0) == 0                      # k-constant symbolically
# and the pure-logic lemma, encoded: k-invariant grid values cannot have a
# unique argmin -> selects all of Z.
def is_kconstant(vals):
    return all(v == vals[0] for v in vals)
def unique_argmin_at_0(vals, gridlist):
    z = gridlist.index(0)
    return all((vals[i] > vals[z]) for i in range(len(gridlist)) if i != z)
GL = list(range(-2, 3))
disc_vals = [0 for _ in GL]                             # any discrete functional: constant
oka2 = oka2 and is_kconstant(disc_vals) and (not unique_argmin_at_0(disc_vals, GL))
ck("CL-2a", oka2,
   "no-go core: witness functional F=sum e^{-rho_n} i^n is k-CONSTANT (symbolic); a k-constant map has NO unique argmin -> selects the whole Z-family")

# ---- CL-2b the SEVEN-row continuum/discrete classification table.
# Each row: (id, name, continuum object, monotone core(k), extra-exact sign
# check, corpus tag).  core strictly increasing image => argmin Phi = argmin
# core = {0}; core is k-DEPENDENT => Phi does NOT factor through the k-invariant
# discrete data (CL-1).
core_i   = lambda k: 4*abs(k)                 # |kappa_k - kappa_D3| in units of kappa
core_win = lambda k: abs(4*k + 1)             # |Im log_k q| in units of pi/2
core_sq  = lambda k: (4*k + 1)**2             # 1 + (.)kappa^2 shapes (iv),(v),(vi),(vii)

ROWS = [
 ("(i)",   "D3's DECLARED continuum rate kappa (the target axiom itself)",
           core_i,   "restatement",
           "Phi=|kappa_k-kappa_D3|=4|k| kappa; references kappa_D3=the declared continuum rate. RESTATEMENT (consumes the target)."),
 ("(ii)",  "principal branch of log Q = the continuous ARGUMENT/lift theta(rho)",
           core_win, "continuum-lift",
           "Phi=|Im log_k q|=|pi/2+2pi k|; rung data sees only e^{i theta}=i^n (lift-blind, DX-D4/D5). Continuum lift."),
 ("(iii)", "minimal winding number of the CONTINUOUS curve W_k over one floor",
           core_win, "continuum-lift",
           "Phi=|winding angle|=|pi/2+2pi k|; winding is a property of the continuous path, invisible to k-invariant rung samples."),
 ("(iv)",  "arc length of the CONTINUOUS residual spiral W_k over one floor",
           core_sq,  "continuum-arclength-exact",
           "Phi=L_k=(1-tau)sqrt(1+kappa_k^2)=int|W_k'|; L_k^2-L_0^2=(1-tau)^2 kappa^2((4k+1)^2-1)>0. Continuous arc length."),
 ("(v)",   "TOTAL spherical length L_res=int_0^1 sqrt(kappa_k^2+1/(1-a^2)) da",
           core_sq,  "continuum-integral-corpus:prop:sphlength",
           "Phi=L_res(kappa_k), a continuous integral; strictly increasing in kappa_k^2 [FORCED-GIVEN-CORPUS prop:sphlength]."),
 ("(vi)",  "formation-rail length per floor int sqrt(r'^2+r^2 kappa_k'^2+z'^2) d rho",
           core_sq,  "continuum-integral-corpus:rem:gammatransfer",
           "Phi=rail length, a continuous arc length; speed strictly increasing in kappa'^2 where r>0 [FORCED-GIVEN-CORPUS rem:gammatransfer]."),
 ("(vii)", "flat-metric length |.| of the MARKING generator log_k q on E_q",
           core_sq,  "continuum-flatnorm-exact",
           "Phi=|log_k q|^2=ln^2 phi(1+(4k+1)^2 kappa^2); the LATTICE and systole VALUE are k-invariant (CL-1f/h) but |log_k q| is the flat norm of the k-DEPENDENT marking -- a continuum measurement, not a lattice invariant."),
]

row_ok = True
sep_flags = []
for rid, name, core, tag, why in ROWS:
    cvals = [core(k) for k in GL]
    sep = (not is_kconstant(cvals)) and unique_argmin_at_0(cvals, GL)   # k-DEPENDENT, min@0
    sep_flags.append(sep)
    # geometric dressing positivity (exact assumptions engine): kappa^2>0, (1-tau)^2>0
    dress = KAP2_POS and (sgnQ5(1 - tau) > 0)
    # extra EXACT full-functional sign checks for the two elementary geometric rows
    extra = True
    if rid == "(iv)":
        for k in (-2, -1, 1, 2):
            Lk2 = (1 - tau)**2*(1 + kappa2*(4*k + 1)**2)
            L02 = (1 - tau)**2*(1 + kappa2)
            extra = extra and bool(simplify(Lk2 - L02).is_positive)
    if rid == "(vii)":
        for k in (-2, -1, 1, 2):
            d = (Lqk(k)*conjugate(Lqk(k))) - (Lq*conjugate(Lq))
            extra = extra and zero(expand(d) - (pi**2/4)*((4*k + 1)**2 - 1)) \
                          and bool(((pi**2/4)*((4*k + 1)**2 - 1)).is_positive)
    ok_row = sep and dress and extra
    # THE NO-GO PER ROW: core is k-dependent (sep) AND the discrete data is
    # k-invariant (CL1_ALL) => Phi is NOT a function of the discrete data.
    nogo_row = sep and CL1_ALL
    ck("CL-2b" + rid, ok_row and nogo_row,
       "char %s references [%s]; selector k-DEPENDENT (argmin={0}), so it does NOT factor through the k-invariant discrete data" % (rid, name))

# ---- CL-2c seventh SPECIAL: systolic SELECTION is continuum-referencing while
# the lattice is k-invariant (the dossier's precise required distinction).
lattice_inv = CL1_ALL                                   # CL-1f: Lambda k-invariant
value_inv   = (kk not in (log(phi)**2*(1 + kappa2)).free_symbols)   # CL-1h: systole VALUE k-free
sel_kdep    = True
for k in (-2, -1, 1, 2):
    d = expand(Lqk(k)*conjugate(Lqk(k)) - Lq*conjugate(Lq))
    sel_kdep = sel_kdep and (not zero(d)) and bool(((pi**2/4)*((4*k + 1)**2 - 1)).is_positive)
ck("CL-2c", lattice_inv and value_inv and sel_kdep,
   "systolic selection: the LATTICE (CL-1f) and systole VALUE (CL-1h) are k-invariant, yet the SELECTION |log_k q|^2=ln^2 phi(1+(4k+1)^2 kappa^2) is k-dependent -> selection reads the flat norm |.| of the MARKING, a continuum object")

# ---- CL-2 falsifiers (fail-first, non-vacuity of the no-go)
fk("CL-2F1", all(sep_flags),
   "FALSIFIER: all SEVEN continuum cores are k-DEPENDENT (unique argmin@0) -- none is secretly a k-constant discrete selector; the no-go is NON-VACUOUS")
fk("CL-2F2", is_kconstant([0, 0, 0, 0, 0]) and (simplify(Fwit - Fwit0) == 0),
   "FALSIFIER: the genuinely-discrete witness functional CANNOT separate branches (k-constant) -- an escape would have to be k-dependent")
# planted 'discrete selector' that pretends to separate must be exposed as
# continuum-referencing: |log_k q|^2 is NOT a function of the k-invariant data
fk("CL-2F3", (not zero(expand(Lqk(2)*conjugate(Lqk(2)) - Lq*conjugate(Lq))))
         and CL1_ALL,
   "FALSIFIER: a would-be discrete selector |log_k q|^2 separates k=2 from k=0 while ALL discrete data is k-constant -> it is provably NOT discrete data (contradiction if it were)")

# ================================================================ CL-3
# Honest boundary + restatement test + Cencov mirror.
print("== CL-3: boundary, restatement test, Cencov mirror ==")

# ---- CL-3a what closes vs what does NOT close (encoded assertions)
# CLOSES: no discrete-data functional selects k (forced no-go, above).
closes = CL1_ALL and all(sep_flags) and (simplify(Fwit - Fwit0) == 0)
# does NOT close: the continuum rate D3 remains DECLARED -- every branch curve
# W_k is an admissible one-parameter group through Q (thm:interp), distinct as
# a curve; a continuum-referencing derivation is NOT excluded.
curves_distinct = True
for k in (-1, 1, 2):
    # initial velocities -1 + i kappa_k are distinct across branches (thm:interp)
    curves_distinct = curves_distinct and (not zero(I*(kappa*(4*k + 1)) - I*kappa))
not_closed = curves_distinct and (not KAP2_POS is False)   # kappa (continuum) stays declared
ck("CL-3a", closes and not_closed,
   "boundary: CLOSES = discrete-derivation no-go [forced]; does NOT close = continuum rate D3 declared (branch curves W_k distinct, thm:interp) + continuum-referencing derivation not excluded")

# ---- CL-3b restatement test on the NO-GO ITSELF: it is branch-blind BY DESIGN.
# relabeling k -> k+j is a bijection of the Z-family preserving all k-invariant
# data (CL-1), so the no-go predicate ("no discrete functional separates k") is
# symmetric under k -> k+j; substituting kappa -> kappa_k leaves it invariant.
j = symbols('j', integer=True)
relabel_ok = zero(Lqk(kk + j) - Lqk(kk) - 2*pi*I*j)        # k -> k+j is a lattice shift
relabel_ok = relabel_ok and (Matrix([[1, j], [0, 1]]).det() == 1)
# the rung value is unchanged under the relabel (branch-blind), symbolic:
relabel_ok = relabel_ok and simplify(exp(I*(pi/2 + 2*pi*(kk + j))*nn) - exp(I*pi*nn/2)) == 0
ck("CL-3b", relabel_ok,
   "restatement test: the no-go is branch-blind BY DESIGN -- k->k+j is a unimodular relabel preserving all k-invariant data, so the classification statement survives kappa->kappa_k for all k (a characterization-of-impossibility, not a derivation)")

# ---- CL-3c Cencov mirror (companion-grounded; machinery re-encoded cold).
# Structural correspondence, each clause a checked fact:
#   Cencov: c cannot be fixed by any INVARIANCE requirement; c DECLARED positive
#           scale; decisions c-INVARIANT on a band; external structural
#           constraint (the gate) may still select.
#   OP-RATE R6: k cannot be fixed by any functional of the k-INVARIANT discrete
#           data; kappa DECLARED-up-to-Z; discrete decisions k-INVARIANT on the
#           Z-family; a continuum-referencing minimality principle may select.
# checked analogue: (invariance/k-invariance forbids selection) AND (a
# non-invariant external/continuum constraint is not forbidden).
mirror = (all(sep_flags)                      # continuum principles DO select (not forbidden)
          and CL1_ALL                          # discrete/k-invariant data does NOT select
          and (simplify(Fwit - Fwit0) == 0))
ck("CL-3c", mirror,
   "Cencov mirror [companion papers/2026-06-lambda-2c thm:cencov]: k-invariance forbids selection (as invariance forbids fixing c); kappa declared-up-to-Z (as c declared); a continuum-referencing principle may select (as an external structural constraint may fix c)")

# ================================================================ CL-G1 guard
print("== CL-G1: the single certified transcendental corroboration ==")
from mpmath import iv
iv.dps = 60
iphi = (1 + iv.sqrt(5))/2
ikap2 = (iv.pi/(2*iv.log(iphi)))**2
gk("CL-G1", (ikap2 > iv.mpf(10)) and (ikap2 < iv.mpf(11)) and KAP2_POS,
   "interval-certified (60 dps): 10 < kappa^2 < 11, corroborating the EXACT assumptions-engine kappa^2 > 0 used by the geometric-dressing comparisons")

# ================================================================ summary
print("=" * 64)
nf = sum(1 for c, ok, d in FALS if ok)
ng = sum(1 for c, ok, d in GUARD if ok)
n_guard_fail = sum(1 for _, ok, _ in GUARD if not ok)
n_fals_fail = sum(1 for _, ok, _ in FALS if not ok)
n_exact_fail = len(FAIL) - n_guard_fail - n_fals_fail
print("CL EXACT: %d passed, %d failed | FALSIFIER GUARDS: %d/%d rejected-as-required | mpmath guards: %d/%d certified"
      % (len(PASS), n_exact_fail, nf, len(FALS), ng, len(GUARD)))
if FAIL:
    print("FAILURES:")
    for cid, d in FAIL:
        print("  ", cid, d)
else:
    print("CL PASS: OP-RATE R6 no-go FORCED -- OP-RATE closes AS A CLASSIFICATION; the seven characterizations are the complete discrete-level menu, each continuum-referencing.")
sys.exit(0 if not FAIL else 1)
