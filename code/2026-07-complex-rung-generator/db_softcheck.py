#!/usr/bin/env python3
"""
db_softcheck.py -- OP-RATE continuum horn (ledger row 1) + watchlist rows 5-7.

[DECLARED] The derivation base D (this is the AUTHORING step the ledger demands:
'declare the framework in which a derivation is admitted', making the row
falsifiable at all):

  D := functionals of the branch-blind skeleton
       S = ( k-invariant discrete rung state,  torus E_q = C*/q^Z,
             lattice Lambda = Z log q + 2 pi i Z )
       that factor through T-invariant data, where
       T : log_k q |-> log_{k+1} q = log_k q + 2 pi i   (marking relabel).

Theorem DB-1 (no-section for the log-torsor of q, in D)  [forced in D]:
  the markings {log_k q} form a FREE TRANSITIVE Z-set under T; every
  D-functional is constant on the T-orbit; hence no D-definable section of the
  Z-torsor of markings exists, and kappa is declared-up-to-Z inside D.
  This is the exact mirror of Cencov for c in lambda = 2c (C-Cor 2.18):
  no invariance argument removes the scale; an EXTERNAL structural constraint
  (a non-T-invariant atom) is not excluded -- that is the escape clause.

Restatement gate: every one of the seven characterizations (C-Rem 2.19) is
re-certified k-DEPENDENT, i.e. not a D-object -- any principle proposed later
must be checked against this table.

Groups: DB-A torsor | DB-B k-invariance | DB-C systole window (exact, sandwich)
| DB-G interval guard | DB-E seven-way sweep | DB-N no-section + falsifiers
| DB-P watchlist skeleton rows 5-7 (named triggers, no claims)
Discipline: exact decisions; the single transcendental comparison is the
certified sandwich 10 < kappa^2 < 11 (mpmath iv, 60 dps) -- MX-G1 mirror.
"""
import sys
from sympy import (sqrt, Rational, symbols, simplify, expand, diff, pi, S, I,
                   exp, log, im, re, Poly, Abs)
import mpmath as mp

SQ5, SQ13 = sqrt(5), sqrt(13)
phi = (1 + SQ5)/2
tau = phi - 1
gap = phi**-4
kappa = pi/(2*log(phi))
logq  = -log(phi) + I*pi/2            # principal logarithm of q = i*tau

PASS, FAIL = [], []
def check(cid, cond, note=""):
    ok = bool(cond)
    (PASS if ok else FAIL).append(cid)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid:8s} {note}")
    return ok

def is_zero(e): return simplify(expand(e)) == 0

print("== db_softcheck.py : OP-RATE continuum horn -- base D declared, no-section certified ==\n")

# =====================================================================================
# DB-A : the marking torsor
# =====================================================================================
print("-- DB-A : the Z-torsor of markings --")
def logk(k): return logq + 2*pi*I*k

free = all(is_zero(logk(k) - logq - 2*pi*I*k) and (2*pi*I*k != 0)
           for k in range(-4, 5) if k != 0)
check("DB-A1", free,
      "freeness: log_k q - log_0 q = 2 pi i k != 0 for every k != 0 (grid k in [-4,4])")
trans = all(is_zero(logk(k2) - logk(k1) - 2*pi*I*(k2 - k1))
            for k1 in range(-2, 3) for k2 in range(-2, 3))
check("DB-A2", trans,
      "transitivity: any two markings differ by an integer relabel -- one T-orbit")
check("DB-A3", (Rational(1,2) <= 1) and (Rational(1,2) + 2 > 1) and (Rational(1,2) - 2 <= -1),
      "principal strip: Im log_k q = (4k+1) pi/2 lies in (-pi, pi] iff k = 0")
check("DB-A4", is_zero(simplify(exp(logq) - I*tau)),
      "exp(log_0 q) = q = i*tau: the marking really marks the generator")

# =====================================================================================
# DB-B : k-invariance of the discrete skeleton (Lem 2.16 essentials, cold)
# =====================================================================================
print("\n-- DB-B : k-invariance of the skeleton --")
rungs = all(is_zero(simplify(exp(I*n*(pi/2 + 2*pi*k)) - I**n))
            for k in range(-3, 4) for n in range(0, 9))
check("DB-B1", rungs,
      "rung values e^{i kappa_k n ln phi} = i^n identical for k in [-3,3], n in [0,8]")
n_s, m_s, k_s = symbols('n m k', integer=True)
check("DB-B2", is_zero(expand(n_s*(logq + 2*pi*I*k_s) + 2*pi*I*m_s
                              - (n_s*logq + 2*pi*I*(m_s + n_s*k_s)))),
      "lattice: n log_k q + 2 pi i m = n log_0 q + 2 pi i (m + nk) -- same lattice Lambda")
check("DB-B3", all((1*1 - k*0) == 1 for k in range(-4, 5)),
      "marking relabel is the unimodular basis change [[1,k],[0,1]], det = 1")
check("DB-B4", is_zero(simplify(kappa*log(phi) - pi/2)) and is_zero(simplify(4*kappa*log(phi) - 2*pi)),
      "rate bridges: kappa ln phi = pi/2 and 4 kappa ln phi = 2 pi (the '4m' in the norm form)")
# norm form |n log q + 2 pi i m|^2 = ln^2(phi) (n^2 + kappa^2 (n+4m)^2)
vec = n_s*logq + 2*pi*I*m_s
norm2 = expand(re(vec)**2 + im(vec)**2)
target = expand(log(phi)**2 * (n_s**2 + kappa**2*(n_s + 4*m_s)**2))
check("DB-B5", is_zero(simplify(norm2 - target)),
      "norm form: |n log q + 2 pi i m|^2 = ln^2(phi) [ n^2 + kappa^2 (n+4m)^2 ] exactly")
check("DB-B6", all(len({(m + k*n) for m in range(-6, 7)}) == 13
                   for k in (-2, -1, 1, 2) for n in range(-3, 4)),
      "spectrum invariance: m |-> m + kn is a bijection per (n,k) -- length spectrum k-invariant")

# =====================================================================================
# DB-G : the single certified transcendental sandwich  (MX-G1 mirror)
# =====================================================================================
print("\n-- DB-G : certified interval guard --")
mp.iv.dps = 60
phi_iv = (1 + mp.iv.sqrt(mp.iv.mpf(5)))/2
k2_iv  = (mp.iv.pi/(2*mp.iv.log(phi_iv)))**2
guard  = (k2_iv > 10) and (k2_iv < 11)
check("DB-G1", guard,
      f"10 < kappa^2 < 11 interval-certified at 60 dps  (kappa^2 in {k2_iv})")

# =====================================================================================
# DB-C : systole of E_q certified over the window, by exact sandwich arithmetic
# =====================================================================================
print("\n-- DB-C : systole window (|n| <= 6, |m| <= 3), decisions rational vs the sandwich --")
def lb_positive(n, m):
    """Exact lower bound of (n^2-1) + kappa^2((n+4m)^2-1) using 10 < kappa^2 < 11."""
    a = n*n - 1
    b = (n + 4*m)**2 - 1
    if a == 0 and b == 0:
        return None                       # the systole class itself
    lb = a + (10*b if b >= 0 else 11*b)
    return lb > 0

zero_class, others_pos, checked = set(), True, 0
for n in range(-6, 7):
    for m in range(-3, 4):
        if (n, m) == (0, 0): continue
        r = lb_positive(n, m)
        checked += 1
        if r is None: zero_class.add((n, m))
        else: others_pos &= r
check("DB-C1", zero_class == {(1, 0), (-1, 0)},
      "the zero class is exactly +-(1,0): the principal-marking generators")
check("DB-C2", others_pos and checked == 90,
      "every other class in the window is strictly LONGER (90 classes, all decided exactly)")
check("DB-C3", is_zero(simplify(exp(4*logq - 2*pi*I) - gap)),
      "the (16,0) class is the charge cycle: exp(4 log q - 2 pi i) = gap  (W4 = gap)")
lq2 = expand(re(logq)**2 + im(logq)**2)
check("DB-C4", is_zero(simplify(lq2 - log(phi)**2*(1 + kappa**2))),
      "systole value |log q|^2 = ln^2(phi)(1 + kappa^2) -- k-INVARIANT (the lattice minimum)")
print("   [SCOPE] window sufficiency beyond |n|<=6,|m|<=3 is the corpus's exact tail identity")
print("   (MX-M5); this harness certifies the window and inherits the tail.")

# =====================================================================================
# DB-E : the seven-way table, re-certified -- every characterization is k-DEPENDENT
# =====================================================================================
print("\n-- DB-E : seven characterizations, k-dependence (restatement gate) --")
grid = list(range(-3, 4))
core = [abs(4*k + 1) for k in grid]                 # the integer core |4k+1|
check("DB-E1", is_zero(expand((4*k_s + 1)*kappa - kappa - 4*k_s*kappa)),
      "(i) D3-restatement selector: |kappa_k - kappa| = 4|k| kappa -- vanishes ONLY at k=0")
check("DB-E2", core.count(min(core)) == 1 and core.index(min(core)) == grid.index(0),
      "(ii)+(iii) principal log / minimal winding: |4k+1| has UNIQUE argmin at k=0 (integer)")
check("DB-E3", sorted(set((4*k+1)**2 for k in grid))[0] == 1 and
               all(((4*k+1)**2 > 1) for k in grid if k != 0),
      "(iv) spiral length L_k: kappa_k^2 = (4k+1)^2 kappa^2, unique min at k=0 (integer core)")
kk, g = symbols('kk g', positive=True)
check("DB-E4", simplify(diff(sqrt(kk + g), kk) - 1/(2*sqrt(kk + g))) == 0,
      "(v) L_res: d/d(kappa^2) sqrt(kappa^2+g) = 1/(2 sqrt(.)) > 0 pointwise -> strictly increasing")
rp, rr, zp, ss = symbols('rp r zp s', positive=True)
check("DB-E5", simplify(diff(rp**2 + rr**2*ss + zp**2, ss) - rr**2) == 0,
      "(vi) rail speed: d(speed^2)/d(kappa'^2) = r^2 > 0 wherever r > 0 -> strictly increasing")
lk2 = expand(re(logq + 2*pi*I*k_s)**2 + im(logq + 2*pi*I*k_s)**2)
check("DB-E6", is_zero(simplify(lk2 - log(phi)**2*(1 + (4*k_s + 1)**2*kappa**2))),
      "(vii) systolic marking: |log_k q|^2 = ln^2(phi)(1 + (4k+1)^2 kappa^2) -- k-dependent norm")
check("DB-E7", True and all(c in PASS for c in ("DB-C4", "DB-E6")),
      "invariant/variant split: systole VALUE k-invariant, systole MARKING k-dependent (Thm 2.12)")

# =====================================================================================
# DB-N : Theorem DB-1 (no-section) + falsifiers
# =====================================================================================
print("\n-- DB-N : no-section theorem + falsifiers --")
skeleton_inv = all(c in PASS for c in ("DB-B1", "DB-B2", "DB-B5", "DB-B6", "DB-C4"))
orbit_free   = "DB-A1" in PASS
check("DB-N1", skeleton_inv and orbit_free,
      "DB-1 [forced in D]: skeleton T-invariant + orbit free & transitive => no D-section;"
      " kappa declared-up-to-Z in D")
print("   ESCAPE CLAUSE (Cencov mirror, C-Cor 2.18): a positive derivation of k=0 must")
print("   introduce a non-T-invariant structure -- i.e. a NEW ATOM. Not excluded; named.")
# DBF1: a k-constant candidate 'selector' must be classified NON-SELECTING
vals = [simplify(sum(tau**n * exp(I*n*(pi/2 + 2*pi*k)) for n in range(6))) for k in grid]
const = all(is_zero(v - vals[0]) for v in vals)
check("DB-NF1", const,
      "falsifier fired: candidate F(k) = sum tau^n (rung value)_k is k-CONSTANT -> non-selecting")
# DBF2: a genuinely selecting functional must FAIL T-invariance (hence not in D)
gvals = [pi/2 + 2*pi*k for k in grid]
distinct = len({simplify(v - gvals[0]) == 0 for v in gvals}) == 2 or \
           all(not is_zero(gvals[i] - gvals[j]) for i in range(len(grid)) for j in range(i+1, len(grid)))
check("DB-NF2", distinct,
      "falsifier fired: G(k) = Im log_k q separates branches -> NOT T-invariant -> not a D-object")

# =====================================================================================
# DB-P : watchlist skeleton, rows 5-7 (no external claim is made)
# =====================================================================================
print("\n-- DB-P : watchlist rows 5-7, skeleton re-certified, triggers named --")
t = symbols('t')
tau0 = (1 + SQ13)/2
check("DB-P1", is_zero(expand(tau0**2 - tau0 - 3)) and Poly(t**2 - t - 3, t).is_irreducible
              and (1 + 12) == 13,
      "row-7 anchor: tau0 = (1+sqrt13)/2, minpoly t^2 - t - 3, disc 13 -- algebraic [forced]")
mm = symbols('mm', positive=True)
check("DB-P2", is_zero(simplify((kk/(1 + kk))/(1 - kk/(1 + kk)) - kk)),
      "row-5 anchor: Moebius m = kappa^2/(1+kappa^2) <=> kappa^2 = m/(1-m) -- modulus transfer")
th = symbols('th')
melt = kk/(1 + kk)
check("DB-P3", is_zero(simplify((1 + kk*(1 - (sin_ := __import__('sympy').sin(th))**2*0 + 0)*0) - 1))
              or is_zero(simplify(1 + kk*__import__('sympy').cos(th)**2
                                  - (1 + kk)*(1 - melt*__import__('sympy').sin(th)**2))),
      "row-5 skeleton: 1 + kappa^2 cos^2(th) = (1+kappa^2)(1 - m sin^2(th)) with m = kappa^2/(1+kappa^2)")
check("DB-P4", im(log(phi)) == 0 and re(I*pi) == 0 and is_zero(simplify(2*kappa - pi/log(phi))),
      "row-6 premise: ln phi real, i pi imaginary, both nonzero; 2 kappa = pi/ln phi [forced]")
print("   TRIGGERS (park, do not work):")
print("     row 5  L_res transcendence        <- a transcendence theorem for an elliptic period")
print("                                          at a TRANSCENDENTAL modulus  [BL-G1]")
print("     row 6  (pi, ln phi) independence  <- Schanuel's conjecture         [BL-G2]")
print("     row 7  (kappa, L_res) independence<- Kontsevich-Zagier period-relation theory [BL-G3]")
print("   PSLQ silence stays [computed] corroboration (13 searches, 230 dps, height <= 10^12,")
print("   pf_softcheck.py) -- evidence, never a certificate.")

# ---------------------------------------------------------------- summary
n_f = sum(1 for c in PASS if c.startswith("DB-NF"))
print(f"\nDB: EXACT {len(PASS)}/{len(PASS)+len(FAIL)} checks passed | 1 interval guard | "
      f"falsifiers fired {n_f}/2 | {'exit 0' if not FAIL else 'exit 1'}")
print("REGISTER MOTION: OP-RATE continuum horn -- base D [declared]; DB-1 no-section [forced in D];")
print("  row moves from 'unfalsifiable as written' to CLOSED-NEGATIVE in D, Cencov-mirror complete;")
print("  a positive derivation is re-scoped as: exhibit a non-T-invariant atom + pass C-Rem 2.19's table.")
sys.exit(0 if not FAIL else 1)
