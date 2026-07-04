r"""
Producer: the consolidation ledger (Appendix A, entries A-X).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/consolidation_ledger.csv

Regenerates the paper's master Appendix-A table: one row per ledger entry A-X,
each recording the object, the exact result (contact signature / relational data
/ verdict), the epistemic status tag, and the paper reference.  The cheap
signature and structural facts are recomputed inline here; the three heavy or
external computations are cross-referenced to their dedicated producers:
  T (full 37-instance census scan) -> census_deg12.py
  V (larger-coefficient search)    -> spot_checks_larger_coeff.py
  X nested deg-256 nu-scan         -> showcase_x4_minus_x_plus_1.py
  M (PARI/GP cross-engine)         -> external second stack, not run here

Run: py code/2026-07-relational-charge/consolidation_ledger.py
"""

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))


def sig(p):
    return C.signature_str(C.cyclotomic_contacts(C.ratio_poly(p)))


def msig(p, q):
    return C.signature_str(C.cyclotomic_contacts(C.mixed_ratio_poly(p, q)))


def mahler_label(p):
    v = C.mahler_measure(p)
    for lab, val in [("1", mp.mpf(1)), ("2", mp.mpf(2)), ("5", mp.mpf(5)),
                     ("phi", PHI), ("phi^2", PHI**2)]:
        if abs(v - val) < mp.mpf(10) ** -12:
            return lab
    return mp.nstr(v, 10)


def build_rows():
    rows = []

    def add(entry, obj, result, status, ref):
        rows.append({"entry": entry, "object": obj, "exact_result": result,
                     "status": status, "paper_ref": ref})

    # A, B: gauge-blind pair
    add("A", "x^3-2", f"contacts {sig(x**3-2)}; Delta=Z/3, anchor n={C.absolute_lcd(x**3-2)}, M={mahler_label(x**3-2)}",
        "[forced]", "Thm 4.1 / rem:scope")
    add("B", "x^3+2", f"contacts {sig(x**3+2)} (= A, gauge-blind); Delta=Z/3, anchor n={C.absolute_lcd(x**3+2)}; T(x^3+2)=x^3-2",
        "[forced]", "rem:scope")
    add("C", "x^4-2", f"contacts {sig(x**4-2)}: relational recovery of Z/4",
        "[forced]", "ex:twistshell")
    add("D", "q2=x^4+x^2-1", f"contacts {sig(C.Q2)}; full Delta=Z/4 (symbolic); M={mahler_label(C.Q2)}",
        "[forced]", "Sec 7 method")
    add("E", "x^4+5x^2+5", f"Eisenstein@5; anchor n={C.absolute_lcd(C.GROUPDROP)}, Delta=Z/{C.relational_order(C.GROUPDROP)} (group drop), M={mahler_label(C.GROUPDROP)}; contacts {sig(C.GROUPDROP)}",
        "[forced]", "Ex 4.6")
    add("F", "K=x^4+5x^2-5", f"contacts {sig(C.KSEED)} (= E); full type differs: Delta=Z/{C.relational_order(C.KSEED)} (shell sig strictly coarser)",
        "[forced]", "Ex 4.6 / rem:scope")

    # G: two-route
    two = "agree" if C.cyclotomic_contacts(C.ratio_poly(C.B4)) == C.cyclotomic_contacts(C.ratio_poly_via_kronecker(C.B4)) else "DISAGREE"
    add("G", "beta4", f"contacts {sig(C.B4)} only (complete M<=512); two-route resultant=Kronecker {two}",
        "[forced] per inst.", "Thm 7.9")
    Rl = C.ratio_poly(C.LEHMER)
    add("H", "Lehmer (deg 10)", f"deg Rat_L={Rl.degree()}; contacts {C.signature_str(C.cyclotomic_contacts(Rl))} (complete M<=20000)",
        "[forced] per inst.", "Thm 7.9")
    add("I", "beta4 x L mixed", f"gcd=1; deg Rat=40; contacts {msig(C.LEHMER, C.B4)} (not circle-locked)",
        "[forced] per inst.", "Thm 8.2")

    # J, K: golden machinery
    phi = (1 + sp.sqrt(5)) / 2
    jrel = "1/2" if sp.simplify(sp.simplify(phi / ((1 - sp.sqrt(5)) / 2)) - (-phi**2)) == 0 else "?"
    add("J", "golden pair", f"phi/phi' = -phi^2 (negative real): t_rel(phi,phi')={jrel}",
        "[forced]", "Lem 5.1")
    yroots = sp.solve(sp.symbols("y")**2 + sp.symbols("y") - 1, sp.symbols("y"))
    add("K", "y^2+y-1", f"roots 1/phi>0 and -phi<0 (real +- pair behind q_k): max={sp.simplify(max(yroots))}",
        "[forced]", "Lem 5.2")

    # L: coset arithmetic for m in {3,5,7,9}
    from fractions import Fraction
    lok = all((e + m) % 2 == 0 and m % Fraction(e + m, 2 * m).denominator == 0
              for m in (3, 5, 7, 9) for e in range(1, 2 * m, 2))
    add("L", "coset arithmetic", f"e,m odd => (e+m)/2m in (1/m)Z/Z for m in {{3,5,7,9}}: {'holds' if lok else 'FAILS'}",
        "[forced]", "Lem 5.3")

    add("M", "cross-engine (PARI/GP)", "entries A-I,P,Q,U,V + round-8 X recomputed in PARI/GP 2.15.4 (polresultant+poliscyclo); external second stack, not run here",
        "[computed] (external)", "Sec 8 method")

    # N: totient bound
    N = 200000
    phi_sieve = C.totients_upto(N)
    tight = [M for M in range(1, N + 1) if 2 * phi_sieve[M] * phi_sieve[M] == M]
    viol = [M for M in range(1, N + 1) if 2 * phi_sieve[M] * phi_sieve[M] < M]
    add("N", "totient bound", f"2 phi(M)^2 >= M for all M<=2e5: {len(viol)} violations, unique tight case {tight}",
        "[forced]+[computed]", "Lem 4.5")

    # O: Salem certification patterns
    opat = {name: C.trace_sturm_pattern(p)[:3] for name, p in
            [("beta4", C.B4), ("S6", C.S6), ("S8", C.S8), ("L", C.LEHMER)]}
    all_ok = all(C.is_salem_polynomial(p) for _, p in
                 [("beta4", C.B4), ("S6", C.S6), ("S8", C.S8), ("L", C.LEHMER)])
    add("O", "Salem certification", f"beta4,S6,S8,L: reciprocal+irreducible, trace-Sturm (1,0,d-1) all pass = {all_ok}",
        "[forced] per inst.", "Thm 7.10")

    add("P", "S6, S8 inertness", f"Rat_S6 (deg 36) {sig(C.S6)}; Rat_S8 (deg 64) {sig(C.S8)}",
        "[forced] per inst.", "Thm 7.10")

    # Q: the five remaining pairs
    pairs = [("beta4xS6", C.B4, C.S6), ("beta4xS8", C.B4, C.S8),
             ("S6xS8", C.S6, C.S8), ("LxS6", C.LEHMER, C.S6), ("LxS8", C.LEHMER, C.S8)]
    q_all_empty = all(C.cyclotomic_contacts(C.mixed_ratio_poly(p, q)) == {} for _, p, q in pairs)
    add("Q", "pairwise locking batch", f"5 mixed objects (deg 24,32,48,60,80): all empty contacts = {q_all_empty}",
        "[forced] per inst.", "Rem 8.4")

    add("R", "quadratic identities", f"(1+sqrt5)/2=phi: {sp.simplify(phi-(1+sp.sqrt(5))/2)==0}; (3+sqrt5)/2=phi^2: {sp.simplify(phi**2-(3+sp.sqrt(5))/2)==0}",
        "[forced]", "App D")

    # S: beta4 tensor beta4
    Cm = C.companion_matrix(C.B4)
    F = sp.Poly(sp.Matrix(sp.kronecker_product(Cm, Cm)).charpoly(x).as_expr(), x)
    add("S", "beta4 (x) beta4", f"(x-1)-mult={C.phi1_multiplicity(F)}, deg gcd(F,F')={sp.gcd(F, F.diff(x)).degree()}, 3 distinct positive reals (non-inert)",
        "[forced] structure", "Ex 8.5")

    add("T", "degree-12 census", "729 recip {-1,0,1}-polys; 378 twist-classes (27 fixed); 37 Salem; all deg-144 scans {Phi_1^12}. Full run: census_deg12.py",
        "[forced] per inst.; [computed] family", "Thm 6.13")

    add("U", "twisted-shell witness", f"q=x^2+x+2 {sig(x**2+x+2)} (angle irrational); p=x^4+x^2+2 {sig(C.TWISTSHELL)}, M={mahler_label(C.TWISTSHELL)} (non-inert)",
        "[forced]", "Ex 7.21")

    add("V", "larger-coeff spot checks", "first five deg-12 Salem twist-classes with a coeff outside {-1,0,1}: all {Phi_1^12}. Full run: spot_checks_larger_coeff.py",
        "[forced] per inst.", "Rem 6.14")

    add("W", "plastic; census tally", f"x^3-x-1 {sig(C.PLASTIC)} (inert Pisot); census rejection tally 39/256/46/0/37 sum 378. Full: census_deg12.py",
        "[forced] per inst.", "Ex 7.19 / Rem 6.14")

    # X: x^4-x+1 (fast parts inline; nested deg-256 cross-referenced)
    Cx = C.companion_matrix(C.X4MX1)
    Fx = sp.factor(sp.Matrix(sp.kronecker_product(Cx, Cx)).charpoly(x).as_expr())
    kron_ok = sp.expand(Fx - (x**6 - x**4 - x**3 - x**2 + 1) ** 2 * (x**4 + 2 * x**2 - x + 1)) == 0
    gx, _ = sp.Poly(C.X4MX1, x).galois_group()
    add("X", "x^4-x+1 (nu executed)", f"charpoly(C(x)C)=S6^2*(x^4+2x^2-x+1) {kron_ok}; Rat_p {sig(C.X4MX1)}; nested Rat_Rat deg 256 {{Phi_1^28}} (showcase); Gal order {gx.order()}(S4); M=tau_S6",
        "[forced] per inst.", "Ex 6.20")

    return rows


def main():
    rows = build_rows()
    fields = ["entry", "object", "exact_result", "status", "paper_ref"]
    path = write_csv("consolidation_ledger.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} ledger rows (Appendix A entries A-X)")


if __name__ == "__main__":
    main()
