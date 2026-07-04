r"""
Producer: the reference-free parity floor (Section 5).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/parity_floor_qk.csv
              data/2026-07-relational-charge/parity_floor_summary.json

Refactors Lemma 5.1 (golden internal relation phi/phi' = -phi^2, so
t_rel(phi,phi') = 1/2), Lemma 5.2 (q_k = x^{2k}+x^k-1: real root of each sign,
angle set (1/2k)Z/Z, Delta = Z/2kZ, anchor n = 2k, M = phi), Lemma 5.3 (the sign
twist T(p)=(-1)^{deg p} p(-x): preserves M and Delta, exchanges anchor sectors
at odd m), Theorem 5.5 (mu_rel(even)=phi attained by q_{m/2}; mu_rel(odd)=mu(odd),
value 2 attained by x^m-2), and Corollary 5.7 (parity bit = evenness of |Delta|).

Statuses per the paper: even attainment at phi and odd attainment at 2 are
[forced]; the identification of the even infimum AS phi is conditional on the
inherited emission-gap hypothesis (EG); this is recorded in the summary JSON.

Run: py code/2026-07-relational-charge/parity_floor.py
"""

from fractions import Fraction

import mpmath as mp
import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv, write_json

x = C.x
PHI = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))


def qk_rows():
    """Lemma 5.2 / Theorem 5.5: the even-floor family q_k."""
    rows = []
    for k in range(1, 8):
        q = C.qk(k)
        ang = set(C.angles_of(q))
        expected = {Fraction(j, 2 * k) for j in range(2 * k)}
        m = C.relational_order(q)
        n = C.absolute_lcd(q)
        M = C.mahler_measure(q)
        diffs = {((a - b) % 1) for a in ang for b in ang}
        rows.append({
            "k": k,
            "q_k": f"x^{2*k}+x^{k}-1",
            "relational_order_m": m,
            "absolute_order_n": n,
            "anchor": "n=2k" if n == 2 * k else f"n={n}",
            "angle_set_is_full_1_over_2k": "yes" if ang == expected else "no",
            "delta_is_Z_2k": "yes" if m == 2 * k else "no",
            "half_in_difference_set": "yes" if Fraction(1, 2) in diffs else "no",
            "mahler": mp.nstr(M, 18),
            "mahler_is_phi": "yes" if abs(M - PHI) < mp.mpf(10) ** -20 else "no",
        })
    return rows


def twist_rows():
    """Lemma 5.3: the sign twist preserves M and Delta; exchanges anchors at
    odd m.  Emit before/after for a spread of objects."""
    rows = []
    objs = [
        ("x^3+2", x**3 + 2),
        ("golden x^2-x-1", x**2 - x - 1),
        ("q2=x^4+x^2-1", C.Q2),
        ("K=x^4+5x^2-5", C.KSEED),
        ("x^5-2", x**5 - 2),
        ("q3=x^6+x^3-1", x**6 + x**3 - 1),
    ]
    for name, p in objs:
        deg = sp.Poly(p, x).degree()
        T = sp.expand((-1) ** deg * p.subs(x, -x))
        rows.append({
            "object": name,
            "twist_T": str(sp.Poly(T, x).as_expr()),
            "m_before": C.relational_order(p),
            "m_after": C.relational_order(T),
            "n_before": C.absolute_lcd(p),
            "n_after": C.absolute_lcd(T),
            "delta_preserved": "yes" if C.relational_order(p) == C.relational_order(T) else "no",
            "mahler_preserved": "yes" if abs(C.mahler_measure(p) - C.mahler_measure(T)) < mp.mpf(10) ** -20 else "no",
            "monic_integral": "yes" if sp.Poly(T, x).is_monic else "no",
        })
    return rows


def build_summary():
    # golden internal relation
    phi = (1 + sp.sqrt(5)) / 2
    phip = (1 - sp.sqrt(5)) / 2
    ratio = sp.simplify(phi / phip)
    golden_ok = (sp.simplify(ratio - (-phi**2)) == 0) and (ratio < 0)
    trel = float((mp.arg(mp.mpf(str(sp.N(ratio, 40)))) / (2 * mp.pi)) % 1)

    # parity bit: 1/2 in (1/m)Z/Z iff m even
    parity = {m: (m % 2 == 0) for m in range(1, 21)}
    parity_bit_matches = all(
        (any(Fraction(j, m) == Fraction(1, 2) for j in range(m)) == (m % 2 == 0))
        for m in range(1, 21)
    )

    # odd attainment at 2 by x^m-2
    odd_attain = {}
    for m in (3, 5, 7):
        p = x**m - 2
        odd_attain[str(m)] = {
            "object": f"x^{m}-2",
            "relational_order": C.relational_order(p),
            "mahler": mp.nstr(C.mahler_measure(p), 12),
        }

    return {
        "golden_relation": {
            "identity": "phi/phi' = -phi^2 (negative real)",
            "t_rel_phi_phip": "1/2",
            "t_rel_numeric": round(trel, 12),
            "verified": bool(golden_ok and abs(trel - 0.5) < 1e-20),
            "status": "[forced] (ledger J)",
        },
        "parity_bit": {
            "statement": "1/2 in Delta(O) iff |Delta(O)| even (Cor 5.7)",
            "verified": bool(parity_bit_matches),
            "half_in_group_by_m": {str(m): v for m, v in parity.items()},
        },
        "theorem_5_5": {
            "even_m": "mu_rel(m) = phi, attained by q_{m/2}",
            "odd_m": "mu_rel(m) = mu(m); value 2 attained by x^m-2",
            "even_attainment_status": "[forced] (q_k gives M=phi)",
            "even_infimum_identification_status": "conditional on (EG), Rem 2.3",
            "odd_identity_status": "[forced] (twist reduction, Lem 5.3)",
            "odd_value_2_status": "[forced] at m=3; [computed]/[plausible] general",
            "odd_attainment_x_m_minus_2": odd_attain,
        },
    }


def main():
    qr = qk_rows()
    fq = ["k", "q_k", "relational_order_m", "absolute_order_n", "anchor",
          "angle_set_is_full_1_over_2k", "delta_is_Z_2k",
          "half_in_difference_set", "mahler", "mahler_is_phi"]
    pq = write_csv("parity_floor_qk.csv", fq, qr, __file__)
    print(f"wrote {pq}  ({len(qr)} q_k rows, Lem 5.2 / Thm 5.5)")

    tr = twist_rows()
    ft = ["object", "twist_T", "m_before", "m_after", "n_before", "n_after",
          "delta_preserved", "mahler_preserved", "monic_integral"]
    pt = write_csv("parity_floor_twist.csv", ft, tr, __file__)
    print(f"wrote {pt}  ({len(tr)} sign-twist rows, Lem 5.3)")

    ps = write_json("parity_floor_summary.json", build_summary(), __file__)
    print(f"wrote {ps}  (golden relation, parity bit, Thm 5.5 statuses)")


if __name__ == "__main__":
    main()
