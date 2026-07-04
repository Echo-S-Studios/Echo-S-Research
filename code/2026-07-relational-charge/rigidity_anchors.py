r"""
Producer: rigidity, anchors, and coherence groups (Section 3-4).

Source paper: papers/2026-07-relational-charge/relational_charge_paper.tex
Produces    : data/2026-07-relational-charge/rigidity_anchors.csv

Refactors Theorem 4.1 (thm:rigidity), Example 4.6 (ex:drop, the group drop),
Lemma 4.4 (lem:oddfull), Example 3.5 (ex:firsttype), and Proposition 3.6
(cyclicity).  For each charge-admissible object it recovers the angle set from
high-precision roots, then emits: absolute charge order n (lcd of the angles),
relational group order m (lcm of pairwise-difference denominators), the anchor
class n in {m, 2m}, whether the relational difference set is exactly the cyclic
group (1/m)Z/Z, and the Mahler measure.  The two anchors (x^m-2 -> n=m;
x^3+2 -> n=2m) and the group drop (x^4+5x^2+5: absolute Z/4, relational Z/2)
are exhibited directly.

Run: py code/2026-07-relational-charge/rigidity_anchors.py
"""

from fractions import Fraction

import sympy as sp

import relcharge_core as C
from relcharge_io import write_csv

x = C.x

CASES = [
    ("x^2-2", x**2 - 2),
    ("x^3-2", x**3 - 2),
    ("x^4-2", x**4 - 2),
    ("x^5-2", x**5 - 2),
    ("x^6-2", x**6 - 2),
    ("x^7-2", x**7 - 2),
    ("x^3+2", x**3 + 2),
    ("x^4+2", x**4 + 2),
    ("x^4+5x^2+5 (group drop)", C.GROUPDROP),
    ("K=x^4+5x^2-5", C.KSEED),
    ("q2=x^4+x^2-1", C.Q2),
    ("golden x^2-x-1", x**2 - x - 1),
    ("q3=x^6+x^3-1", x**6 + x**3 - 1),
]


def build_rows():
    rows = []
    for name, p in CASES:
        ang = C.angles_of(p)
        n = C.absolute_lcd(p)
        m = C.relational_order(p)
        diffs = {((a - b) % 1) for a in ang for b in ang}
        cyclic = diffs == {Fraction(j, m) for j in range(m)}
        if n == m:
            anchor = "n=m"
        elif n == 2 * m:
            anchor = "n=2m"
        else:
            anchor = f"OTHER(n={n},m={m})"
        rows.append({
            "object": name,
            "degree": sp.Poly(p, x).degree(),
            "angles": " ".join(f"{a.numerator}/{a.denominator}" for a in sorted(set(ang))),
            "absolute_order_n": n,
            "relational_order_m": m,
            "anchor": anchor,
            "n_in_m_or_2m": "yes" if n in (m, 2 * m) else "no",
            "delta_is_cyclic_Zm": "yes" if cyclic else "no",
            "odd_n_full": ("n odd => m=n" if (n % 2 == 1 and m == n)
                           else ("n odd MISMATCH" if n % 2 == 1 else "n even")),
            "mahler": mahler_str(p),
        })
    return rows


def mahler_str(p):
    """Mahler measure with the small closed forms that appear here recognised."""
    import mpmath as mp
    v = C.mahler_measure(p)
    phi = mp.mpf(str(sp.N((1 + sp.sqrt(5)) / 2, 45)))
    for label, val in [("1", mp.mpf(1)), ("2", mp.mpf(2)), ("5", mp.mpf(5)),
                       ("phi", phi)]:
        if abs(v - val) < mp.mpf(10) ** -12:
            return label
    return mp.nstr(v, 12)


def main():
    rows = build_rows()
    fields = ["object", "degree", "angles", "absolute_order_n",
              "relational_order_m", "anchor", "n_in_m_or_2m",
              "delta_is_cyclic_Zm", "odd_n_full", "mahler"]
    path = write_csv("rigidity_anchors.csv", fields, rows, __file__)
    print(f"wrote {path}")
    print(f"  {len(rows)} admissible objects (Thm 4.1, Ex 4.6, Lem 4.4, Prop 3.6)")


if __name__ == "__main__":
    main()
