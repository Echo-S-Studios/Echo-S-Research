"""
Producer for Section 7 -- Fixed points and the spectrum-preserving operators,
Proposition 7.1 (prop:fixed) [OA-FX-01, OA-FX-02] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

M(psi^2 A) = M(A) holds iff M(A) = 1 (Kronecker / cyclotomic locus: every
eigenvalue a root of unity or 0).  M(A)^2 = M(A) forces M in {0,1}; as M >= 1
the only Adams-fixed measure is the trivial one.  minpoly / Phi re-present an
object without changing its eigenvalue multiset -- idempotent, spectrum- and
character-preserving.

Emits data/2026-06-operator-algebra/fixed_points.json.

Difference from tests/: verifier asserts; this producer records the M^2=M
solution set, a table of cyclotomic objects (M=1, Adams-fixed) versus the
golden object (M=phi>1, squares instead), and the minpoly reconstruction check.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

from opalg_core import (
    charge,
    charge_multiset_str,
    exact_str,
    golden_seed,
    mahler_exact,
    modulus_mp,
    ms_equal,
    psi,
    seed_from_poly,
    write_json,
    x,
)


def main():
    # fixed-point algebra: M^2 = M  <=>  M in {0,1}; feasible (M>=1) = {1}
    M = sp.Symbol("M")
    sols = sorted(int(s) for s in sp.solve(sp.Eq(M**2, M), M))
    feasible = [s for s in sols if s >= 1]

    # cyclotomic objects: M = 1 and Adams-fixed
    cyclo_polys = [
        ("Phi_4 = x^2+1", x**2 + 1),
        ("Phi_3 = x^2+x+1", x**2 + x + 1),
        ("Phi_8 = x^4+1", x**4 + 1),
        ("Phi_9 = x^6+x^3+1", x**6 + x**3 + 1),
        ("Phi_1 = x-1", x - 1),
        ("Phi_2 = x+1", x + 1),
    ]
    cyclo = []
    mp.mp.dps = 40
    for label, poly in cyclo_polys:
        A = seed_from_poly(poly)
        m = mahler_exact(A)
        m_psi2 = mahler_exact(psi(2, A))
        on_circle = all(abs(modulus_mp(a, 40) - 1) < mp.mpf(10) ** (-30) for a in A)
        cyclo.append(
            {
                "object": label,
                "degree": len(A),
                "M": exact_str(m),
                "M_is_1": bool(sp.simplify(m - 1) == 0),
                "all_eigenvalues_on_unit_circle": bool(on_circle),
                "M_psi2": exact_str(m_psi2),
                "adams_fixed": bool(sp.simplify(m_psi2 - m) == 0),
            }
        )

    # golden: M = phi > 1, NOT Adams-fixed; squares instead
    G = golden_seed()
    mG = mahler_exact(G)
    mG2 = mahler_exact(psi(2, G))
    golden = {
        "object": "golden = x^2-x-1",
        "M": exact_str(mG),
        "M_gt_1": bool(sp.simplify(mG - 1) != 0 and sp.N(mG) > 1),
        "M_psi2": exact_str(mG2),
        "adams_fixed": bool(sp.simplify(mG2 - mG) == 0),
        "squares_instead_M_psi2_eq_M_squared": bool(sp.simplify(mG2 - mG**2) == 0),
    }

    # minpoly / Phi idempotent: rebuild golden from x^2-x-1, spectrum + both
    # characters preserved
    recon = seed_from_poly(x**2 - x - 1)
    idempotent = {
        "object": "golden rebuilt from minpoly x^2-x-1",
        "spectrum_preserved": bool(ms_equal(G, recon)),
        "measure_preserved": bool(sp.simplify(mahler_exact(G) - mahler_exact(recon)) == 0),
        "charge_preserved": bool(charge(G) == charge(recon)),
        "charge_multiset": charge_multiset_str(charge(recon)),
    }

    payload = {
        "fixed_point_equation": {
            "statement": "M(psi^2 A) = M(A)  <=>  M(A)^2 = M(A)",
            "solutions_of_M2_eq_M": sols,
            "feasible_with_M_ge_1": feasible,
            "conclusion": "only Adams-fixed measure is M(A)=1 (cyclotomic / Kronecker locus)",
        },
        "cyclotomic_objects_are_adams_fixed": cyclo,
        "golden_not_adams_fixed": golden,
        "spectrum_preserving_idempotents": idempotent,
    }
    path = write_json("fixed_points.json", "fixed_points.py", payload)
    print(f"wrote {path}")
    print(f"  M^2=M solutions={sols}, feasible={feasible}; "
          f"golden adams_fixed={golden['adams_fixed']} "
          f"(squares={golden['squares_instead_M_psi2_eq_M_squared']})")


if __name__ == "__main__":
    main()
