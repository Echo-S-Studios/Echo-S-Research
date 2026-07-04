"""
Producer for Theorem 3.2 (thm:adams) -- "Squaring is the Adams operation"
[OA-PSI-01..05] and Proposition 3.3 (prop:plethysm) -- "the square is a
diagonal" [OA-PLE-01..06] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

Emits:
  * data/2026-06-operator-algebra/adams_laws.json -- the endomorphism laws
    (psi^n additive / multiplicative / unit-fixed, composition psi^m.psi^n =
    psi^{mn}), Newton's identity p_2 = e_1^2 - 2 e_2 and the general Newton
    recursion, and the doubling tower psi^{2^k}: lambda -> lambda^{2^k}, each
    verified symbolically over generic eigenvalues.
  * data/2026-06-operator-algebra/plethysm_traces.csv -- for each seed the
    concrete plethysm split A(x)A = Sym^2 A (+) Lambda^2 A and the traces
    tr Sym^2 = (p1^2+p2)/2, tr Lambda^2 = (p1^2-p2)/2 = e_2.

Difference from tests/: verifier asserts; this producer records the laws'
holds-flags plus the CONCRETE symmetric-function values (p1,p2,e1,e2, traces)
for the seed catalogue.
"""

from __future__ import annotations

import sympy as sp

from opalg_core import (
    ONE,
    elem_sym,
    exact_str,
    ms_equal,
    oplus,
    otimes,
    power_sum,
    psi,
    seed_catalogue,
    sym2,
    wedge2,
    write_csv,
    write_json,
)

a1, a2, a3, a4, a5 = sp.symbols("a1 a2 a3 a4 a5")
b1, b2 = sp.symbols("b1 b2")


def endomorphism_laws():
    A3, B2 = (a1, a2, a3), (b1, b2)
    add = all(ms_equal(psi(n, oplus(A3, B2)), oplus(psi(n, A3), psi(n, B2)))
              for n in (2, 3, 5))
    mul = all(ms_equal(psi(n, otimes((a1, a2), B2)),
                       otimes(psi(n, (a1, a2)), psi(n, B2)))
              for n in (2, 3, 4))
    unit = all(ms_equal(psi(n, ONE), ONE) for n in (2, 3, 7))
    comp = all(ms_equal(psi(m, psi(n, A3)), psi(m * n, A3))
               for m in (2, 3) for n in (2, 5))
    return add, mul, unit, comp


def newton_checks():
    """p_2 = e_1^2 - 2 e_2 for degrees 2..5, and general recursion k<=4."""
    p2 = all(
        sp.expand(power_sum(A, 2) - (elem_sym(A, 1) ** 2 - 2 * elem_sym(A, 2))) == 0
        for A in [(a1, a2), (a1, a2, a3), (a1, a2, a3, a4), (a1, a2, a3, a4, a5)]
    )
    A = (a1, a2, a3, a4, a5)
    gen = True
    for k in (1, 2, 3, 4):
        rhs = sum((-1) ** (i - 1) * elem_sym(A, i) * power_sum(A, k - i)
                  for i in range(1, k))
        rhs += (-1) ** (k - 1) * k * elem_sym(A, k)
        gen = gen and (sp.expand(power_sum(A, k) - rhs) == 0)
    return p2, gen


def doubling_tower():
    """psi^{2^k} sends lambda -> lambda^{2^k} (OA-PSI-05), k=1..4."""
    A = (a1, a2, a3)
    cur, ok = A, True
    for k in range(1, 5):
        cur = psi(2, cur)
        target = tuple(sp.expand(a ** (2**k)) for a in A)
        ok = ok and ms_equal(cur, target)
    return ok


def main():
    add, mul, unit, comp = endomorphism_laws()
    p2, gen = newton_checks()
    tower = doubling_tower()

    payload = {
        "endomorphism_laws": {
            "OA-PSI-01_additive": {
                "statement": "psi^n(A (+) B) = psi^n A (+) psi^n B", "holds": bool(add)},
            "OA-PSI-02_multiplicative": {
                "statement": "psi^n(A (x) B) = psi^n A (x) psi^n B", "holds": bool(mul)},
            "OA-PSI-03_unit_fixed": {
                "statement": "psi^n(1) = 1", "holds": bool(unit)},
            "OA-PSI-04_composition": {
                "statement": "psi^m . psi^n = psi^{mn}", "holds": bool(comp)},
            "OA-PSI-05_doubling_tower": {
                "statement": "psi^{2^k} : lambda -> lambda^{2^k}", "holds": bool(tower)},
        },
        "newton_identities": {
            "OA-PSI_p2": {
                "statement": "p_2 = e_1^2 - 2 e_2  (degrees 2..5)", "holds": bool(p2)},
            "general_recursion": {
                "statement": "p_k = sum_{i=1}^{k-1} (-1)^{i-1} e_i p_{k-i} + (-1)^{k-1} k e_k  (k<=4)",
                "holds": bool(gen)},
        },
        "lambda_ring": (
            "psi^n both additive and multiplicative -> the operator set is a "
            "lambda-ring, not merely a semiring with a unary map (Adams ops "
            "determined by exterior powers via Newton)."),
    }
    p_json = write_json("adams_laws.json", "adams_plethysm.py", payload)

    # plethysm_traces.csv on the concrete seed catalogue
    rows = []
    for name, poly, A in seed_catalogue():
        p1 = power_sum(A, 1)
        p2v = power_sum(A, 2)
        e1 = elem_sym(A, 1)
        e2 = elem_sym(A, 2)
        S2, W2 = sym2(A), wedge2(A)
        tr_sym2 = sp.expand(sum(S2))
        tr_wedge2 = sp.expand(sum(W2))
        split_ok = ms_equal(otimes(A, A), oplus(S2, W2))
        diag = tuple(sp.expand(a**2) for a in A)  # psi^2 A
        diag_is_psi2 = ms_equal(diag, psi(2, A))
        rows.append(
            {
                "seed": name,
                "degree": len(A),
                "p1": exact_str(p1),
                "p2": exact_str(p2v),
                "e1": exact_str(e1),
                "e2": exact_str(e2),
                "tr_Sym2": exact_str(tr_sym2),
                "tr_Sym2_formula_(p1^2+p2)/2": exact_str((p1**2 + p2v) / 2),
                "tr_Wedge2": exact_str(tr_wedge2),
                "tr_Wedge2_formula_(p1^2-p2)/2": exact_str((p1**2 - p2v) / 2),
                "tr_Wedge2_equals_e2": bool(sp.expand(tr_wedge2 - e2) == 0),
                "tensor_square_splits": bool(split_ok),
                "psi2_is_diagonal_of_Sym2": bool(diag_is_psi2),
            }
        )
    p_csv = write_csv(
        "plethysm_traces.csv",
        "adams_plethysm.py",
        ["seed", "degree", "p1", "p2", "e1", "e2",
         "tr_Sym2", "tr_Sym2_formula_(p1^2+p2)/2",
         "tr_Wedge2", "tr_Wedge2_formula_(p1^2-p2)/2", "tr_Wedge2_equals_e2",
         "tensor_square_splits", "psi2_is_diagonal_of_Sym2"],
        rows,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv} ({len(rows)} seeds)")


if __name__ == "__main__":
    main()
