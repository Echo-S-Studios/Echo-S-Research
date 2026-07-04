"""
Producer for Character II -- the angle charge, Theorem 5.2 (thm:charge)
[OA-CH-01..05] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

The charge q(lambda) = round(2 arg lambda / pi) mod 4 is a finite character
onto Z/4Z:
    charge(A (x) B) = charge(A) + charge(B)   (sumset mod 4)
    charge(psi^2 A) = 2 * charge(A)
    charge(A (+) B) = charge(A) union charge(B)
Real seeds carry {0,2}; the Lorentzian seed K = x^4+5x^2-5 realises the full
group {0,1,2,3}, its imaginary pair +- i beta at +- pi/2.

Emits:
  * data/2026-06-operator-algebra/charge_sumset.csv -- for each seed pair the
    tensor charge multiset vs the sumset-mod-4, and the union-on-(+) check.
  * data/2026-06-operator-algebra/charge_laws.json -- axis charges, per-seed
    charge multisets, the psi^2 doubling law, group realisation, and the
    double-psi^2 collapse to {0}.

Difference from tests/: verifier asserts each law; this producer emits the
concrete charge multisets and the two sides of each transformation law.
"""

from __future__ import annotations

import sympy as sp

from opalg_core import (
    I,
    K_seed,
    charge,
    charge_multiset_str,
    charge_one,
    golden_seed,
    oplus,
    otimes,
    psi,
    seed_catalogue,
    seed_from_poly,
    sumset_mod4,
    write_csv,
    write_json,
    x,
)


def _counter_to_dict(c):
    return {int(k): int(v) for k, v in sorted(c.items())}


def main():
    cat = seed_catalogue()

    # ---- sumset on (x) and union on (+) across pairs -------------------
    rows = []
    for na, _, A in cat:
        cA = charge(A)
        for nb, _, B in cat:
            cB = charge(B)
            c_tensor = charge(otimes(A, B))
            c_sumset = sumset_mod4(cA, cB)
            c_union = charge(oplus(A, B))
            rows.append(
                {
                    "A": na,
                    "B": nb,
                    "charge_A": charge_multiset_str(cA),
                    "charge_B": charge_multiset_str(cB),
                    "charge_A_tensor_B": charge_multiset_str(c_tensor),
                    "sumset_mod4": charge_multiset_str(c_sumset),
                    "sumset_ok": bool(c_tensor == c_sumset),
                    "charge_A_oplus_B": charge_multiset_str(c_union),
                    "union_ok": bool(c_union == (cA + cB)),
                }
            )
    p_csv = write_csv(
        "charge_sumset.csv",
        "charge_character.py",
        ["A", "B", "charge_A", "charge_B", "charge_A_tensor_B", "sumset_mod4",
         "sumset_ok", "charge_A_oplus_B", "union_ok"],
        rows,
    )

    # ---- scalar laws ---------------------------------------------------
    # doubling under psi^2: charge(psi^2 A) = 2*charge(A) (mod 4)
    doubling = []
    for name, _, A in cat:
        cA = charge(A)
        expected = {}
        for q, n in cA.items():
            expected[(2 * q) % 4] = expected.get((2 * q) % 4, 0) + n
        got = _counter_to_dict(charge(psi(2, A)))
        doubling.append(
            {
                "seed": name,
                "charge_A": _counter_to_dict(cA),
                "two_times_charge_A": {int(k): int(v) for k, v in sorted(expected.items())},
                "charge_psi2_A": got,
                "doubling_ok": got == {int(k): int(v) for k, v in sorted(expected.items())},
            }
        )

    # K's imaginary pair at +- pi/2 (charges 1,3)
    Kr = K_seed()
    imag = [r for r in Kr if sp.simplify(sp.re(r)) == 0]
    imag_charges = sorted(charge_one(r) for r in imag)

    # double psi^2 collapses every charge to 0 (4q = 0 mod 4)
    double_psi2 = _counter_to_dict(charge(psi(2, psi(2, K_seed()))))

    payload = {
        "definition": "q(lambda) = round(2*arg(lambda)/pi) mod 4  in Z/4Z",
        "axis_charges": {
            "q(+1)": charge_one(sp.Integer(1)),
            "q(-1)": charge_one(sp.Integer(-1)),
            "q(+i)": charge_one(I),
            "q(-i)": charge_one(-I),
        },
        "per_seed_charge": {
            name: _counter_to_dict(charge(A)) for name, _, A in cat},
        "doubling_under_psi2": doubling,
        "group_realisation": {
            "real_seeds_carry": {
                "golden": sorted(charge(golden_seed()).keys()),
                "x^2-2": sorted(charge(seed_from_poly(x**2 - 2)).keys()),
                "x^2-3": sorted(charge(seed_from_poly(x**2 - 3)).keys()),
            },
            "real_seeds_are_{0,2}": all(
                set(charge(s).keys()) == {0, 2}
                for s in [golden_seed(), seed_from_poly(x**2 - 2),
                          seed_from_poly(x**2 - 3)]),
            "K_realises_full_Z4": sorted(charge(K_seed()).keys()) == [0, 1, 2, 3],
            "K_imaginary_pair_charges": imag_charges,
            "K_imaginary_pair_at_half_pi": imag_charges == [1, 3],
        },
        "double_psi2_collapses_to_0": {
            "charge_psi2_psi2_K": double_psi2,
            "collapses_to_{0}": set(double_psi2.keys()) == {0},
            "note": "psi^2 doubles charge; applied twice, 4q = 0 mod 4",
        },
    }
    p_json = write_json("charge_laws.json", "charge_character.py", payload)
    print(f"wrote {p_csv} ({len(rows)} pairs)")
    print(f"wrote {p_json}")
    all_ok = all(r["sumset_ok"] for r in rows) and all(d["doubling_ok"] for d in doubling)
    print(f"  sumset/doubling all_ok={all_ok}; K charge={sorted(charge(K_seed()).keys())}")


if __name__ == "__main__":
    main()
