"""
Producer: the group-theoretic backbone (Section 3), exact Z/L arithmetic.

Source paper: papers/2026-06-charge-measure-coupling/charge-measure-coupling-whitepaper-v4.tex
Produces    :
  * data/2026-06-charge-measure-coupling/group_structure.json
        Thm 3.1 (cyclicity), Thm 3.4 (CRT + Adams primary projectors)
  * data/2026-06-charge-measure-coupling/lcm_composition.csv
        Thm 3.5 (lcm law) + Thm 3.7 (safe composition <=> coprime), one row per
        ordered pair (m,k) of charge-group orders

Everything here is exact integer computation in Z/L -- no floating point -- so
these are proven-structure results (all tagged [forced] in the paper).

Run: py code/2026-06-charge-measure-coupling/charge_group_structure.py
"""

from math import gcd

import sympy as sp

import cmc_core as core
from cmc_io import write_json, write_csv


# --------------------------------------------------------------------------
# Thm 3.1: every finite subgroup of Q/Z is cyclic
# --------------------------------------------------------------------------
def cyclicity_evidence(Lmax: int = 24):
    """For a spread of charge-lattice orders L and generator sets, confirm the
    generated subgroup of Z/L is exactly the cyclic group <gcd(gens,L)>."""
    gen_sets = [[1], [2, 3], [4, 6], [3, 6, 9]]
    checked, all_cyclic = 0, True
    for L in range(2, Lmax + 1):
        for gens in gen_sets:
            g = core.cyclic_generator_mod(gens, L)
            S = core.subgroup_closure(gens, L)
            expected = {(g * j) % L for j in range(L // g)}
            checked += 1
            if S != expected or len(S) != L // g:
                all_cyclic = False
    # explicit counterexample: Z/p x Z/p is NOT cyclic (exponent p < p^2)
    non_cyclic = {}
    for p in (2, 3, 5):
        max_order = max(
            _order_in_ZpxZp(a, b, p) for a in range(p) for b in range(p))
        non_cyclic[f"Z/{p} x Z/{p}"] = {
            "group_order": p * p,
            "max_element_order (exponent)": max_order,
            "cyclic": max_order == p * p,
        }
    return {
        "checked_subgroups": checked,
        "all_generated_subgroups_cyclic": all_cyclic,
        "L_range": [2, Lmax],
        "noncyclic_counterexamples": non_cyclic,
        "paper_ref": "Thm 3.1",
    }


def _order_in_ZpxZp(a, b, p):
    o = 1
    while (a * o) % p != 0 or (b * o) % p != 0:
        o += 1
    return o


# --------------------------------------------------------------------------
# Thm 3.4: CRT isomorphism + Adams primary projectors psi^{n/p^e}
# --------------------------------------------------------------------------
def crt_adams_evidence(ns=(6, 12, 30, 36, 60)):
    out = {}
    for n in ns:
        fac = sp.factorint(n)
        # CRT bijection Z/n -> prod Z/p^e
        images = {tuple(a % (p**e) for p, e in fac.items()) for a in range(n)}
        projectors = {}
        for p, e in fac.items():
            pe = p**e
            k = n // pe
            img = {(k * a) % n for a in range(n)}
            ker = {a for a in range(n) if (k * a) % n == 0}
            projectors[f"psi^{k}"] = {
                "targets_component": f"Z/{pe}",
                "image_order": len(img),
                "kernel_order": len(ker),
                "image_is_unique_order_pe_subgroup":
                    img == {((n // pe) * j) % n for j in range(pe)},
            }
        out[str(n)] = {
            "factorization": {str(p): int(e) for p, e in fac.items()},
            "crt_bijection": len(images) == n,
            "num_independent_primary_labels (omega(n))": len(fac),
            "adams_projectors": projectors,
        }
    return {"cases": out, "paper_ref": "Thm 3.4"}


# --------------------------------------------------------------------------
# Thm 3.5 + Thm 3.7: lcm law and safe composition <=> coprime
# --------------------------------------------------------------------------
def lcm_composition_rows(max_order: int = 12):
    rows = []
    for m in range(1, max_order + 1):
        for k in range(1, max_order + 1):
            L = core.lcm(m, k)
            # generate <L/m, L/k> in Z/L; order must be L (Thm 3.5)
            S = core.subgroup_closure([L // m, L // k], L)
            kernel_order = (m * k) // L  # = gcd(m,k) (Thm 3.7)
            rows.append({
                "m": m,
                "k": k,
                "gcd": gcd(m, k),
                "lcm": L,
                "tensor_charge_group": f"Z/{L}",
                "generated_order": len(S),
                "full_lattice": "yes" if len(S) == L else "no",
                "kernel_order": kernel_order,
                "lossless_safe": "yes" if kernel_order == 1 else "no",
                "coprime": "yes" if gcd(m, k) == 1 else "no",
            })
    return rows


def main():
    payload = {
        "cyclicity": cyclicity_evidence(),
        "crt_adams": crt_adams_evidence(),
        "lcm_law_note": "tensor charge group = Z/lcm(m,k); kernel order = gcd(m,k); "
                        "safe (lossless) iff coprime -- see lcm_composition.csv",
    }
    p1 = write_json("group_structure.json", payload, __file__)
    rows = lcm_composition_rows()
    fields = ["m", "k", "gcd", "lcm", "tensor_charge_group", "generated_order",
              "full_lattice", "kernel_order", "lossless_safe", "coprime"]
    p2 = write_csv("lcm_composition.csv", fields, rows, __file__)
    print(f"wrote {p1}")
    print(f"wrote {p2}  ({len(rows)} ordered (m,k) pairs)")


if __name__ == "__main__":
    main()
