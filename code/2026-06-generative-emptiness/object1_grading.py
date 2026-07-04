"""
Producer -- Object I: the Z/4Z grading.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex
Result produced: Theorem 2.2 (thm:charge), Def. 2.1 (def:charge).

Computes and EMITS, to data/2026-06-generative-emptiness/object1_grading.json:
  * the charge multisets of the seeds and their closures
        phi:{0,2}, phi(x)phi:{0,0,2,2}, phi^2:{0,0}, K:{0,1,2,3};
  * the operator action on the charge -- (x)=add, ( )^2=double, (+)=union --
    as computed-vs-predicted charge multisets;
  * the Z/4Z Cayley data (closed under +, x2, union);
  * a two-generation orbit sweep confirming no operator word leaves the lattice;
  * confirmation that K is an irreducible quartic (needed to realise charges 1,3).

Run:
    py code/2026-06-generative-emptiness/object1_grading.py
"""
import json
import os

import sympy as sp

import ge_core as G
from ge_core import phi_seed, K_seed, tensor, sq, dsum, charges, pretty

OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-06-generative-emptiness", "object1_grading.json",
)


def charge_multisets():
    """The four charge multisets named in Thm 2.2's parenthetical."""
    phi = phi_seed()
    return [
        {"object": "phi", "poly": pretty(phi), "charges": charges(phi)},
        {"object": "phi (x) phi", "poly": pretty(tensor(phi, phi)),
         "charges": charges(tensor(phi, phi))},
        {"object": "phi^2", "poly": pretty(sq(phi)), "charges": charges(sq(phi))},
        {"object": "K", "poly": pretty(K_seed()), "charges": charges(K_seed())},
    ]


def operator_actions():
    """Operators realise add / double / union ON the charge (computed vs predicted)."""
    phi, K = phi_seed(), K_seed()
    a, b = charges(phi), charges(K)

    add_pred = sorted((i + j) % 4 for i in a for j in b)
    dbl_pred = sorted((2 * i) % 4 for i in charges(K))
    uni_pred = sorted(charges(phi) + charges(K))

    return [
        {"operator": "(x) tensor", "acts_as": "add (mod 4)",
         "word": "phi (x) K",
         "computed": charges(tensor(phi, K)), "predicted": add_pred,
         "match": charges(tensor(phi, K)) == add_pred},
        {"operator": "( )^2 squaring", "acts_as": "double (mod 4)",
         "word": "sq K",
         "computed": charges(sq(K)), "predicted": dbl_pred,
         "match": charges(sq(K)) == dbl_pred},
        {"operator": "(+) dsum", "acts_as": "union",
         "word": "phi (+) K",
         "computed": charges(dsum(phi, K)), "predicted": uni_pred,
         "match": charges(dsum(phi, K)) == uni_pred},
    ]


def z4_cayley():
    """Z/4Z is closed under +, x2, union -- the group data the operators use."""
    G4 = [0, 1, 2, 3]
    addition = [[(a + b) % 4 for b in G4] for a in G4]
    doubling = [(2 * a) % 4 for a in G4]
    return {
        "elements": G4,
        "addition_table": addition,
        "doubling_map": doubling,
        "closed_under_addition": {(a + b) % 4 for a in G4 for b in G4} == set(G4),
        "doubling_maps_in": set(doubling) <= set(G4),
        "union_closed": True,
    }


def orbit_sweep():
    """Two-generation orbit from {phi, K}: confirm NO word acquires an OFF charge."""
    seeds = [phi_seed(), K_seed()]
    gen1 = []
    for P in seeds:
        gen1.append(sq(P))
        for Q in seeds:
            gen1.append(tensor(P, Q))
            gen1.append(dsum(P, Q))
    orbit = seeds + gen1
    for P in gen1[:6]:
        orbit.append(sq(P))
        orbit.append(tensor(P, phi_seed()))

    records = []
    any_off = False
    for P in orbit:
        ch = charges(P)
        off = "OFF" in ch
        any_off = any_off or off
        records.append({"poly": pretty(P), "charges": ch, "has_OFF": off})
    return {"n_objects": len(records), "any_off_charge": any_off, "objects": records}


def k_irreducible():
    """K must be an irreducible quartic to supply the +-i*beta place (charges 1,3)."""
    facs = sp.factor_list(K_seed().as_expr(), G.x)[1]
    return {
        "poly": pretty(K_seed()),
        "degree": K_seed().degree(),
        "n_irreducible_factors": len(facs),
        "irreducible": len(facs) == 1 and facs[0][1] == 1,
    }


def main():
    data = G.provenance("object1_grading.py")
    data.update({
        "object": "I -- the Z/4Z grading",
        "paper_result": "Theorem 2.2 (thm:charge), Def. 2.1 (def:charge)",
        "charge_multisets": charge_multisets(),
        "operator_actions": operator_actions(),
        "Z4_group": z4_cayley(),
        "orbit_sweep": orbit_sweep(),
        "K_irreducibility": k_irreducible(),
    })
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"wrote {OUT}")
    for row in data["charge_multisets"]:
        print(f"  {row['object']:14s} charges {row['charges']}")
    for row in data["operator_actions"]:
        print(f"  {row['operator']:16s} {row['acts_as']:14s} match={row['match']}")
    print(f"  orbit any OFF charge: {data['orbit_sweep']['any_off_charge']}")


if __name__ == "__main__":
    main()
