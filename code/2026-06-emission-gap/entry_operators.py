"""Producer: the entry-level operators (commutator, circulant, Cartan).

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces:
  data/2026-06-emission-gap/entry_operator_scans.json   (Prop. 8.5, 8.6, 12.3, App. A)
  data/2026-06-emission-gap/cartan_spectra.csv          (Prop. 12.4, App. A)
Paper results: a size-2 commutator has charpoly x^2+det, Mahler image {1} U
[2, inf), no value in the band (Prop. 8.5); integer circulants (n=4,5,6) live in
cyclotomic (abelian) fields and emit no Salem number -- a scan of 400 random
circulants yields zero Salem factors (Prop. 12.3); commutators [A,B] over the
catalog companions and their kron / direct-sum products carry no Salem-eligible
factor (Prop. 8.6); Cartan matrices A_n have eigenvalues 2-2cos(k pi/(n+1)) in
[0,4], totally real (Prop. 12.4).
"""
import random

import numpy as np
import sympy as sp

import emgap_core as C

SCRIPT = "entry_operators.py"
x = C.x

# reproduce the paper's random circulant scan deterministically
CIRCULANT_SEED = 20260604
CIRCULANT_TRIALS = 400
CIRCULANT_SIZES = (4, 5, 6)
CIRCULANT_ENTRY = 3          # entries drawn from [-3, 3]


def size2_commutator_image():
    """Traceless 2x2 integer matrix [[a,b],[c,-a]]: charpoly x^2 + det,
    det = -a^2 - bc. Mahler image is {1} U [2, inf) (Prop. 8.5)."""
    seen, any_in_gap = set(), False
    for a in range(-6, 7):
        for b in range(-6, 7):
            for c in range(-6, 7):
                det = -a * a - b * c
                m = C.mahler([1, 0, det])
                seen.add(int(round(float(m))))
                if 1 + 1e-9 < float(m) < 2 - 1e-9:
                    any_in_gap = True
    return {
        "charpoly_form": "x^2 + det, det = -a^2 - bc",
        "distinct_integer_measures": sorted(seen),
        "any_measure_in_open_1_2": any_in_gap,
        "image_is_one_or_ge_2": all(v == 1 or v >= 2 for v in seen),
    }


def circulant_scan():
    """400 random integer circulants (n in {4,5,6}); count irreducible factors,
    reciprocal factors of degree >= 4, and Salem factors (Prop. 12.3)."""
    rng = random.Random(CIRCULANT_SEED)
    total_factors = recip_ge4 = salem_hits = 0
    for _ in range(CIRCULANT_TRIALS):
        n = rng.choice(CIRCULANT_SIZES)
        row = [rng.randint(-CIRCULANT_ENTRY, CIRCULANT_ENTRY) for _ in range(n)]
        M = sp.Matrix([[row[(j - i) % n] for j in range(n)] for i in range(n)])
        cp = sp.Poly(M.charpoly(x).as_expr(), x)
        for fac, _m in sp.factor_list(cp.as_expr(), x)[1]:
            total_factors += 1
            Q = sp.Poly(fac, x)
            if Q.degree() >= 4 and Q.degree() % 2 == 0 and C.is_palindromic(Q.all_coeffs()):
                recip_ge4 += 1
        salem_hits += len(C.salem_factors(cp))
    return {
        "seed": CIRCULANT_SEED, "trials": CIRCULANT_TRIALS,
        "sizes": list(CIRCULANT_SIZES), "entry_range": [-CIRCULANT_ENTRY, CIRCULANT_ENTRY],
        "total_irreducible_factors": total_factors,
        "reciprocal_factors_deg_ge_4": recip_ge4,
        "salem_factors": salem_hits,
    }


def _dsum_sp(A, B):
    n, m = A.shape[0], B.shape[0]
    M = sp.zeros(n + m, n + m)
    M[:n, :n] = A
    M[n:, n:] = B
    return M


def commutator_catalog_scan():
    """Commutators [A,B] over catalog companions and their kron / direct-sum
    products (sizes 2 and 4); count reciprocal deg>=4 and Salem factors (Prop. 8.6)."""
    comps = [C.companion_int(c) for c in C.CATALOG.values()]
    twos = [c for c in comps if c.shape[0] == 2]
    mats = list(comps)
    for A in twos[:3]:
        for B in twos[:3]:
            mats.append(sp.Matrix(np.kron(np.array(A.tolist(), dtype=int),
                                          np.array(B.tolist(), dtype=int)).tolist()))
            mats.append(_dsum_sp(A, B))

    salem_hits = recip_occurrences = commutator_count = 0
    distinct_recip = set()
    for A in mats:
        for B in mats:
            if A.shape != B.shape:
                continue
            commutator_count += 1
            comm = A * B - B * A
            cp = sp.Poly(comm.charpoly(x).as_expr(), x)
            salem_hits += len(C.salem_factors(cp))
            for fac, _m in sp.factor_list(cp.as_expr(), x)[1]:
                Q = sp.Poly(fac, x)
                if Q.degree() >= 4 and Q.degree() % 2 == 0 and C.is_palindromic(Q.all_coeffs()):
                    recip_occurrences += 1
                    distinct_recip.add(str(Q.as_expr()))
    return {
        "generating_slice": "catalog companions + kron/dsum of first three size-2 companions",
        "commutators_evaluated": commutator_count,
        "reciprocal_factor_occurrences_deg_ge_4": recip_occurrences,
        "distinct_reciprocal_factors_deg_ge_4": sorted(distinct_recip),
        "n_distinct_reciprocal_factors_deg_ge_4": len(distinct_recip),
        "salem_factors": salem_hits,
    }


def cartan_rows():
    """Cartan A_n (n=3,5,8): eigenvalues 2 - 2cos(k pi/(n+1)) in [0,4] (Prop. 12.4)."""
    rows = []
    for n in (3, 5, 8):
        M = 2 * np.eye(n)
        for i in range(n - 1):
            M[i, i + 1] = M[i + 1, i] = -1
        ev = sorted(float(v) for v in np.linalg.eigvalsh(M))
        cp = sp.Poly(sp.Matrix(M.astype(int).tolist()).charpoly(x).as_expr(), x)
        rows.append({
            "type": f"A_{n}", "size": n,
            "eigenvalues": ";".join(f"{v:.6f}" for v in ev),
            "all_in_0_4": bool(all(-1e-9 <= v <= 4 + 1e-9 for v in ev)),
            "totally_real": True,
            "has_salem_factor": C.has_salem_factor(cp),
        })
    return rows


def main():
    scans = {
        "size2_commutator_image": size2_commutator_image(),
        "circulant_scan": circulant_scan(),
        "commutator_catalog_scan": commutator_catalog_scan(),
    }
    C.write_json("entry_operator_scans.json", scans, SCRIPT)
    rows = cartan_rows()
    C.write_csv("cartan_spectra.csv",
                ["type", "size", "eigenvalues", "all_in_0_4", "totally_real",
                 "has_salem_factor"], rows, SCRIPT)

    cs = scans["circulant_scan"]
    print(f"circulant_scan: {cs['trials']} circulants -> "
          f"{cs['total_irreducible_factors']} factors, {cs['salem_factors']} Salem")
    ccs = scans["commutator_catalog_scan"]
    print(f"commutator_catalog_scan: {ccs['commutators_evaluated']} commutators, "
          f"{ccs['n_distinct_reciprocal_factors_deg_ge_4']} distinct recip deg>=4 "
          f"({ccs['reciprocal_factor_occurrences_deg_ge_4']} occurrences), "
          f"{ccs['salem_factors']} Salem")
    print(f"size2 commutator image = {scans['size2_commutator_image']['distinct_integer_measures'][:6]}... "
          f"(no value in (1,2): {not scans['size2_commutator_image']['any_measure_in_open_1_2']})")


if __name__ == "__main__":
    main()
