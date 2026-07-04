"""Producer: the degree-12 Salem census (Section 7.1).

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the Section 7.1 census over palindromic degree-12 monic polynomials
with free coefficients c1..c6 in {-1,0,1} (729 vectors).  Under the twist
x -> -x there are 27 twist-fixed vectors, hence Burnside gives (729+27)/2 = 378
twist-classes.  The certification cascade is (CORRECTED 2026-07-04):

    39 (+-1 root) / 256 (trace-Sturm reject) / 46 (reducible) / 37 Salem,

with 302 = 256 + 46 combined non-Salem-non-(+-1) classes.  The Phi_10-factored
class x^12-x^11-x^10-x^9-x^7-x^6-x^5-x^3-x^2-x+1 carries a coincidental (1,0,5)
trace pattern yet is reducible (a cyclotomic factor makes it reducible
regardless), and is counted as reducible -- this is the 257/45 -> 256/46
reclassification.  Every Salem twist-class scans to {Phi_1^12} (Cor. 7.14):
deg Rat = 144, bound 2*12^4 = 41472.

Emits:
  data/2026-07-pisot-residue/salem_census_summary.json
  data/2026-07-pisot-residue/salem_reps.csv           (the 37 Salem twist-class reps)
"""
import itertools

import mpmath as mp
import sympy as sp
from sympy import symbols, Poly, oo, factor_list

from pisot_lib import (rat_object, cyclotomic_scan, hp_roots, scan_json,
                       poly_str, write_json, write_csv)

x, z = symbols('x z')

# how many Salem representatives to scan end-to-end (each deg-144 Rat+factor ~1s)
SCAN_SAMPLE = 5


def build(cvec):
    """Palindromic degree-12 coefficient list c[k] = coeff of x^k."""
    c = [0] * 13
    c[0] = c[12] = 1
    for i in range(1, 7):
        c[i] = cvec[i - 1]
        c[12 - i] = cvec[i - 1]
    return c


def hi2lo(c):
    return [c[k] for k in range(12, -1, -1)]


def twist(v):
    return tuple(((-1) ** i) * v[i - 1] for i in range(1, 7))


def has_pm1(c):
    return sum(c) == 0 or sum(c[k] * (-1) ** k for k in range(13)) == 0


def trace_poly(c):
    """Degree-6 trace polynomial t(z), z = x + 1/x, with p(x) = x^6 t(x+1/x)."""
    m = 6
    Pj = [None] * (m + 1)
    Pj[0] = sp.Integer(2)
    Pj[1] = z
    for j in range(2, m + 1):
        Pj[j] = sp.expand(z * Pj[j - 1] - Pj[j - 2])
    t = sp.Integer(c[m])
    for j in range(1, m + 1):
        t = t + c[m + j] * Pj[j]
    return Poly(sp.expand(t), z)


def survivor(c):
    """Exact trace-Sturm survivor: trace poly has pattern
    (>2 : 1 root, <-2 : 0 roots, in(-2,2) : 5 roots)."""
    t = trace_poly(c)
    if t.degree() != 6:
        return False
    return (t.count_roots(2, oo) == 1
            and t.count_roots(-oo, -2) == 0
            and t.count_roots(-2, 2) == 5)


def is_salem_poly(c):
    """Salem minimal polynomial: irreducible, one root real > 1, one strictly
    inside, remaining 10 on the unit circle."""
    if not Poly(hi2lo(c), x).is_irreducible:
        return False
    with mp.workdps(40):
        eps = mp.mpf(10) ** -12
        rr = hp_roots(hi2lo(c), 40)
        outside = [r for r in rr if abs(r) > 1 + eps]
        oncirc = [r for r in rr if abs(abs(r) - 1) < eps]
        inside = [r for r in rr if abs(r) < 1 - eps]
        if len(outside) != 1 or abs(outside[0].imag) > eps or outside[0].real <= 1:
            return False
        return len(oncirc) == 10 and len(inside) == 1


def orbit_reps():
    """One representative per twist-orbit (378 of them)."""
    seen, reps = set(), []
    for v in itertools.product((-1, 0, 1), repeat=6):
        if v in seen:
            continue
        seen |= {v, twist(v)}
        reps.append(v)
    return reps


def salem_theta(c):
    with mp.workdps(40):
        rr = hp_roots(hi2lo(c), 40)
        return float(max(r.real for r in rr if abs(r) > 1))


def has_cyclotomic_factor(c):
    """Detect a cyclotomic (unit-circle) factor: some irreducible factor whose
    roots all lie on the unit circle -> palindromic with all roots |.|=1."""
    for fac, _ in factor_list(Poly(hi2lo(c), x))[1]:
        rr = hp_roots([int(k) for k in Poly(fac, x).all_coeffs()], 30)
        if rr and all(abs(abs(r) - 1) < mp.mpf(10) ** -8 for r in rr):
            return True, poly_str([int(k) for k in Poly(fac, x).all_coeffs()])
    return False, None


def run_census():
    reps = orbit_reps()
    vecs = list(itertools.product((-1, 0, 1), repeat=6))
    twist_fixed = [v for v in vecs if twist(v) == v]

    pm1 = reject = reducible = salem = 0
    salem_reps = []          # (coeffs_hi2lo, theta)
    reducible_survivors = []
    phi10_class = None
    for v in reps:
        c = build(v)
        if has_pm1(c):
            pm1 += 1
            continue
        ct = build(twist(v))
        sc, sct = survivor(c), survivor(ct)
        if not (sc or sct):
            reject += 1
            continue
        sel = c if sc else ct
        if Poly(hi2lo(sel), x).is_irreducible:
            salem += 1
            salem_reps.append((hi2lo(sel), salem_theta(sel)))
        else:
            reducible += 1
            cyc, fac = has_cyclotomic_factor(sel)
            reducible_survivors.append(hi2lo(sel))
            if cyc and phi10_class is None:
                phi10_class = {"poly": poly_str(hi2lo(sel)),
                               "coeffs_hi_to_lo": hi2lo(sel),
                               "cyclotomic_factor": fac}
    return dict(num_vecs=len(vecs), twist_fixed=len(twist_fixed),
                orbits=len(reps), pm1=pm1, reject=reject, reducible=reducible,
                salem=salem, salem_reps=salem_reps, phi10_class=phi10_class)


def main():
    r = run_census()

    # scan a sample of Salem representatives end-to-end -> {Phi_1^12}
    scans = []
    for coeffs, _ in r["salem_reps"][:SCAN_SAMPLE]:
        R = rat_object(coeffs)
        scan = cyclotomic_scan(R)
        scans.append({"poly": poly_str(coeffs), "deg_Rat": R.degree(),
                      "scan": scan_json(scan), "inert": scan == {1: 12}})

    summary = {
        "description": "Section 7.1 degree-12 Salem census; Burnside 729/27/378, "
                       "corrected cascade 39/256/46/37, 302 combined, Salem scans {Phi_1^12}.",
        "family": {"palindromic_vectors": r["num_vecs"], "twist_fixed": r["twist_fixed"],
                   "burnside_orbits": (r["num_vecs"] + r["twist_fixed"]) // 2,
                   "orbits_enumerated": r["orbits"]},
        "cascade": {
            "pm1_root": r["pm1"],
            "trace_sturm_reject": r["reject"],
            "reducible": r["reducible"],
            "salem": r["salem"],
            "combined_reject_plus_reducible": r["reject"] + r["reducible"],
            "partitions_378": r["pm1"] + r["reject"] + r["reducible"] + r["salem"] == 378,
        },
        "phi10_reducible_class": r["phi10_class"],
        "scan_sizing": {"deg_Rat": 144, "bound_2n4": 2 * 12 ** 4},
        "salem_rep_scans_sample": scans,
        "all_sampled_scans_are_phi1_12": all(s["inert"] for s in scans),
    }
    jpath = write_json("salem_census_summary.json", summary, "salem_census.py")

    rows = [("[" + " ".join(str(int(v)) for v in c) + "]",
             poly_str(c).replace(",", ";"), f"{th:.6f}")
            for (c, th) in r["salem_reps"]]
    cpath = write_csv("salem_reps.csv",
                      ["coeffs_hi_to_lo", "poly", "theta_display"],
                      rows, "salem_census.py")

    print(f"wrote {jpath}")
    print(f"wrote {cpath}")
    print(f"  Burnside 729/27/378 -> orbits={r['orbits']}")
    print(f"  cascade pm1/reject/reducible/salem = "
          f"{r['pm1']}/{r['reject']}/{r['reducible']}/{r['salem']} "
          f"(combined 302 = {r['reject'] + r['reducible']})")
    print(f"  sampled Salem scans all {{Phi_1^12}}: {summary['all_sampled_scans_are_phi1_12']}")


if __name__ == "__main__":
    main()
