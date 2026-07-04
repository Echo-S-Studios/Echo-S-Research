"""Producer: the exhaustive quintic census (Section 6, stage 1).

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
Produces: the Section 6 certification table and pattern split.  Over the box of
monic quintics x^5+ax^4+bx^3+cx^2+dx+e with (a,...,e) in [-2,2]^5 (3125
candidates), the certification cascade rejects, in order:

    e=0 (625) | +-reciprocal (50) | reducible (638) | real pattern (1318) |
    disk certificate (411),

certifying 83 Pisot quintics, whose non-real pattern splits
    real5 = 0 | mixed (one pair) = 16 | two-pair = 67,
with the first two-pair instance in enumeration order x^5-2x^4-2x^3-2x^2-2x-2
(Theorem 6.1).

Pisot classification uses exact irreducibility plus high-precision root moduli;
this is sound because the certified set is irreducible and non-reciprocal, hence
has NO root on the unit circle (a genuine modulus gap to classify against).

Emits:
  data/2026-07-pisot-residue/quintic_census_summary.json   (the reject table + split)
  data/2026-07-pisot-residue/quintic_pisots.csv            (the 83 certified Pisots)
"""
import mpmath as mp
from sympy import symbols, Poly

from pisot_lib import is_pm_reciprocal, hp_roots, poly_str, write_json, write_csv

x = symbols('x')

PATTERN_NAME = {0: "real5", 1: "mixed", 2: "two-pair"}


def run_census(dps=40):
    """The exact stage-1 certification cascade over [-2,2]^5."""
    tol = mp.mpf(10) ** -12
    rej = dict(e0=0, recip=0, reducible=0, realpat=0, disk=0)
    patt = dict(real5=0, mixed=0, twopair=0)
    pisots = []                     # (coeffs, theta, n_real_in, n_pairs, pattern)
    first_twopair = None
    order = 0
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                for d in range(-2, 3):
                    for e in range(-2, 3):
                        coeffs = [1, a, b, c, d, e]
                        if e == 0:
                            rej['e0'] += 1
                            continue
                        if is_pm_reciprocal(coeffs):
                            rej['recip'] += 1
                            continue
                        if not Poly(coeffs, x).is_irreducible:
                            rej['reducible'] += 1
                            continue
                        rts = hp_roots(coeffs, dps)
                        reals = [r.real for r in rts if abs(r.imag) < tol]
                        gt1 = [r for r in reals if r > 1]
                        lem1 = [r for r in reals if r <= -1]
                        if len(gt1) != 1 or len(lem1) > 0:
                            rej['realpat'] += 1
                            continue
                        oncirc = any(abs(abs(r) - 1) < tol for r in rts)
                        outside = [r for r in rts if abs(r) > 1]
                        if oncirc or len(outside) != 1 or abs(outside[0].imag) > tol:
                            rej['disk'] += 1
                            continue
                        # certified Pisot
                        order += 1
                        theta = float(gt1[0])
                        nonreal = [r for r in rts if abs(r.imag) > tol]
                        npairs = len(nonreal) // 2
                        n_real_in = len([r for r in rts
                                         if abs(r.imag) < tol and -1 < r.real < 1])
                        if npairs == 0:
                            patt['real5'] += 1
                        elif npairs == 1:
                            patt['mixed'] += 1
                        elif npairs == 2:
                            patt['twopair'] += 1
                            if first_twopair is None:
                                first_twopair = coeffs
                        pisots.append((order, coeffs, theta, n_real_in, npairs,
                                       PATTERN_NAME[npairs]))
    return rej, patt, pisots, first_twopair


def main():
    rej, patt, pisots, first_twopair = run_census()
    total = sum(rej.values()) + len(pisots)

    summary = {
        "description": "Section 6 stage-1 certification: 3125 monic quintics in [-2,2]^5, "
                       "reject cascade and non-real pattern split.",
        "box": {"coeff_range": [-2, 2], "num_free_coeffs": 5, "size": 5 ** 5},
        "reject_tally": {
            "e_equals_0": rej['e0'],
            "pm_reciprocal": rej['recip'],
            "reducible": rej['reducible'],
            "real_pattern": rej['realpat'],
            "disk_certificate": rej['disk'],
        },
        "certified_pisot": len(pisots),
        "cascade_partitions_box": total == 3125,
        "pattern_split": {
            "real5": patt['real5'],
            "mixed": patt['mixed'],
            "two_pair": patt['twopair'],
            "mixed_plus_two_pair_equals_pisot": patt['mixed'] + patt['twopair'] == len(pisots),
        },
        "first_two_pair_instance": {
            "poly": poly_str(first_twopair),
            "coeffs_hi_to_lo": first_twopair,
        },
    }
    jpath = write_json("quintic_census_summary.json", summary, "quintic_census.py")

    rows = [(o, poly_str(c).replace(",", ";"),
             "[" + " ".join(str(int(v)) for v in c) + "]",
             f"{th:.6f}", nri, npairs, pat)
            for (o, c, th, nri, npairs, pat) in pisots]
    cpath = write_csv(
        "quintic_pisots.csv",
        ["order", "poly", "coeffs_hi_to_lo", "theta_display", "n_real_inside",
         "n_nonreal_pairs", "pattern"],
        rows, "quintic_census.py")

    print(f"wrote {jpath}")
    print(f"wrote {cpath}")
    print(f"  rejects e0/recip/reducible/realpat/disk = "
          f"{rej['e0']}/{rej['recip']}/{rej['reducible']}/{rej['realpat']}/{rej['disk']}")
    print(f"  certified Pisot = {len(pisots)}; partitions box = {total == 3125}")
    print(f"  patterns real5/mixed/two-pair = "
          f"{patt['real5']}/{patt['mixed']}/{patt['twopair']}")
    print(f"  first two-pair instance = {poly_str(first_twopair)}")


if __name__ == "__main__":
    main()
