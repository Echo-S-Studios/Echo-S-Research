"""
Producer: the named numeric constants of *Lehmer's Box*.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Produces the paper's load-bearing constants -- phi (floor), tau = phi^{-1},
mu_S (Smyth / plastic / smallest Pisot), the two floor logarithms log mu_S and
log phi (Cor. 4.9 / eq. floorval), Lehmer's number Mah(L) (Sec. 1.2), the
minimal degree-four Salem beta_4 (Cor. 5.5), the "gap" seed value phi^4
(Def. 2.9), and the off-circle complex place of the Lorentzian quartic (Sec. 5).

Each constant is DERIVED here from its defining polynomial / closed form (not
copied from the paper) and emitted with full provenance.

Outputs:
    data/2026-06-lehmers-box/constants.json
    data/2026-06-lehmers-box/constants.csv
Backs: abstract, Sec. 1.2, Def. 2.9, Cor. 4.9 (eq. floorval), Cor. 5.5, Sec. 5.
"""

from __future__ import annotations

import mpmath as mp

import box_core as C
from box_io import write_csv, write_json

mp.mp.dps = 60
_DP = 40                                    # digits emitted per value


def _s(x):
    """Fixed high-precision decimal string."""
    return mp.nstr(x, _DP)


def compute():
    phi = (1 + mp.sqrt(5)) / 2
    tau = (mp.sqrt(5) - 1) / 2
    muS = mp.findroot(lambda t: t ** 3 - t - 1, mp.mpf("1.3"))
    lehmer = C.mahler_measure(C.LEHMER)
    beta4 = C.dominant_real_root(C.BETA4)
    gap = phi ** 4
    # imaginary root modulus of K = x^4 + 5x^2 - 5  (the Lorentzian complex place)
    imag_place = mp.sqrt((5 + 3 * mp.sqrt(5)) / 2)
    fifth_root = mp.mpf(5) ** mp.mpf("0.25")

    records = [
        dict(key="phi", symbol="varphi",
             description="golden ratio; the floor wall of Lehmer's Box",
             minimal_polynomial="x^2 - x - 1", closed_form="(1+sqrt5)/2",
             value=_s(phi), paper_location="abstract; Def. 2.9; Sec. 1.3",
             role="floor wall"),
        dict(key="tau", symbol="tau",
             description="tau = phi^{-1}, smaller root of x^2 + x - 1",
             minimal_polynomial="x^2 + x - 1", closed_form="(-1+sqrt5)/2",
             value=_s(tau), paper_location="Def. 2.9",
             role="catalog seed (tau)"),
        dict(key="mu_S", symbol="mu_S",
             description="plastic number; smallest Pisot; Smyth non-reciprocal bound",
             minimal_polynomial="x^3 - x - 1", closed_form="real root of x^3-x-1",
             value=_s(muS), paper_location="Sec. 1.2; Lem. 2.7",
             role="Smyth floor mu_S"),
        dict(key="log_mu_S", symbol="log mu_S",
             description="natural log of the Smyth floor (nats)",
             minimal_polynomial="", closed_form="ln(mu_S)",
             value=_s(mp.log(muS)), paper_location="Cor. 4.9 (eq. floorval)",
             role="cost floor (Smyth)"),
        dict(key="log_phi", symbol="log varphi",
             description="natural log of the realised golden floor (nats)",
             minimal_polynomial="", closed_form="ln(phi)",
             value=_s(mp.log(phi)), paper_location="Cor. 4.9 (eq. floorval)",
             role="cost floor (realised)"),
        dict(key="Mah_L", symbol="Mah(L)",
             description="Mahler measure of Lehmer's number L (smallest known)",
             minimal_polynomial="x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1",
             closed_form="dominant root of L (Salem)",
             value=_s(lehmer), paper_location="Sec. 1.2",
             role="excluded emission (Salem, sub-mu_S)"),
        dict(key="beta_4", symbol="beta_4",
             description="minimal degree-four Salem number; beta_4 > phi",
             minimal_polynomial="x^4 - x^3 - x^2 - x + 1",
             closed_form="dominant root (Salem)",
             value=_s(beta4), paper_location="Cor. 5.5",
             role="minimal Salem in K, above the floor"),
        dict(key="gap", symbol="phi^4",
             description="'gap' seed value; roots phi^4, phi^-4 of x^2-7x+1",
             minimal_polynomial="x^2 - 7x + 1", closed_form="phi^4",
             value=_s(gap), paper_location="Def. 2.9; Rem. 3.3",
             role="catalog seed (gap)"),
        dict(key="imag_place_K", symbol="|i beta|",
             description="modulus of the off-circle complex place of Q(5^{1/4})",
             minimal_polynomial="x^4 + 5x^2 - 5 (imag root)",
             closed_form="sqrt((5+3sqrt5)/2)",
             value=_s(imag_place), paper_location="Sec. 5 (closing)",
             role="complex place OFF the unit circle"),
        dict(key="fourth_root_5", symbol="5^{1/4}",
             description="fourth root of 5; other off-circle modulus in Q(5^{1/4})",
             minimal_polynomial="x^4 - 5", closed_form="5^{1/4}",
             value=_s(fifth_root), paper_location="Sec. 5 (closing)",
             role="complex place OFF the unit circle"),
    ]

    ordering = {
        "Mah(L) < mu_S": bool(lehmer < muS),
        "mu_S < phi": bool(muS < phi),
        "phi < beta_4": bool(phi < beta4),
        "0 < log_mu_S < log_phi": bool(0 < mp.log(muS) < mp.log(phi)),
    }
    return records, ordering


def main():
    records, ordering = compute()

    csv_cols = ["key", "symbol", "description", "minimal_polynomial",
                "closed_form", "value", "paper_location", "role"]
    p_csv = write_csv("constants.csv", csv_cols, records, __file__)

    payload = {
        "_description": "Named numeric constants of Lehmer's Box, each derived "
                        "from its defining polynomial or closed form.",
        "_precision_digits": _DP,
        "constants": records,
        "ordering_facts": ordering,
    }
    p_json = write_json("constants.json", payload, __file__)

    print("wrote", p_csv)
    print("wrote", p_json)
    for r in records:
        print(f"  {r['key']:14s} = {r['value']}")
    for k, v in ordering.items():
        print(f"  order: {k:24s} -> {v}")


if __name__ == "__main__":
    main()
