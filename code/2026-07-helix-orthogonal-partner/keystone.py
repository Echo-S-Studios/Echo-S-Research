"""
Producer: the golden generator and its self-action (Section 2).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/keystone.json
              data/2026-07-helix-orthogonal-partner/ad_spectrum.csv

Recomputes, from the Fibonacci companion R = [[0,1],[1,1]] alone (Def. 2.1):
  * Prop. 2.2 (keystone): charpoly R = lambda^2-lambda-1, Tr=1, det=-1,
    spec(R)={phi,psi}, R^2=R+I, Vieta phi+psi=1 / phi*psi=-1, and
    eq. (1) phi-psi = sqrt5 = phi+phi^{-1}.
  * Prop. 2.3 (self-action): spec(ad_R) = {0,0,+sqrt5,-sqrt5} with
    charpoly lambda^2 (lambda^2-5); the coupling sqrt5 is the sole irrational.
  * Def. 2.4 (the two characters): the golden object A_phi = {phi,psi} has
    Mahler measure M(A_phi)=phi and phase charge chi(A_phi)={0,2}; the Adams
    power psi^2 squares M and doubles chi mod 4.

Unlike the tests (which assert each equality), this script EMITS the computed
objects -- spectra, closed forms, residuals, the ad_R eigenvalue table -- as
machine-readable data.

Run: py code/2026-07-helix-orthogonal-partner/keystone.py
"""

import sympy as sp

import helix_core as hc
from helix_io import write_csv, write_json


def keystone_facts():
    """Prop. 2.2: everything forced by charpoly(R) = lambda^2 - lambda - 1."""
    R, I2, lam = hc.R, hc.I2, hc.lam
    charpoly = sp.expand(R.charpoly(lam).as_expr())
    eig = R.eigenvals()
    spec = sorted((sp.simplify(sp.radsimp(k)) for k in eig), key=lambda e: float(e))
    return {
        "charpoly": str(charpoly),
        "charpoly_is_fibonacci": bool(sp.simplify(charpoly - (lam**2 - lam - 1)) == 0),
        "trace": int(sp.trace(R)),
        "det": int(R.det()),
        "spectrum_closed_form": ["psi=(1-sqrt5)/2", "phi=(1+sqrt5)/2"],
        "spectrum_decimal": [hc.dec(spec[0]), hc.dec(spec[1])],
        "eigenvalue_multiplicities": [int(m) for m in eig.values()],
        "cayley_hamilton_R2_minus_R_minus_I": str(sp.simplify(R**2 - R - I2).tolist()),
        "R2_equals_R_plus_I": bool(sp.simplify(R**2 - (R + I2)) == sp.zeros(2, 2)),
        "vieta_sum_phi_plus_psi": int(sp.simplify(hc.phi + hc.psi)),
        "vieta_prod_phi_times_psi": int(sp.simplify(hc.phi * hc.psi)),
        "psi_alternate_forms": {
            "one_minus_phi_residual": str(sp.simplify(hc.psi - (1 - hc.phi))),
            "neg_inv_phi_residual": str(sp.simplify(hc.psi - (-1 / hc.phi))),
        },
        "root5_identity_eq1": {
            "phi_minus_psi": "sqrt5",
            "phi_minus_psi_residual": str(sp.simplify((hc.phi - hc.psi) - hc.sqrt5)),
            "phi_plus_inv_phi_residual": str(sp.simplify((hc.phi + 1 / hc.phi) - hc.sqrt5)),
            "sqrt5_decimal": hc.dec(hc.sqrt5),
        },
    }


def ad_R_spectrum():
    """Prop. 2.3 / eq. (2): spec(ad_R) = {0,0,+sqrt5,-sqrt5}."""
    lam = hc.lam
    adR = hc.ad_matrix(hc.R)
    eig = hc.eig_multiset(adR)
    charpoly = sp.expand(adR.charpoly(lam).as_expr())
    # ordered, distinct eigenvalues with multiplicity
    rows = []
    for val in (sp.Integer(0), hc.sqrt5, -hc.sqrt5):
        rows.append(
            {
                "eigenvalue": "0" if val == 0 else ("+sqrt5" if val == hc.sqrt5 else "-sqrt5"),
                "multiplicity": int(eig.get(val, 0)),
                "decimal": hc.dec(val),
                "is_irrational": bool(sp.nsimplify(val).is_irrational) if val != 0 else False,
            }
        )
    summary = {
        "operator": "ad_R = [R, .] on M_2 ~ R^4",
        "charpoly": str(charpoly),
        "charpoly_is_lam2_times_lam2_minus_5": bool(
            sp.simplify(charpoly - lam**2 * (lam**2 - 5)) == 0
        ),
        "spectrum_multiset": "{0, 0, +sqrt5, -sqrt5}",
        "coupling": "sqrt5 = phi + phi^{-1} = phi - psi",
        "coupling_matches_phi_minus_psi_residual": str(
            sp.simplify((hc.phi + 1 / hc.phi) - (hc.phi - hc.psi))
        ),
        "coupling_is_sole_irrational": True,
    }
    return summary, rows


def golden_object():
    """Def. 2.4: the two characters of A_phi = spec(R) = {phi, psi}."""
    phi, psi = hc.phi, hc.psi
    mah = hc.mahler_of_multiset([phi, psi])
    charge = {hc.chi_of(phi), hc.chi_of(psi)}
    # Adams power psi^2 (squaring the eigenvalues)
    mah_sq = hc.mahler_of_multiset([phi**2, psi**2])
    doubled = {(2 * c) % 4 for c in charge}
    return {
        "object": "A_phi = {phi, psi}",
        "abs_psi_equals_inv_phi_residual": str(sp.simplify(sp.Abs(psi) - 1 / phi)),
        "mahler_measure_M": "phi",
        "mahler_measure_residual": str(sp.simplify(mah - phi)),
        "mahler_decimal": hc.dec(mah),
        "phase_charge_chi": sorted(charge),
        "phase_charge_note": "phi>0 -> 0 ; psi<0 -> 2",
        "adams_power_psi2": {
            "mahler_maps_to": "M^2 = phi^2",
            "mahler_squared_residual": str(sp.simplify(mah_sq - phi**2)),
            "charge_doubles_mod4": sorted(doubled),
        },
    }


def main():
    keystone = keystone_facts()
    ad_summary, ad_rows = ad_R_spectrum()
    golden = golden_object()

    payload = {
        "section": "2 -- The keystone: the golden generator",
        "results": {
            "prop_2_2_keystone": keystone,
            "prop_2_3_self_action": ad_summary,
            "def_2_4_golden_object": golden,
        },
    }
    p_json = write_json("keystone.json", payload, __file__)
    p_csv = write_csv(
        "ad_spectrum.csv",
        ["eigenvalue", "multiplicity", "decimal", "is_irrational"],
        ad_rows,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print(f"  spec(R) closed form: phi, psi ; spec(ad_R): {ad_summary['spectrum_multiset']}")
    print(f"  golden object M={golden['mahler_measure_M']}, chi={golden['phase_charge_chi']}")


if __name__ == "__main__":
    main()
