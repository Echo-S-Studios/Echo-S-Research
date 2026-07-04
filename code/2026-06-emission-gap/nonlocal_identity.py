"""Producer: the non-local identity -- one gap, four domains.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces: data/2026-06-emission-gap/nonlocal_gap_endpoints.csv
Paper results: the four forms of the Emission-Gap Theorem are one "excluded
middle" read across four domains (Sec. 7 table). The gaps share EMPTINESS, not a
number: the height end phi, the entropy end log phi (related by the
Lind-Schmidt-Ward identity h = log M), and the channel end sqrt5 (an ad_R
eigenvalue) are three DISTINCT quantities (Remark 7.2). This producer emits the
four-domain table together with the distinct right endpoints and the h = log M map.
"""
import mpmath as mp

import emgap_core as C

SCRIPT = "nonlocal_identity.py"
PHI = (1 + mp.sqrt(5)) / 2


def rows():
    s5 = mp.sqrt(5)
    logphi = mp.log(PHI)
    # the four-domain table (Sec. 7): domain, gap, right endpoint value, equivalence
    data = [
        {"domain": "height (Mahler)", "gap": "(1, phi)",
         "right_endpoint_symbol": "phi", "right_endpoint_value": C.s(PHI, 16),
         "equivalence_to_form_i": "definitional (empty for integer quadratics)"},
        {"domain": "dynamics (entropy h=log M)", "gap": "(0, log phi)",
         "right_endpoint_symbol": "log phi", "right_endpoint_value": C.s(logphi, 16),
         "equivalence_to_form_i": "LSW: M = e^h (image of the height gap under log)"},
        {"domain": "decision (self-action ad_R)", "gap": "(0, sqrt5)",
         "right_endpoint_symbol": "sqrt5", "right_endpoint_value": C.s(s5, 16),
         "equivalence_to_form_i": "spectrum ad_R = {-sqrt5, 0, sqrt5}; no fourth channel"},
        {"domain": "geometry (trace form)", "gap": "Riemannian <-> Lorentzian at D=0",
         "right_endpoint_symbol": "D=0", "right_endpoint_value": "signature split",
         "equivalence_to_form_i": "Salem => indefinite (2,m-1) trace form"},
    ]
    return data


def endpoint_summary():
    s5, logphi = mp.sqrt(5), mp.log(PHI)
    return {
        "channel_end_sqrt5": C.s(s5, 20),
        "height_end_phi": C.s(PHI, 20),
        "entropy_end_log_phi": C.s(logphi, 20),
        "three_endpoints_distinct": bool(
            abs(s5 - PHI) > mp.mpf("0.5") and abs(PHI - logphi) > 1 and abs(s5 - logphi) > mp.mpf("1.5")),
        "h_equals_log_M": {
            "identity": "h = log M (Lind-Schmidt-Ward)",
            "exp_log_phi_equals_phi": bool(abs(mp.exp(logphi) - PHI) < mp.mpf(10) ** (-30)),
            "golden_seed_entropy": C.s(mp.log(C.mahler(C.CATALOG["phi"])), 16),
        },
        "note": "the correspondence across domains is of emptiness, not a shared endpoint",
    }


def main():
    data = rows()
    C.write_csv("nonlocal_gap_endpoints.csv",
                ["domain", "gap", "right_endpoint_symbol", "right_endpoint_value",
                 "equivalence_to_form_i"], data, SCRIPT)
    # append the endpoint-distinctness summary as a JSON sidecar for provenance
    C.write_json("nonlocal_gap_summary.json", endpoint_summary(), SCRIPT)
    summ = endpoint_summary()
    print("nonlocal_gap_endpoints.csv: four domains, three distinct right endpoints")
    print(f"  channel sqrt5 = {summ['channel_end_sqrt5'][:9]}, "
          f"height phi = {summ['height_end_phi'][:9]}, "
          f"entropy log phi = {summ['entropy_end_log_phi'][:9]}")
    print(f"  distinct: {summ['three_endpoints_distinct']}")


if __name__ == "__main__":
    main()
