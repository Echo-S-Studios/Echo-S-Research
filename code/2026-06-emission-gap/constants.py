"""Producer: the numeric constants of the Emission-Gap paper.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces: data/2026-06-emission-gap/constants.json
Paper results: the constants named in the abstract, Sec. 1, Cor. 6.1, Cor. 10.4,
and Appendix A -- the plastic number mu_S (Smyth's floor), the golden ratio phi
and its powers/log, Lehmer's number, the minimal degree-4 Salem number beta_4,
and the channel endpoint sqrt5. Each constant is rebuilt from its DEFINING
polynomial / expression at 50 dps, then emitted with its closed form.
"""
import mpmath as mp

import emgap_core as C

SCRIPT = "constants.py"


def compute():
    phi = (1 + mp.sqrt(5)) / 2
    mu_S = mp.findroot(lambda z: z**3 - z - 1, mp.mpf("1.324"))
    beta4 = mp.findroot(lambda z: z**4 - z**3 - z**2 - z + 1, mp.mpf("1.72"))
    lehmer = C.mahler(C.LEHMER)                       # Mahler = largest root
    sqrt5 = mp.sqrt(5)

    constants = {
        "mu_S": {
            "value": C.s(mu_S), "float": float(mu_S),
            "definition": "real root of x^3 - x - 1 (plastic number)",
            "paper_location": "Sec. 1 (Smyth bound), App. A",
            "role": "non-reciprocal Mahler floor; right end of the Salem band (1, mu_S)",
        },
        "log_mu_S": {
            "value": C.s(mp.log(mu_S)), "float": float(mp.log(mu_S)),
            "definition": "log of the plastic number, in nats",
            "paper_location": "Cor. 6.1, App. A",
            "role": "cost floor lambda*log M >= log mu_S",
        },
        "phi": {
            "value": C.s(phi), "float": float(phi),
            "definition": "golden ratio (1+sqrt5)/2, root of x^2 - x - 1",
            "paper_location": "Cor. 6.1, App. A",
            "role": "smallest realised catalog Mahler measure",
        },
        "log_phi": {
            "value": C.s(mp.log(phi)), "float": float(mp.log(phi)),
            "definition": "log of the golden ratio, in nats",
            "paper_location": "Cor. 6.1, App. A",
            "role": "realised cost floor (entropy end of the gap)",
        },
        "phi_squared": {
            "value": C.s(phi**2), "float": float(phi**2),
            "definition": "phi^2 = phi + 1 = (3+sqrt5)/2",
            "paper_location": "App. A",
            "role": "reciprocal quadratic Mahler jump / squaring image",
        },
        "phi_fourth": {
            "value": C.s(phi**4), "float": float(phi**4),
            "definition": "phi^4 = (7+3 sqrt5)/2, the 'gap' seed x^2-7x+1 measure",
            "paper_location": "App. A",
            "role": "Mahler measure of the gap seed",
        },
        "lehmer_number": {
            "value": C.s(lehmer), "float": float(lehmer),
            "definition": ("largest root and Mahler measure of the Lehmer polynomial "
                           "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1"),
            "paper_location": "Sec. 1, Lemma 5.1",
            "role": "smallest known Salem number; sits below mu_S inside the Salem band",
        },
        "beta_4": {
            "value": C.s(beta4), "float": float(beta4),
            "definition": "real root of x^4 - x^3 - x^2 - x + 1 (minimal degree-4 Salem number)",
            "paper_location": "Cor. 10.4, App. A",
            "role": "floor for any Salem number housed in K; exceeds phi",
        },
        "sqrt5": {
            "value": C.s(sqrt5), "float": float(sqrt5),
            "definition": "sqrt(disc(x^2-x-1)) = sqrt5; nonzero ad_R eigenvalue",
            "paper_location": "Sec. 7, Lemma 10.1, App. A",
            "role": "channel-gap endpoint (distinct from phi and log phi)",
        },
    }

    relations = {
        "lehmer < mu_S": bool(lehmer < mu_S),
        "mu_S < 2": bool(mu_S < 2),
        "phi < beta_4": bool(phi < beta4),
        "beta_4 exceeds phi": bool(beta4 > phi),
        "phi^2 = phi + 1": bool(abs(phi**2 - (phi + 1)) < mp.mpf(10) ** (-40)),
        "phi^4 = (7+3 sqrt5)/2": bool(abs(phi**4 - (7 + 3 * mp.sqrt(5)) / 2) < mp.mpf(10) ** (-40)),
        "exp(log phi) = phi": bool(abs(mp.exp(mp.log(phi)) - phi) < mp.mpf(10) ** (-40)),
        "salem_band": f"(1, {C.s(mu_S, 12)})",
    }

    return {"constants": constants, "relations": relations}


def main():
    payload = compute()
    path = C.write_json("constants.json", payload, SCRIPT)
    print(f"wrote {path}")
    for name, d in payload["constants"].items():
        print(f"  {name:16s} = {d['float']!r}")


if __name__ == "__main__":
    main()
