"""
Producer for Section 6 -- Self-action: the three powers, Proposition 6.2
(prop:rates) [OA-IT-01..04; GR-table] of
    papers/2026-06-operator-algebra/operator-algebra-whitepaper.tex

From the golden seed (degree 2, M = phi) the three ways an operator acts on
itself are the semiring's three powers:
    additive        (+)^k          : degree 2k , M = phi^k , density (1/2)logphi
    Adams        (psi^2)^{o k}      : degree 2  , M = phi^{2^k}
    multiplicative  (x)^k           : degree 2^k, M = phi^{S_k}
with  S_k = k * C(k-1, floor((k-1)/2))  = 1,2,6,12,30,60,140,...
      S_k ~ 2^k sqrt(k/2pi) , and S_k = the mean-absolute-deviation sum of the
      symmetric binomial.

Emits:
  * data/2026-06-operator-algebra/self_action_rates.csv -- k=1..7: S_k by three
    independent formulas, and each power's degree + Mahler measure (measures
    re-verified on the actual objects for k<=5).
  * data/2026-06-operator-algebra/sk_asymptotic.csv -- S_k / (2^k sqrt(k/2pi))
    -> 1 for growing k.

Difference from tests/: verifier asserts the sequence and measures; this
producer emits the S_k table, the per-power degrees/measures, and the
asymptotic ratio as data.
"""

from __future__ import annotations

from math import comb

import mpmath as mp
import sympy as sp

from opalg_core import (
    golden_seed,
    mahler_exact,
    oplus,
    otimes,
    phi,
    psi,
    write_csv,
)


def Sk_closed(k):
    """Paper's closed form S_k = k * C(k-1, floor((k-1)/2))."""
    return k * comb(k - 1, (k - 1) // 2)


def Sk_from_tensor_spectrum(k):
    """Exponent of phi in M((x)^k golden): word choosing phi j times gives
    phi^{2j-k} with multiplicity C(k,j); only |.|>1 enters M."""
    return sum(comb(k, j) * (2 * j - k) for j in range(k + 1) if 2 * j - k > 0)


def Sk_mad_sum(k):
    """Total absolute-deviation sum of Binomial(k, 1/2)."""
    return sum(sp.Rational(comb(k, j)) * abs(sp.Rational(2 * j - k, 2))
               for j in range(k + 1))


def main():
    G = golden_seed()

    rows = []
    # objects re-verified up to k=5 (degree 2^5 = 32 stays cheap)
    add_obj, adams_obj, tensor_obj = G, G, G
    for k in range(1, 8):
        s_closed = Sk_closed(k)
        s_tensor = Sk_from_tensor_spectrum(k)
        s_mad = int(Sk_mad_sum(k))
        agree = s_closed == s_tensor == s_mad

        verified = None
        if k <= 5:
            # (+)^k, (psi^2)^{o k}, (x)^k of the golden seed, built incrementally
            if k > 1:
                add_obj = oplus(add_obj, G)
                tensor_obj = otimes(tensor_obj, G)
            adams_obj = psi(2, adams_obj)  # one more squaring each step
            add_ok = sp.simplify(mahler_exact(add_obj) - phi**k) == 0 and len(add_obj) == 2 * k
            adams_ok = sp.simplify(mahler_exact(adams_obj) - phi ** (2**k)) == 0 and len(adams_obj) == 2
            tensor_ok = sp.simplify(mahler_exact(tensor_obj) - phi**s_closed) == 0 and len(tensor_obj) == 2**k
            verified = bool(add_ok and adams_ok and tensor_ok)

        rows.append(
            {
                "k": k,
                "S_k_closed_form": s_closed,
                "S_k_tensor_spectrum": s_tensor,
                "S_k_mad_sum": s_mad,
                "S_k_three_formulas_agree": agree,
                "additive_degree": 2 * k,
                "additive_measure": f"phi^{k}",
                "adams_degree": 2,
                "adams_measure": f"phi^(2^{k})=phi^{2**k}",
                "multiplicative_degree": 2**k,
                "multiplicative_measure": f"phi^S_{k}=phi^{s_closed}",
                "measures_verified_on_object": verified,
            }
        )
    p_rates = write_csv(
        "self_action_rates.csv",
        "self_action_rates.py",
        ["k", "S_k_closed_form", "S_k_tensor_spectrum", "S_k_mad_sum",
         "S_k_three_formulas_agree", "additive_degree", "additive_measure",
         "adams_degree", "adams_measure", "multiplicative_degree",
         "multiplicative_measure", "measures_verified_on_object"],
        rows,
    )

    # asymptotic ratio S_k / (2^k sqrt(k/2pi)) -> 1
    mp.mp.dps = 50
    arows = []
    for k in (5, 10, 50, 100, 200, 500, 1000):
        approx = mp.mpf(2) ** k * mp.sqrt(mp.mpf(k) / (2 * mp.pi))
        ratio = mp.mpf(Sk_closed(k)) / approx
        arows.append(
            {
                "k": k,
                "S_k": str(Sk_closed(k)),
                "asymptotic_2^k_sqrt(k/2pi)": mp.nstr(approx, 12),
                "ratio": mp.nstr(ratio, 18),
                "abs_ratio_minus_1": mp.nstr(abs(ratio - 1), 6),
            }
        )
    p_asym = write_csv(
        "sk_asymptotic.csv",
        "self_action_rates.py",
        ["k", "S_k", "asymptotic_2^k_sqrt(k/2pi)", "ratio", "abs_ratio_minus_1"],
        arows,
    )

    seq = [r["S_k_closed_form"] for r in rows]
    print(f"wrote {p_rates} (S_k = {seq})")
    print(f"wrote {p_asym} (ratio@k=1000 = {arows[-1]['ratio']})")


if __name__ == "__main__":
    main()
