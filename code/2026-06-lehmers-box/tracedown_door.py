"""
Producer: the trace-down face and the one door.  Sections 6-7.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Emits:
  * tracedown_flip.csv -- Def. 6.1 / Lem. 6.2: the trace-down T of Lehmer and
        beta_4 and its flip-straddle pattern (one root > 2, the rest in (-2,2)),
        the Salem criterion; plus the cyclotomic lattice rho(mu_4) = {2,0,-2}
        (Prop. 6.3);
  * guard_ladder.csv   -- Prop. 7.5: the closure guard's three verdicts (FORCED /
        FORCED_ABOVE_FLOOR / INVALID_CLOSURE) on the catalog, beta_4, and Lehmer,
        with the exact Q(sqrt5) test element m_beta(phi) = a + b sqrt5 and its
        sign;
  * door_summary.json  -- Lem. 7.4 self-action spectrum {0, +/-sqrt5}; Prop. 7.1
        circulant emits no Salem; Lem. 7.2 (Shoda) commutators are traceless and
        a traceless carrier of Lehmer exists; Lem. 6.2 flip discriminant t^2-4.

Backs: Def. 6.1, Lem. 6.2, Prop. 6.3, Prop. 7.1, Lem. 7.2, Lem. 7.4, Prop. 7.5.
"""

from __future__ import annotations

import random

import mpmath as mp
import numpy as np
import sympy as sp

import box_core as C
from box_io import write_csv, write_json

mp.mp.dps = 50
_x = sp.symbols('x')


def tracedown_rows():
    """Def. 6.1 / Lem. 6.2: trace-down and flip-straddle pattern for the two
    Salem instances."""
    rows = []
    for name, coeffs in [("beta_4", C.BETA4), ("Lehmer_L", C.LEHMER)]:
        m = (len(coeffs) - 1) // 2
        T = C.tracedown(coeffs)
        tot_real, above2, inside, _ = C.tracedown_root_pattern(coeffs)
        rows.append(dict(
            name=name,
            degree=len(coeffs) - 1,
            tracedown_T=sp.sstr(sp.Poly(T, sp.symbols('t')).as_expr()),
            totally_real=bool(tot_real),
            roots_above_2=above2,
            roots_in_open_minus2_2=inside,
            expected_inside=m - 1,
            is_flip_straddle_salem=bool(tot_real and above2 == 1 and inside == m - 1)))
    return rows


def guard_rows():
    """Prop. 7.5 guard ladder over catalog + beta_4 + Lehmer."""
    items = list(C.CATALOG.items()) + [("beta_4", C.BETA4), ("Lehmer_L", C.LEHMER)]
    rows = []
    for name, coeffs in items:
        v = C.validate(coeffs)
        ab = v["m_beta_phi_ab"]
        rows.append(dict(
            name=name,
            minimal_polynomial=sp.sstr(sp.Poly(coeffs, _x).as_expr()),
            has_salem_factor=bool(v["salem_factor"] is not None),
            m_beta_phi=(f"{ab[0]} + ({ab[1]})*sqrt5" if ab else ""),
            sign_m_beta_phi=("" if v["sign"] is None else v["sign"]),
            verdict=v["verdict"]))
    return rows


def selfaction_facts():
    """Lem. 7.4: spec(ad_R) = difference set of eigenvalues of R; golden seed
    gives {0, +sqrt5, -sqrt5}."""
    ev = C.selfaction_spectrum(C.CATALOG["phi"])
    s5 = float(mp.sqrt(5))
    return {
        "seed": "R = companion(x^2 - x - 1) (golden)",
        "spectrum_sorted": [round(v, 12) for v in ev],
        "expected": [round(-s5, 12), 0.0, 0.0, round(s5, 12)],
        "difference_set": "{0, +sqrt5, -sqrt5}",
        "matches": bool(abs(ev[0] + s5) < 1e-9 and abs(ev[-1] - s5) < 1e-9
                        and sum(1 for v in ev if abs(v) < 1e-9) == 2),
    }


def rho_lattice_facts():
    """Prop. 6.3: rho(z) = z + 1/z sends mu_4 = {1,i,-1,-i} to {2,0,-2}."""
    rho = {str(z): sp.simplify(z + 1 / z)
           for z in (sp.Integer(1), sp.I, sp.Integer(-1), -sp.I)}
    return {
        "rho": "rho(z) = z + 1/z",
        "images": {k: str(v) for k, v in rho.items()},
        "image_set": sorted(str(v) for v in set(rho.values())),
        "cyclotomic_lattice": "{2, 0, -2}",
    }


def flip_discriminant_facts():
    """Lem. 6.2: for x^2 - t x + 1, D = t^2 - 4;  D<0 on circle, D>0 off, D=0 at
    t = +/-2 (the lambda = 2c flip)."""
    samples = []
    for tv in ("-3", "-2", "-1.7", "0", "1.3", "2", "2.3"):
        t = mp.mpf(tv)
        D = t ** 2 - 4
        samples.append(dict(
            t=tv, D=mp.nstr(D, 8),
            location=("off-circle" if D > 0 else
                      ("boundary" if D == 0 else "on-circle"))))
    return {"discriminant": "D = t^2 - 4", "samples": samples}


def circulant_facts(n_trials=8, seed=3):
    """Prop. 7.1: integer circulant eigenvalues are sum_k c_k omega^{jk}, and no
    integer circulant char-poly factor is a Salem minimal polynomial."""
    random.seed(seed)
    all_match = True
    any_salem = False
    trials = []
    for _ in range(n_trials):
        n = random.randint(3, 6)
        c = [random.randint(-3, 3) for _ in range(n)]
        circ = np.array([[c[(j - k) % n] for k in range(n)] for j in range(n)],
                        dtype=complex)
        ev_np = list(np.linalg.eigvals(circ))
        w = mp.e ** (2j * mp.pi / n)
        ev_f = [complex(sum(c[k] * w ** (j * k) for k in range(n)))
                for j in range(n)]
        remaining = ev_np[:]
        matched = True
        for b in ev_f:
            j = min(range(len(remaining)), key=lambda idx: abs(remaining[idx] - b))
            if abs(remaining[j] - b) >= 1e-9:
                matched = False
                break
            remaining.pop(j)
        all_match = all_match and matched
        cp = sp.Matrix([[c[(j - k) % n] for k in range(n)]
                        for j in range(n)]).charpoly(_x)
        salem_here = False
        for fac, _mult in sp.factor_list(cp.as_expr())[1]:
            fc = [int(cc) for cc in sp.Poly(fac, _x).all_coeffs()]
            if fc and fc[0] < 0:
                fc = [-v for v in fc]
            if C.is_salem(fc):
                salem_here = True
        any_salem = any_salem or salem_here
        trials.append(dict(n=n, circulant_row=c, salem_factor=salem_here))
    return {
        "eigenvalue_formula": "lambda_j = sum_k c_k omega^{jk}, omega = e^{2 pi i/n}",
        "n_trials": n_trials,
        "all_eigenvalue_formulas_match": bool(all_match),
        "any_salem_factor_found": bool(any_salem),
        "trials": trials,
    }


def shoda_facts():
    """Lem. 7.2 (Shoda), necessary direction: every commutator is traceless; a
    traceless integer carrier of Lehmer exists (the one door)."""
    random.seed(5)
    np.random.seed(5)  # deterministic: drift CI regenerates this exactly
    max_tr = 0.0
    for _ in range(30):
        n = random.randint(2, 5)
        X = np.random.randint(-3, 4, (n, n))
        Y = np.random.randint(-3, 4, (n, n))
        max_tr = max(max_tr, abs(float(np.trace(X @ Y - Y @ X))))
    # traceless carrier: companion(L) (trace -1) direct-sum [1]
    companion = sp.zeros(10, 10)
    for i in range(9):
        companion[i + 1, i] = 1
    for i, coeff in enumerate(C.LEHMER[1:][::-1]):
        companion[i, 9] = -coeff
    carrier = sp.zeros(11, 11)
    carrier[:10, :10] = companion
    carrier[10, 10] = 1
    cp = carrier.charpoly(_x).as_expr()
    _, rem = sp.div(cp, sp.Poly(C.LEHMER, _x).as_expr(), _x)
    return {
        "necessary_direction": "tr[X,Y] = 0 for all commutators",
        "max_abs_trace_over_random_commutators": mp.nstr(mp.mpf(max_tr), 4),
        "traceless_lehmer_carrier": "companion(L) (+) [1], size 11, trace 0",
        "carrier_trace": int(carrier.trace()),
        "lehmer_divides_carrier_charpoly": bool(sp.simplify(rem) == 0),
        "note": "Shoda's theorem makes any traceless matrix a commutator, so a "
                "FREE commutator is the one door out of the box.",
    }


def main():
    td_rows = tracedown_rows()
    gd_rows = guard_rows()

    p1 = write_csv("tracedown_flip.csv",
                   ["name", "degree", "tracedown_T", "totally_real",
                    "roots_above_2", "roots_in_open_minus2_2", "expected_inside",
                    "is_flip_straddle_salem"], td_rows, __file__)
    p2 = write_csv("guard_ladder.csv",
                   ["name", "minimal_polynomial", "has_salem_factor",
                    "m_beta_phi", "sign_m_beta_phi", "verdict"], gd_rows, __file__)
    payload = {
        "_description": "The trace-down face (Section 6) and the one door "
                        "(Section 7) of Lehmer's Box.",
        "rho_cyclotomic_lattice": rho_lattice_facts(),
        "flip_discriminant": flip_discriminant_facts(),
        "selfaction_difference_spectrum": selfaction_facts(),
        "circulant_no_salem": circulant_facts(),
        "shoda_commutator_door": shoda_facts(),
    }
    p3 = write_json("door_summary.json", payload, __file__)

    print("wrote", p1)
    print("wrote", p2)
    print("wrote", p3)
    for r in gd_rows:
        print(f"  guard {r['name']:10s} -> {r['verdict']:20s} "
              f"m_beta(phi)= {r['m_beta_phi']}")


if __name__ == "__main__":
    main()
