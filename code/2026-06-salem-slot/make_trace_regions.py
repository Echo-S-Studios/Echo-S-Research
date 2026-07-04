"""
Producer: the downstairs trace structure of Section 2 (+ the lift/interval
structure of Lemma 4.4 and Section 'entry').

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces (Section 2 "There is no fourth channel ...", Lemma 4.4):
  * data/trace_regions.json
      - flip D(t)=t^2-4 as the discriminant of x^2 - t x + 1        (Def 2.1)
      - the trace-down of Lehmer, and its round-trip lift back to Lehmer
        (Def 2.1: P(x)=x^m T(x+1/x))
      - the three regions cut by the flip t=+-2, with samples        (Lem 2.2)
      - "no fourth channel": (-2,2) and (2,oo) are disjoint          (Cor 2.3)
      - Lehmer trace-down totally-real split (1 root >2, m-1 inside)  (Cor 2.3)
      - lift on the captured lattice {-2,0,2} = roots of unity        (Sec entry)
      - interval endpoints: t=2 |-> x=1 ; t=sqrt5 |-> {phi,1/phi}     (Lem 4.4)

Run:  py code/2026-06-salem-slot/make_trace_regions.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

x, t = sc.x, sc.t

# Lehmer's polynomial (degree 10) and its stated trace-down (degree 5).
LEHMER = x ** 10 + x ** 9 - x ** 7 - x ** 6 - x ** 5 - x ** 4 - x ** 3 + x + 1


def three_region_samples():
    """Lem 2.2: sample theta in each of the three pre-images and report t."""
    mp = sc.mp
    mp.mp.dps = 40
    rows = []
    # on-circle: theta = e^{i psi} -> t = 2 cos psi in [-2,2]
    for k in (3, 7, 12):
        psi = mp.pi * k / 21
        theta = mp.e ** (1j * psi)
        tv = theta + 1 / theta
        rows.append({
            "region": "captured (|t|<2)", "theta": f"e^(i*{k}pi/21)",
            "t_real": sc.mp.nstr(tv.real, 12), "t_imag": sc.mp.nstr(tv.imag, 3),
            "two_cos_psi": sc.mp.nstr(2 * mp.cos(psi), 12),
        })
    # grow: theta > 1 -> t > 2
    for theta in (mp.mpf('1.0001'), mp.mpf('1.5'), mp.mpf('10')):
        rows.append({
            "region": "grow (t>2)", "theta": sc.mp.nstr(theta, 8),
            "t_real": sc.mp.nstr(theta + 1 / theta, 12), "t_imag": "0",
            "two_cos_psi": "-",
        })
    # decay: theta < -1 -> t < -2
    for theta in (mp.mpf('-1.0001'), mp.mpf('-3')):
        rows.append({
            "region": "decay (t<-2)", "theta": sc.mp.nstr(theta, 8),
            "t_real": sc.mp.nstr(theta + 1 / theta, 12), "t_imag": "0",
            "two_cos_psi": "-",
        })
    return rows


def lehmer_split():
    """Cor 2.3: Lehmer trace-down is totally real, 1 root >2, m-1 in (-2,2)."""
    T = sc.trace_down(LEHMER).as_expr()
    roots = [complex(r) for r in sp.nroots(sp.Poly(T, t), n=40)]
    reals = sorted(r.real for r in roots)
    return {
        "trace_down_T": str(sp.Poly(T, t).as_expr()),
        "all_real": all(abs(r.imag) < 1e-25 for r in roots),
        "roots": [round(r, 10) for r in reals],
        "count_past_2": sum(1 for r in reals if r > 2),
        "count_inside_pm2": sum(1 for r in reals if -2 < r < 2),
        "tau0_dominant": round(max(reals), 10),
    }


def lattice_lift():
    """Sec entry: L on the captured lattice {-2,0,2} gives roots of unity."""
    out = {}
    for label, T_of_t, expect in (
        ("t-2", t - 2, (x - 1) ** 2),
        ("t (=0)", t, x ** 2 + 1),
        ("t+2", t + 2, (x + 1) ** 2),
    ):
        L = sc.lift(T_of_t).as_expr()
        out[label] = {"lift": str(sp.expand(L)), "equals": str(sp.expand(expect)),
                      "match": sp.simplify(L - sp.expand(expect)) == 0}
    return out


def interval_endpoints():
    """Lem 4.4: t=2 -> x=1 ; t=sqrt5 -> {phi, 1/phi}."""
    # t = 2 endpoint
    p2 = sp.expand((x + 1 / x - 2) * x)  # x^2 - 2x + 1
    # t = sqrt5 endpoint
    P = x ** 2 - sp.sqrt(5) * x + 1
    roots = sp.solve(P, x)
    prod = sp.simplify(roots[0] * roots[1])
    ssum = sp.simplify(roots[0] + roots[1])
    phi_is_root = any(sp.simplify(r - sc.phi_sym) == 0 for r in roots)
    return {
        "t=2": {"lift": str(sp.expand(p2)), "roots": [str(r) for r in sp.solve(sp.Eq(p2, 0), x)]},
        "t=sqrt5": {
            "lift": "x^2 - sqrt5*x + 1",
            "root_product": str(prod), "root_sum": str(ssum),
            "golden_pair": phi_is_root and prod == 1 and sp.simplify(ssum - sp.sqrt(5)) == 0,
        },
    }


def main():
    T = sc.trace_down(LEHMER)
    roundtrip = sc.lift(T).as_expr()
    payload = {
        "_description": "Downstairs trace structure (Section 2) and lift/interval structure (Lemma 4.4).",
        "flip_discriminant": {
            "claim": "D(t) = discriminant_x(x^2 - t x + 1) = t^2 - 4",
            "computed": str(sc.flip_discriminant()),
            "match": sp.simplify(sc.flip_discriminant() - (t ** 2 - 4)) == 0,
        },
        "trace_down_lehmer": {
            "P_lehmer": str(sp.expand(LEHMER)),
            "trace_down_T": str(T.as_expr()),
            "roundtrip_lift_equals_lehmer": sp.simplify(roundtrip - LEHMER) == 0,
        },
        "three_regions": three_region_samples(),
        "no_fourth_channel": {
            "claim": "(-2,2) intersect (2,oo) is empty; 'on-circle and expanding' is a contradiction",
            "intersection": str(sp.Interval.open(-2, 2).intersect(sp.Interval.open(2, sp.oo))),
            "disjoint": sp.Interval.open(-2, 2).intersect(sp.Interval.open(2, sp.oo)) == sp.EmptySet,
        },
        "lehmer_totally_real_split": lehmer_split(),
        "lattice_lift_roots_of_unity": lattice_lift(),
        "interval_endpoints": interval_endpoints(),
    }
    path = io.write_json("trace_regions.json", payload, __file__)
    print("wrote", path)


if __name__ == "__main__":
    main()
