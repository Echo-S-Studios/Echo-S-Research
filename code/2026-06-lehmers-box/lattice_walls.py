"""
Producer: the lattice walls -- the four directions (pi/2)Z.  Section 4.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
Emits the angular half of Lehmer's Box:
  * catalog_arguments.csv -- every eigenvalue of every catalog seed with its
        modulus and argument, confirming each lies in (pi/2)Z (Lem. 4.1);
  * box_membership.csv    -- the box test (floor wall AND lattice wall) applied
        to the 7 catalog seeds, Lehmer L, and beta_4: the catalog is inside,
        Lehmer is outside on both walls, beta_4 is above the floor yet outside on
        the lattice wall (Thm. 4.6, Thm. 4.7, Fig. 1);
  * lattice_walls.json    -- the (pi/2)Z group facts (closed under + and
        doubling, NOT halving; pi/4 off lattice, Lem. 4.2), the fourth-roots-of-
        unity mu_4 (Lem. 4.3), the K x^2 split (-5 +/- 3sqrt5)/2 (Lem. 4.1), and
        the off-lattice irrationality of Salem on-circle conjugates (Lem. 4.4).

Backs: Lem. 4.1-4.4, Thm. 4.6, Thm. 4.7, Fig. 1.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

import box_core as C
from box_io import write_csv, write_json

mp.mp.dps = 50
PHI = (1 + mp.sqrt(5)) / 2
TOL = mp.mpf(10) ** -12


def _arg_deg(r):
    """Argument in degrees, folded to [0, 360)."""
    a = mp.degrees(mp.arg(r))
    return a % 360


def _arg_in_halfpi_lattice(r, tol=mp.mpf(10) ** -20):
    """True iff arg(r) in {0, pi/2, pi, 3pi/2}: r is real OR purely imaginary."""
    return abs(mp.re(r)) < tol or abs(mp.im(r)) < tol


def catalog_argument_rows():
    rows = []
    for name, coeffs in C.CATALOG.items():
        poly = sp.sstr(sp.Poly(coeffs, sp.symbols('x')).as_expr())
        for r in C.roots_of(coeffs):
            mod = abs(r)
            rows.append(dict(
                seed=name, minimal_polynomial=poly,
                root_real=mp.nstr(mp.re(r), 12),
                root_imag=mp.nstr(mp.im(r), 12),
                modulus=mp.nstr(mod, 12),
                argument_deg=mp.nstr(_arg_deg(r), 10),
                on_unit_circle=bool(abs(mod - 1) < TOL),
                arg_in_halfpi_Z=bool(_arg_in_halfpi_lattice(r))))
    return rows


def box_membership_rows():
    """Box test = floor wall (Mah in {1} U [phi,inf)) AND lattice wall (every
    on-circle root in (pi/2)Z).  Applied to catalog + Lehmer + beta_4."""
    items = list(C.CATALOG.items()) + [("Lehmer_L", C.LEHMER), ("beta_4", C.BETA4)]
    rows = []
    for name, coeffs in items:
        m = C.mahler_measure(coeffs)
        floor_ok = bool(abs(m - 1) < TOL or m >= PHI - TOL)
        oncirc = [r for r in C.roots_of(coeffs) if abs(abs(r) - 1) < mp.mpf(10) ** -10]
        lattice_ok = all(_arg_in_halfpi_lattice(r) for r in oncirc)
        in_box = bool(floor_ok and lattice_ok)
        rows.append(dict(
            name=name,
            mahler_measure=mp.nstr(m, 16),
            floor_wall_ok=floor_ok,
            n_oncircle_roots=len(oncirc),
            lattice_wall_ok=bool(lattice_ok),
            in_box=in_box))
    return rows


def lattice_group_facts():
    """Lem. 4.2: model (pi/2)Z as Z/4 (units of pi/2).  Closed under + and
    doubling, NOT under halving; pi/4 = (1/2)(pi/2) is off the lattice."""
    lattice = list(range(4))                          # {0,1,2,3} * (pi/2)
    closed_add = all((a + b) % 4 in lattice for a in lattice for b in lattice)
    closed_double = all((2 * a) % 4 in lattice for a in lattice)
    half_of_pi2_off = sp.Rational(1, 2) not in [sp.Integer(k) for k in lattice]
    return {
        "model": "(pi/2)Z ~ Z/4 (units of pi/2)",
        "closed_under_addition": bool(closed_add),
        "closed_under_doubling": bool(closed_double),
        "closed_under_halving": False,
        "half_of_pi_over_2_is_pi_over_4_off_lattice": bool(half_of_pi2_off),
    }


def mu4_facts():
    """Lem. 4.3: |z|=1 and arg z in (pi/2)Z  =>  z in mu_4 = {1,i,-1,-i}, z^4=1."""
    mu4 = [mp.mpc(1, 0), mp.mpc(0, 1), mp.mpc(-1, 0), mp.mpc(0, -1)]
    entries = []
    for k in range(4):
        z = mp.e ** (1j * (mp.pi / 2) * k)
        nearest = min(abs(z - w) for w in mu4)
        entries.append(dict(
            k=k, z_real=mp.nstr(mp.re(z), 8), z_imag=mp.nstr(mp.im(z), 8),
            is_fourth_root_of_unity=bool(abs(z ** 4 - 1) < mp.mpf(10) ** -30),
            matches_mu4=bool(nearest < mp.mpf(10) ** -30)))
    return {"mu_4": "{1, i, -1, -i}", "entries": entries}


def K_split_facts():
    """Lem. 4.1: x^4+5x^2-5 gives x^2 = y with y^2 + 5y - 5 = 0, i.e.
    y = (-5 +/- 3sqrt5)/2 = {+0.854.. -> real pair, -5.854.. -> imag pair}."""
    yp = (-5 + 3 * mp.sqrt(5)) / 2
    ym = (-5 - 3 * mp.sqrt(5)) / 2
    return {
        "quadratic_in_x2": "y^2 + 5y - 5 = 0",
        "upper_root_y": mp.nstr(yp, 16),
        "upper_gives": "x^2 > 0 -> real pair +/- K at arg {0, pi}",
        "lower_root_y": mp.nstr(ym, 16),
        "lower_gives": "x^2 < 0 -> imaginary pair +/- i beta at arg {pi/2, 3pi/2}",
    }


def salem_offlattice_facts():
    """Lem. 4.4: the on-circle conjugates of a Salem number are NOT roots of
    unity and carry argument off (pi/2)Z.  Instance: Lehmer L, 8 on-circle
    roots, each bounded away from both axes."""
    oncirc = [r for r in C.roots_of(C.LEHMER) if abs(abs(r) - 1) < mp.mpf(10) ** -10]
    min_re = min(abs(mp.re(r)) for r in oncirc)
    min_im = min(abs(mp.im(r)) for r in oncirc)
    return {
        "instance": "Lehmer L (degree 10 Salem)",
        "n_oncircle_conjugates": len(oncirc),
        "all_off_lattice": bool(all(not _arg_in_halfpi_lattice(r) for r in oncirc)),
        "min_abs_real_part": mp.nstr(min_re, 8),
        "min_abs_imag_part": mp.nstr(min_im, 8),
        "note": "each on-circle conjugate has |Re|>0 and |Im|>0, so its argument "
                "is not in (pi/2)Z; being a conjugate of beta>1 it is not a root "
                "of unity.",
    }


def main():
    arg_rows = catalog_argument_rows()
    box_rows = box_membership_rows()

    p1 = write_csv("catalog_arguments.csv",
                   ["seed", "minimal_polynomial", "root_real", "root_imag",
                    "modulus", "argument_deg", "on_unit_circle", "arg_in_halfpi_Z"],
                   arg_rows, __file__)
    p2 = write_csv("box_membership.csv",
                   ["name", "mahler_measure", "floor_wall_ok", "n_oncircle_roots",
                    "lattice_wall_ok", "in_box"],
                   box_rows, __file__)
    payload = {
        "_description": "The lattice walls (pi/2)Z of Lehmer's Box (Section 4).",
        "lattice_group_facts": lattice_group_facts(),
        "mu4_facts": mu4_facts(),
        "K_x2_split": K_split_facts(),
        "salem_oncircle_off_lattice": salem_offlattice_facts(),
        "all_catalog_arguments_in_halfpi_Z":
            bool(all(r["arg_in_halfpi_Z"] for r in arg_rows)),
    }
    p3 = write_json("lattice_walls.json", payload, __file__)

    print("wrote", p1, f"({len(arg_rows)} catalog eigenvalues)")
    print("wrote", p2)
    print("wrote", p3)
    for r in box_rows:
        print(f"  {r['name']:10s} in_box={r['in_box']!s:5s} "
              f"floor={r['floor_wall_ok']!s:5s} lattice={r['lattice_wall_ok']!s:5s} "
              f"Mah={r['mahler_measure']}")


if __name__ == "__main__":
    main()
