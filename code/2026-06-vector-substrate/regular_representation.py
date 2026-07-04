"""
PRODUCER: regular representation & spectral picture.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex, Sec. 2.1-2.2
  - Prop. 2.2 (regular representation), Ex. 2.3 (rho of phi, sqrt5, 3)
  - Prop. 2.6 / Ex. 2.7 (M diagonalises rho: spec rho(x) = conjugates)
  - Ex. 2.8 (quartic spectrum = 4 conjugates)
  - Rem. 2.5 (Mahler measure as a spectral invariant of the companion matrix)

Emits:
  data/regular_representation.json  -- companion/regular-rep matrices, field
                                        trace/norm, charpoly vs minpoly, and the
                                        Minkowski diagonalisation M rho M^{-1}.
  data/mahler_spectrum.csv          -- eigenvalue magnitudes and the spectral
                                        Mahler measure prod_{|lambda|>1}|lambda|.
"""
import mpmath as mp
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "regular_representation.py"


def regrep_examples():
    """Ex. 2.3 + Prop. 2.2: rho(phi), rho(sqrt5), rho(3), and the n/d
    multiplicity made visible via charpoly vs minpoly."""
    C = vc.companion_from_poly(x**2 - x - 1)          # golden field generator phi
    out = {}

    # rho(phi) = companion(x^2-x-1)
    out["rho_phi"] = {
        "coord": [0, 1],
        "matrix": vc.mat_to_list(C),
        "charpoly": vc.sval(sp.expand(C.charpoly(x).as_expr())),
        "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(C))),
        "field_trace": vc.sval(C.trace()),           # Tr(phi) = 1
        "field_norm": vc.sval(C.det()),              # N(phi) = -1
    }

    # rho(sqrt5) = 2 rho(phi) - I ; sqrt5 = 2 phi - 1 -> coords (-1,2)
    Rs = vc.rho((-1, 2), C)
    out["rho_sqrt5"] = {
        "coord": [-1, 2],
        "matrix": vc.mat_to_list(Rs),
        "charpoly": vc.sval(sp.expand(Rs.charpoly(x).as_expr())),
        "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(Rs))),
        "field_trace": vc.sval(vc.field_trace((-1, 2), C)),   # 0
        "field_norm": vc.sval(vc.field_norm((-1, 2), C)),     # -5
    }

    # rho(3) = 3I is derogatory: charpoly (x-3)^2 but minpoly x-3
    R3 = vc.rho((3, 0), C)
    out["rho_3_derogatory"] = {
        "coord": [3, 0],
        "matrix": vc.mat_to_list(R3),
        "charpoly": vc.sval(sp.expand(R3.charpoly(x).as_expr())),
        "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(R3))),
        "multiplicity_n_over_d": 2,
    }

    # Prop. 2.2 Eq.(2): x=sqrt2 inside Q(sqrt2+sqrt3) (n=4,d=2):
    # charpoly=(x^2-2)^2, minpoly=x^2-2.  sqrt2 = (theta^3-9theta)/2.
    C4 = vc.companion_from_poly(x**4 - 10 * x**2 + 1)
    s2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
    Rs2 = vc.rho(s2, C4)
    out["multiplicity_sqrt2_in_deg4"] = {
        "coord": [vc.sval(v) for v in s2],
        "charpoly": vc.sval(sp.expand(Rs2.charpoly(x).as_expr())),
        "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(Rs2))),
        "n": 4,
        "d": 2,
    }
    return out


def spectral_diagonalisation():
    """Prop. 2.6 / Ex. 2.7: the embedding matrix M diagonalises the whole regular
    representation, M rho(x) M^{-1} = diag(sigma_k(x)); shown for phi and sqrt5
    in Q(sqrt5) with the same M."""
    phi = (1 + sp.sqrt(5)) / 2
    phic = (1 - sp.sqrt(5)) / 2
    C = vc.companion_from_poly(x**2 - x - 1)
    M = sp.Matrix([[1, phi], [1, phic]])
    D_phi = sp.simplify(M * C * M.inv())
    D_s5 = sp.simplify(M * vc.rho((-1, 2), C) * M.inv())
    return {
        "field": "Q(sqrt5)",
        "M_rows_are_embeddings": vc.mat_to_list(M),
        "diag_rho_phi": vc.mat_to_list(D_phi),      # diag(phi, phi')
        "diag_rho_sqrt5": vc.mat_to_list(D_s5),     # diag(sqrt5, -sqrt5)
        "charpoly_rho_phi": vc.sval(sp.expand(C.charpoly(x).as_expr())),
    }


def quartic_spectrum():
    """Ex. 2.8: rho(sqrt2+sqrt3)=C(x^4-10x^2+1) has eigenvalues the four
    conjugates {+-sqrt2+-sqrt3}, two outside the unit circle."""
    mp.mp.dps = 40
    roots = mp.polyroots([1, 0, -10, 0, 1], extraprec=200)
    mags = sorted(float(abs(r)) for r in roots)
    return {
        "poly": "x^4 - 10*x^2 + 1",
        "conjugates_closed_form": ["sqrt2+sqrt3", "sqrt2-sqrt3", "-sqrt2+sqrt3", "-sqrt2-sqrt3"],
        "eigenvalue_magnitudes": mags,
        "num_outside_unit_circle": int(sum(1 for r in roots if abs(r) > 1)),
        "mahler_measure_closed_form": "5 + 2*sqrt(6)",
        "mahler_measure_float": float(mp.mpf(5) + 2 * mp.sqrt(6)),
    }


def mahler_spectrum_rows():
    """Rem. 2.5: the Mahler measure read off the companion spectrum,
    Mah = prod_{|lambda|>1} |lambda|, for a spread of seeds."""
    fields = [
        ("phi", x**2 - x - 1, "phi"),
        ("plastic (Smyth mu_S)", x**3 - x - 1, "mu_S = 1.324717957..."),
        ("sqrt2+sqrt3", x**4 - 10 * x**2 + 1, "5 + 2*sqrt(6)"),
        ("2sqrt6", x**2 - 24, "24"),
        ("sqrt7", x**2 - 7, "7"),
    ]
    rows = []
    for name, poly, closed in fields:
        mp.mp.dps = 45
        p = sp.Poly(poly, x)
        coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
        roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
        n_out = sum(1 for r in roots if abs(r) > 1)
        mah = vc.mahler_measure_mp(poly, dps=45)
        rows.append([
            name,
            sp.sstr(poly),
            p.degree(),
            n_out,
            mp.nstr(mah, 15),
            closed,
        ])
    return rows


def main():
    payload = {
        "regular_representation": regrep_examples(),
        "spectral_diagonalisation": spectral_diagonalisation(),
        "quartic_spectrum": quartic_spectrum(),
    }
    p1 = vc.write_json("regular_representation.json", payload, SCRIPT)

    rows = mahler_spectrum_rows()
    p2 = vc.write_csv(
        "mahler_spectrum.csv",
        ["seed", "min_poly", "degree", "num_eig_outside_unit_circle",
         "mahler_measure", "closed_form"],
        rows,
        SCRIPT,
    )
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print("  rho(phi) trace/norm =",
          payload["regular_representation"]["rho_phi"]["field_trace"],
          payload["regular_representation"]["rho_phi"]["field_norm"])
    print("  quartic: #eig outside unit circle =",
          payload["quartic_spectrum"]["num_outside_unit_circle"],
          "Mahler = 5+2sqrt6")


if __name__ == "__main__":
    main()
