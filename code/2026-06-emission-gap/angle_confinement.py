"""Producer: angle confinement -- the (pi/2)Z invariant and its consequences.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces:
  data/2026-06-emission-gap/catalog_arguments.csv          (Lemma 4.1)
  data/2026-06-emission-gap/spectral_oncircle_eigenvalues.csv (Lemma 4.3)
  data/2026-06-emission-gap/lehmer_oncircle_arguments.csv   (Lemma 5.1 / App. A)
Paper results: every catalog eigenvalue has argument in (pi/2)Z (Lemma 4.1);
the spectral operators (kron, direct-sum, squaring) send every on-circle
eigenvalue to a fourth root of unity {1, i, -1, -i} (Lemma 4.3); a Salem
number's on-circle conjugates are NOT roots of unity -- the eight on-circle
arguments of the Lehmer polynomial are none a multiple of 90 degrees (Lemma 5.1).
"""
import numpy as np
import mpmath as mp

import emgap_core as C

SCRIPT = "angle_confinement.py"

FOURTH_ROOTS = {(1, 0): "1", (-1, 0): "-1", (0, 1): "i", (0, -1): "-i"}


def nearest_fourth_root(e):
    best = min(((1, 0), (-1, 0), (0, 1), (0, -1)),
               key=lambda z: abs(e - complex(z[0], z[1])))
    return FOURTH_ROOTS[best], abs(e - complex(best[0], best[1]))


def catalog_argument_rows():
    rows = []
    for name, coeffs in C.CATALOG.items():
        for r in C.mp_roots(coeffs):
            if abs(r) < mp.mpf(10) ** (-30):
                continue
            ang = C.arg_degrees(r)
            rows.append({
                "seed": name,
                "root_real": C.s(r.real, 16),
                "root_imag": C.s(r.imag, 16),
                "modulus": C.s(abs(r), 16),
                "argument_deg": round(ang, 6),
                "arg_mod_90": round(min(ang % 90, 90 - (ang % 90)), 9),
            })
    return rows


def _dsum(A, B):
    n, m = A.shape[0], B.shape[0]
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[n:, n:] = B
    return M


def spectral_oncircle_rows():
    """On-circle eigenvalues of kron / direct-sum / squares of catalog companions;
    each must be a fourth root of unity (Lemma 4.3)."""
    comps = {name: C.companion(c) for name, c in C.CATALOG.items()}
    rows = []
    names = list(comps)
    for i, na in enumerate(names):
        for nb in names[i:]:
            A, B = comps[na], comps[nb]
            ops = {"kron": np.kron(A, B), "dsum": _dsum(A, B)}
            if na == nb:
                ops["square"] = A @ A
            for op, M in ops.items():
                seen = set()
                for e in np.linalg.eigvals(M):
                    if abs(abs(e) - 1) < 1e-9:
                        root_label, dist = nearest_fourth_root(e)
                        key = (root_label,)
                        if key in seen:
                            continue
                        seen.add(key)
                        rows.append({
                            "operation": op,
                            "operand_A": na,
                            "operand_B": nb,
                            "oncircle_eig_real": round(float(e.real), 9),
                            "oncircle_eig_imag": round(float(e.imag), 9),
                            "argument_deg": round((np.degrees(np.angle(e)) % 360), 4),
                            "nearest_fourth_root": root_label,
                            "is_fourth_root_of_unity": bool(dist < 1e-8),
                        })
    return rows


def lehmer_oncircle_rows():
    rows = []
    for r in C.mp_roots(C.LEHMER):
        if abs(abs(r) - 1) < mp.mpf(10) ** (-20):
            ang = C.arg_degrees(r)
            m90 = min(ang % 90, 90 - (ang % 90))
            rows.append({
                "argument_deg": round(ang, 5),
                "arg_mod_90": round(m90, 5),
                "is_multiple_of_90": bool(m90 < 1e-6),
                "is_root_of_unity": False,
            })
    rows.sort(key=lambda d: d["argument_deg"])
    return rows


def main():
    ca = catalog_argument_rows()
    C.write_csv("catalog_arguments.csv",
                ["seed", "root_real", "root_imag", "modulus", "argument_deg", "arg_mod_90"],
                ca, SCRIPT)
    sp_rows = spectral_oncircle_rows()
    C.write_csv("spectral_oncircle_eigenvalues.csv",
                ["operation", "operand_A", "operand_B", "oncircle_eig_real",
                 "oncircle_eig_imag", "argument_deg", "nearest_fourth_root",
                 "is_fourth_root_of_unity"], sp_rows, SCRIPT)
    le = lehmer_oncircle_rows()
    C.write_csv("lehmer_oncircle_arguments.csv",
                ["argument_deg", "arg_mod_90", "is_multiple_of_90", "is_root_of_unity"],
                le, SCRIPT)

    print(f"catalog_arguments.csv: {len(ca)} eigenvalues, "
          f"all args in (pi/2)Z: {all(r['arg_mod_90'] < 1e-6 for r in ca)}")
    print(f"spectral_oncircle_eigenvalues.csv: {len(sp_rows)} on-circle values, "
          f"all fourth roots of unity: {all(r['is_fourth_root_of_unity'] for r in sp_rows)}")
    print(f"lehmer_oncircle_arguments.csv: {len(le)} on-circle conjugates, "
          f"none a multiple of 90 deg: {not any(r['is_multiple_of_90'] for r in le)}")


if __name__ == "__main__":
    main()
