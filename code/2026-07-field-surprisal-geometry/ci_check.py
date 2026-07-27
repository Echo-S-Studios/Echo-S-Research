#!/usr/bin/env python3
"""ci_check.py -- the per-push CI fast lane for the Field Surprisal Geometry harnesses.

Runs the load-bearing, quick, exact-arithmetic harnesses (each fail-first: it exits
non-zero on its first failed assertion) that establish the paper's forced results:
the five foundational layers (Z, Fisher = Var(log M), dual flatness, C = beta^2 I,
the embedding obstruction), the k=2 ruled classification together with its conceptual
master window identity  P = Z^2 * sum_{|s|=4} q4(s) w_s  (Sylvester + Cauchy-Binet),
the temperature dichotomy, k-necessity (2<=k<=5), and the compositum coupling.
Exits 0 iff every one passes.

The exhaustive / expensive lanes are shipped in harnesses/ and reproduced IN FULL
offline via run_all.sh (422 checks across all 20 harnesses, ~15 min); they are held
out of the per-push CI lane only for time, and their results are covered here by the
proved master identity and the rank-<=1 window criterion:
  * field_surprisal_classification.py  (~5000-sample non-indicator exclusion, now
                                        superseded by the master-identity proof)
  * t1_core.py                         (explicit ruled-system enumeration, redundant
                                        with the conceptual proof in t1_windowproof)
  * t3_suspension.py                   (iterated indicator joins, k=3,4,5)
  * t4_kwindows.py / t4b_census_fast.py / t5_catalog_census.py / t9_landscape.py
                                        (the higher-k / catalog / landscape censuses)
See MANIFEST.md for the full harness -> checks -> paper-section map.

    python3 ci_check.py     # ~2.5 min; exits 0 iff every fast-lane harness passes
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HARN = Path(__file__).resolve().parent / "harnesses"

# Dependency order matters once: t1_engine writes P_coeffs.pkl into the run cwd, and
# t1_windowproof reads it for its full 44/44 (43/43 standalone if absent).
FAST_LANE = [
    "field_surprisal_v2.py",       # 12/12  five forced layers: Z, Fisher=Var, dual-flat, C=b^2 I
    "suspension_theorem.py",       #  9/9   indicator-suspension: constant curvature 1/4
    "field_surprisal_tier2.py",    # 15/15  multi-statistic landscape, product geometry, charge grading
    "t1_reduction.py",             # 17/17  Z^2 Dq = P(w); four-point rank-one; collision lattice
    "t1_engine.py",                #  4/4   462-monomial dictionary of P (writes P_coeffs.pkl)
    "t1_branches.py",              # 16/16  branch survival across all 8 families
    "t1_windowproof.py",           # 44/44  the MASTER window identity (conceptual proof)
    "t3c_partC_exact.py",          # 16/16  geodesic-threshold exact certificates
    "t2_temperature.py",           # 12/12  forced Gibbs form; C = beta^2 I anchors
    "t6_selection.py",             # 18/18  the temperature dichotomy
    "t7_knecessity.py",            # 28/28  k-necessity 2<=k<=5 via the squarefree collision
    "t8_compositum.py",            # 13/13  compositum coupling (log 2, -6 log phi, rank 5)
    "t10_coupled.py",              #  9/9   the coupled family curves (Gauss obstruction != 0)
]


def main() -> int:
    # Force UTF-8 in the children so the harnesses' Unicode output (the integral sign,
    # phi, the blackboard letters, ...) never crashes on a legacy-codepage stdout
    # (e.g. Windows cp1252). CI runners are already UTF-8; this only helps local runs.
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    failures = []
    # A shared scratch cwd keeps the generated P_coeffs.pkl / P_dict.pkl out of the repo
    # tree while preserving the t1_engine -> t1_windowproof handoff.
    with tempfile.TemporaryDirectory() as cwd:
        print(f"{'harness':<34}result")
        print("-" * 70)
        for name in FAST_LANE:
            r = subprocess.run(
                [sys.executable, str(HARN / name)],
                cwd=cwd, env=env, capture_output=True, text=True,
                encoding="utf-8", errors="replace",
            )
            summary = ""
            for line in reversed((r.stdout or "").splitlines()):
                if "PASSED" in line or "passed" in line:
                    summary = line.strip()
                    break
            ok = r.returncode == 0
            print(f"{name:<34}{'OK   ' if ok else 'FAIL '}{summary}", flush=True)
            if not ok:
                failures.append(name)
                sys.stderr.write(f"\n----- {name} failed (exit {r.returncode}) -----\n")
                sys.stderr.write((r.stdout or "")[-2000:])
                sys.stderr.write((r.stderr or "")[-2000:])
    print("-" * 70)
    if failures:
        print(f"FAST-LANE FAILURES: {failures}")
        return 1
    print(f"ALL {len(FAST_LANE)} FAST-LANE HARNESSES PASSED "
          f"(full 422-check suite: run_all.sh)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
