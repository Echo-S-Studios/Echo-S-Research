#!/usr/bin/env python3
"""Pre-flight: run the EXACT SAME structure + drift validation locally that CI runs,
with a human-readable pass/fail per contribution. Run this BEFORE pushing so
problems are caught before a badge ever goes red.

    py scripts/check.py                 # check everything (structure + full drift)
    py scripts/check.py <shortname>     # structure for all, drift only for one contribution (fast)

Exit code 0 = green (safe to push); non-zero = fix the reported items first.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_structure  # noqa: E402
import run_drift  # noqa: E402


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    bar = "=" * 68
    print(bar)
    print("Echo S Research  —  pre-flight  (same checks CI runs)")
    if only:
        print(f"(drift scoped to contribution: {only})")
    print(bar)

    print("\n[1/2] STRUCTURE — the universal contract for every contribution\n")
    s_ok, s_lines = check_structure.run()
    print("\n".join(s_lines))

    print("\n[2/2] DRIFT — each contribution's declared checks (no-repro ones are skipped)\n")
    d_ok, d_lines = run_drift.run(only)
    print("\n".join(d_lines))

    print("\n" + bar)
    if s_ok and d_ok:
        print("PRE-FLIGHT PASSED  —  safe to push; both CI workflows should be green.")
        return 0
    print("PRE-FLIGHT FAILED  —  fix these before pushing:")
    if not s_ok:
        print("  * structure problems  (see section [1/2] above)")
    if not d_ok:
        print("  * drift problems      (see section [2/2] above)")
    print(bar)
    return 1


if __name__ == "__main__":
    sys.exit(main())
