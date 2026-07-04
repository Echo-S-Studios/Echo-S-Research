#!/usr/bin/env python3
"""Run every DATA producer under code/ to regenerate data/ in place.

Used by the `drift` workflow, which snapshots data/ first and then diffs it against
the regenerated result. Each producer is run as an isolated subprocess (so sibling
imports resolve per-folder and a failure in one is contained).

Skips:
  * make_figures.py            (figure builders need XeLaTeX; figures aren't data)
  * modules without an `if __name__ == "__main__"` guard (shared helper libraries)

Exits 1 if any producer errors. Run:  py scripts/regen_all.py
"""

import glob
import os
import subprocess
import sys

# Self-contained (no third-party imports) so the drift job needs only the
# producers' own deps (sympy/mpmath/numpy) — not the registry's pyyaml.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SKIP_BASENAMES = {"make_figures.py"}


def is_runnable(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return "__main__" in fh.read()
    except Exception:
        return False


def main():
    scripts = sorted(glob.glob(os.path.join(REPO_ROOT, "code", "*", "*.py")))
    ran, failed = 0, []
    for s in scripts:
        if os.path.basename(s) in SKIP_BASENAMES or not is_runnable(s):
            continue
        rel = os.path.relpath(s, REPO_ROOT).replace("\\", "/")
        r = subprocess.run([sys.executable, s], cwd=REPO_ROOT, capture_output=True, text=True)
        ran += 1
        if r.returncode != 0:
            failed.append((rel, (r.stderr or r.stdout or "")[-1000:]))
            print(f"FAIL {rel}")
        else:
            print(f"ok   {rel}")

    print(f"\nran {ran} producers, {len(failed)} failed")
    for rel, err in failed:
        print(f"\n=== {rel} ===\n{err}")
    if failed:
        print("\n::error::one or more producers failed to run — see errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
