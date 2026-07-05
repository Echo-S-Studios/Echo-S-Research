#!/usr/bin/env python3
"""Run all data producers under code/<shortname>/ (as isolated subprocesses).

A small helper referenced by the math contributions' manifests
(`run: "{PY} scripts/run_producers.py <shortname>"`). Skips make_figures.py
(figures need a TeX toolchain) and helper modules without a `__main__` guard.
Self-contained (stdlib only), so it needs no extra deps.

    py scripts/run_producers.py <shortname>
"""

import glob
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_producers.py <shortname>")
    sn = sys.argv[1]
    code_dir = os.path.join(REPO_ROOT, "code", sn)
    if not os.path.isdir(code_dir):
        raise SystemExit(f"ERROR: no code/{sn}/ directory")
    scripts = sorted(glob.glob(os.path.join(code_dir, "*.py")))
    ran = 0
    for s in scripts:
        if os.path.basename(s) == "make_figures.py":
            continue
        if "__main__" not in open(s, encoding="utf-8").read():
            continue
        r = subprocess.run([sys.executable, s], cwd=REPO_ROOT, capture_output=True, text=True)
        ran += 1
        if r.returncode != 0:
            rel = os.path.relpath(s, REPO_ROOT).replace("\\", "/")
            sys.stderr.write(f"producer failed: {rel}\n{(r.stderr or r.stdout)[-1200:]}\n")
            sys.exit(1)
    print(f"ran {ran} producers for code/{sn}/")


if __name__ == "__main__":
    main()
