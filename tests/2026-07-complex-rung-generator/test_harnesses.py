"""Drift check for the complex-rung-generator contribution.

The paper's verification "pipes" are the softcheck harnesses in
code/2026-07-complex-rung-generator/. Each is a self-contained exact-arithmetic
verifier (sympy over Q/Q(sqrt5) extended by i; mpmath only for separately-counted
interval guards / numeric displays) that prints a pass-line and exits 0 iff every
one of its checks passes. This test runs each harness as an isolated subprocess and
asserts exit 0 — so a single `verify: pytest ...` in contribution.yml re-certifies
the whole verification layer on every push.

Discovery is by glob, so adding a harness to code/ is picked up automatically.
"""
import glob
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))            # tests/<sn>/ -> repo root
CODE = os.path.join(REPO, "code", "2026-07-complex-rung-generator")

HARNESSES = sorted(os.path.basename(p) for p in glob.glob(os.path.join(CODE, "*.py")))


def test_harnesses_present():
    # The 15 shipped softcheck harnesses + the decimal audit.
    assert len(HARNESSES) >= 16, f"expected >=16 pipes, found {HARNESSES}"


@pytest.mark.parametrize("harness", HARNESSES)
def test_harness_exits_zero(harness):
    path = os.path.join(CODE, harness)
    r = subprocess.run(
        [sys.executable, path], cwd=REPO,
        capture_output=True, text=True, timeout=600,
    )
    tail = (r.stderr or r.stdout or "")[-2000:]
    assert r.returncode == 0, f"{harness} exited {r.returncode}\n{tail}"
