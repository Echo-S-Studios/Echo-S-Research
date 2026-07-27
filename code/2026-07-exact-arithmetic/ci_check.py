#!/usr/bin/env python3
"""ci_check.py -- per-push CI wrapper for the Exact Arithmetic evidence bundle.

Runs the shipped reproduction harness  run_all.py --full  (verify SHA256SUMS, then run
every lesson / census / revision / synthesis / certificate script; exit 0 iff all
checksums match and every SUMMARY reports no failures). We wrap it only to force a UTF-8
child environment so the lessons' Unicode output (mu_S, phi, the blackboard letters, ...)
does not crash a legacy-codepage stdout on a local Windows pre-flight; CI runners are
already UTF-8. This wrapper is NOT part of the checksummed bundle (SHA256SUMS) and changes
nothing about what run_all.py verifies -- it only fixes the child stdout encoding.

    python3 ci_check.py     # exits with run_all.py --full's status (0 == ALL GREEN)
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
sys.exit(
    subprocess.run(
        [sys.executable, str(HERE / "run_all.py"), "--full"], env=env
    ).returncode
)
