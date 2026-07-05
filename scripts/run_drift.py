#!/usr/bin/env python3
"""Manifest-driven drift check: for every contribution, run its declared checks.

Domain-agnostic. For each papers/<sn>/contribution.yml:
  * `checks: none`  -> SKIP, logged as "no-repro (reason)" (never a failure);
  * otherwise, for each check:
      - snapshot its `produces` files/dirs,
      - run `run` (fail if non-zero exit),
      - diff each `produces` against the snapshot, EOL-agnostic (fail on any change),
      - run `verify` (fail if non-zero exit).

Commands are arbitrary shell (any language/tool). `{PY}` is substituted with this
interpreter, so the same manifest works on a contributor's machine and in CI.

`run()` returns (ok, lines) for reuse by scripts/check.py. Run standalone:
    py scripts/run_drift.py
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _manifest import contribution_shortnames, load_manifest, manifest_path  # noqa: E402
from _registry import REPO_ROOT  # noqa: E402


def _subst(cmd):
    # Quote so a path with spaces (e.g. C:\Program Files\...) works on cmd and sh.
    return cmd.replace("{PY}", '"' + sys.executable + '"')


def _shell(cmd):
    r = subprocess.run(_subst(cmd), shell=True, cwd=REPO_ROOT, capture_output=True, text=True)
    return r.returncode, (r.stderr or r.stdout or "")


def _norm(b):
    return b.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def _files_under(path_rel):
    """Yield (rel_to_repo, abspath) for a produces entry that is a file or a dir."""
    ap = os.path.join(REPO_ROOT, path_rel)
    if os.path.isdir(ap):
        for root, _dirs, files in os.walk(ap):
            for f in files:
                fp = os.path.join(root, f)
                yield os.path.relpath(fp, REPO_ROOT).replace("\\", "/"), fp
    else:
        yield path_rel.replace("\\", "/"), ap


def _snapshot(produces):
    snap = {}
    for p in produces:
        for rel, ap in _files_under(p):
            snap[rel] = _norm(open(ap, "rb").read()) if os.path.exists(ap) else None
    return snap


def _diff_after(produces, snap):
    """Return list of changed/added/removed paths (EOL-agnostic)."""
    changed = []
    now = {}
    for p in produces:
        for rel, ap in _files_under(p):
            now[rel] = _norm(open(ap, "rb").read()) if os.path.exists(ap) else None
    for rel in sorted(set(snap) | set(now)):
        before, after = snap.get(rel), now.get(rel)
        if before is None and after is not None:
            changed.append(rel + " (new / not committed)")
        elif before is not None and after is None:
            changed.append(rel + " (removed by run)")
        elif before != after:
            changed.append(rel + " (content changed)")
    return changed


def run(only=None):
    """Run drift checks. If `only` is a shortname, check just that contribution."""
    lines = []
    fails = []
    n_checks = n_skip = 0
    targets = [only] if only else contribution_shortnames()
    for sn in targets:
        if not os.path.exists(manifest_path(sn)):
            fails.append(f"{sn}: missing contribution.yml")
            continue
        m = load_manifest(sn)
        checks = m.get("checks")
        if checks == "none":
            n_skip += 1
            lines.append(f"  skip {sn}: no-repro ({m.get('reason', '')[:70]})")
            continue
        if not isinstance(checks, list):
            fails.append(f"{sn}: invalid `checks` (run structure check)")
            continue
        for c in checks:
            name = c.get("name", "?")
            produces = c.get("produces") or []
            snap = _snapshot(produces) if produces else {}
            run_cmd, verify_cmd = c.get("run"), c.get("verify")
            n_checks += 1
            failed = False
            if run_cmd:
                rc, out = _shell(run_cmd)
                if rc != 0:
                    fails.append(f"{sn} / {name}: `run` failed (exit {rc})\n{out[-1000:]}")
                    failed = True
            if not failed and produces:
                drift = _diff_after(produces, snap)
                if drift:
                    fails.append(f"{sn} / {name}: DRIFT — regenerated output differs from committed:\n    "
                                 + "\n    ".join(drift[:20]))
                    failed = True
            if not failed and verify_cmd:
                rc, out = _shell(verify_cmd)
                if rc != 0:
                    fails.append(f"{sn} / {name}: `verify` failed (exit {rc})\n{out[-1000:]}")
                    failed = True
            if not failed:
                lines.append(f"  ok   {sn} / {name}")

    ok = not fails
    if fails:
        lines.append("")
        lines.append(f"DRIFT FAILED - {len(fails)} problem(s):")
        for f in fails:
            lines.append(f"  FAIL {f}")
    else:
        lines.append(f"\nDRIFT OK - {n_checks} checks passed, {n_skip} no-repro skipped.")
    return ok, lines


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    ok, lines = run(only)
    print("\n".join(lines))
    if not ok:
        print("::error title=Drift detected::A contribution's declared checks failed — see above.")
    sys.exit(0 if ok else 1)
