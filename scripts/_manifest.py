"""The per-contribution manifest: load, schema-validate, and interpret it.

Every papers/<shortname>/ owns a `contribution.yml` — the domain-neutral contract.
This module is the single place that knows the schema, shared by check_structure.py,
run_drift.py, check.py and new_contribution.py.

Schema (see docs/ANTI_DRIFT.md for the friendly version):

    shortname: <string>          # must equal the folder name
    member: <github-handle>      # must be in .github/members.yml
    title: <string>
    domain: <free-form string>   # math | physics-theory | music-theory | bio | ...
    artifacts:                   # non-empty list of files/dirs this contribution owns
      - papers/<shortname>/<doc>
      - code/<shortname>         # dir entries are fine
      - ...
    checks:                      # EITHER a list of checks, OR the literal `none`
      - name: <string>
        run: "<shell command>"   # optional; regenerates outputs. {PY} = the python interpreter.
        produces: [<path>, ...]  # optional (requires run); diffed vs the committed copy
        verify: "<shell command>"# optional; must exit 0 (tests / lint / schema / notebook)
    # --- OR, for a pure-theory / no-repro contribution: ---
    # checks: none
    # reason: "<why there is nothing to reproduce>"   # required when checks: none
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _registry import REPO_ROOT  # noqa: E402

import yaml  # noqa: E402

MANIFEST_NAME = "contribution.yml"
FOLDER_RE = re.compile(r"^\d{4}-\d{2}-[a-z0-9]+(-[a-z0-9]+)*$")
OUTPUT_TREES = ("tests", "code", "data", "figures")


def contribution_shortnames():
    """Every date-named papers/<sn> dir (skips _TEMPLATE and dotfiles)."""
    d = os.path.join(REPO_ROOT, "papers")
    return sorted(
        n for n in os.listdir(d)
        if os.path.isdir(os.path.join(d, n)) and not n.startswith(("_", "."))
    )


def all_tree_shortnames():
    """Union of <sn> folders across papers/ and the four output trees (for orphan checks)."""
    names = set(contribution_shortnames())
    for tree in OUTPUT_TREES:
        td = os.path.join(REPO_ROOT, tree)
        if os.path.isdir(td):
            for n in os.listdir(td):
                if os.path.isdir(os.path.join(td, n)) and not n.startswith(("_", ".")):
                    names.add(n)
    return names


def manifest_path(sn):
    return os.path.join(REPO_ROOT, "papers", sn, MANIFEST_NAME)


def load_manifest(sn):
    with open(manifest_path(sn), encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def declared_trees(m):
    """Which output trees {tests,code,data,figures} the manifest declares (via artifacts)."""
    out = set()
    for a in (m.get("artifacts") or []):
        if not isinstance(a, str):
            continue
        parts = a.replace("\\", "/").strip("/").split("/")
        if len(parts) >= 2 and parts[0] in OUTPUT_TREES:
            out.add(parts[0])
    return out


def is_no_repro(m):
    return isinstance(m, dict) and m.get("checks") == "none"


def validate_manifest(sn, m, members):
    """Return a list of human-readable error strings ([] = valid)."""
    errs = []
    if not isinstance(m, dict):
        return [f"{sn}: {MANIFEST_NAME} is not a YAML mapping"]

    def req_str(field):
        v = m.get(field)
        if not isinstance(v, str) or not v.strip():
            errs.append(f"{sn}: `{field}` must be a non-empty string")
            return None
        return v

    scn = req_str("shortname")
    if scn is not None and scn != sn:
        errs.append(f"{sn}: `shortname` ({scn!r}) must equal the folder name '{sn}'")
    member = req_str("member")
    if member is not None and member not in members:
        errs.append(f"{sn}: `member` '{member}' is not a handle in .github/members.yml")
    req_str("title")
    req_str("domain")

    arts = m.get("artifacts")
    if not isinstance(arts, list) or not arts:
        errs.append(f"{sn}: `artifacts` must be a non-empty list of paths")
    else:
        for a in arts:
            if not isinstance(a, str):
                errs.append(f"{sn}: artifact entries must be strings, got {a!r}")
            elif not os.path.exists(os.path.join(REPO_ROOT, a)):
                errs.append(f"{sn}: declared artifact does not exist: {a}")

    checks = m.get("checks")
    if checks == "none":
        if not (isinstance(m.get("reason"), str) and m["reason"].strip()):
            errs.append(f"{sn}: `checks: none` requires a non-empty `reason` (why nothing reproduces)")
    elif isinstance(checks, list):
        if not checks:
            errs.append(f"{sn}: `checks` is an empty list — use `checks: none` + `reason` for no-repro")
        seen = set()
        for i, c in enumerate(checks):
            if not isinstance(c, dict):
                errs.append(f"{sn}: check #{i} is not a mapping")
                continue
            nm = c.get("name")
            if not isinstance(nm, str) or not nm.strip():
                errs.append(f"{sn}: check #{i} needs a non-empty `name`")
                nm = f"#{i}"
            elif nm in seen:
                errs.append(f"{sn}: duplicate check name {nm!r}")
            else:
                seen.add(nm)
            has_run = isinstance(c.get("run"), str) and c["run"].strip()
            has_verify = isinstance(c.get("verify"), str) and c["verify"].strip()
            if not (has_run or has_verify):
                errs.append(f"{sn}: check {nm!r} needs at least one of `run` or `verify`")
            prod = c.get("produces")
            if prod is not None:
                if not (isinstance(prod, list) and all(isinstance(p, str) for p in prod)):
                    errs.append(f"{sn}: check {nm!r} `produces` must be a list of paths")
                elif not has_run:
                    errs.append(f"{sn}: check {nm!r} declares `produces` but no `run` to regenerate them")
    else:
        errs.append(f"{sn}: `checks` must be a list of checks OR the literal `none`")

    return errs
