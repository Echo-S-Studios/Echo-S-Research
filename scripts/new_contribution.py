#!/usr/bin/env python3
"""Scaffold a new contribution — domain-agnostic.

Usage:
  py scripts/new_contribution.py --member <handle> --shortname <name> --date YYYY-MM \
        --domain <domain> [--title "..."] [--trees code,data,tests,figures]

Always creates papers/<date>-<shortname>/ with a paper stub, paper.cff (attribution
pre-filled from .github/members.yml), and a contribution.yml manifest (the contract).
By DEFAULT the manifest is `checks: none` (pure-theory ready) — a compliant, green
contribution with nothing to reproduce yet. Pass --trees to also scaffold any of the
code/ data/ tests/ figures/ folders (and have them declared in the manifest); then
edit the manifest's `checks` to add reproduction. See docs/ANTI_DRIFT.md.

Refuses an unregistered handle or a folder that already exists.
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _registry import REPO_ROOT, load_members  # noqa: E402

TEMPLATE = os.path.join(REPO_ROOT, "papers", "_TEMPLATE")
DATE_RE = re.compile(r"^\d{4}-\d{2}$")
SHORT_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
VALID_TREES = ("code", "data", "tests", "figures")

# Stub files scaffolded for each optional tree (template-relative -> dest builder).
TREE_FILES = {
    "code": [("code/producer.py", "code/{f}/producer.py"), ("code/README.md", "code/{f}/README.md")],
    "data": [("data/README.md", "data/{f}/README.md")],
    "tests": [("tests/test_example.py", "tests/{f}/test_example.py"), ("tests/NOTES.md", "tests/{f}/NOTES.md")],
    "figures": [("figures/README.md", "figures/{f}/README.md")],
}


def render(text, subs):
    for token, value in subs.items():
        text = text.replace(token, value)
    return text


def build_manifest(folder, shortname, handle, title, domain, trees):
    """Build the contribution.yml text (checks: none default, artifacts per trees)."""
    arts = [f"papers/{folder}/{shortname}.tex", f"papers/{folder}/paper.cff"]
    arts += [f"{t}/{folder}" for t in VALID_TREES if t in trees]
    lines = [
        f"# Contribution contract for papers/{folder}/  (see ../../docs/ANTI_DRIFT.md)",
        f"shortname: {folder}",
        f"member: {handle}",
        'title: "' + title.replace('"', '\\"') + '"',
        f"domain: {domain}",
        "artifacts:",
    ]
    lines += [f"  - {a}" for a in arts]
    lines += [
        "",
        "# `checks` is EITHER the literal `none` (pure theory/prose) with a `reason`, OR a",
        "# list of named checks (run/produces/verify). Replace this once you have something",
        "# reproducible — see docs/ANTI_DRIFT.md for worked examples across domains.",
        "checks: none",
        'reason: "No automated reproduction yet. If this stays pure prose/theory, say so here; '
        'otherwise replace this block with a `checks:` list (see docs/ANTI_DRIFT.md)."',
        "",
    ]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Scaffold a new Echo S Research contribution (any domain).")
    ap.add_argument("--member", required=True, help="your GitHub handle (must be in .github/members.yml)")
    ap.add_argument("--shortname", required=True, help="short hyphenated slug, e.g. spectral-gap")
    ap.add_argument("--date", required=True, help="YYYY-MM")
    ap.add_argument("--domain", required=True,
                    help="free-form domain: math, physics-theory, music-theory, bio, metacybernetics, other, ...")
    ap.add_argument("--title", default="TITLE OF YOUR CONTRIBUTION", help="title (optional)")
    ap.add_argument("--trees", default="",
                    help="comma-separated subset of code,data,tests,figures to also scaffold (default: none)")
    args = ap.parse_args()

    members = load_members()
    if args.member not in members:
        raise SystemExit(
            f"ERROR: '{args.member}' is not a registered member.\n"
            f"       Valid handles (.github/members.yml): {', '.join(sorted(members))}\n"
            f"       To onboard a new member, add ONE line to .github/members.yml first, then re-run."
        )
    if not DATE_RE.match(args.date):
        raise SystemExit(f"ERROR: --date must be YYYY-MM (e.g. 2026-08), got '{args.date}'")
    if not SHORT_RE.match(args.shortname):
        raise SystemExit(f"ERROR: --shortname must be lowercase letters/digits with single hyphens, "
                         f"got '{args.shortname}'")
    if not args.domain.strip():
        raise SystemExit("ERROR: --domain must be a non-empty string")
    trees = [t.strip() for t in args.trees.split(",") if t.strip()]
    bad = [t for t in trees if t not in VALID_TREES]
    if bad:
        raise SystemExit(f"ERROR: --trees may only contain {', '.join(VALID_TREES)}; got {bad}")

    display = members[args.member]
    folder = f"{args.date}-{args.shortname}"
    subs = {
        "__DISPLAY_NAME__": display, "__HANDLE__": args.member, "__SHORTNAME__": args.shortname,
        "__FOLDER__": folder, "__DATE__": args.date, "__DATE_RELEASED__": f"{args.date}-01",
        "__TITLE__": args.title, "__DOMAIN__": args.domain,
    }

    for tree in ("papers", *VALID_TREES):
        if os.path.exists(os.path.join(REPO_ROOT, tree, folder)):
            raise SystemExit(f"ERROR: {tree}/{folder} already exists — pick a different --date/--shortname.")

    # Papers folder: doc stub + paper.cff (from template), plus the built manifest.
    jobs = [("paper.cff", f"papers/{folder}/paper.cff"), ("paper.tex", f"papers/{folder}/{args.shortname}.tex")]
    for t in trees:
        jobs += [(src, dst.format(f=folder)) for src, dst in TREE_FILES[t]]

    print(f"Scaffolding '{folder}' for {display} (@{args.member})  [domain: {args.domain}]:")
    for src_rel, dst_rel in jobs:
        with open(os.path.join(TEMPLATE, src_rel), encoding="utf-8") as fh:
            content = render(fh.read(), subs)
        dst = os.path.join(REPO_ROOT, dst_rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"  created {dst_rel}")

    man_rel = f"papers/{folder}/contribution.yml"
    with open(os.path.join(REPO_ROOT, man_rel), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(build_manifest(folder, args.shortname, args.member, args.title, args.domain, trees))
    print(f"  created {man_rel}")

    print(
        "\nDone. Next steps:\n"
        f"  1. Write your paper in papers/{folder}/ (the {args.shortname}.tex stub, or your own .md/.pdf —\n"
        f"     then update `artifacts` in the manifest).\n"
        f"  2. Decide reproduction: keep `checks: none` (pure theory) with an honest reason, OR add a\n"
        f"     `checks:` list. Read docs/ANTI_DRIFT.md — it has a checklist and worked examples.\n"
        f"  3. Run  py scripts/check.py {folder}  before you push (same checks CI runs).\n"
        f"  4. Commit + open a PR. Reference example: papers/2026-06-salem-slot/."
    )


if __name__ == "__main__":
    main()
