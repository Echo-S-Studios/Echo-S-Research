#!/usr/bin/env python3
"""Scaffold a new contribution across all five parallel trees.

Usage:
  py scripts/new_contribution.py --member <github_handle> --shortname <name> --date YYYY-MM [--title "..."]

Validates <github_handle> against .github/members.yml, then generates (mirroring
the layout of every existing paper, e.g. papers/2026-06-salem-slot/):

  papers/<date>-<shortname>/   paper.cff (attribution pre-filled), <shortname>.tex
  tests/<date>-<shortname>/    test_example.py, NOTES.md
  code/<date>-<shortname>/     producer.py, README.md
  data/<date>-<shortname>/     README.md
  figures/<date>-<shortname>/  README.md

Refuses to run for an unregistered handle or if any target folder already exists.
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

# (template-relative source, destination-relative-to-repo-root builder)
FILE_MAP = [
    ("paper.cff", lambda f, s: f"papers/{f}/paper.cff"),
    ("paper.tex", lambda f, s: f"papers/{f}/{s}.tex"),
    ("tests/test_example.py", lambda f, s: f"tests/{f}/test_example.py"),
    ("tests/NOTES.md", lambda f, s: f"tests/{f}/NOTES.md"),
    ("code/producer.py", lambda f, s: f"code/{f}/producer.py"),
    ("code/README.md", lambda f, s: f"code/{f}/README.md"),
    ("data/README.md", lambda f, s: f"data/{f}/README.md"),
    ("figures/README.md", lambda f, s: f"figures/{f}/README.md"),
]


def render(text, subs):
    for token, value in subs.items():
        text = text.replace(token, value)
    return text


def main():
    ap = argparse.ArgumentParser(description="Scaffold a new Echo S Research contribution.")
    ap.add_argument("--member", required=True, help="your GitHub handle (must be in .github/members.yml)")
    ap.add_argument("--shortname", required=True, help="short hyphenated slug, e.g. spectral-gap")
    ap.add_argument("--date", required=True, help="YYYY-MM (year and month of the paper)")
    ap.add_argument("--title", default="TITLE OF YOUR PAPER", help="paper title (optional)")
    args = ap.parse_args()

    members = load_members()
    if args.member not in members:
        valid = ", ".join(sorted(members))
        raise SystemExit(
            f"ERROR: '{args.member}' is not a registered member.\n"
            f"       Valid handles (.github/members.yml): {valid}\n"
            f"       To onboard a new member, add ONE line to .github/members.yml first, then re-run."
        )
    if not DATE_RE.match(args.date):
        raise SystemExit(f"ERROR: --date must be YYYY-MM (e.g. 2026-08), got '{args.date}'")
    if not SHORT_RE.match(args.shortname):
        raise SystemExit(
            f"ERROR: --shortname must be lowercase letters/digits with single hyphens "
            f"(e.g. spectral-gap), got '{args.shortname}'"
        )

    display = members[args.member]
    folder = f"{args.date}-{args.shortname}"
    subs = {
        "__DISPLAY_NAME__": display,
        "__HANDLE__": args.member,
        "__SHORTNAME__": args.shortname,
        "__FOLDER__": folder,
        "__DATE__": args.date,
        "__DATE_RELEASED__": f"{args.date}-01",
        "__TITLE__": args.title,
    }

    # Refuse to clobber an existing contribution.
    for tree in ("papers", "tests", "code", "data", "figures"):
        d = os.path.join(REPO_ROOT, tree, folder)
        if os.path.exists(d):
            raise SystemExit(f"ERROR: {tree}/{folder} already exists — pick a different --date/--shortname.")

    print(f"Scaffolding '{folder}' for {display} (@{args.member}):")
    for src_rel, dst_builder in FILE_MAP:
        src = os.path.join(TEMPLATE, src_rel)
        with open(src, encoding="utf-8") as fh:
            content = render(fh.read(), subs)
        dst_rel = dst_builder(folder, args.shortname)
        dst = os.path.join(REPO_ROOT, dst_rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"  created {dst_rel}")

    print(
        "\nDone. Next steps:\n"
        f"  1. Write papers/{folder}/{args.shortname}.tex and compile its PDF (XeLaTeX/LuaLaTeX).\n"
        f"  2. Fill code/{folder}/ producers (emit to data/{folder}/) and tests/{folder}/ verifiers.\n"
        f"  3. Add any new citations to papers/references.bib.\n"
        f"  4. Commit and open a PR — the drift and structure workflows will check it.\n"
        f"  Reference example to copy: papers/2026-06-salem-slot/ (+ its tests/code/data/figures)."
    )


if __name__ == "__main__":
    main()
