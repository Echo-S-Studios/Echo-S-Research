#!/usr/bin/env python3
"""Validate repository invariants. Exits 1 with human-readable messages on any
violation (it reports ALL problems, not just the first).

Checks:
  1. Every papers/<YYYY-MM-shortname>/ has its required files (paper.cff + a .tex).
  2. Every paper.cff author `alias` is a github handle in .github/members.yml
     (missing alias or unknown handle = fail); same for CITATION.cff members.
  3. The four parallel trees tests/ code/ data/ figures/ each have a matching
     <shortname> folder for every paper (gaps or orphans = fail).
  4. Every contribution folder name follows YYYY-MM-shortname.
  5. papers/references.bib parses (and is non-empty).
  6. .zenodo.json is valid JSON; every paper.cff / CITATION.cff parses as YAML.

Deep CFF-schema validation (cffconvert) of CITATION.cff and each paper.cff is done
by the workflow as a separate step. Run locally with: py scripts/check_structure.py
"""

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _registry import REPO_ROOT, load_members  # noqa: E402

import yaml  # noqa: E402

FOLDER_RE = re.compile(r"^\d{4}-\d{2}-[a-z0-9]+(-[a-z0-9]+)*$")
TREES = ("papers", "tests", "code", "data", "figures")
FAILS = []
NOTES = []


def fail(msg):
    FAILS.append(msg)


def note(msg):
    NOTES.append(msg)


def list_folders(tree):
    """Contribution folders in a tree (dirs only; skip _-prefixed like _TEMPLATE)."""
    d = os.path.join(REPO_ROOT, tree)
    if not os.path.isdir(d):
        return {}
    return {
        name: os.path.join(d, name)
        for name in os.listdir(d)
        if os.path.isdir(os.path.join(d, name))
        and not name.startswith("_")
        and not name.startswith(".")
    }


def main():
    members = load_members()
    print(f"registry: {len(members)} members ({', '.join(sorted(members))})")

    papers = list_folders("papers")
    print(f"papers: {len(papers)} contribution folders")

    # 4. naming, across all trees
    for tree in TREES:
        for name in list_folders(tree):
            if not FOLDER_RE.match(name):
                fail(f"[naming] {tree}/{name} does not match YYYY-MM-shortname")

    # 1 + 2 + 6. required files, attribution, cff-parse
    for name, path in sorted(papers.items()):
        cff = os.path.join(path, "paper.cff")
        files = os.listdir(path)
        if not any(f.endswith(".tex") for f in files):
            fail(f"[required] papers/{name}/ has no .tex source")
        if not any(f.endswith(".pdf") for f in files):
            note(f"papers/{name}/ has no compiled .pdf yet (recommended before release)")
        if not os.path.exists(cff):
            fail(f"[required] papers/{name}/ is missing paper.cff")
            continue
        try:
            doc = yaml.safe_load(open(cff, encoding="utf-8"))
        except Exception as e:
            fail(f"[cff] papers/{name}/paper.cff does not parse as YAML: {e}")
            continue
        authors = (doc or {}).get("authors") or []
        if not authors:
            fail(f"[attribution] papers/{name}/paper.cff has no authors")
        for a in authors:
            a = a or {}
            alias = a.get("alias")
            if not alias:
                fail(f"[attribution] papers/{name}/paper.cff author {a.get('name', '?')!r} "
                     f"has no `alias` (GitHub handle)")
            elif alias not in members:
                fail(f"[attribution] papers/{name}/paper.cff alias '{alias}' is not in "
                     f".github/members.yml (register the handle or fix the alias)")

    # CITATION.cff member aliases
    cit_path = os.path.join(REPO_ROOT, "CITATION.cff")
    try:
        cit = yaml.safe_load(open(cit_path, encoding="utf-8")) or {}
        for a in cit.get("authors") or []:
            alias = (a or {}).get("alias")
            if alias and alias not in members:
                fail(f"[attribution] CITATION.cff alias '{alias}' is not in .github/members.yml")
    except Exception as e:
        fail(f"[cff] CITATION.cff does not parse as YAML: {e}")

    # 3. parallel-tree parity
    paper_set = set(papers)
    for tree in ("tests", "code", "data", "figures"):
        s = set(list_folders(tree))
        for g in sorted(paper_set - s):
            fail(f"[parity] {tree}/ is missing a folder for paper '{g}'")
        for o in sorted(s - paper_set):
            fail(f"[parity] {tree}/{o} has no matching papers/{o}/ (orphan)")

    # 5. references.bib parses
    bib = os.path.join(REPO_ROOT, "papers", "references.bib")
    try:
        import bibtexparser
        with open(bib, encoding="utf-8") as fh:
            db = bibtexparser.load(fh)
        if not db.entries:
            fail("[bib] papers/references.bib parsed but has no entries")
        else:
            print(f"references.bib: {len(db.entries)} entries")
    except Exception as e:
        fail(f"[bib] papers/references.bib does not parse: {e}")

    # 6. .zenodo.json valid JSON
    zen = os.path.join(REPO_ROOT, ".zenodo.json")
    try:
        json.load(open(zen, encoding="utf-8"))
    except Exception as e:
        fail(f"[zenodo] .zenodo.json is not valid JSON: {e}")

    for n in NOTES:
        print(f"  note: {n}")

    if FAILS:
        print(f"\nSTRUCTURE CHECK FAILED - {len(FAILS)} problem(s):", file=sys.stderr)
        for f in FAILS:
            print(f"  FAIL: {f}", file=sys.stderr)
        sys.exit(1)
    print(f"\nSTRUCTURE OK - {len(papers)} papers, all invariants hold.")


if __name__ == "__main__":
    main()
