#!/usr/bin/env python3
"""Validate the universal contribution contract for EVERY contribution, regardless
of domain. Manifest-driven and domain-agnostic (no assumption of math/LaTeX/pytest).

For each papers/<shortname>/ it checks:
  * contribution.yml present and schema-valid (see scripts/_manifest.py);
  * `member` is a handle in .github/members.yml (attribution present);
  * every declared artifact exists;
  * the parallel trees tests/ code/ data/ figures/ exist EXACTLY where the manifest
    declares artifacts of that kind (a prose-only contribution isn't forced to have
    code/, and an undeclared tree folder is flagged as an orphan);
  * folder naming YYYY-MM-shortname;
  * if a paper.cff is present, its author aliases are registered.
Plus repo-level: references.bib parses, .zenodo.json is valid JSON, CITATION.cff
member aliases are registered.

`run()` returns (ok, lines) so scripts/check.py can reuse it. Run standalone:
    py scripts/check_structure.py
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _manifest import (  # noqa: E402
    FOLDER_RE, MANIFEST_NAME, OUTPUT_TREES, all_tree_shortnames,
    contribution_shortnames, declared_trees, load_manifest, manifest_path,
    validate_manifest,
)
from _registry import REPO_ROOT, load_members  # noqa: E402

import yaml  # noqa: E402


def run():
    lines = []
    fails = []

    def fail(msg):
        fails.append(msg)

    members = load_members()
    lines.append(f"registry: {len(members)} members ({', '.join(sorted(members))})")

    every_sn = sorted(all_tree_shortnames())

    # Naming (every folder in any tree).
    for sn in every_sn:
        if not FOLDER_RE.match(sn):
            fail(f"[naming] '{sn}' is not YYYY-MM-shortname")

    contributions = set(contribution_shortnames())
    for sn in every_sn:
        if not os.path.exists(manifest_path(sn)):
            if sn in contributions:
                fail(f"[manifest] papers/{sn}/ is missing {MANIFEST_NAME}")
            else:
                fail(f"[orphan] a tree folder '{sn}/' exists but there is no papers/{sn}/ (no contribution)")
            continue
        try:
            m = load_manifest(sn)
        except Exception as e:
            fail(f"[manifest] papers/{sn}/{MANIFEST_NAME} does not parse: {e}")
            continue

        errs = validate_manifest(sn, m, members)

        # Parallel-tree parity: declared-in-manifest  <=>  folder-exists.
        decl = declared_trees(m)
        for tree in OUTPUT_TREES:
            exists = os.path.isdir(os.path.join(REPO_ROOT, tree, sn))
            if tree in decl and not exists:
                errs.append(f"{sn}: manifest declares {tree}/ artifacts but {tree}/{sn}/ does not exist")
            if exists and tree not in decl:
                errs.append(f"{sn}: {tree}/{sn}/ exists but the manifest declares no {tree}/ artifacts (orphan)")

        # paper.cff (optional) — validate aliases if it's there.
        cff = os.path.join(REPO_ROOT, "papers", sn, "paper.cff")
        if os.path.exists(cff):
            try:
                doc = yaml.safe_load(open(cff, encoding="utf-8")) or {}
                for a in doc.get("authors") or []:
                    al = (a or {}).get("alias")
                    if al and al not in members:
                        errs.append(f"{sn}: paper.cff alias '{al}' is not in .github/members.yml")
            except Exception as e:
                errs.append(f"{sn}: paper.cff does not parse: {e}")

        if errs:
            for e in errs:
                fail(f"[contract] {e}")
        else:
            tag = m.get("domain", "?")
            extra = f", no-repro: {m.get('reason','')[:50]}" if m.get("checks") == "none" else ""
            lines.append(f"  ok  {sn}  [{tag}{extra}]")

    # Repo-level metadata.
    try:
        import bibtexparser
        with open(os.path.join(REPO_ROOT, "papers", "references.bib"), encoding="utf-8") as fh:
            db = bibtexparser.load(fh)
        if not db.entries:
            fail("[bib] papers/references.bib parsed but has no entries")
        else:
            lines.append(f"references.bib: {len(db.entries)} entries")
    except Exception as e:
        fail(f"[bib] papers/references.bib does not parse: {e}")

    try:
        json.load(open(os.path.join(REPO_ROOT, ".zenodo.json"), encoding="utf-8"))
    except Exception as e:
        fail(f"[zenodo] .zenodo.json is not valid JSON: {e}")

    try:
        cit = yaml.safe_load(open(os.path.join(REPO_ROOT, "CITATION.cff"), encoding="utf-8")) or {}
        for a in cit.get("authors") or []:
            al = (a or {}).get("alias")
            if al and al not in members:
                fail(f"[attribution] CITATION.cff alias '{al}' is not in .github/members.yml")
    except Exception as e:
        fail(f"[cff] CITATION.cff does not parse: {e}")

    ok = not fails
    if fails:
        lines.append("")
        lines.append(f"STRUCTURE FAILED - {len(fails)} problem(s):")
        lines.extend(f"  FAIL {f}" for f in fails)
    else:
        lines.append(f"\nSTRUCTURE OK - {len(contributions)} contributions, all invariants hold.")
    return ok, lines


if __name__ == "__main__":
    ok, lines = run()
    print("\n".join(lines))
    sys.exit(0 if ok else 1)
