#!/usr/bin/env python3
"""Assemble the GitHub Pages site into _site/ from the papers/ tree.

Manifest-driven, like everything else in this repo: each papers/YYYY-MM-shortname/
contributes one card, built from its contribution.yml (title, member) and paper.cff
(date-released). Adding a paper by the normal contract makes it appear on the site
with NO edit to this script or the workflow.

The site serves the COMMITTED PDFs — the same sealed builds Zenodo archives on each
release — rather than recompiling the .tex in CI. In this repo the committed PDF is
the canonical artifact (see docs/ANTI_DRIFT.md); the site must match the archive.
PDFs are copied to _site/papers/<shortname>.pdf so served URLs are stable slugs even
when the committed filename carries a version suffix or a space.

Optional per-paper card copy lives in site/blurbs.json:
    { "<shortname>": { "series": "...", "description": "..." }, ... }
A paper without a blurb still gets a card (title, member, date, page count).

Inputs are UTF-8 (titles contain λ, Z₅, 𝔄, φ). Everything is read and written with
an explicit encoding so the build is identical on Windows and the Ubuntu runner.

Run from anywhere:  py site/build_site.py   (deps: pyyaml, pypdf)
"""

import html
import json
import re
import shutil
import sys
from datetime import date
from pathlib import Path

import yaml
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
OUT = ROOT / "_site"
REPO_URL = "https://github.com/Echo-S-Studios/Echo-S-Research"
FOLDER_RE = re.compile(r"^\d{4}-\d{2}-[a-z0-9][a-z0-9-]*$")

MONTHS = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# A card heading shows the paper's main title; the descriptive subtitle after the
# first colon or spaced em-dash goes only in the aria-label / catalog. This keeps
# long academic titles legible as cards without discarding the full title.
SUBTITLE_RE = re.compile(r":\s|\s—\s")


def display_title(title: str) -> str:
    return SUBTITLE_RE.split(title, maxsplit=1)[0].strip()


def fail(msg: str) -> None:
    print(f"::error title=site build::{msg}")
    sys.exit(1)


def load_yaml(path: Path):
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def pick_artifact(manifest: dict, sn: str, suffix: str) -> str | None:
    """The papers/<sn>/ artifact with this suffix, as declared in the manifest."""
    prefix = f"papers/{sn}/"
    hits = [a for a in manifest.get("artifacts", [])
            if isinstance(a, str) and a.startswith(prefix) and a.lower().endswith(suffix)]
    return hits[0] if hits else None


def collect_papers():
    papers = []
    for folder in sorted((ROOT / "papers").iterdir()):
        sn = folder.name
        if not folder.is_dir() or not FOLDER_RE.match(sn):
            continue  # _TEMPLATE, loose files
        manifest = load_yaml(folder / "contribution.yml")

        pdf_rel = pick_artifact(manifest, sn, ".pdf")
        if pdf_rel is None:  # manifest declares none — fall back to the folder itself
            pdfs = sorted(folder.glob("*.pdf"))
            if len(pdfs) != 1:
                fail(f"papers/{sn}: expected exactly one PDF (manifest declares none, "
                     f"folder has {len(pdfs)})")
            pdf_rel = f"papers/{sn}/{pdfs[0].name}"
        pdf = ROOT / pdf_rel
        if not pdf.is_file():
            fail(f"papers/{sn}: declared PDF missing on disk: {pdf_rel}")

        tex_rel = pick_artifact(manifest, sn, ".tex")

        released = ""
        cff = folder / "paper.cff"
        if cff.is_file():
            released = str(load_yaml(cff).get("date-released", "") or "")

        try:
            pages = len(PdfReader(pdf).pages)
        except Exception:
            pages = None  # unreadable page tree — card just omits the count

        papers.append({
            "shortname": sn,
            "title": str(manifest.get("title", sn)),
            "member": str(manifest.get("member", "")),
            "date": released,
            "pages": pages,
            "src_pdf": pdf,
            "file": f"papers/{sn}.pdf",
            "tex": f"{REPO_URL}/blob/main/{tex_rel}" if tex_rel else None,
        })
    if not papers:
        fail("no papers/YYYY-MM-*/ contributions found")
    # Newest first; folder name (YYYY-MM-...) breaks ties when date-released is absent.
    papers.sort(key=lambda p: (p["date"], p["shortname"]), reverse=True)
    return papers


def month_label(released: str, sn: str) -> str:
    m = re.match(r"^(\d{4})-(\d{2})", released or sn)
    return f"{MONTHS[int(m.group(2))]} {m.group(1)}" if m else ""


def render_card(p: dict) -> str:
    full = html.escape(p["title"])
    heading = html.escape(display_title(p["title"]))
    series = html.escape(p["series"]) if p.get("series") else "Echo S Research"
    member = html.escape(p["member"])
    desc = html.escape(p["description"]) if p.get("description") else ""
    pp = f"{p['pages']} pp" if p["pages"] else ""
    when = month_label(p["date"], p["shortname"])
    meta = " · ".join(x for x in (pp, when) if x)
    aria_pp = f", {p['pages']} pages" if p["pages"] else ""
    lines = [
        f'      <a class="card" href="{html.escape(p["file"])}" target="_blank" rel="noopener"',
        f'         aria-label="{full} — PDF{aria_pp}, opens in a new tab">',
        f'        <div class="toprow"><span class="tag pdf">PDF ↗</span><span class="pp">{meta}</span></div>',
        f'        <p class="series">{series} · {member}</p>',
        f'        <h3>{heading}</h3>',
    ]
    if desc:
        lines.append(f'        <p>{desc}</p>')
    lines += [
        '        <span class="go">Read PDF <span class="arr">→</span></span>',
        '      </a>',
    ]
    return "\n".join(lines)


def main() -> None:
    blurbs_path = SITE / "blurbs.json"
    blurbs = json.loads(blurbs_path.read_text(encoding="utf-8")) if blurbs_path.is_file() else {}

    papers = collect_papers()
    for p in papers:
        # blurbs.json may override the manifest title (correct Unicode rendering) and
        # supplies the series label + card description.
        p.update({k: v for k, v in blurbs.get(p["shortname"], {}).items()
                  if k in ("title", "series", "description")})

    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "papers").mkdir(parents=True)

    for p in papers:
        shutil.copyfile(p["src_pdf"], OUT / p["file"])

    catalog = [{k: p[k] for k in
                ("shortname", "title", "member", "date", "pages", "file", "tex")}
               | {"series": p.get("series"), "description": p.get("description")}
               for p in papers]
    (OUT / "papers" / "catalog.json").write_text(
        json.dumps(catalog, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    template = (SITE / "index.html").read_text(encoding="utf-8")
    cards = "\n\n".join(render_card(p) for p in papers)
    out_html = (template
                .replace("<!--@CARDS@-->", cards)
                .replace("@COUNT@", str(len(papers)))
                .replace("@UPDATED@", date.today().isoformat()))
    if out_html == template:
        fail("site/index.html has no <!--@CARDS@--> placeholder — nothing was injected")
    (OUT / "index.html").write_text(out_html, encoding="utf-8")

    print(f"built _site: {len(papers)} papers, "
          f"{sum(1 for p in papers if p.get('description'))} with blurbs")


if __name__ == "__main__":
    main()
