"""
Producer: the two inline TikZ figures of the paper, rendered to PNG + PDF.

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * figures/2026-06-salem-slot/figure1.png / figure1.pdf
        -> Figure 1 (label fig:trace, in Section 3): the trace-line
           trifurcation, a would-be Salem redirected by trace to tau0 in GROW.
  * figures/2026-06-salem-slot/figure2.png / figure2.pdf
        -> Figure 2 (label fig:edge, in Section 4): the grow root beta(t) with
           the square-root branch point at the flip t=2.
  (Figures use the article-class flat counter, so they are Figure 1 and 2 in
   source order; only theorems are section-numbered in this paper.)

The paper's own preamble is reused verbatim so every macro/color/tikzlibrary
(\trace, growc, capc, openc, arrows.meta, ...) resolves. Each tikzpicture is
extracted from the source and previewed on a tight page.

Run:  py code/2026-06-salem-slot/make_figures.py
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import subprocess

import fitz  # PyMuPDF

XELATEX = os.environ.get("XELATEX") or shutil.which("xelatex") or "xelatex"
SHORT = "2026-06-salem-slot"
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEX = glob.glob(os.path.join(ROOT, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(ROOT, "figures", SHORT)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(TEX, encoding="utf-8").read()
    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)
    if not tikz:
        print("no tikzpicture environments found")
        return
    print(f"found {len(tikz)} tikzpicture environment(s)")

    doc = (preamble
           + "\n\\usepackage[active,tightpage]{preview}\n"
           + "\\PreviewEnvironment{tikzpicture}\n"
           + "\\begin{document}\n"
           + "\n\n".join(tikz)
           + "\n\\end{document}\n")
    open(os.path.join(OUT, "_figbuild.tex"), "w", encoding="utf-8").write(doc)

    for _ in range(2):
        r = subprocess.run(
            [XELATEX, "-interaction=nonstopmode", "_figbuild.tex"],
            cwd=OUT, capture_output=True, text=True, timeout=240,
        )
    pdf = os.path.join(OUT, "_figbuild.pdf")
    if not os.path.exists(pdf):
        print("XELATEX did not produce a PDF; tail of log:")
        print(r.stdout[-2500:])
        return

    d = fitz.open(pdf)
    made = []
    for i, page in enumerate(d):
        png = os.path.join(OUT, f"figure{i + 1}.png")
        page.get_pixmap(dpi=200).save(png)
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i + 1}.pdf"))
        one.close()
        made.append(f"figure{i + 1}")
    d.close()

    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)
    print("wrote:", ", ".join(f"{m}.png/{m}.pdf" for m in made))


if __name__ == "__main__":
    main()
