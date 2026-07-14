r"""
Producer: the inline TikZ figures of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

The paper has three \begin{tikzpicture} environments:
  figure1  -> Fig 1 (label fig:twofaces)  : the two faces of residual return
  figure2  -> Fig 2 (label fig:nondisjoint): the non-disjoint degree-4 witness
  figure3  -> Fig 3 (label fig:loop)       : the automatic detector-driven loop

This script extracts the paper's own preamble + each tikzpicture, compiles them
with XeLaTeX under the preview package (so all of the paper's macros/colours
resolve), and rasterises each page to figureN.png + figureN.pdf.

Run:  py code/2026-06-residual-return-learning/make_figures.py
"""
import glob
import os
import re
import shutil
import subprocess

import fitz

XELATEX = os.environ.get("XELATEX") or shutil.which("xelatex") or "xelatex"
SHORT = "2026-06-residual-return-learning"
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEX = glob.glob(os.path.join(REPO, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(REPO, "figures", SHORT)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(TEX, encoding="utf-8").read()
    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)
    if not tikz:
        print("No tikzpicture environments found.")
        return 0

    doc = (preamble
           + "\n\\usepackage[active,tightpage]{preview}\n"
           + "\\PreviewEnvironment{tikzpicture}\n"
           + "\\begin{document}\n"
           + "\n\n".join(tikz)
           + "\n\\end{document}\n")
    open(os.path.join(OUT, "_figbuild.tex"), "w", encoding="utf-8").write(doc)

    for _ in range(2):
        subprocess.run([XELATEX, "-interaction=nonstopmode", "_figbuild.tex"],
                       cwd=OUT, capture_output=True, text=True, timeout=240)

    pdf_path = os.path.join(OUT, "_figbuild.pdf")
    if not os.path.exists(pdf_path):
        log = os.path.join(OUT, "_figbuild.log")
        tail = open(log, encoding="utf-8", errors="ignore").read()[-2000:] \
            if os.path.exists(log) else "(no log)"
        raise SystemExit(f"XeLaTeX did not produce a PDF. Log tail:\n{tail}")

    d = fitz.open(pdf_path)
    n = len(d)
    for i, page in enumerate(d):
        page.get_pixmap(dpi=200).save(os.path.join(OUT, f"figure{i+1}.png"))
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i+1}.pdf"))
        one.close()
    d.close()

    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)
    print(f"[make_figures.py] {len(tikz)} tikzpictures -> {n} figure pages in {OUT}")
    return n


if __name__ == "__main__":
    main()
