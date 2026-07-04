"""
Producer: reconstruct the inline TikZ figures of *Lehmer's Box*.

Source paper: papers/2026-06-lehmers-box/lehmers_box.tex.
The paper renders its figures inline with TikZ.  This script extracts every
\\begin{tikzpicture}...\\end{tikzpicture} block, wraps them in the paper's OWN
preamble (so every custom macro and colour -- \\forced, forcedc, boxfill,
forbidfill, ... -- resolves identically), compiles with XeLaTeX under the
`preview` package (tight bounding box), and rasterises each page to
figureN.png / figureN.pdf.

Figure 1 (label fig:box) is the height-angle rendering of Lehmer's Box: the
green floor wall at log(phi), the four blue lattice posts at (pi/2)Z, the
forbidden strip (1,phi), phi and beta_4 on the theta=0 post, and Lehmer L
outside on both walls.

Outputs (in figures/2026-06-lehmers-box/):
    figure1.png, figure1.pdf   (paper Figure 1, label fig:box)
"""

from __future__ import annotations

import glob
import os
import re
import subprocess

import fitz

XELATEX = r"C:\Users\acead\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
SHORT = "2026-06-lehmers-box"
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEX = glob.glob(os.path.join(ROOT, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(ROOT, "figures", SHORT)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(TEX, encoding="utf-8").read()
    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)

    if not tikz:
        print("No tikzpicture environments found in", TEX)
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

    pdf = os.path.join(OUT, "_figbuild.pdf")
    if not os.path.exists(pdf):
        print("ERROR: XeLaTeX did not produce", pdf)
        return 0

    d = fitz.open(pdf)
    n = 0
    for i, page in enumerate(d):
        page.get_pixmap(dpi=200).save(os.path.join(OUT, f"figure{i + 1}.png"))
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i + 1}.pdf"))
        one.close()
        n += 1
    d.close()

    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)

    print(f"rendered {n} figure(s) from {len(tikz)} tikzpicture block(s) to {OUT}")
    return n


if __name__ == "__main__":
    main()
