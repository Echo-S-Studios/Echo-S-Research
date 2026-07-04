"""
Producer -- figures.

Source: papers/2026-06-generative-emptiness/generative_emptiness.tex

Extracts every tikzpicture from the paper, wraps it in the paper's own preamble
(so all macros/colors -- forcedc, growc, capc, openc, etc. -- resolve), compiles
with XeLaTeX via the preview package, and renders each to figureN.png + figureN.pdf
in figures/2026-06-generative-emptiness/.

The paper contains one tikzpicture: Figure 1 (label fig:charge), "The Z/4Z angle
charge" (Sec 3).

Run:
    py code/2026-06-generative-emptiness/make_figures.py
"""
import glob
import os
import re
import subprocess

import fitz

XELATEX = r"C:\Users\acead\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
SHORT = "2026-06-generative-emptiness"
ROOT = r"C:\Users\acead\projects\Echo-S-Research"
TEX = glob.glob(os.path.join(ROOT, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(ROOT, "figures", SHORT)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(TEX, encoding="utf-8").read()
    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)

    if not tikz:
        print("No tikzpicture found; nothing to render.")
        return 0

    doc = (preamble
           + "\n\\usepackage[active,tightpage]{preview}\n"
           + "\\PreviewEnvironment{tikzpicture}\n"
           + "\\begin{document}\n"
           + "\n\n".join(tikz)
           + "\n\\end{document}\n")
    open(os.path.join(OUT, "_figbuild.tex"), "w", encoding="utf-8").write(doc)

    for _ in range(2):
        subprocess.run(
            [XELATEX, "-interaction=nonstopmode", "_figbuild.tex"],
            cwd=OUT, capture_output=True, text=True, timeout=240,
        )

    pdf_path = os.path.join(OUT, "_figbuild.pdf")
    if not os.path.exists(pdf_path):
        log = os.path.join(OUT, "_figbuild.log")
        msg = open(log, encoding="utf-8", errors="ignore").read()[-2000:] \
            if os.path.exists(log) else "(no log)"
        raise RuntimeError("XeLaTeX did not produce a PDF. Log tail:\n" + msg)

    d = fitz.open(pdf_path)
    for i, page in enumerate(d):
        page.get_pixmap(dpi=200).save(os.path.join(OUT, f"figure{i + 1}.png"))
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i + 1}.pdf"))
        one.close()
    n = d.page_count
    d.close()

    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)

    print(f"rendered {n} figure(s) to {OUT}")
    return n


if __name__ == "__main__":
    main()
