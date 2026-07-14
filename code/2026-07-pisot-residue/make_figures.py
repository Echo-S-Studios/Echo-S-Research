"""Producer: figure reconstruction for the Pisot Cross-Shell Residue paper.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex

The paper's body contains inline figures only if its .tex has a
\\begin{tikzpicture} environment.  This script scans the source: if any tikz
pictures are present it rebuilds each (reusing the paper's own preamble so all
macros/colors resolve) into figures/2026-07-pisot-residue/figureN.{png,pdf};
otherwise it writes a README recording that the paper has no figures.

This paper has NO tikzpicture / figure / includegraphics environments, so this
script writes the "no figures" README.

Emits (this paper): figures/2026-07-pisot-residue/README.md
"""
import glob
import os
import re
import shutil
import subprocess

XELATEX = os.environ.get("XELATEX") or shutil.which("xelatex") or "xelatex"
SHORT = "2026-07-pisot-residue"
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEX = glob.glob(os.path.join(ROOT, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(ROOT, "figures", SHORT)


def main():
    os.makedirs(OUT, exist_ok=True)
    src = open(TEX, encoding="utf-8").read()
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)

    if not tikz:
        readme = os.path.join(OUT, "README.md")
        with open(readme, "w", encoding="utf-8") as fh:
            fh.write(f"No figures appear in papers/{SHORT}/{os.path.basename(TEX)}.\n\n")
            fh.write("The paper is a computational whitepaper: all results are reported as "
                     "inline tables and prose. It contains no `\\begin{tikzpicture}`, "
                     "`\\begin{figure}`, or `\\includegraphics` environments, so there is "
                     "nothing to reconstruct. The machine-readable data behind the paper's "
                     "tables is produced under `data/" + SHORT + "/` by the scripts in "
                     "`code/" + SHORT + "/`.\n")
        print(f"no tikzpicture found; wrote {readme}")
        return

    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    doc = (preamble
           + "\n\\usepackage[active,tightpage]{preview}\n"
           + "\\PreviewEnvironment{tikzpicture}\n\\begin{document}\n"
           + "\n\n".join(tikz) + "\n\\end{document}\n")
    open(os.path.join(OUT, "_figbuild.tex"), "w", encoding="utf-8").write(doc)
    for _ in range(2):
        subprocess.run([XELATEX, "-interaction=nonstopmode", "_figbuild.tex"],
                       cwd=OUT, capture_output=True, text=True, timeout=240)
    import fitz
    d = fitz.open(os.path.join(OUT, "_figbuild.pdf"))
    for i, page in enumerate(d):
        page.get_pixmap(dpi=200).save(os.path.join(OUT, f"figure{i+1}.png"))
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i+1}.pdf"))
        one.close()
    d.close()
    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)
    print(f"reconstructed {len(tikz)} figure(s) into {OUT}")


if __name__ == "__main__":
    main()
