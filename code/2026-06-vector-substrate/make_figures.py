"""
PRODUCER: rebuild the paper's four inline tikz/pgfplots figures, and regenerate
the two data-driven figures from first principles with matplotlib.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex
  Figure 1 (fig:triptych)  -- three views of K glued by the trace-form Gram
  Figure 2 (fig:proj)      -- projector/residual in the Minkowski plane of Q(sqrt5)
  Figure 3 (fig:mahler)    -- Mahler picture for the plastic number (companion spectrum)
  Figure 4 (fig:threshold) -- growth decision in the (cost, gain) plane (pgfplots)

Reuses the paper's own preamble so every macro/color resolves.  Figures 3 and 4
are additionally regenerated from computed data as figure3_regen.png /
figure4_regen.png.

Run: py code/2026-06-vector-substrate/make_figures.py
"""
import glob
import os
import re
import subprocess

import fitz  # pymupdf

XELATEX = r"C:\Users\acead\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
SHORT = "2026-06-vector-substrate"
ROOT = r"C:\Users\acead\projects\Echo-S-Research"
TEX = glob.glob(os.path.join(ROOT, "papers", SHORT, "*.tex"))[0]
OUT = os.path.join(ROOT, "figures", SHORT)
os.makedirs(OUT, exist_ok=True)


def build_tikz_figures():
    """Extract every tikzpicture and typeset each on its own cropped page."""
    src = open(TEX, encoding="utf-8").read()
    preamble = re.search(r"(.*?)\\begin\{document\}", src, re.DOTALL).group(1)
    tikz = re.findall(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", src, re.DOTALL)
    if not tikz:
        print("no tikzpicture environments found")
        return 0
    doc = (
        preamble
        + "\n\\usepackage[active,tightpage]{preview}\n"
        + "\\PreviewEnvironment{tikzpicture}\n"
        + "\\begin{document}\n"
        + "\n\n".join(tikz)
        + "\n\\end{document}\n"
    )
    open(os.path.join(OUT, "_figbuild.tex"), "w", encoding="utf-8").write(doc)
    for _ in range(2):
        subprocess.run(
            [XELATEX, "-interaction=nonstopmode", "_figbuild.tex"],
            cwd=OUT, capture_output=True, text=True, timeout=240,
        )
    pdf_path = os.path.join(OUT, "_figbuild.pdf")
    if not os.path.exists(pdf_path):
        log = os.path.join(OUT, "_figbuild.log")
        tail = ""
        if os.path.exists(log):
            tail = "".join(open(log, encoding="utf-8", errors="ignore").readlines()[-40:])
        print("XELATEX FAILED; log tail:\n" + tail)
        return 0
    d = fitz.open(pdf_path)
    made = 0
    for i, page in enumerate(d):
        page.get_pixmap(dpi=200).save(os.path.join(OUT, f"figure{i + 1}.png"))
        one = fitz.open()
        one.insert_pdf(d, from_page=i, to_page=i)
        one.save(os.path.join(OUT, f"figure{i + 1}.pdf"))
        one.close()
        made += 1
    d.close()
    for f in glob.glob(os.path.join(OUT, "_figbuild.*")):
        os.remove(f)
    print(f"built {made} tikz figures -> figure1..{made}")
    return made


def regen_mahler():
    """Figure 3 regen: eigenvalues of the companion of x^3-x-1 (the plastic
    number) with the unit circle; Mahler = product of |eigenvalues| > 1 = mu_S."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    roots = np.roots([1, 0, -1, -1])           # x^3 - x - 1
    outside = [r for r in roots if abs(r) > 1]
    mahler = float(np.prod([abs(r) for r in outside]))

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="black", lw=1.2)
    ax.fill(np.cos(th), np.sin(th), color="0.92", zorder=0)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.axvline(0, color="0.6", lw=0.6)
    for r in roots:
        inside = abs(r) <= 1
        ax.plot(r.real, r.imag, "o", color=("0.15" if inside else "crimson"), ms=7)
    ax.annotate(r"$\lambda_1=\mu_S=%.4f$" % mahler, (outside[0].real, 0),
                textcoords="offset points", xytext=(8, 10))
    ax.set_xlabel(r"$\Re$")
    ax.set_ylabel(r"$\Im$")
    ax.set_title(r"Mahler picture: plastic number $x^3-x-1$"
                 "\n" r"$\mathrm{M}=\prod_{|\lambda|>1}|\lambda|=%.6f$" % mahler)
    ax.set_aspect("equal")
    ax.set_xlim(-1.6, 2.0)
    ax.set_ylim(-1.3, 1.4)
    fig.tight_layout()
    path = os.path.join(OUT, "figure3_regen.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"regenerated {path}  (mu_S = {mahler:.6f})")


def regen_threshold():
    """Figure 4 regen: growth decision in the (cost, gain) plane (log-log), from
    the produced threshold data: 2sqrt6 (96 vs 6.356), sqrt7 (56 vs 3.892), the
    tiny STOP (0.1), the gain=cost boundary, and the two floors."""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data_csv = os.path.join(ROOT, "data", SHORT, "threshold_decisions.csv")
    pts = []
    if os.path.exists(data_csv):
        with open(data_csv, encoding="utf-8") as f:
            rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
        header = rows[0]
        for r in rows[1:]:
            rec = dict(zip(header, r))
            try:
                gain = float(eval(rec["gain"])) if rec["gain"] not in ("-", "") else None
                cost = float(rec["cost_lambda_logM"]) if rec["cost_lambda_logM"] not in ("-", "") else None
            except Exception:
                gain, cost = None, None
            pts.append((rec["seed"], cost, gain, rec["decision"]))

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    xs = np.logspace(np.log10(0.35), np.log10(16), 50)
    ax.plot(xs, xs, "k--", lw=1, label=r"gain = cost")
    ax.axhline(0.562, color="tab:blue", lw=1.2, label="constant floor 0.56")
    ax.axhline(4.5, color="tab:red", ls=":", lw=1.4, label="deg. floor n=8 (4.5)")
    for name, cost, gain, decision in pts:
        if gain is None:
            continue
        if cost is None:                       # tiny STOP: illustrative x
            cost = 6.356
        marker = "o" if decision == "GROW" else "x"
        color = "black" if decision == "GROW" else "red"
        ax.plot(cost, gain, marker, color=color, ms=8, mew=2)
        ax.annotate(f"{name} ({decision})", (cost, gain),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.35, 16)
    ax.set_ylim(0.04, 240)
    ax.set_xlabel(r"cost $= \lambda \log \mathrm{M}(\theta)$")
    ax.set_ylabel(r"gain $= \|r\|_G^2$")
    ax.set_title("Growth decision in the (cost, gain) plane")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, which="both", color="0.9")
    fig.tight_layout()
    path = os.path.join(OUT, "figure4_regen.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"regenerated {path}")


def main():
    n = build_tikz_figures()
    regen_mahler()
    regen_threshold()
    print(f"done: {n} tikz figures + 2 matplotlib regenerations")


if __name__ == "__main__":
    main()
