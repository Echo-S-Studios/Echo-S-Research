# figures/__FOLDER__/

Figures for **papers/__FOLDER__/**.

If your paper has inline TikZ/pgfplots, add a `code/__FOLDER__/make_figures.py`
that reconstructs each one (reuse the paper's own preamble + the `preview`
package + XeLaTeX, then render to PNG/PDF with `pymupdf` — copy the reference
example [`code/2026-06-salem-slot/make_figures.py`](../../code/2026-06-salem-slot)),
and map each `figureN.*` here to its paper Figure number.

If your paper has no figures, keep this note saying so.
