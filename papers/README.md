# Papers

One subfolder per paper. Name each subfolder `YYYY-MM-shortname`, using the year and month of the draft plus a short, hyphenated slug — for example:

```
papers/
  2026-07-spectral-gap/
    spectral-gap.pdf     # compiled paper
    spectral-gap.tex     # LaTeX source
    references.bib       # bibliography
```

Each paper folder should contain, at minimum:

- the compiled **PDF**, and
- the complete **LaTeX source** (`.tex`, `.bib`, and any custom class/style files) needed to rebuild it.

Shared artifacts that support the paper — larger datasets, reusable figures, and code — live in the top-level [`data/`](../data), [`figures/`](../figures), and [`code/`](../code) folders. Reference them from the paper as needed.

Intermediate LaTeX build files (`.aux`, `.log`, `.out`, …) are ignored by [`.gitignore`](../.gitignore); commit only the sources and the final PDF.
