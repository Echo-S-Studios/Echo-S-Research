# Papers

One subfolder per paper, named `YYYY-MM-shortname` (year and month of the paper
plus a short, hyphenated slug) — for example:

```
papers/
  2026-07-spectral-gap/
    spectral-gap.tex     # LaTeX source (\author set to your attribution)
    spectral-gap.pdf     # compiled paper
    paper.cff            # per-paper citation metadata (author block)
```

Each paper folder must contain:

- the complete **LaTeX source** (`.tex` and any custom class/style files),
- the compiled **PDF**, and
- **`paper.cff`** — per-paper citation metadata recording the author(s).

Copy [`_TEMPLATE/`](_TEMPLATE) to start a new paper. Authorship is per-member —
see **[`../docs/ATTRIBUTION.md`](../docs/ATTRIBUTION.md)** for the members table and
the step-by-step workflow.

Citations for **all** papers live in one consolidated bibliography,
[`references.bib`](references.bib) — do not keep a separate `.bib` per paper.

Shared artifacts that support a paper — larger datasets, reusable figures, and
code — live in the top-level [`data/`](../data), [`figures/`](../figures), and
[`code/`](../code) folders. Reference them from the paper as needed.

Intermediate LaTeX build files (`.aux`, `.log`, `.out`, …) are ignored by
[`.gitignore`](../.gitignore); commit only the sources, the final PDF, and
`paper.cff`.
