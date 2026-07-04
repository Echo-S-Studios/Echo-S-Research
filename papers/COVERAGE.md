# Coverage report — papers/ intake

Generated 2026-07-03 while populating `papers/` from local research work on this
machine. **Originals were only ever copied, never moved or deleted.**

_Updated 2026-07-03: the five originally source-less PDFs received their LaTeX
sources (uploaded to `Downloads\Zenodo\`); all 14 papers are now matched
`.tex` + `.pdf` pairs, and `references.bib` grew from 42 to 51 entries._

## Search scope

Recursively searched `C:\Users\acead\projects`, `C:\Users\acead\Downloads`, and
`C:\Users\acead\Documents` for `.tex` and `.pdf` files, skipping
build/venv/node_modules/.git and obvious non-papers. No `.bib` or `.bbl` files
exist anywhere in the source material.

## Counts

| Metric | Count |
|--------|-------|
| `.tex` files found | 18 (14 distinct papers + 4 duplicate/older-version copies) |
| `.pdf` files found | 23 |
| Distinct papers archived | 14 |
| — with matched `.tex` + `.pdf` | 14 |
| — PDF only (no LaTeX source) | 0 |
| `.tex` with no compiled PDF | 0 |
| PDFs excluded (non-paper / out of scope) | 4 |
| Superseded/duplicate PDFs not copied | 5 |
| Consolidated references in `references.bib` | 51 |

---

## Matched `.tex` + `.pdf` pairs (14)

Each was placed in its `papers/<folder>/` with both source and compiled PDF.

| Folder | `.tex` original | `.pdf` original |
|--------|-----------------|-----------------|
| `2026-06-charge-measure-coupling` | `Downloads\Foundational Pipes\charge-measure-coupling-whitepaper-v4.tex` | `Downloads\Foundational Pipes\papers\charge-measure-coupling-whitepaper-v4.pdf` |
| `2026-06-emission-gap` | `Downloads\Foundational Pipes\lambda2c-emissiongap-verification\lambda2c-emissiongap-verification\papers\emission_gap_paper.tex` | `Downloads\Foundational Pipes\papers\emission_gap_paper.pdf` |
| `2026-06-lambda-2c` | `Downloads\Foundational Pipes\lambda2c-emissiongap-verification\lambda2c-emissiongap-verification\papers\lambda_2c_paper.tex` | `Downloads\Foundational Pipes\papers\lambda_2c_paper.pdf` |
| `2026-06-lehmers-box` | `Downloads\Foundational Pipes\lehmers_box.tex` | `Downloads\Foundational Pipes\papers\lehmers box.pdf` |
| `2026-06-z5-no-salem-dichotomy` | `Downloads\Foundational Pipes\Z5-no-salem-dichotomy-whitepaper.tex` | `Downloads\Foundational Pipes\papers\Z5-no-salem-dichotomy-whitepaper.pdf` |
| `2026-06-residual-return-learning` | `Downloads\Foundational Pipes\residual-return-verification-v2\residual-return-verification\L00M\paper\residual_return_learning.tex` | (same folder) `residual_return_learning.pdf` |
| `2026-06-vector-substrate` | `Downloads\Foundational Pipes\residual-return-verification-v2\residual-return-verification\L00M\paper\vector_substrate.tex` | (same folder) `vector_substrate.pdf` |
| `2026-07-pisot-residue` | `Downloads\Zenodo\pisot_residue_whitepaper.tex` | `Downloads\Zenodo\pisot_residue_whitepaper.v13.pdf` |
| `2026-07-relational-charge` | `Downloads\Zenodo\relational_charge_paper.tex` | `Downloads\Zenodo\relational_charge_paper.l.pdf` |
| `2026-06-generative-emptiness` | `Downloads\Zenodo\generative_emptiness.tex` † | `Downloads\Foundational Pipes\papers\generative_emptiness.pdf` |
| `2026-06-operator-algebra` | `Downloads\Zenodo\operator-algebra-whitepaper.tex` † | `Downloads\Foundational Pipes\papers\operator-algebra-whitepaper.pdf` |
| `2026-06-salem-slot` | `Downloads\Zenodo\salem_slot.tex` † | `Downloads\Foundational Pipes\papers\salem_slot.pdf` |
| `2026-07-emission-algebra-primer` | `Downloads\Zenodo\emission_algebra_primer.tex` † | `Downloads\Foundational Pipes\papers\emission_algebra_primer.pdf` |
| `2026-07-helix-orthogonal-partner` | `Downloads\Zenodo\helix_orthogonal_partner.tex` † | `Downloads\Foundational Pipes\papers\helix_orthogonal_partner.pdf` |

† LaTeX source uploaded 2026-07-03, after the initial intake; the PDF was already archived as source-less. Titles in the `.tex` exactly match the archived PDFs.

Titles (from `\title{}`):
- *The Charge–Measure Coupling on a Spectral Semiring*
- *The Emission–Gap Theorem*
- *The Exchange Rate λ = 2c: A Conformal Identity, Its Gate, and Its Flip*
- *Lehmer's Box*
- *The Z₅ Case of the No–Salem Dichotomy*
- *Residual Return: Exact Learning Dynamics and Language over the Vector Substrate*
- *The Vector Substrate: Number Fields as Exact Learning Geometry*
- *The Pisot Cross-Shell Residue*
- *Relational Charge on the Spectral Semiring*
- *The Generative Content of a Conserved Emptiness*
- *The Operator Algebra of the Emission Semiring*
- *The Occupant of the Salem Slot*
- *The Emission Algebra 𝔄*
- *The Dissolved Helix and Its Orthogonal Partner*

---

## Missing LaTeX — PDFs archived with NO source (0)

**None — resolved.** All five originally source-less PDFs (emission-algebra-primer,
generative-emptiness, helix-orthogonal-partner, operator-algebra, salem-slot)
received their LaTeX sources on 2026-07-03. The sources were copied into the
matching folders and the papers moved to the matched-pairs table above.

## `.tex` with no compiled PDF

**None.** All 14 LaTeX sources have a corresponding compiled PDF.

---

## Ambiguous / judgment calls

1. **Version selection for duplicated papers.** `residual_return_learning` and
   `vector_substrate` each existed twice: an older copy in
   `projects\L00M\paper\` (2026-06-22/24) and a newer, larger copy in the
   `residual-return-verification-v2` bundle (2026-06-28). **Chose the newer v2
   copies.** The `projects\L00M` copies were left in place, not archived.
2. **`relational_charge_paper` PDF version.** Four PDFs existed:
   `relational_charge_paper.pdf` (07-01, in two locations), `.k.pdf` (07-02),
   and `.l.pdf` (07-02, largest). **Chose `.l.pdf`** as the newest revision
   (alphabetical suffix + size). Guess — if `.k` or the base is actually the
   intended final, swap it in.
3. **`pisot_residue_whitepaper`** had one PDF, `…​.v13.pdf` — treated `v13` as the
   compiled version of the `.tex`.
4. **`lehmers box.pdf`** (with a space) was matched to `lehmers_box.tex` (with an
   underscore) by title/content — near-certain, but the names differ.
5. **`cognitive-compiler-formalized.pdf` — dropped on request.** It was initially
   archived as a genuine authored paper (abstract, appendices, exact-arithmetic
   verification record, by J. Turner & C. Snodgrass), but it is systems /
   AI-orchestration research from the `Clif\` DARPA-proposal working folder rather
   than the pure-math spectral-semiring line, so it was **removed to keep the
   archive math-focused.** See the excluded-PDFs table below; the original in
   `Downloads\Clif\` is untouched if it should be added back.
6. **`.tex` and `.pdf` in different folders.** For several papers the source and
   the compiled PDF lived in different directories (e.g. charge-measure `.tex` in
   the `Foundational Pipes` root but its PDF in `Foundational Pipes\papers\`;
   relational-charge `.tex` in `Zenodo\` with a duplicate PDF also in
   `Foundational Pipes\papers\`). Matched by filename stem + title.
7. **Folder dates.** `YYYY-MM` uses the paper's `\date{}` where present, else the
   file's modification month. Several `.tex` had an empty `\date{}`; those used
   the file date (all June 2026).
8. **Late-arriving sources (2026-07-03).** The five LaTeX sources uploaded to
   `Downloads\Zenodo\` were matched to already-archived PDFs by identical filename
   stem and exact title match. The same upload also included an older
   `charge-measure-coupling-whitepaper-v2.tex` and a re-uploaded
   `Z5-no-salem-dichotomy-whitepaper.tex`; **both were ignored** because those
   papers already have their canonical sources archived (v4 and the Foundational
   Pipes copy, respectively).

---

## Excluded PDFs (found but not archived)

| Original path | Reason |
|---------------|--------|
| `Downloads\DMV Drivers License _ State ID _ CDL Renewal.pdf` | Personal/administrative document — not research. |
| `Downloads\Clif\DARPA_26BZ_R3.pdf` | External government document: *DARPA DoW 2026 SBIR Broad Agency Announcement* (40 pp solicitation), not the team's own work. |
| `Downloads\Clif\KT_MAESTRO_Jun_30_2026_1221_sonnet_4_6.pdf` | 2-page AI chat export (Claude Sonnet 4.6) summarizing *Lehmer's Box*; a derivative session artifact, not a paper. |
| `Downloads\Clif\cognitive-compiler-formalized.pdf` | Genuine authored research paper (13 pp), but systems / AI-orchestration rather than the pure-math line; dropped on request to keep the archive math-focused. |

## Superseded / duplicate PDFs (not copied)

- `projects\L00M\paper\residual_return_learning.pdf` — older than the v2-bundle copy archived.
- `projects\L00M\paper\vector_substrate.pdf` — older than the v2-bundle copy archived.
- `Downloads\Zenodo\relational_charge_paper.pdf` (07-01) — superseded by `.l.pdf`.
- `Downloads\Zenodo\relational_charge_paper.k.pdf` (07-02) — superseded by `.l.pdf`.
- `Downloads\Foundational Pipes\papers\relational_charge_paper.pdf` (07-01) — duplicate of the Zenodo base version.

## Related assets seen but out of scope (not `.tex`/`.pdf`)

Companion non-paper files in `Foundational Pipes\papers\` and `Clif\` were left
alone: `emission-algebra-compendium.html`, `lehmers_box_instrument.html`,
`lehmers-problem-an-introduction.md`, `Cognitive_Compiler_Rewrite_Draft.md`.

**Archives not extracted** (may contain further sources — not searched):
`Downloads\Zenodo\files.zip`, `files-fixed.zip`;
`Downloads\Foundational Pipes\Zipped-Tarred-Pipelines\*.tar`, `Bulk\*.zip`;
`C:\Users\acead\projects\L00M.tar`;
`C:\Users\acead\residual-return-verification.tar.gz`.

---

## `references.bib`

51 unique entries, harvested and de-duplicated from the inline
`\begin{thebibliography}` environments of the `.tex` sources (no `.bib`/`.bbl`
existed). 42 came from the initial intake; 9 more were added from
`operator-algebra-whitepaper.tex` when its source arrived on 2026-07-03. Keys are
preserved verbatim so existing `\cite{}` calls resolve — note that operator-algebra
uses lowercase keys, so `lehmer` and `smyth` are kept as separate entries
duplicating `Lehmer` and `Smyth` (same works, different `\cite` keys).

Entries flagged as incomplete or internal: `CMC`, `relcharge`, `EG`,
`VectorSubstrate`, `corpus`, `companions`, `repo` (shared/missing DOIs or
unpublished companion material), plus `Rissanen` and `Yuzvinskii` (source bundled
a second work into one bibitem, preserved in a `note`) and `Bertin` (source author
list was "et al."). Validated: parses with 51/51 unique keys.
