# Errata register

Adjudications of errata reported against the corpus. Each entry records what
was alleged, what an audit of this repository found, and where the erratum
actually applies. Corrections whose targets exist in this lineage are fixed
in the papers directly (see `tests/README.md` for the six adjudicated-and-
corrected flags, and the git history for the fix commits); this file tracks
reports whose targets do **not** exist here.

An external review batch (2026-07) raised four items, D1 and E1–E3. D1 — a
stale "twenty-three cross-engine signatures" count in the relational-charge
Methodology replication list — was confirmed in-repo and corrected (the bullet
now reads "all twenty-six cross-engine items of ledger entry M (twenty-four
signatures, two factorization checks)", matching the paper's own Second-engine
paragraph and Appendix A entry M). E1–E3 were adjudicated as follows.

## E1 — degree-2 census 6→7; prediction P4 12→13; prediction P3 2661→2662

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

E1 alleged three wrong numerals in the pisot/relational lineage: (a) a
degree-2 census count printed 6 that should be 7; (b) a prediction "P4" count
12 that should be 13; (c) a prediction "P3" count 2661 that should be 2662.
Audit findings:

- Neither `2661` nor `2662` appears in any `papers/**/*.tex`, in the committed
  deposit builds (`pisot_residue_whitepaper` v13, `relational_charge_paper`
  rev. l), or in `tests/`.
- Both P3 and P4 in this lineage (`pisot_residue_whitepaper.tex` §"Predictions",
  `relational_charge_paper.tex` open-problems list) are count-free statements.
- The only degree-2 census count anywhere in the archive
  (`lehmers_box.tex`, subfield-signature census, the (2,2,0) column) already
  reads 7.

E1 therefore targets an external Zenodo deposit revision that was never
committed to this repository. Action: none in-repo; if the affected deposit
version is identified, file the erratum against that DOI version record.

## E2 — degree-7 Pisot polynomial display vs. its c-vector

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

No degree-7 polynomial display and no "c-vector" notation exists in any
committed `.tex`, PDF, test, or data file. Every `x^7` occurrence is an
interior term of a higher-degree polynomial (Lehmer's degree-10 polynomial;
the degree-12 reducible twist-class). Degree 7 appears only as prose future
work (pisot-residue item N4). If this erratum is real it targets an external
Zenodo deposit revision only.

## E3 — "both-straddle" count 163 vs. 160

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

The token "both-straddle" appears nowhere in the repository; every "straddle"
in this lineage is a boolean predicate (the Salem flip-straddle; K-formation
y-roots straddling zero), never a counted category. 163 occurs only as an RGB
color component in paper preambles; no 160/163 tally pair exists in any
`.tex`, committed PDF, test, or data file. If this erratum is real it targets
an external Zenodo deposit revision only.
