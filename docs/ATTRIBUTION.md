# Attribution

How authorship is recorded in this archive. **Read this before adding a paper.**

## Team members

Members are attributed by GitHub username for now (ORCIDs may be added later).

| GitHub username | Display name | Author string (use in papers) |
|-----------------|--------------|-------------------------------|
| [`AceTheDactyl`](https://github.com/AceTheDactyl) | AceTheDactyl | `AceTheDactyl (@AceTheDactyl)` |
| [`eMKa7-tY`](https://github.com/eMKa7-tY) | mK@ | `mK@ (@eMKa7-tY)` |
| [`sKiDaGgAbAtEe`](https://github.com/sKiDaGgAbAtEe) | sKiDaGgAbAtEe | `sKiDaGgAbAtEe (@sKiDaGgAbAtEe)` |

The umbrella organization credit is always **Echo S Studios Research Developments**.

> **The machine-readable registry is [`.github/members.yml`](../.github/members.yml) — the single source of truth.** The table above is the human-readable view of it. The `alias:` (GitHub handle) you put in your `paper.cff` is the self-identification, and it is validated against that registry by the `structure` CI workflow; [`scripts/new_contribution.py`](../scripts/new_contribution.py) also refuses to scaffold for a handle that isn't listed. Onboarding a new teammate is a **one-line edit** to `members.yml` — nothing else changes.

## The umbrella rule

Two layers of credit, and they do not conflict:

- **Repository / release level — the organization.** The umbrella creator for the
  whole archive is **Echo S Studios Research Developments**. It is the single
  `creators` entry in the repo-level [`.zenodo.json`](../.zenodo.json), so every
  Zenodo deposition (minted per GitHub Release) is credited to the org. **Do not
  add individual people to `.zenodo.json`.**
- **Paper level — the individual.** Each paper records its own author(s) in that
  paper's `paper.cff` and in the LaTeX `\author{}`. This is where individual
  members get credit.

The repo-level [`CITATION.cff`](../CITATION.cff) lists the org first (umbrella)
followed by all members, giving a single "cite the whole archive" entry. Because
`.zenodo.json` is present it takes precedence for the Zenodo deposition, so listing
members in `CITATION.cff` does **not** change the org-only Zenodo credit — the two
files are consistent, not conflicting.

## How to attribute your paper

When you add a paper, follow these four steps.

**a. Scaffold the folder** — run the generator. It validates your handle against
[`.github/members.yml`](../.github/members.yml) and creates `papers/YYYY-MM-shortname/`
with a paper stub, `paper.cff` (attribution pre-filled), and a `contribution.yml`
manifest. Add `--trees code,data,tests,figures` for any compute folders you want:

```bash
py scripts/new_contribution.py --member <your-github-handle> --shortname my-new-result --date 2026-08 --domain <your-field>
```

Attribution is carried by the manifest's `member` field (validated against the
registry) and, if present, your `paper.cff` `alias`. See
[`CONTRIBUTING.md`](../CONTRIBUTING.md) and [`ANTI_DRIFT.md`](ANTI_DRIFT.md) for the
end-to-end workflow, and copy [`papers/2026-06-salem-slot/`](../papers/2026-06-salem-slot)
as the reference example.

**b. Set the LaTeX author** — in your `.tex`, set `\author{}` to your author string
from the table above (org affiliation on the second line):

```latex
\author{mK@ (@eMKa7-tY)\\ Echo S Studios Research Developments}
```

**c. Fill `paper.cff`** — put your **display name** in `name:` and your **GitHub
username** in `alias:`, and set `title:` and `date-released:`:

```yaml
authors:
  - name: "mK@"
    alias: "eMKa7-tY"
```

For a co-authored paper, add one `- name: / alias:` block per author.

**d. Add citations** — put any new references into the shared
[`papers/references.bib`](../papers/references.bib) — the single consolidated
bibliography for the whole archive. Keep BibTeX keys stable so existing `\cite{}`
calls keep resolving. Do **not** start a separate per-paper `.bib`.

Then commit and open a PR. Keep your `paper.cff` valid CFF 1.2.0 — check with:

```bash
cffconvert --validate -i papers/2026-08-my-new-result/paper.cff
```

## Required files per paper

Each `papers/YYYY-MM-shortname/` folder must contain:

- the LaTeX source (`.tex`) — with `\author{}` set to your author string,
- the compiled paper (`.pdf`),
- `paper.cff` — per-paper citation metadata with the author block.

See [`papers/_TEMPLATE/`](../papers/_TEMPLATE/) for a ready-to-copy starter.
