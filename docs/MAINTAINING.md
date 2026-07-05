# Maintaining Echo S Research

Two GitHub Actions keep the archive internally consistent across **every domain**.
Both run on every `push` and `pull_request`. Everything is driven by the
per-contribution manifest (`papers/<sn>/contribution.yml`), the member registry
([`.github/members.yml`](../.github/members.yml)), and folder conventions — so
onboarding a member or adding a contribution needs **no workflow edit**.

The status badges at the top of the [README](../README.md) are the at-a-glance
signal: two green badges = consistent repo; a red badge = something drifted.
Contributors keep themselves green with [`docs/ANTI_DRIFT.md`](ANTI_DRIFT.md) and
the `scripts/check.py` pre-flight; this page is what the maintainers see.

## The contract

Each `papers/<sn>/contribution.yml` declares `shortname`, `member`, `title`,
`domain` (free text), `artifacts` (files/dirs it owns), and `checks` — either a
list of named checks (`run` / `produces` / `verify`) or the literal `none` with a
`reason`. The schema lives in one place:
[`scripts/_manifest.py`](../scripts/_manifest.py). Both workflows are thin wrappers
over that contract, so the *contract* is what you maintain, not the YAML.

## `structure` — [`.github/workflows/structure.yml`](../.github/workflows/structure.yml)

Runs [`scripts/check_structure.py`](../scripts/check_structure.py) (reports **all**
problems at once) plus `cffconvert`. Failure prefixes:

| Prefix | Meaning |
|--------|---------|
| `[manifest]` | `contribution.yml` missing or doesn't parse |
| `[contract]` | a manifest field is wrong: unknown `member`, missing `title`/`domain`, a declared artifact doesn't exist, a bad `checks` block, or a tree folder present/absent that disagrees with what the manifest declares |
| `[orphan]` | a `tests/ code/ data/ figures/` folder exists with no matching `papers/<sn>/` (no manifest) |
| `[naming]` | a folder isn't `YYYY-MM-shortname` |
| `[attribution]` | a `CITATION.cff`/`paper.cff` alias isn't a registered handle |
| `[bib]` / `[zenodo]` / `[cff]` | `references.bib` won't parse, `.zenodo.json` isn't valid JSON, or a `.cff` isn't valid CFF 1.2.0 |

Note it's **domain-neutral**: a prose contribution with `checks: none` and no
`code/` folder passes. Reproduce locally: `py scripts/check_structure.py`.

## `drift` — [`.github/workflows/drift.yml`](../.github/workflows/drift.yml)

Runs [`scripts/run_drift.py`](../scripts/run_drift.py): for each contribution it
executes the manifest's declared `checks` — runs each `run`, diffs every
`produces` against the committed copy (line-ending-agnostic), and runs each
`verify` (must exit 0). `checks: none` contributions are **skipped and logged**
(`skip <sn>: no-repro (reason)`), never failed.

The base image pins Python 3.12 + `sympy`/`mpmath`/`numpy`/`pytest` to the versions
that generated the committed math data (so a library upgrade can't masquerade as
drift). Contributions needing other tools install them inside their own check
commands.

**A red drift badge means, for the named contribution:**

1. **Stale output** — a producer changed but its committed `produces` file wasn't
   regenerated. Fix: rerun the producer and commit the output.
2. **Nondeterministic producer** — unseeded RNG, timestamps, unordered `set`,
   signed-zero (`-0.0` vs `0.0`), or an unpinned dependency. Fix per
   [ANTI_DRIFT.md §c](ANTI_DRIFT.md).
3. **A failing `verify`** — that contribution's test/lint/schema command returned
   non-zero.

The failure line names the contribution and check. Reproduce locally:
`py scripts/run_drift.py <sn>` (or the whole pre-flight, `py scripts/check.py`).

## When you change things

- **New teammate:** one-line edit to [`.github/members.yml`](../.github/members.yml).
- **New contribution (any domain):**
  `py scripts/new_contribution.py --member <handle> --shortname <name> --date YYYY-MM --domain <field>`
  (add `--trees …` for compute folders). Fill in the manifest; both workflows validate it.
- **Never** edit a workflow to accommodate a member or contribution — change the
  registry, the manifest, or the folder convention instead.
