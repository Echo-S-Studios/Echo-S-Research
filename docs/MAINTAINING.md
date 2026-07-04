# Maintaining Echo S Research

Two GitHub Actions keep the archive internally consistent. Both run on every
`push` and `pull_request`, and both are **registry-driven** (from
[`.github/members.yml`](../.github/members.yml)) and **convention-driven** (folder
globs) — onboarding a member or adding a paper needs **no workflow edit**.

The status badges at the top of the [README](../README.md) are the at-a-glance
signal: two green badges = consistent repo; a red badge = something went haywire.

## `structure` — [`.github/workflows/structure.yml`](../.github/workflows/structure.yml)

Validates repository invariants via [`scripts/check_structure.py`](../scripts/check_structure.py)
plus `cffconvert`. It reports **all** problems at once. It fails if:

| Check | Failure message prefix | Meaning |
|-------|------------------------|---------|
| Required files | `[required]` | a `papers/<folder>/` is missing `paper.cff` or a `.tex` |
| Attribution | `[attribution]` | a `paper.cff` (or `CITATION.cff`) author has no `alias`, or an `alias` that isn't a handle in `.github/members.yml` |
| Parallel-tree parity | `[parity]` | `tests/`, `code/`, `data/`, or `figures/` is missing a folder for some paper, or has an orphan folder with no matching paper |
| Naming | `[naming]` | a contribution folder isn't named `YYYY-MM-shortname` |
| Bibliography | `[bib]` | `papers/references.bib` doesn't parse |
| Metadata | `[zenodo]` / `[cff]` | `.zenodo.json` isn't valid JSON, or a `.cff` doesn't parse / isn't valid CFF 1.2.0 |

**How to read a failure:** open the failed run, read the `✗ [...]` lines — each
names the exact file and the fix. Reproduce locally with:

```bash
py scripts/check_structure.py
cffconvert --validate -i CITATION.cff
```

Most common causes: forgot to add your handle to `members.yml`; created only some
of the five folders (use the generator, which makes all of them); a typo in a
folder name.

## `drift` — [`.github/workflows/drift.yml`](../.github/workflows/drift.yml)

Guards the "code produces data" contract. It snapshots the committed `data/`,
re-runs **every producer** under `code/` via
[`scripts/regen_all.py`](../scripts/regen_all.py), and diffs the regenerated
`data/` against the snapshot (line-ending-agnostic). If anything differs it fails
with `Data drift detected` and prints the diff. It then runs the full `pytest`
suite across `tests/`.

Figure builders (`make_figures.py`) and helper libraries (modules without a
`__main__` guard) are skipped — the check is about **data**, and figures need TeX.

**A red drift badge means one of two things:**

1. **Stale data.** You changed a producer but didn't commit the regenerated data.
   Fix: `py scripts/regen_all.py`, then commit the updated `data/`.
2. **Nondeterministic producer.** A producer emits different output each run
   (unseeded RNG, timestamps, unordered `set` iteration). Fix: make it
   deterministic — seed RNG (`random.seed`, `np.random.seed`), sort before
   emitting, don't write wall-clock time.

Reproduce locally:

```bash
py scripts/regen_all.py
git diff -- data          # any output here = drift
py -m pytest tests -q
```

## When you change things

- **New teammate:** one-line edit to [`.github/members.yml`](../.github/members.yml).
  Both workflows pick it up automatically.
- **New paper:** `py scripts/new_contribution.py --member <handle> --shortname <name> --date YYYY-MM`
  creates all five folders in the right shape; fill them in and both workflows validate them.
- **Never** edit a workflow to accommodate a new member or paper — if you feel the
  need to, the registry or the folder convention is the thing to change instead.
