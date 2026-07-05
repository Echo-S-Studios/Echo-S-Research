# Contributing to Echo S Research

Every contribution — **any field**: math, engineering, physics theory, music
theory, bio, metacybernetics, or pure prose — follows one small, domain-neutral
contract: a **`contribution.yml` manifest** at the root of your
`papers/<folder>/`. Satisfy the manifest and both CI checks stay green. You do
**not** need to do math, use LaTeX, or write code.

**Start here:** [`docs/ANTI_DRIFT.md`](docs/ANTI_DRIFT.md) — a plain checklist with
worked examples for several domains. This file is the short version.

**Reference examples:**
- Computational (math): [`papers/2026-06-salem-slot/`](papers/2026-06-salem-slot)
  and its `code/ data/ tests/ figures/` folders + its
  [`contribution.yml`](papers/2026-06-salem-slot/contribution.yml).
- Pure prose: use `checks: none` with a `reason` (see ANTI_DRIFT.md, example 3).

## 0. Be in the registry

Your GitHub handle must be in [`.github/members.yml`](.github/members.yml) — the
single source of truth. New teammate? Add **one line** (or uncomment a reserved
slot). Nothing else needs changing; everything is registry-driven.

## 1. Scaffold your contribution

```bash
py scripts/new_contribution.py --member <your-handle> --shortname my-idea \
    --date 2026-08 --domain <your-field> [--trees code,data,tests,figures]
```

This validates your handle, then creates `papers/<folder>/` with a paper stub,
`paper.cff` (attribution pre-filled), and a **`contribution.yml`** manifest.
`--domain` is free text. Add `--trees` only for the output folders you actually
want (a prose-only piece needs none). The generator refuses unregistered handles
and existing folders.

## 2. Fill it in

- **Your paper** goes in `papers/<folder>/` — the `.tex` stub, or your own
  `.md`/`.pdf`. List every file you keep under `artifacts` in the manifest.
- **`contribution.yml`** is the contract. Set `title`/`domain`, and choose:
  - **`checks: none`** + a `reason` — for pure theory/prose (passes cleanly), or
  - a **`checks:`** list — one entry per thing that should keep reproducing:
    `run` regenerates an output, `produces` diffs it against the committed copy,
    `verify` runs a command that must exit 0 (a test, linter, schema, notebook…).
    `{PY}` = the Python interpreter (portable). See ANTI_DRIFT.md for the shapes.
- **Keep reproducible outputs deterministic** — seed RNG, round floats, sort
  before writing (ANTI_DRIFT.md §c). That's the whole trick to a green `drift`.
- **Citations** go in the one shared
  [`papers/references.bib`](papers/references.bib); attribution rules are in
  [`docs/ATTRIBUTION.md`](docs/ATTRIBUTION.md) (org credit stays org-only;
  per-paper credit lives in your `paper.cff` / manifest `member`).

## 3. Pre-flight, then open a PR

Run the **exact same checks CI runs**, locally, before pushing:

```bash
py scripts/check.py <your-folder>     # your contribution (fast)
py scripts/check.py                   # everything
```

Green here → green in CI. Then commit and open a PR. Two GitHub Actions run
(details in [`docs/MAINTAINING.md`](docs/MAINTAINING.md)):

- **structure** — manifest present & valid; your handle is registered; declared
  artifacts exist; you have exactly the `code/ data/ tests/ figures/` folders your
  manifest declares; naming is `YYYY-MM-shortname`; repo metadata is valid.
- **drift** — runs each contribution's declared `checks` and diffs `produces`
  against committed; `checks: none` contributions are skipped, not failed.

Two green badges at the top of the [README](README.md) mean the repo is consistent.
