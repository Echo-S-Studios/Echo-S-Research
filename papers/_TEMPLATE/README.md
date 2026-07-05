# _TEMPLATE — contribution skeleton

This folder is the skeleton the generator instantiates. Its files contain
`__TOKEN__` placeholders (`__DISPLAY_NAME__`, `__HANDLE__`, `__SHORTNAME__`,
`__FOLDER__`, `__DATE__`, `__TITLE__`, `__DOMAIN__`) that
[`scripts/new_contribution.py`](../../scripts/new_contribution.py) substitutes.
You normally do **not** copy this by hand — run the generator:

```bash
py scripts/new_contribution.py --member <your-handle> --shortname my-result \
    --date 2026-08 --domain <your-field> [--trees code,data,tests,figures]
```

It validates your handle against [`.github/members.yml`](../../.github/members.yml)
and always creates `papers/<folder>/` with a paper stub, `paper.cff`, and a
**`contribution.yml`** manifest (the contract — see
[`../../docs/ANTI_DRIFT.md`](../../docs/ANTI_DRIFT.md)). It scaffolds the optional
`code/ data/ tests/ figures/` folders **only** for the `--trees` you ask for, so a
prose-only contribution stays minimal.

Skeleton → destination mapping:

| Template file | Generated as | When |
|---------------|--------------|------|
| `paper.cff` | `papers/<folder>/paper.cff` | always |
| `paper.tex` | `papers/<folder>/<shortname>.tex` | always |
| `contribution.yml` | `papers/<folder>/contribution.yml` | always (built from your flags; `checks: none` default) |
| `code/producer.py`, `code/README.md` | `code/<folder>/` | `--trees code` |
| `data/README.md` | `data/<folder>/` | `--trees data` |
| `tests/test_example.py`, `tests/NOTES.md` | `tests/<folder>/` | `--trees tests` |
| `figures/README.md` | `figures/<folder>/` | `--trees figures` |

Full workflow: [`../../CONTRIBUTING.md`](../../CONTRIBUTING.md). Reference example:
[`../2026-06-salem-slot/`](../2026-06-salem-slot). This `_TEMPLATE` folder is
skipped by the CI checks (it is a skeleton, not a contribution).
