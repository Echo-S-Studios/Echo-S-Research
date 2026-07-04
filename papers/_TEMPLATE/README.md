# _TEMPLATE — contribution skeleton

This folder is the **cross-folder skeleton** that the generator instantiates. Its
files contain `__TOKEN__` placeholders (`__DISPLAY_NAME__`, `__HANDLE__`,
`__SHORTNAME__`, `__FOLDER__`, `__DATE__`, `__TITLE__`) that
[`scripts/new_contribution.py`](../../scripts/new_contribution.py) substitutes.
You normally do **not** copy this folder by hand — run the generator:

```bash
py scripts/new_contribution.py --member <your-github-handle> --shortname my-result --date 2026-08
```

It validates your handle against [`.github/members.yml`](../../.github/members.yml)
and creates the **five parallel folders**, matching the layout of every existing
paper (reference example: [`papers/2026-06-salem-slot/`](../2026-06-salem-slot)):

```
papers/YYYY-MM-shortname/    paper.cff (attribution pre-filled), <shortname>.tex
tests/YYYY-MM-shortname/     test_example.py, NOTES.md
code/YYYY-MM-shortname/      producer.py, README.md
data/YYYY-MM-shortname/      README.md   (producers write outputs here)
figures/YYYY-MM-shortname/   README.md
```

Skeleton → destination mapping:

| Template file | Generated as |
|---------------|--------------|
| `paper.cff` | `papers/<folder>/paper.cff` |
| `paper.tex` | `papers/<folder>/<shortname>.tex` |
| `tests/test_example.py`, `tests/NOTES.md` | `tests/<folder>/` |
| `code/producer.py`, `code/README.md` | `code/<folder>/` |
| `data/README.md` | `data/<folder>/` |
| `figures/README.md` | `figures/<folder>/` |

Full workflow: [`../../CONTRIBUTING.md`](../../CONTRIBUTING.md). This `_TEMPLATE`
folder is skipped by the CI structure checks (it is a skeleton, not a paper).
