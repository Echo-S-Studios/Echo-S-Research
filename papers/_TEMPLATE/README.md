# _TEMPLATE — starter for a new paper

Copy this whole folder to start a new paper. **Do not edit the template in place.**

```bash
cp -r papers/_TEMPLATE papers/2026-08-my-shortname             # macOS / Linux
```
```powershell
Copy-Item -Recurse papers\_TEMPLATE papers\2026-08-my-shortname   # Windows
```

## Folder naming: `YYYY-MM-shortname`

- `YYYY-MM` — year and month of the paper (its `\date`, or when it was drafted).
- `shortname` — a short, hyphenated slug, e.g. `spectral-gap`.

Example: `papers/2026-08-spectral-gap/`.

## Required files

| File | What to do |
|------|------------|
| `<shortname>.tex` | LaTeX source. Set `\author{}` to your author string — see [../../docs/ATTRIBUTION.md](../../docs/ATTRIBUTION.md). |
| `<shortname>.pdf` | The compiled paper. |
| `paper.cff` | Fill in `title`, `authors` (`name` = display name, `alias` = GitHub username), and `date-released`. |

Citations go in the shared [`../references.bib`](../references.bib), **not** here.

Full workflow: [../../docs/ATTRIBUTION.md](../../docs/ATTRIBUTION.md).
