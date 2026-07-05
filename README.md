# Echo S Research

[![structure](https://github.com/Echo-S-Studios/Echo-S-Research/actions/workflows/structure.yml/badge.svg)](https://github.com/Echo-S-Studios/Echo-S-Research/actions/workflows/structure.yml)
[![drift](https://github.com/Echo-S-Studios/Echo-S-Research/actions/workflows/drift.yml/badge.svg)](https://github.com/Echo-S-Studios/Echo-S-Research/actions/workflows/drift.yml)

Research papers and their computational products from **Echo S Studios Research Developments**.

> A red **structure** or **drift** badge means the repo went inconsistent — see
> [`docs/MAINTAINING.md`](docs/MAINTAINING.md) for what each workflow checks and how
> to read a failure. Contributors in **any field** (math, physics theory, music
> theory, bio, metacybernetics, prose, …) keep the repo green by following one small
> per-contribution manifest: start with [`docs/ANTI_DRIFT.md`](docs/ANTI_DRIFT.md)
> and [`CONTRIBUTING.md`](CONTRIBUTING.md).

This repository is the group's permanent archive for mathematics research. Each paper is stored alongside the data, figures, and code needed to reproduce its results, and every tagged release is mirrored to [Zenodo](https://zenodo.org/) with a citable DOI.

## Purpose

We publish mathematics research that comes paired with computational products — proofs, experiments, datasets, and the programs that generate them. Keeping each paper and its supporting artifacts together in one versioned archive means every published result stays reproducible and independently citable.

## Repository layout

| Folder     | Contents |
|------------|----------|
| `papers/`  | One subfolder per paper — the PDF plus its LaTeX source. |
| `data/`    | Datasets and machine-readable inputs/outputs supporting the papers. |
| `figures/` | Figures, plots, and images used in the papers. |
| `code/`    | Scripts, notebooks, and programs that produce the results. |

Each folder has its own `README.md` describing local conventions.

## How the Zenodo pipeline works

This repository is enabled in Zenodo, so archiving is automatic:

1. A maintainer publishes a **GitHub Release** (from a git tag).
2. Zenodo receives the release webhook, downloads a snapshot of the repository at that tag, and creates a new deposition.
3. Zenodo mints a **version DOI** for that release, plus a **concept DOI** that always resolves to the latest version.
4. Deposition metadata (title, authors, license, keywords, description) is read from [`.zenodo.json`](.zenodo.json).

No secrets or tokens live in this repository — the GitHub ↔ Zenodo link is configured once in the Zenodo account settings.

The included GitHub Actions workflow ([`.github/workflows/release.yml`](.github/workflows/release.yml)) validates `.zenodo.json` and `CITATION.cff` on every release — and on pull requests that touch them — so a malformed metadata file can't quietly break a deposition.

## How to add a paper

1. Copy [`papers/_TEMPLATE/`](papers/_TEMPLATE) to a new folder named `YYYY-MM-shortname` (year and month of the paper, plus a short slug), e.g. `papers/2026-07-spectral-gap`.
2. Add the full **LaTeX source** and compiled **PDF**, and fill in **`paper.cff`** (title, author, date). Set the LaTeX `\author{}` and the `paper.cff` author to your own attribution, and add any new citations to the consolidated [`papers/references.bib`](papers/references.bib).
3. Put supporting artifacts in the shared top-level folders: datasets in `data/`, figures in `figures/`, code in `code/`. Reference them from the paper or its folder as needed.
4. Open a pull request. CI will confirm the metadata files still validate before you merge.

Authorship is per-member and umbrella-org — follow **[`docs/ATTRIBUTION.md`](docs/ATTRIBUTION.md)** (members table + step-by-step workflow).

## How to cut a release

1. Make sure `main` is up to date and the metadata files (`.zenodo.json`, `CITATION.cff`) reflect what you're archiving.
2. Tag the commit and push the tag:
   ```bash
   git tag -a v1.0.0 -m "Describe this release"
   git push origin v1.0.0
   ```
3. On GitHub, go to **Releases → Draft a new release**, choose the tag, add a title and notes, and click **Publish release**.
4. Zenodo archives the release automatically and assigns a DOI. Once it appears, copy the DOI badge from Zenodo into this README or the relevant paper folder.

> Use [semantic versioning](https://semver.org/) for tags (`vMAJOR.MINOR.PATCH`). Cut a new release whenever you add or substantially revise an archived paper.

## License

Released under the [Creative Commons Attribution 4.0 International](LICENSE) (CC-BY-4.0) license — you may share and adapt the material with attribution. See [`CITATION.cff`](CITATION.cff) for how to cite this archive.
