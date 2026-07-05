# Anti-drift: keep the repo green (any field)

"Drift" just means the repo has quietly gone **inconsistent** — a file says one
thing, the code that made it says another; an author isn't on the team list; a
folder is half-there. Two GitHub checks catch that automatically:

- **structure** — is your contribution shaped like a contribution?
- **drift** — does whatever you *said* is reproducible still reproduce?

You do **not** need to do math, use LaTeX, or write Python to pass. Everything is
driven by one small file you own: **`papers/<your-folder>/contribution.yml`**. If
you write pure prose, you say so in that file and both checks stay green. This page
is the whole game plan.

---

## (a) The minimal checklist to add a compliant contribution

1. **Be on the team list.** Your GitHub handle must be in
   [`.github/members.yml`](../.github/members.yml). New teammate? Add one line.

2. **Scaffold your folder** (this writes the manifest for you):

   ```bash
   py scripts/new_contribution.py --member <your-handle> --shortname my-idea \
       --date 2026-08 --domain <your-field>
   ```

   `--domain` is free text: `math`, `engineering`, `physics-theory`,
   `music-theory`, `bio`, `metacybernetics`, `other` — whatever fits. Add
   `--trees code,data,tests,figures` if you want those folders scaffolded too
   (leave it off for a prose-only piece).

3. **Write your paper** in `papers/<folder>/` (the `.tex` stub, or delete it and
   drop in a `.md`/`.pdf` — just make sure every file you keep is listed under
   `artifacts` in the manifest).

4. **Decide reproduction** — edit `checks` in your `contribution.yml`:
   - Nothing to run (pure theory/prose)? Leave **`checks: none`** and write an
     honest **`reason`**. ✅ This passes.
   - Have a script/notebook/test that should keep working? List it as a check
     (see the examples below).

5. **Pre-flight before pushing** — run the exact checks CI runs, locally:

   ```bash
   py scripts/check.py <your-folder>     # or: py scripts/check.py  (everything)
   ```

   Green here = green in CI. Fix anything it prints, then commit and open a PR.

That's it. The rules you're satisfying: manifest present & valid, your handle
known, declared files exist, you only have the `code/ data/ tests/ figures/`
folders you actually use, folder named `YYYY-MM-shortname`.

---

## The manifest, in one screen

```yaml
shortname: 2026-08-my-idea      # must equal the folder name
member: your-handle             # must be in .github/members.yml
title: "My Idea"
domain: physics-theory          # free text
artifacts:                      # every file/dir you own — all must exist
  - papers/2026-08-my-idea/my-idea.tex
  - papers/2026-08-my-idea/paper.cff

checks: none                    # EITHER this + a reason (pure prose/theory) ...
reason: "Conceptual paper; nothing to compute."
# checks:                       # ... OR a list of named checks:
#   - name: <label>
#     run: "<command>"          # regenerates outputs (optional)
#     produces: [<file-or-dir>] # diffed against the committed copy (needs `run`)
#     verify: "<command>"       # must exit 0 — tests / lint / schema / notebook (optional)
```

`{PY}` in any command is replaced with the Python interpreter, so the same manifest
works on your machine and on CI. For other tools, just write the command
(`node analyze.js`, `Rscript check.R`, `bash run.sh`) — and if a dependency isn't
already on the runner, install it inside the `run` command (e.g.
`run: "pip install scipy && {PY} sim.py"`).

---

## (b) Three worked examples

### 1. Physics theory — a numerical simulation check

Folder `2026-08-toy-cosmology`, scaffolded with `--domain physics-theory --trees code,data`.
`code/2026-08-toy-cosmology/simulate.py` integrates a model and writes a small,
rounded results table to `data/2026-08-toy-cosmology/spectrum.csv`.

```yaml
shortname: 2026-08-toy-cosmology
member: eMKa7-tY
title: "A Toy Cosmology and Its Power Spectrum"
domain: physics-theory
artifacts:
  - papers/2026-08-toy-cosmology/toy-cosmology.pdf
  - papers/2026-08-toy-cosmology/paper.cff
  - code/2026-08-toy-cosmology
  - data/2026-08-toy-cosmology
checks:
  - name: simulation
    run: "{PY} code/2026-08-toy-cosmology/simulate.py"
    produces:
      - data/2026-08-toy-cosmology/spectrum.csv
```

`simulate.py` seeds its RNG and rounds outputs, so the CSV is identical every run
→ the diff stays clean.

### 2. Music theory — an analysis script plus a data check

Folder `2026-08-cadence-stats`, scaffolded with `--domain music-theory --trees code,data,tests`.
A script tallies cadence types across a corpus into JSON; a small test asserts the
tallies satisfy a claimed invariant.

```yaml
shortname: 2026-08-cadence-stats
member: sKiDaGgAbAtEe
title: "Cadence Frequencies in the Corpus"
domain: music-theory
artifacts:
  - papers/2026-08-cadence-stats/cadence-stats.md
  - papers/2026-08-cadence-stats/paper.cff
  - code/2026-08-cadence-stats
  - data/2026-08-cadence-stats
  - tests/2026-08-cadence-stats
checks:
  - name: analysis
    run: "{PY} code/2026-08-cadence-stats/tally.py"
    produces:
      - data/2026-08-cadence-stats/cadences.json
  - name: invariant
    verify: "{PY} -m pytest tests/2026-08-cadence-stats -q"
```

### 3. Pure prose (bio / metacybernetics) — no reproduction

Folder `2026-08-metacybernetics`, scaffolded with `--domain metacybernetics` (no
`--trees`). It's a conceptual essay — there is nothing to run, and that is fine:

```yaml
shortname: 2026-08-metacybernetics
member: eMKa7-tY
title: "Metacybernetics: A Framework for Living Regulation"
domain: metacybernetics
artifacts:
  - papers/2026-08-metacybernetics/metacybernetics.md
  - papers/2026-08-metacybernetics/paper.cff
checks: none
reason: "Conceptual/theoretical essay; no data, code, or computation to reproduce."
```

**structure** passes (a valid manifest with attribution and existing files);
**drift** logs `skip 2026-08-metacybernetics: no-repro (Conceptual/theoretical …)`
and never fails it. No `code/`, `data/`, `tests/`, or `figures/` folders are
required or wanted.

---

## (c) How to declare checks so they're cheap to satisfy

The only way a `run`+`produces` check fails (other than a real bug) is if your
output isn't **bit-for-bit reproducible**. Keep it that way:

- **Emit small, deterministic outputs.** A rounded CSV/JSON of the numbers that
  matter beats a huge dump. Round floats (`round(x, 6)`); sort collections before
  writing; and normalize signed zero (`x + 0.0`) — `-0.0` vs `0.0` differs between
  machines.
- **Pin every seed.** `random.seed(0)`, `np.random.seed(0)`, and the equivalent in
  any language. No timestamps, no `set()` iteration order, no wall-clock, no
  machine paths in the output.
- **Keep producers fast.** CI reruns them on every push. Seconds, not minutes.
  Sample or cap expensive sweeps and note it.
- **Don't diff binaries.** PDFs, PNGs, and audio rarely reproduce byte-for-byte.
  Put them in `artifacts` (so they're tracked) but don't list them under
  `produces`. Diff the *numbers* that generated them instead.
- **Pin your dependencies.** If your check needs a library, install a pinned
  version in the `run` command so an upstream release can't change your output.
- **When in doubt, `checks: none`.** An honest no-repro with a reason is always
  better than a flaky check. You can add reproduction later.

Then run `py scripts/check.py <your-folder>` and push with confidence. See
[`docs/MAINTAINING.md`](MAINTAINING.md) for what the maintainers see, and
[`CONTRIBUTING.md`](../CONTRIBUTING.md) for the end-to-end flow.
