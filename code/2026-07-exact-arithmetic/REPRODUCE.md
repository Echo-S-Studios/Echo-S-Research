# REPRODUCE

Fifteen files: fourteen check scripts, one harness. Everything is stdlib + SymPy + mpmath; no
network, no data files, no fixtures. Each script is self-generating and prints a `SUMMARY` line.

## Environment

| component | pinned | known-compatible |
|---|---|---|
| CPython | 3.12.3 | 3.12 – 3.13 |
| SymPy | 1.14.0 | 1.13 – 1.14 |
| mpmath | 1.3.0 | 1.3 |

```bash
pip install 'sympy>=1.13,<1.15' 'mpmath~=1.3'
```

The floor is not decorative. SymPy changed the shape of `Poly.intervals(all=True)` from a flat
list to a `(reals, complexes)` pair, which silently broke the census's disk certificate.
`lesson6_census.py` now tolerates both shapes — but nothing older than SymPy 1.12 has been
exercised, so don't assume it.

## Run it

```bash
python3 run_all.py --list      # environment, per-script runtime, what each script decides
python3 run_all.py --verify    # checksums only, against SHA256SUMS
python3 run_all.py             # quick: whole suite, census on a slice     ~50-90 s
python3 run_all.py --full      # everything, including the published census  ~3.5 min
```

Exit status is 0 only if every checksum matches, every script exits 0, and every `SUMMARY`
reports no failures.

## What `--quick` does and does not claim

Quick mode runs `lesson6_census.py --quick`, which exercises every code path over the `a = -1`
slice of the coefficient box in about 12 seconds instead of about two minutes. It says so at the
top of its own output and **does not assert** the published tallies
(625 / 50 / 638 / 1318 / 411 → 83 Pisot, patterns 0 / 16 / 67), nor re-derive the D₅ erratum.
Only `--full` does that. A smoke run that quietly reported the same badge as a full run would be
the exact defect this bundle exists to prevent.

## Expected results

| script | checks | quick | full |
|---|---|---|---|
| `lesson1_checks.py` | 16/16 | 1 s | 1 s |
| `lesson2_checks.py` | 32/32 | 2 s | 2 s |
| `lesson3_checks.py` | 33/33 | 2 s | 2 s |
| `lesson4_checks.py` | 25/25 | 2 s | 2 s |
| `lesson5_checks.py` | 23/23 | 1 s | 1 s |
| `lesson6_checks.py` | 20/20 | 3 s | 3 s |
| `lesson6_census.py` | 6/6 | 12 s | 130 s |
| `lesson7_checks.py` | 39/39 | 2 s | 2 s |
| `lesson8_checks.py` | 20/20 | 12 s | 12 s |
| `lesson9_checks.py` | 32/32 | 3 s | 3 s |
| `lesson10_checks.py` | 26/26 | 3 s | 3 s |
| `revision_checks.py` | 13/13 | 1 s | 1 s |
| `synthesis_checks.py` | 21/21 | 1 s | 1 s |
|  `certificates.py` | 22 executable + 1 cited | 45 s | 45 s |

354 numbered checks across the thirteen lessons. The census (6), the two revision passes (13, 21)
and the certificates (22 + 1) are separate artifacts and are **not** folded into that number.

## Three repairs made to this bundle after the fourth audit

Recorded because a bundle that hides its own repairs is not an audit trail.

1. **`lesson6_checks.py` carried a dead census section** — superseded during authoring by
   `lesson6_census.py`, never removed, and it crashed on load under current SymPy. Excised. The
   file is checks A–C, which are the 20 the Lesson 6 badge names.
2. **`revision_checks.py` hard-coded `/home/claude`** and could not run outside the authoring
   container. Now resolves siblings via `Path(__file__).resolve().parent`.
3. **`lesson6_census.py` check D3 asserted the source note's claim** — `Rat° irreducible on all
   67` — and therefore *failed* when run, because the truth is 66. **That failure was the D₅
   finding**, and reporting the artifact as "6/6" had hidden a finding inside a passing badge. D3
   now asserts the corrected counts (squarefree 67/67, irreducible 66/67, exception named
   `x^5 - x^3 - 2x^2 - 2x - 1`), so the finding lives in the check and the 6/6 is honest.

## On `certificates.py`

It reports **24 executable certificates + 2 cited theorems**, and the two kinds are counted
separately on purpose. Cramér–Rao is imported, not executed: `C15` verifies the *applicability*
conditions for this family (fixed 7-atom support independent of β, score mean zero — symbolically,
in seven free costs), and `T1` cites the bound itself. `T2` cites Rao (1945) for the square-root
isometry onto the radius-2 sphere — classical ambient background, not corpus content. An earlier printing counted a tautology
(`1/(nI) - 1/(nI) == 0`) as a passing check; that entry is withdrawn.

 `C22`/`C23` certify the narrowest unconditional form of the exclusion: the forced band is
(1, μ_S), not (1, φ), and Lehmer's number lies strictly inside it — so the headline conclusion
touches no plausible tag.

`certificates.py` is a **selected headline recheck** of Lessons 11–13 plus every third- and
fourth-pass repair — 22 checks against the 88 in the original authoring scripts, which are
unrecoverable. Stronger per check, weaker in coverage. It does not supersede them.
